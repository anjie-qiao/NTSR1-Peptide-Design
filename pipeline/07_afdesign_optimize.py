#!/usr/bin/env python3
"""Create cumulative one-, two-, and three-substitution AfDesign variants."""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import math
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from pepmlm_pipeline_common import (
    CANONICAL_AAS,
    atomic_write_csv,
    atomic_write_json,
    atomic_write_text,
    canonical_json_hash,
    sequence_properties,
    sha256_file,
)


AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}

REQUIRED_INPUT_FIELDS = {
    "candidate_id", "sequence", "length", "prediction_status", "structure_pass",
    "structure_pdb",
}

OUTPUT_FIELDS = [
    "candidate_id", "parent_candidate_id", "parent_sequence", "sequence", "length",
    "mutation_count", "mutation_1", "mutation_2", "mutation_3",
    "afdesign_parent_loss", "afdesign_round_1_loss", "afdesign_round_2_loss",
    "afdesign_round_3_loss", "afdesign_final_loss", "afdesign_loss_delta",
    "sequence_filter_pass", "sequence_filter_reasons", "net_charge_ph7_4",
    "hydrophobic_fraction", "max_hydrophobic_run", "max_homopolymer_run",
    "normalized_entropy", "cysteine_count", "source_structure_pdb",
]

MUTATION_PATTERN = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def resolve_path(project_root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def load_config(config_path: Path) -> tuple[dict, Path]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    required = {
        "input_csv", "output_dir", "afdesign_data_dir", "expected_input_count",
        "receptor_chain_id", "peptide_chain_id", "mutation_rounds",
        "proposals_per_round", "use_multimer", "num_recycles", "num_models",
        "model_names", "dropout", "base_seed", "loss_weights", "hard_filters",
    }
    missing = required - set(config)
    if missing:
        raise ValueError(f"Missing config keys: {sorted(missing)}")
    if int(config["mutation_rounds"]) != 3:
        raise ValueError("mutation_rounds must be exactly 3 for this fixed protocol")
    if int(config["proposals_per_round"]) < 1:
        raise ValueError("proposals_per_round must be positive")
    if int(config["num_recycles"]) < 0 or int(config["num_models"]) < 1:
        raise ValueError("num_recycles must be non-negative and num_models must be positive")
    if not isinstance(config["model_names"], list) or not config["model_names"]:
        raise ValueError("model_names must be a non-empty list")
    if int(config["num_models"]) > len(config["model_names"]):
        raise ValueError("num_models cannot exceed the number of loaded model_names")
    if config["receptor_chain_id"] == config["peptide_chain_id"]:
        raise ValueError("Receptor and peptide chain IDs must differ")
    expected_weights = {"dgram_cce", "pae", "plddt", "i_pae", "con", "i_con"}
    if set(config["loss_weights"]) != expected_weights:
        raise ValueError(f"loss_weights must contain exactly: {sorted(expected_weights)}")
    for name, value in config["loss_weights"].items():
        if float(value) < 0 or not math.isfinite(float(value)):
            raise ValueError(f"loss_weights.{name} must be finite and non-negative")
    expected_filters = {
        "min_net_charge", "max_net_charge", "max_hydrophobic_fraction",
        "max_hydrophobic_run", "max_homopolymer_run", "min_normalized_entropy",
        "max_cysteine_count",
    }
    if set(config["hard_filters"]) != expected_filters:
        raise ValueError(f"hard_filters must contain exactly: {sorted(expected_filters)}")
    return config, config_path.resolve().parent.parent


def _truthy(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def read_parents(path: Path, project_root: Path) -> list[dict]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_INPUT_FIELDS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Input CSV is missing fields: {sorted(missing)}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise ValueError("Input CSV is empty")

    seen_ids: set[str] = set()
    seen_sequences: set[str] = set()
    for row in rows:
        candidate_id = row["candidate_id"].strip()
        sequence = row["sequence"].strip().upper()
        if not candidate_id or candidate_id in seen_ids:
            raise ValueError(f"Missing or duplicate candidate ID: {candidate_id!r}")
        if not sequence or sequence in seen_sequences or set(sequence) - set(CANONICAL_AAS):
            raise ValueError(f"Invalid or duplicate parent sequence: {candidate_id}")
        if int(row["length"]) != len(sequence):
            raise ValueError(f"Length mismatch for {candidate_id}")
        if row["prediction_status"].strip() != "ok" or not _truthy(row["structure_pass"]):
            raise ValueError(f"Input contains a non-passing structure: {candidate_id}")
        pdb_path = resolve_path(project_root, row["structure_pdb"])
        if not pdb_path.is_file():
            raise FileNotFoundError(f"Structure PDB does not exist for {candidate_id}: {pdb_path}")
        row["candidate_id"] = candidate_id
        row["sequence"] = sequence
        row["_length"] = len(sequence)
        row["_pdb_path"] = pdb_path
        seen_ids.add(candidate_id)
        seen_sequences.add(sequence)
    return rows


def read_pdb_chain_sequence(path: Path, chain_id: str) -> str:
    residues: list[tuple[tuple[str, str], str]] = []
    observed: set[tuple[str, str]] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  ") or len(line) < 27 or line[21].strip() != chain_id:
                continue
            altloc = line[16]
            if altloc not in {" ", "A"}:
                continue
            residue_key = (line[22:26].strip(), line[26].strip())
            if residue_key in observed:
                continue
            residue_name = line[17:20].strip().upper()
            if residue_name not in AA3_TO_AA1:
                raise ValueError(f"Unsupported residue {residue_name!r} in chain {chain_id} of {path}")
            residues.append((residue_key, AA3_TO_AA1[residue_name]))
            observed.add(residue_key)
    if not residues:
        raise ValueError(f"No protein residues found for chain {chain_id!r} in {path}")
    return "".join(aa for _, aa in residues)


def validate_parent_structures(rows: list[dict], receptor_chain: str, peptide_chain: str) -> None:
    for index, row in enumerate(rows, start=1):
        peptide_sequence = read_pdb_chain_sequence(row["_pdb_path"], peptide_chain)
        if peptide_sequence != row["sequence"]:
            raise ValueError(
                f"Peptide sequence mismatch for {row['candidate_id']}: "
                f"CSV={row['sequence']} PDB={peptide_sequence}"
            )
        read_pdb_chain_sequence(row["_pdb_path"], receptor_chain)
        if index % 100 == 0 or index == len(rows):
            print(f"Validated parent PDB chains: {index}/{len(rows)}", flush=True)


def hard_filter_failures(sequence: str, thresholds: dict) -> tuple[dict, list[str]]:
    properties = sequence_properties(sequence)
    failures = []
    if properties["net_charge_ph7_4"] < float(thresholds["min_net_charge"]):
        failures.append("net_charge_too_low")
    if properties["net_charge_ph7_4"] > float(thresholds["max_net_charge"]):
        failures.append("net_charge_too_high")
    if properties["hydrophobic_fraction"] > float(thresholds["max_hydrophobic_fraction"]):
        failures.append("hydrophobic_fraction_too_high")
    if properties["max_hydrophobic_run"] > int(thresholds["max_hydrophobic_run"]):
        failures.append("hydrophobic_run_too_long")
    if properties["max_homopolymer_run"] > int(thresholds["max_homopolymer_run"]):
        failures.append("homopolymer_run_too_long")
    if properties["normalized_entropy"] < float(thresholds["min_normalized_entropy"]):
        failures.append("sequence_complexity_too_low")
    if properties["cysteine_count"] > int(thresholds["max_cysteine_count"]):
        failures.append("too_many_cysteines")
    return properties, failures


def load_afdesign_runtime():
    try:
        from colabdesign import clear_mem, mk_afdesign_model
        from colabdesign.af.alphafold.common import residue_constants
    except ImportError as exc:
        raise RuntimeError(
            "ColabDesign/AfDesign is not installed in this environment. "
            "Install it in a separate afdesign conda environment."
        ) from exc
    alphabet = "".join(residue_constants.restypes)
    if set(alphabet) != set(CANONICAL_AAS) or len(alphabet) != 20:
        raise RuntimeError(f"Unexpected AfDesign amino-acid alphabet: {alphabet}")
    return mk_afdesign_model, clear_mem, alphabet


def validate_afdesign_parameters(data_dir: Path, model_names: list[str]) -> None:
    for model_name in model_names:
        alternatives = [
            data_dir / "params" / f"params_{model_name}.npz",
            data_dir / f"params_{model_name}.npz",
            data_dir / "params" / f"{model_name}.npz",
            data_dir / f"{model_name}.npz",
        ]
        if not any(path.is_file() for path in alternatives):
            expected = ", ".join(str(path) for path in alternatives[:2])
            raise FileNotFoundError(
                f"AfDesign weights for {model_name!r} were not found; expected e.g. {expected}"
            )


def scalar_loss(aux: dict) -> float:
    value = aux.get("loss", aux.get("log", {}).get("loss"))
    if value is None:
        raise RuntimeError("AfDesign output did not contain a loss")
    loss = float(np.asarray(value))
    if not math.isfinite(loss):
        raise RuntimeError(f"AfDesign returned a non-finite loss: {loss}")
    return loss


def gradient_matrix(aux: dict, peptide_length: int) -> np.ndarray:
    try:
        gradient = np.asarray(aux["grad"]["seq"], dtype=float)
    except (KeyError, TypeError) as exc:
        raise RuntimeError("AfDesign output did not contain sequence gradients") from exc
    while gradient.ndim > 2:
        gradient = gradient.mean(axis=0)
    if gradient.shape != (peptide_length, 20):
        raise RuntimeError(
            f"Unexpected AfDesign sequence-gradient shape {gradient.shape}; "
            f"expected ({peptide_length}, 20)"
        )
    if not np.isfinite(gradient).all():
        raise RuntimeError("AfDesign returned non-finite sequence gradients")
    return gradient


def propose_mutants(
    sequence: str,
    gradient: np.ndarray,
    alphabet: str,
    used_positions: set[int],
    proposal_count: int,
) -> list[tuple[float, int, str, str]]:
    aa_to_index = {aa: index for index, aa in enumerate(alphabet)}
    proposals: list[tuple[float, int, str, str]] = []
    for position, current_aa in enumerate(sequence):
        if position in used_positions:
            continue
        current_index = aa_to_index[current_aa]
        for replacement_index, replacement_aa in enumerate(alphabet):
            if replacement_aa == current_aa:
                continue
            estimated_delta = gradient[position, replacement_index] - gradient[position, current_index]
            mutant = sequence[:position] + replacement_aa + sequence[position + 1:]
            proposals.append((float(estimated_delta), position, replacement_aa, mutant))
    proposals.sort(key=lambda item: (item[0], item[1], item[2]))
    if len(proposals) < proposal_count:
        raise RuntimeError(
            f"Only {len(proposals)} substitutions remain, fewer than proposals_per_round={proposal_count}"
        )
    return proposals[:proposal_count]


def predict_loss(model, sequence: str, config: dict) -> float:
    aux = model.predict(
        seq=sequence,
        num_models=int(config["num_models"]),
        num_recycles=int(config["num_recycles"]),
        sample_models=False,
        dropout=bool(config["dropout"]),
        hard=True,
        soft=False,
        return_aux=True,
        verbose=False,
    )
    return scalar_loss(aux)


def configure_model(model, config: dict) -> None:
    weights = {name: float(value) for name, value in config["loss_weights"].items()}
    model.set_weights(weights, set_defaults=True)
    model.set_opt(
        num_recycles=int(config["num_recycles"]),
        num_models=int(config["num_models"]),
        sample_models=False,
        dropout=bool(config["dropout"]),
        soft=1.0,
        hard=0.0,
        temp=1.0,
        set_defaults=True,
    )
    model.set_opt("i_con", num=1, num_pos=5, set_defaults=True)


def optimize_parent(model, row: dict, config: dict, alphabet: str, seed: int) -> dict:
    parent_sequence = row["sequence"]
    current_sequence = parent_sequence
    used_positions: set[int] = set()
    mutations: list[str] = []
    round_losses: list[float] = []

    model.restart(seed=seed, seq=current_sequence)
    parent_loss = predict_loss(model, current_sequence, config)

    for round_index in range(1, 4):
        model.restart(seed=seed + round_index, seq=current_sequence)
        model.set_opt(
            soft=1.0,
            hard=0.0,
            temp=1.0,
            dropout=bool(config["dropout"]),
        )
        aux = model.run(
            num_recycles=int(config["num_recycles"]),
            num_models=int(config["num_models"]),
            sample_models=False,
            backprop=True,
            return_aux=True,
        )
        gradient = gradient_matrix(aux, len(current_sequence))
        proposals = propose_mutants(
            current_sequence,
            gradient,
            alphabet,
            used_positions,
            int(config["proposals_per_round"]),
        )
        evaluated = []
        for estimated_delta, position, replacement_aa, mutant in proposals:
            loss = predict_loss(model, mutant, config)
            evaluated.append((loss, estimated_delta, position, replacement_aa, mutant))
        loss, _, position, replacement_aa, chosen_sequence = min(
            evaluated, key=lambda item: (item[0], item[1], item[2], item[3])
        )
        original_aa = current_sequence[position]
        mutations.append(f"{original_aa}{position + 1}{replacement_aa}")
        used_positions.add(position)
        current_sequence = chosen_sequence
        round_losses.append(loss)

    if len(mutations) != 3 or len(used_positions) != 3:
        raise RuntimeError("The fixed protocol did not produce exactly three distinct substitutions")
    hamming = sum(a != b for a, b in zip(parent_sequence, current_sequence))
    if hamming != 3:
        raise RuntimeError(f"Expected Hamming distance 3, observed {hamming}")

    properties, failures = hard_filter_failures(current_sequence, config["hard_filters"])
    return {
        "candidate_id": f"{row['candidate_id']}__AFD3",
        "parent_candidate_id": row["candidate_id"],
        "parent_sequence": parent_sequence,
        "sequence": current_sequence,
        "length": len(current_sequence),
        "mutation_count": 3,
        "mutation_1": mutations[0],
        "mutation_2": mutations[1],
        "mutation_3": mutations[2],
        "afdesign_parent_loss": parent_loss,
        "afdesign_round_1_loss": round_losses[0],
        "afdesign_round_2_loss": round_losses[1],
        "afdesign_round_3_loss": round_losses[2],
        "afdesign_final_loss": round_losses[2],
        "afdesign_loss_delta": round_losses[2] - parent_loss,
        "sequence_filter_pass": int(not failures),
        "sequence_filter_reasons": ";".join(failures),
        **properties,
        "source_structure_pdb": str(row["_pdb_path"]),
    }


def expand_cumulative_variants(result: dict, thresholds: dict) -> list[dict]:
    """Reconstruct and score the cumulative M1, M2 and M3 variants.

    Restart shards retain the compact final result and the three ordered
    substitutions. Reconstructing intermediate sequences here keeps completed
    shards reusable and avoids repeating any GPU optimization.
    """
    parent_sequence = str(result["parent_sequence"])
    parent_id = str(result["parent_candidate_id"])
    current_sequence = parent_sequence
    used_positions: set[int] = set()
    variants: list[dict] = []

    parent_loss = float(result["afdesign_parent_loss"])
    mutations = [str(result[f"mutation_{index}"]) for index in range(1, 4)]
    round_losses = [float(result[f"afdesign_round_{index}_loss"]) for index in range(1, 4)]

    for round_index, (mutation, round_loss) in enumerate(
        zip(mutations, round_losses), start=1
    ):
        match = MUTATION_PATTERN.fullmatch(mutation)
        if match is None:
            raise ValueError(f"Invalid mutation notation for {parent_id}: {mutation!r}")
        original_aa, position_text, replacement_aa = match.groups()
        position = int(position_text) - 1
        if position < 0 or position >= len(current_sequence):
            raise ValueError(f"Mutation position is out of range for {parent_id}: {mutation}")
        if position in used_positions:
            raise ValueError(f"Mutation position was reused for {parent_id}: {mutation}")
        if current_sequence[position] != original_aa:
            raise ValueError(
                f"Mutation source residue mismatch for {parent_id}: {mutation}; "
                f"observed {current_sequence[position]}"
            )
        if replacement_aa == original_aa or replacement_aa not in CANONICAL_AAS:
            raise ValueError(f"Invalid replacement residue for {parent_id}: {mutation}")

        current_sequence = (
            current_sequence[:position] + replacement_aa + current_sequence[position + 1:]
        )
        used_positions.add(position)
        properties, failures = hard_filter_failures(current_sequence, thresholds)

        variants.append({
            "candidate_id": f"{parent_id}__AFD{round_index}",
            "parent_candidate_id": parent_id,
            "parent_sequence": parent_sequence,
            "sequence": current_sequence,
            "length": len(current_sequence),
            "mutation_count": round_index,
            "mutation_1": mutations[0],
            "mutation_2": mutations[1] if round_index >= 2 else "",
            "mutation_3": mutations[2] if round_index >= 3 else "",
            "afdesign_parent_loss": parent_loss,
            "afdesign_round_1_loss": round_losses[0],
            "afdesign_round_2_loss": round_losses[1] if round_index >= 2 else "",
            "afdesign_round_3_loss": round_losses[2] if round_index >= 3 else "",
            "afdesign_final_loss": round_loss,
            "afdesign_loss_delta": round_loss - parent_loss,
            "sequence_filter_pass": int(not failures),
            "sequence_filter_reasons": ";".join(failures),
            **properties,
            "source_structure_pdb": str(result["source_structure_pdb"]),
        })

    if current_sequence != str(result["sequence"]):
        raise ValueError(
            f"Reconstructed three-substitution sequence does not match shard for {parent_id}"
        )
    return variants


def _serialize(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.8f}"
    if isinstance(value, (int, np.integer)):
        return int(value)
    return value


def mutation_summary(row: dict) -> str:
    return ",".join(
        str(row[f"mutation_{index}"])
        for index in range(1, 4)
        if row[f"mutation_{index}"]
    )


def shard_path(shard_dir: Path, candidate_id: str) -> Path:
    return shard_dir / f"{candidate_id}.json"


def candidate_seed(base_seed: int, candidate_id: str) -> int:
    digest = hashlib.sha256(candidate_id.encode("utf-8")).digest()
    return (int(base_seed) + int.from_bytes(digest[:4], "big")) % (2 ** 31 - 1)


def load_valid_shard(path: Path, parent: dict, fingerprint: str) -> dict | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        result = payload["result"]
        if payload["config_fingerprint"] != fingerprint:
            return None
        if payload["parent_candidate_id"] != parent["candidate_id"]:
            return None
        if payload["parent_sequence"] != parent["sequence"]:
            return None
        if result["mutation_count"] != 3:
            return None
        if sum(a != b for a, b in zip(parent["sequence"], result["sequence"])) != 3:
            return None
        return result
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        return None


def write_shard(path: Path, parent: dict, fingerprint: str, result: dict) -> None:
    atomic_write_json(path, {
        "schema_version": 1,
        "stage": "07_afdesign_optimize",
        "config_fingerprint": fingerprint,
        "parent_candidate_id": parent["candidate_id"],
        "parent_sequence": parent["sequence"],
        "result": {key: _serialize(value) for key, value in result.items()},
    })


def run(config_path: Path, plan_only: bool = False) -> Path | None:
    config, project_root = load_config(config_path)
    input_path = resolve_path(project_root, config["input_csv"])
    output_dir = resolve_path(project_root, config["output_dir"])
    data_dir = resolve_path(project_root, config["afdesign_data_dir"])
    parents = read_parents(input_path, project_root)
    expected = config["expected_input_count"]
    if expected is not None and len(parents) != int(expected):
        raise ValueError(f"Expected {expected} structure-pass parents, observed {len(parents)}")
    validate_parent_structures(
        parents, str(config["receptor_chain_id"]), str(config["peptide_chain_id"])
    )
    plan = {
        "input_candidates": len(parents),
        "mutation_rounds": 3,
        "proposals_per_round": int(config["proposals_per_round"]),
        "variants_per_parent": 3,
        "maximum_new_variants": len(parents) * 3,
        "afdesign_data_dir": str(data_dir),
        "length_counts": {
            str(length): sum(row["_length"] == length for row in parents)
            for length in sorted({row["_length"] for row in parents})
        },
    }
    if plan_only:
        print(json.dumps(plan, indent=2, ensure_ascii=False))
        return None
    if not data_dir.is_dir():
        raise FileNotFoundError(f"AfDesign parameter directory does not exist: {data_dir}")
    validate_afdesign_parameters(data_dir, [str(name) for name in config["model_names"]])

    output_dir.mkdir(parents=True, exist_ok=True)
    shard_dir = output_dir / "shards"
    shard_dir.mkdir(parents=True, exist_ok=True)
    fingerprint = canonical_json_hash(config)
    mk_afdesign_model, clear_mem, alphabet = load_afdesign_runtime()

    results_by_parent: dict[str, dict] = {}
    pending_by_length: dict[int, list[dict]] = defaultdict(list)
    for parent in parents:
        cached = load_valid_shard(shard_path(shard_dir, parent["candidate_id"]), parent, fingerprint)
        if cached is None:
            pending_by_length[parent["_length"]].append(parent)
        else:
            results_by_parent[parent["candidate_id"]] = cached

    completed = len(results_by_parent)
    if completed:
        print(f"Resuming with {completed}/{len(parents)} valid candidate shards", flush=True)
    errors: list[dict] = []
    for length in sorted(pending_by_length):
        group = pending_by_length[length]
        clear_mem()
        gc.collect()
        model = mk_afdesign_model(
            protocol="binder",
            use_multimer=bool(config["use_multimer"]),
            num_recycles=int(config["num_recycles"]),
            recycle_mode="last",
            num_models=int(config["num_models"]),
            sample_models=False,
            model_names=[str(name) for name in config["model_names"]],
            data_dir=str(data_dir),
        )
        print(f"Optimizing length {length}: {len(group)} parents", flush=True)
        for parent in group:
            try:
                seed = candidate_seed(int(config["base_seed"]), parent["candidate_id"])
                model.prep_inputs(
                    pdb_filename=str(parent["_pdb_path"]),
                    target_chain=str(config["receptor_chain_id"]),
                    binder_chain=str(config["peptide_chain_id"]),
                    use_binder_template=False,
                    ignore_missing=True,
                    seed=seed,
                )
                configure_model(model, config)
                result = optimize_parent(
                    model,
                    parent,
                    config,
                    alphabet,
                    seed,
                )
                write_shard(
                    shard_path(shard_dir, parent["candidate_id"]), parent, fingerprint, result
                )
                results_by_parent[parent["candidate_id"]] = {
                    key: _serialize(value) for key, value in result.items()
                }
            except Exception as exc:
                errors.append({
                    "parent_candidate_id": parent["candidate_id"],
                    "error": f"{type(exc).__name__}: {exc}",
                })
                print(f"ERROR {parent['candidate_id']}: {type(exc).__name__}: {exc}", file=sys.stderr)
            completed += 1
            print(f"Stage 07 progress: {completed}/{len(parents)}", flush=True)
        del model
        clear_mem()
        gc.collect()

    ordered_parent_results = [
        results_by_parent[parent["candidate_id"]]
        for parent in parents
        if parent["candidate_id"] in results_by_parent
    ]
    ordered_results = [
        variant
        for result in ordered_parent_results
        for variant in expand_cumulative_variants(result, config["hard_filters"])
    ]
    all_path = output_dir / "afdesign_variants_all.csv"
    pass_path = output_dir / "afdesign_variants_pass.csv"
    fasta_path = output_dir / "afdesign_variants_pass.fasta"
    error_path = output_dir / "afdesign_errors.csv"
    atomic_write_csv(all_path, ordered_results, OUTPUT_FIELDS)
    pass_rows = [row for row in ordered_results if _truthy(row["sequence_filter_pass"])]
    atomic_write_csv(pass_path, pass_rows, OUTPUT_FIELDS)
    atomic_write_text(
        fasta_path,
        "".join(
            f">{row['candidate_id']} parent={row['parent_candidate_id']} mutations="
            f"{mutation_summary(row)}\n{row['sequence']}\n"
            for row in pass_rows
        ),
    )
    atomic_write_csv(error_path, errors, ["parent_candidate_id", "error"])

    manifest = {
        "schema_version": 1,
        "stage": "07_afdesign_optimize",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": str(config_path.resolve()),
        "config_fingerprint": fingerprint,
        "inputs": {
            "structure_pass_csv": str(input_path),
            "structure_pass_csv_sha256": sha256_file(input_path),
            "input_candidates": len(parents),
        },
        "protocol": {
            "mutation_rounds": 3,
            "variants_per_parent": 3,
            "proposals_per_round": int(config["proposals_per_round"]),
            "distinct_mutation_positions": True,
            "early_stopping": False,
            "keep_intermediate_variants": True,
            "loss_weights": config["loss_weights"],
        },
        "counts": {
            "input_parents": len(parents),
            "optimized_parents": len(ordered_parent_results),
            "generated_variants": len(ordered_results),
            "sequence_filter_pass_variants": len(pass_rows),
            "sequence_filter_fail_variants": len(ordered_results) - len(pass_rows),
            "by_mutation_count": {
                str(mutation_count): {
                    "generated": sum(
                        int(row["mutation_count"]) == mutation_count for row in ordered_results
                    ),
                    "sequence_filter_pass": sum(
                        int(row["mutation_count"]) == mutation_count
                        and _truthy(row["sequence_filter_pass"])
                        for row in ordered_results
                    ),
                }
                for mutation_count in (1, 2, 3)
            },
            "errors": len(errors),
        },
        "outputs": {
            "all_variants_csv": str(all_path),
            "sequence_pass_csv": str(pass_path),
            "sequence_pass_fasta": str(fasta_path),
            "errors_csv": str(error_path),
        },
    }
    atomic_write_json(output_dir / "afdesign_manifest.json", manifest)
    if errors:
        raise RuntimeError(
            f"AfDesign failed for {len(errors)} parents; rerun after reviewing {error_path}"
        )
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument(
        "--plan-only", action="store_true",
        help="Validate inputs and print the workload without importing AfDesign",
    )
    args = parser.parse_args()
    try:
        output_dir = run(args.config.resolve(), plan_only=args.plan_only)
    except Exception as exc:
        print(f"07_afdesign_optimize FAILED: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    if output_dir is not None:
        print(f"07_afdesign_optimize PASS: {output_dir}")


if __name__ == "__main__":
    main()
