#!/usr/bin/env python3
"""Prepare Boltz jobs for all sequence-pass cumulative AfDesign variants."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shlex
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from pepmlm_pipeline_common import (
    CANONICAL_AAS,
    atomic_write_csv,
    atomic_write_json,
    atomic_write_text,
    canonical_json_hash,
    read_single_fasta,
    sha256_file,
)


REQUIRED_VARIANT_FIELDS = {
    "candidate_id", "parent_candidate_id", "parent_sequence", "sequence", "length",
    "mutation_count", "mutation_1", "mutation_2", "mutation_3",
    "afdesign_parent_loss", "afdesign_round_1_loss", "afdesign_round_2_loss",
    "afdesign_round_3_loss", "afdesign_final_loss", "afdesign_loss_delta",
    "sequence_filter_pass", "sequence_filter_reasons", "source_structure_pdb",
}

TASK_FIELDS = [
    "job_id", "candidate_id", "parent_candidate_id", "parent_sequence", "sequence",
    "length", "mutation_count", "mutation_1", "mutation_2", "mutation_3",
    "afdesign_parent_loss", "afdesign_final_loss", "afdesign_loss_delta",
    "sequence_filter_pass", "sequence_filter_reasons", "source_structure_pdb",
    "boltz_yaml", "receptor_chain_id", "peptide_chain_id",
]

SAFE_ID_PATTERN = re.compile(r"^[A-Za-z0-9._-]+$")
MUTATION_PATTERN = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def resolve_path(project_root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def _truthy(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def load_config(config_path: Path) -> tuple[dict, Path]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    required = {
        "variants_csv", "afdesign_manifest", "target_fasta", "target_manifest",
        "screening_masks", "output_dir", "expected_variant_count",
        "receptor_chain_id", "peptide_chain_id",
        "receptor_msa", "peptide_msa", "boltz_prediction",
    }
    missing = required - set(config)
    if missing:
        raise ValueError(f"Missing config keys: {sorted(missing)}")
    if config["receptor_chain_id"] == config["peptide_chain_id"]:
        raise ValueError("Receptor and peptide chain IDs must differ")
    if int(config["expected_variant_count"]) < 1:
        raise ValueError("Expected variant count must be positive")
    if not isinstance(config["boltz_prediction"], dict):
        raise ValueError("boltz_prediction must be an object")
    return config, config_path.resolve().parent.parent


def parse_finite(row: dict, field: str, candidate_id: str) -> float:
    try:
        value = float(row[field])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {field} for {candidate_id}") from exc
    if not math.isfinite(value):
        raise ValueError(f"Non-finite {field} for {candidate_id}")
    return value


def read_variants(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_VARIANT_FIELDS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"AfDesign variants CSV is missing fields: {sorted(missing)}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise ValueError("AfDesign variants CSV is empty")

    seen_ids: set[str] = set()
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        candidate_id = row["candidate_id"].strip()
        parent_id = row["parent_candidate_id"].strip()
        parent_sequence = row["parent_sequence"].strip().upper()
        sequence = row["sequence"].strip().upper()
        if not candidate_id or candidate_id in seen_ids or not SAFE_ID_PATTERN.fullmatch(candidate_id):
            raise ValueError(f"Missing, duplicate or unsafe candidate ID: {candidate_id!r}")
        if not parent_id or not SAFE_ID_PATTERN.fullmatch(parent_id):
            raise ValueError(f"Missing or unsafe parent candidate ID for {candidate_id}")
        if not parent_sequence or set(parent_sequence) - set(CANONICAL_AAS):
            raise ValueError(f"Invalid parent sequence for {candidate_id}")
        if not sequence or set(sequence) - set(CANONICAL_AAS):
            raise ValueError(f"Invalid optimized sequence for {candidate_id}")
        try:
            length = int(row["length"])
            mutation_count = int(row["mutation_count"])
        except ValueError as exc:
            raise ValueError(f"Invalid length or mutation count for {candidate_id}") from exc
        if length != len(sequence) or len(parent_sequence) != len(sequence):
            raise ValueError(f"Sequence length mismatch for {candidate_id}")
        if mutation_count not in {1, 2, 3}:
            raise ValueError(f"Mutation count must be 1, 2 or 3 for {candidate_id}")
        if candidate_id != f"{parent_id}__AFD{mutation_count}":
            raise ValueError(f"Candidate ID does not match its parent/round: {candidate_id}")
        hamming = sum(a != b for a, b in zip(parent_sequence, sequence))
        if hamming != mutation_count:
            raise ValueError(
                f"Expected Hamming distance {mutation_count} for {candidate_id}, observed {hamming}"
            )
        mutations = [row[f"mutation_{index}"].strip() for index in range(1, 4)]
        if any(not mutations[index] for index in range(mutation_count)):
            raise ValueError(f"Missing cumulative mutation annotation for {candidate_id}")
        if any(mutations[index] for index in range(mutation_count, 3)):
            raise ValueError(f"Unexpected later-round mutation annotation for {candidate_id}")

        parent_loss = parse_finite(row, "afdesign_parent_loss", candidate_id)
        final_loss = parse_finite(row, "afdesign_final_loss", candidate_id)
        loss_delta = parse_finite(row, "afdesign_loss_delta", candidate_id)
        round_loss = parse_finite(row, f"afdesign_round_{mutation_count}_loss", candidate_id)
        if not math.isclose(final_loss, round_loss, rel_tol=1e-7, abs_tol=1e-7):
            raise ValueError(f"Final loss does not match round loss for {candidate_id}")
        if not math.isclose(loss_delta, final_loss - parent_loss, rel_tol=1e-6, abs_tol=1e-6):
            raise ValueError(f"AfDesign loss delta is inconsistent for {candidate_id}")

        row["candidate_id"] = candidate_id
        row["parent_candidate_id"] = parent_id
        row["parent_sequence"] = parent_sequence
        row["sequence"] = sequence
        row["_length"] = length
        row["_mutation_count"] = mutation_count
        row["_sequence_filter_pass"] = _truthy(row["sequence_filter_pass"])
        if not row["_sequence_filter_pass"]:
            raise ValueError(
                f"Sequence-filter failure found in the Stage 07 pass CSV: {candidate_id}"
            )
        seen_ids.add(candidate_id)
        grouped[parent_id].append(row)

    for parent_id, family in grouped.items():
        ordered = sorted(family, key=lambda item: item["_mutation_count"])
        mutation_counts = [row["_mutation_count"] for row in ordered]
        if len(mutation_counts) != len(set(mutation_counts)):
            raise ValueError(f"Parent {parent_id} has a duplicated optimization round")
        if len({row["parent_sequence"] for row in ordered}) != 1:
            raise ValueError(f"Parent sequence is inconsistent within family {parent_id}")
        for row in ordered:
            reconstructed = row["parent_sequence"]
            used_positions: set[int] = set()
            for round_index in range(1, row["_mutation_count"] + 1):
                mutation = row[f"mutation_{round_index}"]
                match = MUTATION_PATTERN.fullmatch(mutation)
                if match is None:
                    raise ValueError(
                        f"Invalid mutation notation for {row['candidate_id']}: {mutation}"
                    )
                source_aa, position_text, replacement_aa = match.groups()
                position = int(position_text) - 1
                if position < 0 or position >= len(reconstructed) or position in used_positions:
                    raise ValueError(f"Invalid or reused mutation position for {row['candidate_id']}")
                if reconstructed[position] != source_aa or replacement_aa == source_aa:
                    raise ValueError(
                        f"Mutation does not match its source sequence for {row['candidate_id']}"
                    )
                reconstructed = (
                    reconstructed[:position] + replacement_aa + reconstructed[position + 1:]
                )
                used_positions.add(position)
            if reconstructed != row["sequence"]:
                raise ValueError(f"Mutation annotations do not reconstruct {row['candidate_id']}")
    return rows


def yaml_scalar(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def boltz_yaml(
    receptor_sequence: str,
    peptide_sequence: str,
    receptor_chain_id: str,
    peptide_chain_id: str,
    receptor_msa: str,
    peptide_msa: str,
) -> str:
    return (
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        f"      id: {yaml_scalar(receptor_chain_id)}\n"
        f"      sequence: {yaml_scalar(receptor_sequence)}\n"
        f"      msa: {yaml_scalar(receptor_msa)}\n"
        "  - protein:\n"
        f"      id: {yaml_scalar(peptide_chain_id)}\n"
        f"      sequence: {yaml_scalar(peptide_sequence)}\n"
        f"      msa: {yaml_scalar(peptide_msa)}\n"
    )


def relative_posix(path: Path, base: Path) -> str:
    return path.relative_to(base).as_posix()


def validate_yaml_directory(path: Path, expected_names: set[str]) -> None:
    observed = {item.name for item in path.glob("*.yaml")}
    unexpected = sorted(observed - expected_names)
    if unexpected:
        raise RuntimeError(
            f"Stale YAML files exist in {path}; move or remove them before rerunning: {unexpected[:5]}"
        )


def build_boltz_flags(boltz_config: dict) -> str:
    flags = (
        "  --accelerator gpu \\\n"
        f"  --devices {int(boltz_config.get('devices', 1))} \\\n"
        f"  --recycling_steps {int(boltz_config.get('recycling_steps', 3))} \\\n"
        f"  --sampling_steps {int(boltz_config.get('sampling_steps', 200))} \\\n"
        f"  --diffusion_samples {int(boltz_config.get('diffusion_samples', 1))} \\\n"
        f"  --output_format {str(boltz_config.get('output_format', 'pdb'))}"
    )
    if bool(boltz_config.get("write_full_pae", True)):
        flags += " \\\n  --write_full_pae"
    if bool(boltz_config.get("use_potentials", True)):
        flags += " \\\n  --use_potentials"
    if bool(boltz_config.get("no_kernels", False)):
        flags += " \\\n  --no_kernels"
    if boltz_config.get("max_msa_seqs") is not None:
        flags += f" \\\n  --max_msa_seqs {int(boltz_config['max_msa_seqs'])}"
    if boltz_config.get("seed") is not None:
        flags += f" \\\n  --seed {int(boltz_config['seed'])}"
    return flags


def run(config_path: Path) -> Path:
    config, project_root = load_config(config_path)
    variants_path = resolve_path(project_root, config["variants_csv"])
    afdesign_manifest_path = resolve_path(project_root, config["afdesign_manifest"])
    target_path = resolve_path(project_root, config["target_fasta"])
    target_manifest_path = resolve_path(project_root, config["target_manifest"])
    masks_path = resolve_path(project_root, config["screening_masks"])
    output_dir = resolve_path(project_root, config["output_dir"])

    rows = read_variants(variants_path)
    afdesign_manifest = json.loads(afdesign_manifest_path.read_text(encoding="utf-8"))
    target_header, receptor_sequence = read_single_fasta(target_path)
    target_manifest = json.loads(target_manifest_path.read_text(encoding="utf-8"))
    masks = json.loads(masks_path.read_text(encoding="utf-8"))

    expected_variants = int(config["expected_variant_count"])
    observed_parents = len({row["parent_candidate_id"] for row in rows})
    if len(rows) != expected_variants:
        raise ValueError(f"Expected {expected_variants} sequence-pass variants; observed {len(rows)}")
    manifest_counts = afdesign_manifest.get("counts", {})
    if int(manifest_counts.get("sequence_filter_pass_variants", -1)) != len(rows):
        raise ValueError("AfDesign manifest sequence-pass count does not match the CSV")
    if int(manifest_counts.get("errors", -1)) != 0:
        raise ValueError("AfDesign manifest reports optimization errors")
    if len(receptor_sequence) != int(target_manifest["target"]["sequence_length"]):
        raise ValueError("Target FASTA length does not match target manifest")

    receptor_msa = str(config["receptor_msa"])
    peptide_msa = str(config["peptide_msa"])
    if receptor_msa != "empty" and not resolve_path(project_root, receptor_msa).is_file():
        raise FileNotFoundError(f"Receptor MSA does not exist: {receptor_msa}")
    if peptide_msa != "empty" and not resolve_path(project_root, peptide_msa).is_file():
        raise FileNotFoundError(f"Peptide MSA does not exist: {peptide_msa}")

    output_dir.mkdir(parents=True, exist_ok=True)
    yaml_dir = output_dir / "boltz_inputs_afdesign"
    yaml_dir.mkdir(parents=True, exist_ok=True)
    expected_yamls = {f"{row['candidate_id']}.yaml" for row in rows}
    validate_yaml_directory(yaml_dir, expected_yamls)

    receptor_chain = str(config["receptor_chain_id"])
    peptide_chain = str(config["peptide_chain_id"])
    task_rows = []
    for row in sorted(
        rows,
        key=lambda item: (
            item["_length"], item["parent_candidate_id"], item["_mutation_count"]
        ),
    ):
        yaml_path = yaml_dir / f"{row['candidate_id']}.yaml"
        atomic_write_text(
            yaml_path,
            boltz_yaml(
                receptor_sequence, row["sequence"], receptor_chain, peptide_chain,
                receptor_msa, peptide_msa,
            ),
        )
        task_rows.append({
            "job_id": f"NTSR1__{row['candidate_id']}",
            "candidate_id": row["candidate_id"],
            "parent_candidate_id": row["parent_candidate_id"],
            "parent_sequence": row["parent_sequence"],
            "sequence": row["sequence"],
            "length": row["_length"],
            "mutation_count": row["_mutation_count"],
            "mutation_1": row["mutation_1"],
            "mutation_2": row["mutation_2"],
            "mutation_3": row["mutation_3"],
            "afdesign_parent_loss": row["afdesign_parent_loss"],
            "afdesign_final_loss": row["afdesign_final_loss"],
            "afdesign_loss_delta": row["afdesign_loss_delta"],
            "sequence_filter_pass": row["sequence_filter_pass"],
            "sequence_filter_reasons": row["sequence_filter_reasons"],
            "source_structure_pdb": row["source_structure_pdb"],
            "boltz_yaml": relative_posix(yaml_path, project_root),
            "receptor_chain_id": receptor_chain,
            "peptide_chain_id": peptide_chain,
        })
    tasks_path = output_dir / "structure_tasks_afdesign.csv"
    atomic_write_csv(tasks_path, task_rows, TASK_FIELDS)

    receptor_context = {
        "target_header": target_header,
        "target_accession": target_manifest["target"]["accession"],
        "target_length": len(receptor_sequence),
        "target_fasta": str(target_path),
        "target_fasta_sha256": sha256_file(target_path),
        "receptor_chain_id": receptor_chain,
        "peptide_chain_id": peptide_chain,
        "orthosteric_pocket_union": masks["orthosteric_pocket_union"],
        "orthosteric_pocket_intersection": masks["orthosteric_pocket_intersection"],
        "putative_extracellular_glycan_sites": masks["putative_extracellular_glycan_sites"],
        "do_not_dock_through_missing_segments": masks["do_not_dock_through_missing_segments"],
        "small_molecule_grids_are_peptide_docking_boxes": masks[
            "small_molecule_grids_are_peptide_docking_boxes"
        ],
    }
    receptor_context_path = output_dir / "receptor_context.json"
    atomic_write_json(receptor_context_path, receptor_context)

    preflight_script = (
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        "command -v boltz >/dev/null || { echo 'ERROR: boltz is not installed' >&2; exit 1; }\n"
        "python - <<'PY'\n"
        "import torch\n"
        "print('torch:', torch.__version__)\n"
        "print('torch CUDA:', torch.version.cuda)\n"
        "print('CUDA available:', torch.cuda.is_available())\n"
        "if not torch.cuda.is_available():\n"
        "    raise SystemExit('ERROR: CUDA is unavailable; refusing to run Boltz on CPU')\n"
        "print('GPU:', torch.cuda.get_device_name(0))\n"
        "PY\n"
        "boltz predict --help >/dev/null\n"
        "echo 'Boltz/GPU preflight PASS'\n"
    )
    preflight_path = output_dir / "run_boltz_preflight.sh"
    atomic_write_text(preflight_path, preflight_script)

    boltz_flags = build_boltz_flags(config["boltz_prediction"])
    run_script = (
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        "export BOLTZ_CACHE=\"${BOLTZ_CACHE:-$PWD/models/boltz_cache}\"\n\n"
        "bash work/08_prepare_optimized_structure_jobs/run_boltz_preflight.sh\n\n"
        "boltz predict work/08_prepare_optimized_structure_jobs/boltz_inputs_afdesign \\\n"
        "  --out_dir work/08_prepare_optimized_structure_jobs/boltz_predictions_afdesign \\\n"
        f"{boltz_flags}\n"
    )
    run_script = run_script.replace("work/08_prepare_optimized_structure_jobs/", shlex.quote(output_dir.as_posix()) + "/")
    run_path = output_dir / "run_boltz_afdesign.sh"
    atomic_write_text(run_path, run_script)

    mutation_counts = Counter(row["_mutation_count"] for row in rows)
    pass_count = sum(row["_sequence_filter_pass"] for row in rows)
    unique_sequences = len({row["sequence"] for row in rows})
    manifest = {
        "schema_version": 1,
        "stage": "08_prepare_optimized_structure_jobs",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": str(config_path.resolve()),
        "config_fingerprint": canonical_json_hash(config),
        "inputs": {
            "afdesign_variants_csv": str(variants_path),
            "afdesign_variants_csv_sha256": sha256_file(variants_path),
            "afdesign_manifest": str(afdesign_manifest_path),
            "afdesign_manifest_sha256": sha256_file(afdesign_manifest_path),
            "target_fasta": str(target_path),
            "target_manifest": str(target_manifest_path),
            "screening_masks": str(masks_path),
        },
        "policy": {
            "input_scope": "all Stage 07 sequence-filter-pass AFD variants",
            "stage07_sequence_filter_applied_upstream": True,
            "additional_sequence_filter_applied": False,
            "afdesign_loss_filter_applied": False,
            "sequence_deduplication_applied": False,
            "structure_filter_applied": False,
            "one_boltz_job_per_variant_row": True,
        },
        "method": {
            "predictor_input": "Boltz YAML v1",
            "receptor_input": "full-length O88319 sequence",
            "receptor_msa": receptor_msa,
            "peptide_msa": peptide_msa,
            "templates": "none",
            "parent_structure_used_as_prediction_template": False,
            "orthosteric_pocket_used_as_positive_constraint": False,
            "boltz_prediction": config["boltz_prediction"],
        },
        "counts": {
            "parent_peptides": observed_parents,
            "input_variant_rows": len(rows),
            "boltz_tasks": len(task_rows),
            "unique_sequences_annotation_only": unique_sequences,
            "duplicate_sequence_rows_retained": len(rows) - unique_sequences,
            "sequence_filter_pass_annotation": pass_count,
            "stage07_sequence_filter_fail_rows_excluded": int(
                manifest_counts.get("sequence_filter_fail_variants", 0)
            ),
            "by_mutation_count": {str(key): mutation_counts[key] for key in (1, 2, 3)},
        },
        "outputs": {
            "tasks_csv": str(tasks_path),
            "boltz_inputs": str(yaml_dir),
            "receptor_context": str(receptor_context_path),
            "boltz_preflight_script": str(preflight_path),
            "boltz_run_script": str(run_path),
            "boltz_predictions": str(output_dir / "boltz_predictions_afdesign"),
        },
    }
    manifest_path = output_dir / "structure_job_manifest.json"
    atomic_write_json(manifest_path, manifest)

    report = (
        "# Stage 08 AfDesign-variant structure jobs\n\n"
        f"- Parent peptides: {observed_parents}\n"
        f"- Variant rows and Boltz jobs: {len(task_rows)}\n"
        f"- AFD1/AFD2/AFD3: {mutation_counts[1]}/{mutation_counts[2]}/{mutation_counts[3]}\n"
        f"- All {pass_count} Stage 07 sequence-filter-pass variants are retained.\n"
        f"- Stage 07 sequence-filter failures excluded: {manifest_counts.get('sequence_filter_fail_variants', 0)}\n"
        "- No additional sequence, AfDesign-loss, diversity, deduplication or structure filter is applied.\n"
        "- Every variant row receives an independent Boltz input and prediction job.\n"
        "- The parent PDB is metadata only and is not used as a Boltz template.\n"
        "- Structural thresholds will be selected only after all predictions are audited.\n"
    )
    atomic_write_text(output_dir / "structure_prediction_readiness.md", report)
    return output_dir


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    try:
        output_dir = run(args.config.resolve())
    except Exception as exc:
        print(f"08_prepare_optimized_structure_jobs FAILED: {exc}", file=sys.stderr)
        return 1
    print(f"08_prepare_optimized_structure_jobs PASS: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
