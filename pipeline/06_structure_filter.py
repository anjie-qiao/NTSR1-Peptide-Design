#!/usr/bin/env python3
"""Audit and threshold-filter Boltz NTSR1–peptide complex predictions."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from pepmlm_pipeline_common import (
    atomic_write_csv,
    atomic_write_json,
    atomic_write_text,
    canonical_json_hash,
    sha256_file,
)


METRIC_FIELDS = [
    "candidate_id", "sequence", "length", "sampling_profile",
    "pseudo_perplexity", "ppl_percentile_within_length", "cluster_id",
    "prediction_status", "prediction_error", "confidence_score", "ptm", "iptm",
    "complex_plddt", "complex_iplddt", "complex_pde", "complex_ipde",
    "receptor_ptm", "peptide_ptm", "pair_iptm_a_to_b", "pair_iptm_b_to_a",
    "pair_iptm_mean", "pair_iptm_min", "receptor_mean_plddt",
    "peptide_mean_plddt", "peptide_q10_plddt", "interface_pae_mean",
    "interface_pae_q90", "minimum_interface_distance", "peptide_contact_residues",
    "peptide_contact_fraction", "receptor_contact_residues", "interface_residue_pairs",
    "pocket_receptor_contact_residues", "pocket_peptide_contact_residues",
    "pocket_peptide_contact_fraction", "severe_clash_atom_pairs",
    "selection_applied", "structure_pass", "filter_failure_reasons",
    "structure_pdb", "confidence_json", "pae_npz", "plddt_npz",
]

SUMMARY_METRICS = [
    "confidence_score", "ptm", "iptm", "complex_plddt", "complex_iplddt",
    "complex_pde", "complex_ipde", "receptor_ptm", "peptide_ptm",
    "pair_iptm_mean", "pair_iptm_min", "receptor_mean_plddt",
    "peptide_mean_plddt", "peptide_q10_plddt", "interface_pae_mean",
    "interface_pae_q90", "minimum_interface_distance", "peptide_contact_residues",
    "peptide_contact_fraction", "receptor_contact_residues", "interface_residue_pairs",
    "pocket_receptor_contact_residues", "pocket_peptide_contact_residues",
    "pocket_peptide_contact_fraction", "severe_clash_atom_pairs",
]

THRESHOLD_RULES = {
    "min_pair_iptm_mean": ("pair_iptm_mean", "min"),
    "min_pair_iptm_min": ("pair_iptm_min", "min"),
    "min_peptide_ptm": ("peptide_ptm", "min"),
    "min_peptide_mean_plddt": ("peptide_mean_plddt", "min"),
    "min_peptide_q10_plddt": ("peptide_q10_plddt", "min"),
    "max_interface_pae_mean": ("interface_pae_mean", "max"),
    "max_interface_pae_q90": ("interface_pae_q90", "max"),
    "min_peptide_contact_fraction": ("peptide_contact_fraction", "min"),
    "min_interface_residue_pairs": ("interface_residue_pairs", "min"),
    "max_severe_clash_atom_pairs": ("severe_clash_atom_pairs", "max"),
}


def resolve_path(project_root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def load_config(config_path: Path) -> tuple[dict, Path]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    required = {
        "tasks_csv", "boltz_output_dir", "receptor_context", "output_dir",
        "model_index", "mode", "geometry", "thresholds",
    }
    missing = required - set(config)
    if missing:
        raise ValueError(f"Missing config keys: {sorted(missing)}")
    if config["mode"] not in {"audit", "thresholds"}:
        raise ValueError("mode must be 'audit' or 'thresholds'")
    geometry = config["geometry"]
    for key in ("contact_cutoff_angstrom", "pocket_cutoff_angstrom", "clash_cutoff_angstrom"):
        if key not in geometry or float(geometry[key]) <= 0:
            raise ValueError(f"geometry.{key} must be positive")
    unknown = set(config["thresholds"]) - set(THRESHOLD_RULES)
    if unknown:
        raise ValueError(f"Unknown threshold keys: {sorted(unknown)}")
    if config["mode"] == "thresholds":
        active = [value for value in config["thresholds"].values() if value is not None]
        if not active:
            raise ValueError("thresholds mode requires at least one active threshold")
    return config, config_path.resolve().parent.parent


def read_tasks(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "candidate_id", "sequence", "length", "sampling_profile",
            "pseudo_perplexity", "ppl_percentile_within_length", "cluster_id",
        }
        if not reader.fieldnames or required - set(reader.fieldnames):
            raise ValueError(f"Task CSV is missing fields: {sorted(required - set(reader.fieldnames or []))}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise ValueError("Task CSV is empty")
    seen = set()
    for row in rows:
        candidate_id = row["candidate_id"].strip()
        sequence = row["sequence"].strip().upper()
        if not candidate_id or candidate_id in seen:
            raise ValueError(f"Missing or duplicate candidate ID: {candidate_id!r}")
        if int(row["length"]) != len(sequence):
            raise ValueError(f"Length mismatch for {candidate_id}")
        row["candidate_id"] = candidate_id
        row["sequence"] = sequence
        seen.add(candidate_id)
    return rows


def _index_files(root: Path, pattern: str, id_pattern: re.Pattern[str]) -> dict[str, Path]:
    index: dict[str, Path] = {}
    for path in root.rglob(pattern):
        match = id_pattern.fullmatch(path.name)
        if not match:
            continue
        candidate_id = match.group("candidate")
        if candidate_id in index:
            raise RuntimeError(f"Duplicate prediction artifact for {candidate_id}: {index[candidate_id]} and {path}")
        index[candidate_id] = path.resolve()
    return index


def index_prediction_files(root: Path, model_index: int) -> dict[str, dict[str, Path]]:
    escaped = re.escape(str(model_index))
    specs = {
        "confidence_json": (
            f"confidence_*_model_{model_index}.json",
            re.compile(rf"confidence_(?P<candidate>.+)_model_{escaped}\.json"),
        ),
        "pae_npz": (
            f"pae_*_model_{model_index}.npz",
            re.compile(rf"pae_(?P<candidate>.+)_model_{escaped}\.npz"),
        ),
        "plddt_npz": (
            f"plddt_*_model_{model_index}.npz",
            re.compile(rf"plddt_(?P<candidate>.+)_model_{escaped}\.npz"),
        ),
        "structure_pdb": (
            f"*_model_{model_index}.pdb",
            re.compile(rf"(?P<candidate>.+)_model_{escaped}\.pdb"),
        ),
    }
    by_kind = {kind: _index_files(root, glob, regex) for kind, (glob, regex) in specs.items()}
    candidates = set().union(*(set(index) for index in by_kind.values()))
    return {
        candidate_id: {
            kind: index[candidate_id]
            for kind, index in by_kind.items()
            if candidate_id in index
        }
        for candidate_id in candidates
    }


def _load_npz_array(path: Path, preferred_key: str) -> np.ndarray:
    with np.load(path) as archive:
        if preferred_key in archive.files:
            value = archive[preferred_key]
        elif len(archive.files) == 1:
            value = archive[archive.files[0]]
        else:
            raise ValueError(f"Cannot identify {preferred_key!r} array in {path}; keys={archive.files}")
        return np.asarray(value, dtype=float)


def _normalize_plddt(values: np.ndarray) -> np.ndarray:
    if values.size and np.nanmax(values) > 1.5:
        values = values / 100.0
    return values


def confidence_metrics(confidence: dict) -> dict[str, float]:
    pair = confidence.get("pair_chains_iptm", {})
    chains = confidence.get("chains_ptm", {})
    ab = float(pair.get("0", {}).get("1", math.nan))
    ba = float(pair.get("1", {}).get("0", math.nan))
    return {
        "confidence_score": float(confidence.get("confidence_score", math.nan)),
        "ptm": float(confidence.get("ptm", math.nan)),
        "iptm": float(confidence.get("iptm", math.nan)),
        "complex_plddt": float(confidence.get("complex_plddt", math.nan)),
        "complex_iplddt": float(confidence.get("complex_iplddt", math.nan)),
        "complex_pde": float(confidence.get("complex_pde", math.nan)),
        "complex_ipde": float(confidence.get("complex_ipde", math.nan)),
        "receptor_ptm": float(chains.get("0", math.nan)),
        "peptide_ptm": float(chains.get("1", math.nan)),
        "pair_iptm_a_to_b": ab,
        "pair_iptm_b_to_a": ba,
        "pair_iptm_mean": float(np.nanmean([ab, ba])),
        "pair_iptm_min": float(np.nanmin([ab, ba])),
    }


def token_metrics(
    pae_path: Path,
    plddt_path: Path,
    receptor_length: int,
    peptide_length: int,
) -> dict[str, float]:
    total_length = receptor_length + peptide_length
    plddt = _normalize_plddt(_load_npz_array(plddt_path, "plddt").reshape(-1))
    pae = _load_npz_array(pae_path, "pae")
    if plddt.shape[0] < total_length:
        raise ValueError(f"pLDDT has {plddt.shape[0]} tokens; expected at least {total_length}")
    if pae.ndim != 2 or pae.shape[0] < total_length or pae.shape[1] < total_length:
        raise ValueError(f"PAE shape {pae.shape}; expected at least ({total_length}, {total_length})")
    receptor = plddt[:receptor_length]
    peptide = plddt[receptor_length:total_length]
    ab = pae[:receptor_length, receptor_length:total_length].reshape(-1)
    ba = pae[receptor_length:total_length, :receptor_length].reshape(-1)
    interface = np.concatenate([ab, ba])
    return {
        "receptor_mean_plddt": float(np.mean(receptor)),
        "peptide_mean_plddt": float(np.mean(peptide)),
        "peptide_q10_plddt": float(np.quantile(peptide, 0.10)),
        "interface_pae_mean": float(np.mean(interface)),
        "interface_pae_q90": float(np.quantile(interface, 0.90)),
    }


def parse_pdb_atoms(path: Path) -> dict[str, dict[tuple[int, str], list[tuple[float, float, float]]]]:
    chains: dict[str, dict[tuple[int, str], list[tuple[float, float, float]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
                continue
            altloc = line[16]
            if altloc not in {" ", "A"}:
                continue
            atom_name = line[12:16].strip()
            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            if not element:
                element = atom_name.lstrip("0123456789")[:1].upper()
            if element == "H":
                continue
            chain = line[21].strip()
            if not chain:
                continue
            try:
                residue_number = int(line[22:26])
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            insertion = line[26].strip()
            chains[chain][(residue_number, insertion)].append(xyz)
    return chains


def geometry_metrics(
    pdb_path: Path,
    receptor_chain: str,
    peptide_chain: str,
    peptide_length: int,
    pocket_residue_numbers: set[int],
    contact_cutoff: float,
    pocket_cutoff: float,
    clash_cutoff: float,
) -> dict[str, float | int]:
    chains = parse_pdb_atoms(pdb_path)
    if receptor_chain not in chains or peptide_chain not in chains:
        raise ValueError(
            f"Expected chains {receptor_chain!r}/{peptide_chain!r}; observed={sorted(chains)}"
        )
    receptor = chains[receptor_chain]
    peptide = chains[peptide_chain]
    query_cutoff = max(contact_cutoff, pocket_cutoff, clash_cutoff)
    cell_size = query_cutoff
    grid: dict[tuple[int, int, int], list[tuple[tuple[int, str], tuple[float, float, float]]]] = defaultdict(list)
    for residue_key, atoms in receptor.items():
        for xyz in atoms:
            cell = tuple(math.floor(value / cell_size) for value in xyz)
            grid[cell].append((residue_key, xyz))

    contact_pairs: set[tuple[tuple[int, str], tuple[int, str]]] = set()
    pocket_pairs: set[tuple[tuple[int, str], tuple[int, str]]] = set()
    clash_pairs: set[tuple[tuple[int, str], int, tuple[int, str], tuple[float, float, float]]] = set()
    minimum_distance = math.inf
    contact_sq = contact_cutoff ** 2
    pocket_sq = pocket_cutoff ** 2
    clash_sq = clash_cutoff ** 2
    for peptide_residue, atoms in peptide.items():
        for atom_index, xyz in enumerate(atoms):
            base = tuple(math.floor(value / cell_size) for value in xyz)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for receptor_residue, receptor_xyz in grid.get(
                            (base[0] + dx, base[1] + dy, base[2] + dz), []
                        ):
                            distance_sq = sum((xyz[i] - receptor_xyz[i]) ** 2 for i in range(3))
                            if distance_sq < minimum_distance ** 2:
                                minimum_distance = math.sqrt(distance_sq)
                            if distance_sq <= contact_sq:
                                contact_pairs.add((peptide_residue, receptor_residue))
                            if receptor_residue[0] in pocket_residue_numbers and distance_sq <= pocket_sq:
                                pocket_pairs.add((peptide_residue, receptor_residue))
                            if distance_sq <= clash_sq:
                                clash_pairs.add((peptide_residue, atom_index, receptor_residue, receptor_xyz))

    peptide_contacts = {pair[0] for pair in contact_pairs}
    receptor_contacts = {pair[1] for pair in contact_pairs}
    pocket_peptide_contacts = {pair[0] for pair in pocket_pairs}
    pocket_receptor_contacts = {pair[1] for pair in pocket_pairs}
    return {
        "minimum_interface_distance": minimum_distance,
        "peptide_contact_residues": len(peptide_contacts),
        "peptide_contact_fraction": len(peptide_contacts) / peptide_length,
        "receptor_contact_residues": len(receptor_contacts),
        "interface_residue_pairs": len(contact_pairs),
        "pocket_receptor_contact_residues": len(pocket_receptor_contacts),
        "pocket_peptide_contact_residues": len(pocket_peptide_contacts),
        "pocket_peptide_contact_fraction": len(pocket_peptide_contacts) / peptide_length,
        "severe_clash_atom_pairs": len(clash_pairs),
    }


def threshold_failures(metrics: dict, thresholds: dict) -> list[str]:
    failures = []
    for threshold_name, threshold_value in thresholds.items():
        if threshold_value is None:
            continue
        metric_name, direction = THRESHOLD_RULES[threshold_name]
        value = float(metrics[metric_name])
        if not math.isfinite(value):
            failures.append(f"{metric_name}_missing")
        elif direction == "min" and value < float(threshold_value):
            failures.append(threshold_name)
        elif direction == "max" and value > float(threshold_value):
            failures.append(threshold_name)
    return failures


def _serialize(value: object) -> object:
    if isinstance(value, float):
        return "" if not math.isfinite(value) else f"{value:.8f}"
    return value


def summarize(rows: list[dict]) -> dict:
    quantiles = [0.0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1.0]
    summary = {}
    valid_rows = [row for row in rows if row["prediction_status"] == "ok"]
    for metric in SUMMARY_METRICS:
        values = np.asarray(
            [float(row[metric]) for row in valid_rows if math.isfinite(float(row[metric]))],
            dtype=float,
        )
        if not values.size:
            summary[metric] = {"count": 0, "quantiles": {}}
            continue
        summary[metric] = {
            "count": int(values.size),
            "mean": float(values.mean()),
            "quantiles": {
                f"q{int(q * 100):02d}": float(np.quantile(values, q)) for q in quantiles
            },
        }
    return summary


def run(config_path: Path) -> Path:
    config, project_root = load_config(config_path)
    tasks_path = resolve_path(project_root, config["tasks_csv"])
    boltz_root = resolve_path(project_root, config["boltz_output_dir"])
    receptor_context_path = resolve_path(project_root, config["receptor_context"])
    output_dir = resolve_path(project_root, config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    if not boltz_root.exists():
        raise FileNotFoundError(f"Boltz output directory does not exist: {boltz_root}")

    tasks = read_tasks(tasks_path)
    receptor_context = json.loads(receptor_context_path.read_text(encoding="utf-8"))
    receptor_length = int(receptor_context["target_length"])
    receptor_chain = str(receptor_context["receptor_chain_id"])
    peptide_chain = str(receptor_context["peptide_chain_id"])
    pocket_numbers = {
        int(re.search(r"\d+", value).group())
        for value in receptor_context.get("orthosteric_pocket_union", [])
        if re.search(r"\d+", value)
    }
    artifacts = index_prediction_files(boltz_root, int(config["model_index"]))
    geometry = config["geometry"]
    selection_applied = config["mode"] == "thresholds"

    rows: list[dict] = []
    for task in tasks:
        candidate_id = task["candidate_id"]
        row = {
            "candidate_id": candidate_id,
            "sequence": task["sequence"],
            "length": int(task["length"]),
            "sampling_profile": task["sampling_profile"],
            "pseudo_perplexity": task["pseudo_perplexity"],
            "ppl_percentile_within_length": task["ppl_percentile_within_length"],
            "cluster_id": task["cluster_id"],
            "prediction_status": "error",
            "prediction_error": "",
            "selection_applied": int(selection_applied),
            "structure_pass": "",
            "filter_failure_reasons": "",
        }
        files = artifacts.get(candidate_id, {})
        for kind in ("structure_pdb", "confidence_json", "pae_npz", "plddt_npz"):
            row[kind] = str(files[kind]) if kind in files else ""
        missing = [kind for kind in ("structure_pdb", "confidence_json", "pae_npz", "plddt_npz") if kind not in files]
        try:
            if missing:
                raise FileNotFoundError(f"missing artifacts: {','.join(missing)}")
            confidence = json.loads(files["confidence_json"].read_text(encoding="utf-8"))
            row.update(confidence_metrics(confidence))
            row.update(token_metrics(
                files["pae_npz"], files["plddt_npz"], receptor_length, int(task["length"])
            ))
            row.update(geometry_metrics(
                files["structure_pdb"], receptor_chain, peptide_chain, int(task["length"]),
                pocket_numbers,
                float(geometry["contact_cutoff_angstrom"]),
                float(geometry["pocket_cutoff_angstrom"]),
                float(geometry["clash_cutoff_angstrom"]),
            ))
            row["prediction_status"] = "ok"
            if selection_applied:
                failures = threshold_failures(row, config["thresholds"])
                row["structure_pass"] = int(not failures)
                row["filter_failure_reasons"] = ";".join(failures)
        except Exception as exc:  # keep the audit complete and report per-candidate errors
            row["prediction_error"] = f"{type(exc).__name__}: {exc}"
            if selection_applied:
                row["structure_pass"] = 0
                row["filter_failure_reasons"] = "prediction_artifact_error"
        for field in SUMMARY_METRICS:
            row.setdefault(field, math.nan)
        rows.append(row)

    serialized_rows = [{field: _serialize(row.get(field, "")) for field in METRIC_FIELDS} for row in rows]
    metrics_path = output_dir / "structure_metrics_all.csv"
    atomic_write_csv(metrics_path, serialized_rows, METRIC_FIELDS)

    summary = {
        "schema_version": 1,
        "stage": "06_structure_filter",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "mode": config["mode"],
        "counts": {
            "input": len(rows),
            "prediction_ok": sum(row["prediction_status"] == "ok" for row in rows),
            "prediction_error": sum(row["prediction_status"] != "ok" for row in rows),
            "structure_pass": (
                sum(row.get("structure_pass") == 1 for row in rows) if selection_applied else None
            ),
        },
        "metric_distributions": summarize(rows),
    }
    summary_path = output_dir / "structure_metrics_summary.json"
    atomic_write_json(summary_path, summary)

    outputs = {
        "structure_metrics_all": str(metrics_path),
        "structure_metrics_summary": str(summary_path),
    }
    if selection_applied:
        pass_rows = [row for row in serialized_rows if row["structure_pass"] == 1]
        pass_csv = output_dir / "structure_pass.csv"
        pass_fasta = output_dir / "structure_pass.fasta"
        atomic_write_csv(pass_csv, pass_rows, METRIC_FIELDS)
        atomic_write_text(
            pass_fasta,
            "".join(f">{row['candidate_id']}\n{row['sequence']}\n" for row in pass_rows),
        )
        outputs.update({"structure_pass_csv": str(pass_csv), "structure_pass_fasta": str(pass_fasta)})

    manifest = {
        "schema_version": 1,
        "stage": "06_structure_filter",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": str(config_path.resolve()),
        "config_fingerprint": canonical_json_hash(config),
        "inputs": {
            "tasks_csv": str(tasks_path),
            "tasks_csv_sha256": sha256_file(tasks_path),
            "boltz_output_dir": str(boltz_root),
            "receptor_context": str(receptor_context_path),
            "model_index": int(config["model_index"]),
        },
        "mode": config["mode"],
        "thresholds": config["thresholds"],
        "geometry": geometry,
        "counts": summary["counts"],
        "outputs": outputs,
    }
    atomic_write_json(output_dir / "structure_filter_manifest.json", manifest)
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    try:
        output_dir = run(args.config.resolve())
    except Exception as exc:
        print(f"06_structure_filter FAILED: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    print(f"06_structure_filter PASS: {output_dir}")


if __name__ == "__main__":
    main()
