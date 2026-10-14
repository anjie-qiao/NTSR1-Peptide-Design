#!/usr/bin/env python3
"""Audit all AfDesign-variant Boltz structures and compare them with parents."""

from __future__ import annotations

import argparse
import csv
import importlib
import json
import math
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from pepmlm_pipeline_common import (
    atomic_write_csv,
    atomic_write_json,
    canonical_json_hash,
    sha256_file,
)


stage06 = importlib.import_module("06_structure_filter")

TASK_METADATA_FIELDS = [
    "candidate_id", "parent_candidate_id", "parent_sequence", "sequence", "length",
    "mutation_count", "mutation_1", "mutation_2", "mutation_3",
    "afdesign_parent_loss", "afdesign_final_loss", "afdesign_loss_delta",
    "sequence_filter_pass", "sequence_filter_reasons", "source_structure_pdb",
]

STRUCTURE_METRIC_FIELDS = [
    "confidence_score", "ptm", "iptm", "complex_plddt", "complex_iplddt",
    "complex_pde", "complex_ipde", "receptor_ptm", "peptide_ptm",
    "pair_iptm_a_to_b", "pair_iptm_b_to_a", "pair_iptm_mean", "pair_iptm_min",
    "receptor_mean_plddt", "peptide_mean_plddt", "peptide_q10_plddt",
    "interface_pae_mean", "interface_pae_q90", "minimum_interface_distance",
    "peptide_contact_residues", "peptide_contact_fraction",
    "receptor_contact_residues", "interface_residue_pairs",
    "pocket_receptor_contact_residues", "pocket_peptide_contact_residues",
    "pocket_peptide_contact_fraction", "severe_clash_atom_pairs",
]

AUDIT_FIELDS = (
    TASK_METADATA_FIELDS
    + ["prediction_status", "prediction_error"]
    + STRUCTURE_METRIC_FIELDS
    + ["structure_pdb", "confidence_json", "pae_npz", "plddt_npz"]
)

COMPARISON_DIRECTIONS = {
    "confidence_score": "higher",
    "pair_iptm_mean": "higher",
    "pair_iptm_min": "higher",
    "peptide_ptm": "higher",
    "peptide_mean_plddt": "higher",
    "peptide_q10_plddt": "higher",
    "interface_pae_mean": "lower",
    "interface_pae_q90": "lower",
    "peptide_contact_fraction": "higher",
    "interface_residue_pairs": "higher",
    "pocket_peptide_contact_fraction": "higher",
    "severe_clash_atom_pairs": "lower",
}

COMPARISON_FIELDS = TASK_METADATA_FIELDS + [
    "prediction_status", "metrics_improved_count", "metrics_worsened_count",
    "metrics_tied_count",
]
for _metric in COMPARISON_DIRECTIONS:
    COMPARISON_FIELDS.extend([
        f"parent_{_metric}", f"optimized_{_metric}", f"delta_{_metric}",
        f"improved_{_metric}",
    ])


def resolve_path(project_root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def load_config(config_path: Path) -> tuple[dict, Path]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    required = {
        "tasks_csv", "structure_job_manifest", "boltz_output_dir", "receptor_context",
        "parent_metrics_csv", "output_dir", "expected_task_count",
        "expected_parent_metrics_count", "model_index", "require_complete", "geometry",
    }
    missing = required - set(config)
    if missing:
        raise ValueError(f"Missing config keys: {sorted(missing)}")
    if int(config["expected_task_count"]) < 1 or int(config["expected_parent_metrics_count"]) < 1:
        raise ValueError("Expected task and parent-metric counts must be positive")
    for key in ("contact_cutoff_angstrom", "pocket_cutoff_angstrom", "clash_cutoff_angstrom"):
        if key not in config["geometry"] or float(config["geometry"][key]) <= 0:
            raise ValueError(f"geometry.{key} must be positive")
    return config, config_path.resolve().parent.parent


def _truthy(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def read_tasks(path: Path, expected_count: int) -> list[dict]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = set(TASK_METADATA_FIELDS) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Stage 08 tasks CSV is missing fields: {sorted(missing)}")
        rows = [dict(row) for row in reader]
    if len(rows) != expected_count:
        raise ValueError(f"Expected {expected_count} Stage 08 tasks; observed {len(rows)}")

    seen: set[str] = set()
    for row in rows:
        candidate_id = row["candidate_id"].strip()
        parent_id = row["parent_candidate_id"].strip()
        sequence = row["sequence"].strip().upper()
        parent_sequence = row["parent_sequence"].strip().upper()
        if not candidate_id or candidate_id in seen or not parent_id:
            raise ValueError(f"Missing or duplicate task identity: {candidate_id!r}")
        try:
            length = int(row["length"])
            mutation_count = int(row["mutation_count"])
        except ValueError as exc:
            raise ValueError(f"Invalid length/mutation count for {candidate_id}") from exc
        if length != len(sequence) or len(parent_sequence) != length:
            raise ValueError(f"Sequence length mismatch for {candidate_id}")
        if mutation_count not in {1, 2, 3}:
            raise ValueError(f"Unexpected mutation count for {candidate_id}")
        if candidate_id != f"{parent_id}__AFD{mutation_count}":
            raise ValueError(f"Candidate ID does not encode its parent/round: {candidate_id}")
        if not _truthy(row["sequence_filter_pass"]):
            raise ValueError(f"Task did not pass the Stage 07 sequence filter: {candidate_id}")
        row["candidate_id"] = candidate_id
        row["parent_candidate_id"] = parent_id
        row["sequence"] = sequence
        row["parent_sequence"] = parent_sequence
        row["_length"] = length
        row["_mutation_count"] = mutation_count
        seen.add(candidate_id)
    return rows


def parse_metric(row: dict, field: str, candidate_id: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Missing/invalid parent metric {field} for {candidate_id}") from exc
    if not math.isfinite(value):
        raise ValueError(f"Non-finite parent metric {field} for {candidate_id}")
    return value


def read_parent_metrics(path: Path, expected_count: int) -> dict[str, dict]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"candidate_id", "sequence", "prediction_status"} | set(COMPARISON_DIRECTIONS)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Parent metrics CSV is missing fields: {sorted(missing)}")
        rows = [dict(row) for row in reader]
    if len(rows) != expected_count:
        raise ValueError(f"Expected {expected_count} parent metric rows; observed {len(rows)}")
    index: dict[str, dict] = {}
    for row in rows:
        candidate_id = row["candidate_id"].strip()
        if not candidate_id or candidate_id in index:
            raise ValueError(f"Missing or duplicate parent metric ID: {candidate_id!r}")
        if row["prediction_status"].strip() != "ok":
            raise ValueError(f"Parent prediction is not valid: {candidate_id}")
        for metric in COMPARISON_DIRECTIONS:
            parse_metric(row, metric, candidate_id)
        index[candidate_id] = row
    return index


def serialize(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        return "" if not math.isfinite(float(value)) else f"{float(value):.8f}"
    if isinstance(value, (int, np.integer)):
        return int(value)
    return value


def audit_predictions(
    tasks: list[dict],
    artifacts: dict[str, dict[str, Path]],
    receptor_context: dict,
    geometry: dict,
) -> list[dict]:
    receptor_length = int(receptor_context["target_length"])
    receptor_chain = str(receptor_context["receptor_chain_id"])
    peptide_chain = str(receptor_context["peptide_chain_id"])
    pocket_numbers = {
        int(re.search(r"\d+", value).group())
        for value in receptor_context.get("orthosteric_pocket_union", [])
        if re.search(r"\d+", value)
    }
    rows: list[dict] = []
    required_artifacts = ("structure_pdb", "confidence_json", "pae_npz", "plddt_npz")
    for task in tasks:
        candidate_id = task["candidate_id"]
        row = {field: task.get(field, "") for field in TASK_METADATA_FIELDS}
        row.update({"prediction_status": "error", "prediction_error": ""})
        files = artifacts.get(candidate_id, {})
        for kind in required_artifacts:
            row[kind] = str(files[kind]) if kind in files else ""
        missing = [kind for kind in required_artifacts if kind not in files]
        try:
            if missing:
                raise FileNotFoundError(f"missing artifacts: {','.join(missing)}")
            confidence = json.loads(files["confidence_json"].read_text(encoding="utf-8"))
            row.update(stage06.confidence_metrics(confidence))
            row.update(stage06.token_metrics(
                files["pae_npz"], files["plddt_npz"], receptor_length, task["_length"]
            ))
            row.update(stage06.geometry_metrics(
                files["structure_pdb"], receptor_chain, peptide_chain, task["_length"],
                pocket_numbers,
                float(geometry["contact_cutoff_angstrom"]),
                float(geometry["pocket_cutoff_angstrom"]),
                float(geometry["clash_cutoff_angstrom"]),
            ))
            row["prediction_status"] = "ok"
        except Exception as exc:
            row["prediction_error"] = f"{type(exc).__name__}: {exc}"
        for field in STRUCTURE_METRIC_FIELDS:
            row.setdefault(field, math.nan)
        rows.append(row)
    return rows


def compare_with_parents(audit_rows: list[dict], parent_index: dict[str, dict]) -> list[dict]:
    comparisons: list[dict] = []
    for row in audit_rows:
        parent_id = row["parent_candidate_id"]
        if parent_id not in parent_index:
            raise ValueError(f"Parent metrics not found for {parent_id}")
        parent = parent_index[parent_id]
        if parent["sequence"].strip().upper() != row["parent_sequence"]:
            raise ValueError(f"Parent sequence mismatch for {parent_id}")
        comparison = {field: row.get(field, "") for field in TASK_METADATA_FIELDS}
        comparison["prediction_status"] = row["prediction_status"]
        improved = worsened = tied = 0
        for metric, direction in COMPARISON_DIRECTIONS.items():
            parent_value = parse_metric(parent, metric, parent_id)
            optimized_value = float(row[metric])
            comparison[f"parent_{metric}"] = parent_value
            comparison[f"optimized_{metric}"] = optimized_value
            if row["prediction_status"] != "ok" or not math.isfinite(optimized_value):
                comparison[f"delta_{metric}"] = math.nan
                comparison[f"improved_{metric}"] = ""
                continue
            delta = optimized_value - parent_value
            is_tied = math.isclose(delta, 0.0, rel_tol=1e-9, abs_tol=1e-9)
            is_improved = (delta > 0) if direction == "higher" else (delta < 0)
            comparison[f"delta_{metric}"] = delta
            comparison[f"improved_{metric}"] = int(is_improved and not is_tied)
            if is_tied:
                tied += 1
            elif is_improved:
                improved += 1
            else:
                worsened += 1
        comparison["metrics_improved_count"] = improved
        comparison["metrics_worsened_count"] = worsened
        comparison["metrics_tied_count"] = tied
        comparisons.append(comparison)
    return comparisons


def summarize_by_round(rows: list[dict]) -> dict:
    return {
        str(round_index): stage06.summarize([
            row for row in rows if int(row["mutation_count"]) == round_index
        ])
        for round_index in (1, 2, 3)
    }


def summarize_paired(comparisons: list[dict]) -> dict:
    result = {}
    for round_index in (1, 2, 3):
        group = [
            row for row in comparisons
            if int(row["mutation_count"]) == round_index and row["prediction_status"] == "ok"
        ]
        metric_summary = {}
        for metric in COMPARISON_DIRECTIONS:
            deltas = np.asarray([
                float(row[f"delta_{metric}"])
                for row in group
                if math.isfinite(float(row[f"delta_{metric}"]))
            ])
            improved = sum(int(row[f"improved_{metric}"]) == 1 for row in group)
            tied = sum(
                math.isclose(float(row[f"delta_{metric}"]), 0.0, rel_tol=1e-9, abs_tol=1e-9)
                for row in group
                if math.isfinite(float(row[f"delta_{metric}"]))
            )
            metric_summary[metric] = {
                "direction": COMPARISON_DIRECTIONS[metric],
                "count": int(deltas.size),
                "optimized_better": improved,
                "optimized_better_fraction": float(improved / deltas.size) if deltas.size else None,
                "tied": tied,
                "parent_better": int(deltas.size) - improved - tied,
                "mean_raw_delta": float(deltas.mean()) if deltas.size else None,
                "median_raw_delta": float(np.median(deltas)) if deltas.size else None,
            }
        result[str(round_index)] = {"prediction_ok": len(group), "metrics": metric_summary}
    return result


def run(config_path: Path) -> Path:
    config, project_root = load_config(config_path)
    tasks_path = resolve_path(project_root, config["tasks_csv"])
    job_manifest_path = resolve_path(project_root, config["structure_job_manifest"])
    boltz_root = resolve_path(project_root, config["boltz_output_dir"])
    receptor_context_path = resolve_path(project_root, config["receptor_context"])
    parent_metrics_path = resolve_path(project_root, config["parent_metrics_csv"])
    output_dir = resolve_path(project_root, config["output_dir"])
    if not boltz_root.is_dir():
        raise FileNotFoundError(f"Boltz prediction directory does not exist: {boltz_root}")

    expected_tasks = int(config["expected_task_count"])
    tasks = read_tasks(tasks_path, expected_tasks)
    job_manifest = json.loads(job_manifest_path.read_text(encoding="utf-8"))
    if int(job_manifest.get("counts", {}).get("boltz_tasks", -1)) != expected_tasks:
        raise ValueError("Stage 08 manifest task count does not match the Stage 09 config")
    parent_index = read_parent_metrics(
        parent_metrics_path, int(config["expected_parent_metrics_count"])
    )
    receptor_context = json.loads(receptor_context_path.read_text(encoding="utf-8"))
    artifacts = stage06.index_prediction_files(boltz_root, int(config["model_index"]))
    unexpected_artifacts = sorted(set(artifacts) - {row["candidate_id"] for row in tasks})
    if unexpected_artifacts:
        raise ValueError(f"Unexpected candidate artifacts in Boltz output: {unexpected_artifacts[:5]}")

    audit_rows = audit_predictions(tasks, artifacts, receptor_context, config["geometry"])
    comparisons = compare_with_parents(audit_rows, parent_index)
    output_dir.mkdir(parents=True, exist_ok=True)
    audit_path = output_dir / "optimized_structure_metrics_all.csv"
    comparison_path = output_dir / "parent_variant_comparison_all.csv"
    atomic_write_csv(
        audit_path,
        [{field: serialize(row.get(field, "")) for field in AUDIT_FIELDS} for row in audit_rows],
        AUDIT_FIELDS,
    )
    atomic_write_csv(
        comparison_path,
        [{field: serialize(row.get(field, "")) for field in COMPARISON_FIELDS} for row in comparisons],
        COMPARISON_FIELDS,
    )

    prediction_ok = sum(row["prediction_status"] == "ok" for row in audit_rows)
    summary = {
        "schema_version": 1,
        "stage": "09_audit_optimized_structures",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selection_applied": False,
        "counts": {
            "input_tasks": len(tasks),
            "prediction_ok": prediction_ok,
            "prediction_error": len(tasks) - prediction_ok,
            "parent_families_represented": len({row["parent_candidate_id"] for row in tasks}),
            "by_mutation_count": dict(Counter(int(row["mutation_count"]) for row in tasks)),
        },
        "overall_metric_distributions": stage06.summarize(audit_rows),
        "metric_distributions_by_mutation_count": summarize_by_round(audit_rows),
        "paired_parent_comparisons_by_mutation_count": summarize_paired(comparisons),
    }
    summary_path = output_dir / "structure_audit_summary.json"
    atomic_write_json(summary_path, summary)

    manifest = {
        "schema_version": 1,
        "stage": "09_audit_optimized_structures",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": str(config_path.resolve()),
        "config_fingerprint": canonical_json_hash(config),
        "inputs": {
            "tasks_csv": str(tasks_path),
            "tasks_csv_sha256": sha256_file(tasks_path),
            "structure_job_manifest": str(job_manifest_path),
            "structure_job_manifest_sha256": sha256_file(job_manifest_path),
            "boltz_output_dir": str(boltz_root),
            "receptor_context": str(receptor_context_path),
            "parent_metrics_csv": str(parent_metrics_path),
            "parent_metrics_csv_sha256": sha256_file(parent_metrics_path),
            "model_index": int(config["model_index"]),
        },
        "policy": {
            "mode": "audit_only",
            "thresholds_applied": False,
            "ranking_applied": False,
            "family_winner_selected": False,
            "comparison_delta_definition": "optimized_minus_parent",
        },
        "geometry": config["geometry"],
        "counts": summary["counts"],
        "outputs": {
            "optimized_metrics_csv": str(audit_path),
            "parent_variant_comparison_csv": str(comparison_path),
            "audit_summary": str(summary_path),
        },
    }
    atomic_write_json(output_dir / "structure_audit_manifest.json", manifest)
    if bool(config["require_complete"]) and prediction_ok != len(tasks):
        raise RuntimeError(
            f"Structure audit found {len(tasks) - prediction_ok} incomplete/invalid predictions; "
            f"review {audit_path}"
        )
    return output_dir


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    try:
        output_dir = run(args.config.resolve())
    except Exception as exc:
        print(f"09_audit_optimized_structures FAILED: {exc}", file=sys.stderr)
        return 1
    print(f"09_audit_optimized_structures PASS: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
