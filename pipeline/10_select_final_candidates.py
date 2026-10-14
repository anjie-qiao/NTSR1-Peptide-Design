#!/usr/bin/env python3
"""Select final NTSR1 peptide candidates from parent and AfDesign families."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from pepmlm_pipeline_common import (
    atomic_write_csv,
    atomic_write_json,
    atomic_write_text,
    canonical_json_hash,
    sha256_file,
)


STRUCTURE_METRICS = [
    "pair_iptm_mean",
    "pair_iptm_min",
    "peptide_ptm",
    "peptide_mean_plddt",
    "peptide_q10_plddt",
    "interface_pae_mean",
    "interface_pae_q90",
    "peptide_contact_fraction",
    "interface_residue_pairs",
    "pocket_peptide_contact_fraction",
    "severe_clash_atom_pairs",
]

CORE_DIRECTIONS = {
    "pair_iptm_mean": "higher",
    "pair_iptm_min": "higher",
    "interface_pae_mean": "lower",
    "interface_pae_q90": "lower",
}

ABSOLUTE_RULES = {
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

MEANINGFUL_RULES = {
    "min_delta_pair_iptm_mean": ("pair_iptm_mean", "min"),
    "min_delta_pair_iptm_min": ("pair_iptm_min", "min"),
    "max_delta_interface_pae_mean": ("interface_pae_mean", "max"),
    "max_delta_interface_pae_q90": ("interface_pae_q90", "max"),
}

GUARDRAIL_RULES = {
    "min_delta_peptide_mean_plddt": ("peptide_mean_plddt", "min"),
    "min_delta_peptide_q10_plddt": ("peptide_q10_plddt", "min"),
    "min_delta_peptide_contact_fraction": ("peptide_contact_fraction", "min"),
    "min_delta_pocket_peptide_contact_fraction": (
        "pocket_peptide_contact_fraction", "min"
    ),
}

BASE_OUTPUT_FIELDS = [
    "candidate_id",
    "parent_candidate_id",
    "sequence",
    "parent_sequence",
    "length",
    "source_type",
    "mutation_count",
    "mutation_1",
    "mutation_2",
    "mutation_3",
    "sampling_profile",
    "parent_pseudo_perplexity",
    "parent_ppl_percentile_within_length",
    "parent_cluster_id",
    "afdesign_parent_loss",
    "afdesign_final_loss",
    "afdesign_loss_delta",
    "structure_pdb",
    "absolute_structure_pass",
    "absolute_structure_failure_reasons",
    "core_improvements_exact",
    "core_improvements_meaningful",
    "guardrail_pass",
    "guardrail_failure_reasons",
    "replacement_eligible",
    "paired_improvement_score",
    "family_selected",
    "family_selection_reason",
    "final_composite_score",
    "final_rank",
    "diversity_cluster_id",
    "diversity_cluster_size",
    "diversity_cluster_representative_id",
    "identity_to_diversity_representative",
    "is_diversity_cluster_representative",
    "final_selected",
    "final_selection_reason",
]

OUTPUT_FIELDS = (
    BASE_OUTPUT_FIELDS
    + STRUCTURE_METRICS
    + [f"delta_{metric}" for metric in STRUCTURE_METRICS]
    + [f"rank_percentile_{metric}" for metric in STRUCTURE_METRICS]
)


def resolve_path(project_root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def _require_keys(value: dict, required: set[str], label: str) -> None:
    missing = required - set(value)
    if missing:
        raise ValueError(f"{label} is missing keys: {sorted(missing)}")


def load_config(config_path: Path) -> tuple[dict, Path]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    _require_keys(
        config,
        {
            "parent_metrics_csv",
            "optimized_metrics_csv",
            "parent_variant_comparison_csv",
            "stage09_manifest",
            "output_dir",
            "expected_parent_count",
            "expected_variant_count",
            "final_count",
            "absolute_structure_thresholds",
            "paired_improvement",
            "ranking_metrics",
            "diversity",
        },
        "config",
    )
    if int(config["expected_parent_count"]) < 1 or int(config["expected_variant_count"]) < 1:
        raise ValueError("Expected row counts must be positive")
    if not 1 <= int(config["final_count"]) <= int(config["expected_parent_count"]):
        raise ValueError("final_count must be between 1 and expected_parent_count")

    absolute = config["absolute_structure_thresholds"]
    unknown = set(absolute) - set(ABSOLUTE_RULES)
    if unknown:
        raise ValueError(f"Unknown absolute thresholds: {sorted(unknown)}")
    if any(value is None for value in absolute.values()):
        raise ValueError("Stage 10 requires explicit, non-null absolute thresholds")

    paired = config["paired_improvement"]
    _require_keys(
        paired,
        {
            "min_exact_core_improvements",
            "min_meaningful_core_improvements",
            "meaningful_deltas",
            "guardrails",
            "paired_score_scales",
        },
        "paired_improvement",
    )
    if not 1 <= int(paired["min_exact_core_improvements"]) <= len(CORE_DIRECTIONS):
        raise ValueError("min_exact_core_improvements must be in [1, 4]")
    if not 1 <= int(paired["min_meaningful_core_improvements"]) <= len(CORE_DIRECTIONS):
        raise ValueError("min_meaningful_core_improvements must be in [1, 4]")
    if set(paired["meaningful_deltas"]) != set(MEANINGFUL_RULES):
        raise ValueError("meaningful_deltas must define all four supported core rules")
    if set(paired["guardrails"]) != set(GUARDRAIL_RULES):
        raise ValueError("guardrails must define all supported paired guardrails")
    if set(paired["paired_score_scales"]) != (
        set(CORE_DIRECTIONS) | set(GUARDRAIL_RULES[key][0] for key in GUARDRAIL_RULES)
    ):
        raise ValueError("paired_score_scales must define all core and guardrail metrics")
    if any(float(value) <= 0 for value in paired["paired_score_scales"].values()):
        raise ValueError("All paired_score_scales must be positive")

    ranking = config["ranking_metrics"]
    if not ranking:
        raise ValueError("ranking_metrics cannot be empty")
    names = [entry.get("metric") for entry in ranking]
    if len(names) != len(set(names)) or set(names) - set(STRUCTURE_METRICS):
        raise ValueError("ranking_metrics contain duplicates or unsupported metrics")
    for entry in ranking:
        if entry.get("direction") not in {"higher", "lower"}:
            raise ValueError("Every ranking direction must be higher or lower")
        if float(entry.get("weight", 0)) <= 0:
            raise ValueError("Every ranking weight must be positive")

    diversity = config["diversity"]
    _require_keys(
        diversity,
        {"method", "same_length_identity_threshold", "fill_to_final_count"},
        "diversity",
    )
    if diversity["method"] != "same_length_hamming":
        raise ValueError("Only same_length_hamming diversity is supported")
    if not 0 < float(diversity["same_length_identity_threshold"]) <= 1:
        raise ValueError("same_length_identity_threshold must be in (0, 1]")
    return config, config_path.resolve().parent.parent


def read_csv(path: Path, required: set[str], expected_count: int, label: str) -> list[dict]:
    if not path.is_file():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{label} is missing fields: {sorted(missing)}")
        rows = [dict(row) for row in reader]
    if len(rows) != expected_count:
        raise ValueError(f"Expected {expected_count} {label} rows; observed {len(rows)}")
    return rows


def parse_float(row: dict, field: str, candidate_id: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Missing/invalid {field} for {candidate_id}") from exc
    if not math.isfinite(value):
        raise ValueError(f"Non-finite {field} for {candidate_id}")
    return value


def parse_int(row: dict, field: str, candidate_id: str) -> int:
    try:
        return int(row[field])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Missing/invalid {field} for {candidate_id}") from exc


def index_rows(rows: list[dict], label: str) -> dict[str, dict]:
    index: dict[str, dict] = {}
    for row in rows:
        candidate_id = row["candidate_id"].strip()
        if not candidate_id or candidate_id in index:
            raise ValueError(f"Missing or duplicate candidate_id in {label}: {candidate_id!r}")
        index[candidate_id] = row
    return index


def absolute_failures(candidate: dict, thresholds: dict) -> list[str]:
    failures = []
    for rule, value in thresholds.items():
        metric, direction = ABSOLUTE_RULES[rule]
        observed = float(candidate[metric])
        if direction == "min" and observed < float(value):
            failures.append(rule)
        elif direction == "max" and observed > float(value):
            failures.append(rule)
    return failures


def exact_core_improvements(candidate: dict) -> tuple[int, dict[str, int]]:
    flags = {}
    for metric, direction in CORE_DIRECTIONS.items():
        delta = float(candidate[f"delta_{metric}"])
        flags[metric] = int(delta > 0) if direction == "higher" else int(delta < 0)
    return sum(flags.values()), flags


def meaningful_core_improvements(candidate: dict, thresholds: dict) -> tuple[int, dict[str, int]]:
    flags = {}
    for rule, threshold in thresholds.items():
        metric, direction = MEANINGFUL_RULES[rule]
        delta = float(candidate[f"delta_{metric}"])
        flags[metric] = (
            int(delta >= float(threshold)) if direction == "min" else int(delta <= float(threshold))
        )
    return sum(flags.values()), flags


def guardrail_failures(candidate: dict, guardrails: dict) -> list[str]:
    failures = []
    for rule, threshold in guardrails.items():
        metric, direction = GUARDRAIL_RULES[rule]
        delta = float(candidate[f"delta_{metric}"])
        if direction == "min" and delta < float(threshold):
            failures.append(rule)
        elif direction == "max" and delta > float(threshold):
            failures.append(rule)
    return failures


def paired_improvement_score(candidate: dict, scales: dict) -> float:
    score = 0.0
    for metric, direction in CORE_DIRECTIONS.items():
        delta = float(candidate[f"delta_{metric}"])
        score += delta / float(scales[metric]) if direction == "higher" else -delta / float(scales[metric])
    guard_metrics = sorted({GUARDRAIL_RULES[key][0] for key in GUARDRAIL_RULES})
    for metric in guard_metrics:
        score += float(candidate[f"delta_{metric}"]) / float(scales[metric])
    return score


def parent_candidate(row: dict) -> dict:
    candidate_id = row["candidate_id"].strip()
    sequence = row["sequence"].strip().upper()
    if row.get("prediction_status", "").strip() != "ok":
        raise ValueError(f"Parent structure is not valid: {candidate_id}")
    if parse_int(row, "length", candidate_id) != len(sequence):
        raise ValueError(f"Parent sequence length mismatch: {candidate_id}")
    candidate = {
        "candidate_id": candidate_id,
        "parent_candidate_id": candidate_id,
        "sequence": sequence,
        "parent_sequence": sequence,
        "length": len(sequence),
        "source_type": "parent",
        "mutation_count": 0,
        "mutation_1": "",
        "mutation_2": "",
        "mutation_3": "",
        "sampling_profile": row.get("sampling_profile", ""),
        "parent_pseudo_perplexity": row.get("pseudo_perplexity", ""),
        "parent_ppl_percentile_within_length": row.get("ppl_percentile_within_length", ""),
        "parent_cluster_id": row.get("cluster_id", ""),
        "afdesign_parent_loss": "",
        "afdesign_final_loss": "",
        "afdesign_loss_delta": "",
        "structure_pdb": row.get("structure_pdb", ""),
        "core_improvements_exact": "",
        "core_improvements_meaningful": "",
        "guardrail_pass": "",
        "guardrail_failure_reasons": "",
        "replacement_eligible": 0,
        "paired_improvement_score": "",
        "family_selected": 0,
        "family_selection_reason": "",
        "final_composite_score": "",
        "final_rank": "",
        "diversity_cluster_id": "",
        "diversity_cluster_size": "",
        "diversity_cluster_representative_id": "",
        "identity_to_diversity_representative": "",
        "is_diversity_cluster_representative": "",
        "final_selected": 0,
        "final_selection_reason": "",
    }
    for metric in STRUCTURE_METRICS:
        candidate[metric] = parse_float(row, metric, candidate_id)
        candidate[f"delta_{metric}"] = ""
        candidate[f"rank_percentile_{metric}"] = ""
    return candidate


def optimized_candidate(
    metric_row: dict,
    comparison_row: dict,
    parent: dict,
) -> dict:
    candidate_id = metric_row["candidate_id"].strip()
    parent_id = metric_row["parent_candidate_id"].strip()
    sequence = metric_row["sequence"].strip().upper()
    if candidate_id != comparison_row["candidate_id"].strip():
        raise ValueError(f"Stage 09 optimized/comparison ID mismatch: {candidate_id}")
    if parent_id != comparison_row["parent_candidate_id"].strip() or parent_id != parent["candidate_id"]:
        raise ValueError(f"Parent ID mismatch for optimized candidate: {candidate_id}")
    if metric_row.get("prediction_status", "").strip() != "ok":
        raise ValueError(f"Optimized structure is not valid: {candidate_id}")
    mutation_count = parse_int(metric_row, "mutation_count", candidate_id)
    if mutation_count not in {1, 2, 3} or candidate_id != f"{parent_id}__AFD{mutation_count}":
        raise ValueError(f"Invalid AfDesign round encoding: {candidate_id}")
    if parse_int(metric_row, "length", candidate_id) != len(sequence) or len(sequence) != parent["length"]:
        raise ValueError(f"Optimized sequence length mismatch: {candidate_id}")
    if metric_row.get("parent_sequence", "").strip().upper() != parent["sequence"]:
        raise ValueError(f"Optimized parent sequence mismatch: {candidate_id}")

    candidate = {
        "candidate_id": candidate_id,
        "parent_candidate_id": parent_id,
        "sequence": sequence,
        "parent_sequence": parent["sequence"],
        "length": len(sequence),
        "source_type": "afdesign",
        "mutation_count": mutation_count,
        "mutation_1": metric_row.get("mutation_1", ""),
        "mutation_2": metric_row.get("mutation_2", ""),
        "mutation_3": metric_row.get("mutation_3", ""),
        "sampling_profile": parent["sampling_profile"],
        "parent_pseudo_perplexity": parent["parent_pseudo_perplexity"],
        "parent_ppl_percentile_within_length": parent["parent_ppl_percentile_within_length"],
        "parent_cluster_id": parent["parent_cluster_id"],
        "afdesign_parent_loss": metric_row.get("afdesign_parent_loss", ""),
        "afdesign_final_loss": metric_row.get("afdesign_final_loss", ""),
        "afdesign_loss_delta": metric_row.get("afdesign_loss_delta", ""),
        "structure_pdb": metric_row.get("structure_pdb", ""),
        "family_selected": 0,
        "family_selection_reason": "",
        "final_composite_score": "",
        "final_rank": "",
        "diversity_cluster_id": "",
        "diversity_cluster_size": "",
        "diversity_cluster_representative_id": "",
        "identity_to_diversity_representative": "",
        "is_diversity_cluster_representative": "",
        "final_selected": 0,
        "final_selection_reason": "",
    }
    for metric in STRUCTURE_METRICS:
        optimized_value = parse_float(metric_row, metric, candidate_id)
        delta = parse_float(comparison_row, f"delta_{metric}", candidate_id)
        expected_delta = optimized_value - float(parent[metric])
        if not math.isclose(delta, expected_delta, rel_tol=1e-6, abs_tol=2e-6):
            raise ValueError(f"Parent delta mismatch for {candidate_id}: {metric}")
        candidate[metric] = optimized_value
        candidate[f"delta_{metric}"] = delta
        candidate[f"rank_percentile_{metric}"] = ""
    return candidate


def average_rank_percentiles(rows: list[dict], metric: str, direction: str) -> dict[str, float]:
    ordered = sorted(rows, key=lambda row: (float(row[metric]), row["candidate_id"]))
    n = len(ordered)
    result: dict[str, float] = {}
    index = 0
    while index < n:
        end = index + 1
        value = float(ordered[index][metric])
        while end < n and math.isclose(float(ordered[end][metric]), value, rel_tol=0, abs_tol=1e-12):
            end += 1
        average_zero_based_rank = (index + end - 1) / 2
        percentile = average_zero_based_rank / max(n - 1, 1)
        if direction == "lower":
            percentile = 1.0 - percentile
        for position in range(index, end):
            result[ordered[position]["candidate_id"]] = percentile
        index = end
    return result


def score_family_representatives(rows: list[dict], ranking_metrics: list[dict]) -> None:
    total_weight = sum(float(entry["weight"]) for entry in ranking_metrics)
    percentile_maps = {
        entry["metric"]: average_rank_percentiles(rows, entry["metric"], entry["direction"])
        for entry in ranking_metrics
    }
    for row in rows:
        score = 0.0
        for entry in ranking_metrics:
            metric = entry["metric"]
            percentile = percentile_maps[metric][row["candidate_id"]]
            row[f"rank_percentile_{metric}"] = percentile
            score += float(entry["weight"]) * percentile
        row["final_composite_score"] = score / total_weight
    rows.sort(key=lambda row: (-float(row["final_composite_score"]), row["candidate_id"]))
    for rank, row in enumerate(rows, start=1):
        row["final_rank"] = rank


def same_length_identity(sequence_a: str, sequence_b: str) -> float | None:
    if len(sequence_a) != len(sequence_b):
        return None
    return sum(a == b for a, b in zip(sequence_a, sequence_b)) / len(sequence_a)


def cluster_by_identity(rows: list[dict], threshold: float) -> list[dict]:
    clusters: list[dict] = []
    for row in rows:
        assignment = None
        identity_to_representative = None
        for cluster in clusters:
            identity = same_length_identity(row["sequence"], cluster["representative"]["sequence"])
            if identity is not None and identity >= threshold:
                assignment = cluster
                identity_to_representative = identity
                break
        if assignment is None:
            assignment = {"representative": row, "members": []}
            clusters.append(assignment)
            identity_to_representative = 1.0
        assignment["members"].append(row)
        row["identity_to_diversity_representative"] = identity_to_representative
    for index, cluster in enumerate(clusters, start=1):
        cluster_id = f"F10-C{index:04d}"
        representative = cluster["representative"]
        for row in cluster["members"]:
            row["diversity_cluster_id"] = cluster_id
            row["diversity_cluster_size"] = len(cluster["members"])
            row["diversity_cluster_representative_id"] = representative["candidate_id"]
            row["is_diversity_cluster_representative"] = int(row is representative)
    return clusters


def serialize(value: object) -> object:
    if isinstance(value, float):
        return "" if not math.isfinite(value) else f"{value:.8f}"
    return value


def output_rows(rows: list[dict]) -> list[dict]:
    return [
        {field: serialize(row.get(field, "")) for field in OUTPUT_FIELDS}
        for row in rows
    ]


def distribution(rows: list[dict], field: str) -> dict[str, int]:
    return dict(sorted(Counter(str(row[field]) for row in rows).items()))


def run(config_path: Path) -> Path:
    config, project_root = load_config(config_path)
    parent_path = resolve_path(project_root, config["parent_metrics_csv"])
    optimized_path = resolve_path(project_root, config["optimized_metrics_csv"])
    comparison_path = resolve_path(project_root, config["parent_variant_comparison_csv"])
    stage09_manifest_path = resolve_path(project_root, config["stage09_manifest"])
    output_dir = resolve_path(project_root, config["output_dir"])

    parent_required = {
        "candidate_id", "sequence", "length", "prediction_status", "structure_pdb"
    } | set(STRUCTURE_METRICS)
    optimized_required = {
        "candidate_id", "parent_candidate_id", "parent_sequence", "sequence", "length",
        "mutation_count", "prediction_status", "structure_pdb",
    } | set(STRUCTURE_METRICS)
    comparison_required = {
        "candidate_id", "parent_candidate_id"
    } | {f"delta_{metric}" for metric in STRUCTURE_METRICS}

    expected_parents = int(config["expected_parent_count"])
    expected_variants = int(config["expected_variant_count"])
    parent_rows = read_csv(parent_path, parent_required, expected_parents, "parent metrics")
    optimized_rows = read_csv(optimized_path, optimized_required, expected_variants, "optimized metrics")
    comparison_rows = read_csv(
        comparison_path, comparison_required, expected_variants, "parent-variant comparisons"
    )
    if not stage09_manifest_path.is_file():
        raise FileNotFoundError(f"Stage 09 manifest does not exist: {stage09_manifest_path}")
    stage09_manifest = json.loads(stage09_manifest_path.read_text(encoding="utf-8"))
    manifest_tasks = int(stage09_manifest.get("counts", {}).get("input_tasks", -1))
    manifest_ok = int(stage09_manifest.get("counts", {}).get("prediction_ok", -1))
    if manifest_tasks != expected_variants or manifest_ok != expected_variants:
        raise ValueError("Stage 09 manifest does not confirm a complete optimized-structure audit")

    parent_source_index = index_rows(parent_rows, "parent metrics")
    optimized_index = index_rows(optimized_rows, "optimized metrics")
    comparison_index = index_rows(comparison_rows, "parent-variant comparisons")
    if set(optimized_index) != set(comparison_index):
        raise ValueError("Optimized metrics and comparison candidate IDs do not match")

    parents = {candidate_id: parent_candidate(row) for candidate_id, row in parent_source_index.items()}
    absolute_thresholds = config["absolute_structure_thresholds"]
    for parent in parents.values():
        failures = absolute_failures(parent, absolute_thresholds)
        parent["absolute_structure_pass"] = int(not failures)
        parent["absolute_structure_failure_reasons"] = ";".join(failures)
    parent_failures = [row["candidate_id"] for row in parents.values() if not row["absolute_structure_pass"]]
    if parent_failures:
        raise ValueError(
            "Configured absolute thresholds reject parent structure-pass rows; "
            f"first failures={parent_failures[:5]}"
        )

    variants: list[dict] = []
    paired = config["paired_improvement"]
    for candidate_id in sorted(optimized_index):
        metric_row = optimized_index[candidate_id]
        parent_id = metric_row["parent_candidate_id"].strip()
        if parent_id not in parents:
            raise ValueError(f"Optimized variant references an unknown parent: {candidate_id}")
        candidate = optimized_candidate(
            metric_row, comparison_index[candidate_id], parents[parent_id]
        )
        absolute = absolute_failures(candidate, absolute_thresholds)
        exact_count, _ = exact_core_improvements(candidate)
        meaningful_count, _ = meaningful_core_improvements(
            candidate, paired["meaningful_deltas"]
        )
        guard_failures = guardrail_failures(candidate, paired["guardrails"])
        candidate["absolute_structure_pass"] = int(not absolute)
        candidate["absolute_structure_failure_reasons"] = ";".join(absolute)
        candidate["core_improvements_exact"] = exact_count
        candidate["core_improvements_meaningful"] = meaningful_count
        candidate["guardrail_pass"] = int(not guard_failures)
        candidate["guardrail_failure_reasons"] = ";".join(guard_failures)
        candidate["paired_improvement_score"] = paired_improvement_score(
            candidate, paired["paired_score_scales"]
        )
        candidate["replacement_eligible"] = int(
            not absolute
            and exact_count >= int(paired["min_exact_core_improvements"])
            and meaningful_count >= int(paired["min_meaningful_core_improvements"])
            and not guard_failures
        )
        variants.append(candidate)

    variants_by_parent: dict[str, list[dict]] = defaultdict(list)
    for variant in variants:
        variants_by_parent[variant["parent_candidate_id"]].append(variant)

    family_representatives: list[dict] = []
    for parent_id in sorted(parents):
        parent = parents[parent_id]
        eligible = [row for row in variants_by_parent.get(parent_id, []) if row["replacement_eligible"]]
        if not eligible:
            parent["family_selected"] = 1
            parent["family_selection_reason"] = "parent_retained_no_eligible_optimized_variant"
            family_representatives.append(parent)
            for variant in variants_by_parent.get(parent_id, []):
                variant["family_selection_reason"] = "optimized_variant_not_replacement_eligible"
            continue
        eligible.sort(key=lambda row: (
            -int(row["core_improvements_meaningful"]),
            -int(row["core_improvements_exact"]),
            -float(row["paired_improvement_score"]),
            row["candidate_id"],
        ))
        winner = eligible[0]
        winner["family_selected"] = 1
        winner["family_selection_reason"] = "optimized_variant_replaced_parent"
        parent["family_selection_reason"] = "parent_replaced_by_eligible_optimized_variant"
        family_representatives.append(winner)
        for variant in variants_by_parent.get(parent_id, []):
            if variant is winner:
                continue
            variant["family_selection_reason"] = (
                "eligible_optimized_variant_not_best_in_family"
                if variant["replacement_eligible"]
                else "optimized_variant_not_replacement_eligible"
            )

    if len(family_representatives) != expected_parents:
        raise RuntimeError("Family selection did not produce exactly one representative per parent")
    score_family_representatives(family_representatives, config["ranking_metrics"])
    clusters = cluster_by_identity(
        family_representatives,
        float(config["diversity"]["same_length_identity_threshold"]),
    )
    cluster_representatives = [cluster["representative"] for cluster in clusters]
    cluster_representatives.sort(key=lambda row: int(row["final_rank"]))
    final_count = int(config["final_count"])
    final_rows = cluster_representatives[:final_count]
    if len(final_rows) < final_count:
        if not bool(config["diversity"]["fill_to_final_count"]):
            raise RuntimeError(
                f"Diversity clustering produced only {len(final_rows)} representatives; "
                f"requested {final_count}"
            )
        selected_ids = {row["candidate_id"] for row in final_rows}
        fillers = [row for row in family_representatives if row["candidate_id"] not in selected_ids]
        final_rows.extend(fillers[: final_count - len(final_rows)])

    final_ids = {row["candidate_id"] for row in final_rows}
    cluster_rep_ids = {row["candidate_id"] for row in cluster_representatives}
    for row in family_representatives:
        row["final_selected"] = int(row["candidate_id"] in final_ids)
        if row["final_selected"] and row["candidate_id"] in cluster_rep_ids:
            row["final_selection_reason"] = "top_ranked_diversity_cluster_representative"
        elif row["final_selected"]:
            row["final_selection_reason"] = "score_fill_after_diversity_clustering"
        elif row["candidate_id"] not in cluster_rep_ids:
            row["final_selection_reason"] = "sequence_redundant_with_higher_ranked_representative"
        else:
            row["final_selection_reason"] = "below_final_rank_cutoff"
    final_rows.sort(key=lambda row: int(row["final_rank"]))

    output_dir.mkdir(parents=True, exist_ok=True)
    all_candidates = list(parents.values()) + variants
    all_candidates.sort(key=lambda row: (
        row["parent_candidate_id"], int(row["mutation_count"]), row["candidate_id"]
    ))
    eligible_variants = [row for row in variants if row["replacement_eligible"]]
    eligible_variants.sort(key=lambda row: (
        row["parent_candidate_id"], int(row["mutation_count"]), row["candidate_id"]
    ))

    evaluation_path = output_dir / "candidate_evaluation_all.csv"
    replacement_path = output_dir / "optimized_replacement_candidates.csv"
    family_path = output_dir / "family_representatives.csv"
    final_path = output_dir / "final_candidates_100.csv"
    final_fasta_path = output_dir / "final_candidates_100.fasta"
    atomic_write_csv(evaluation_path, output_rows(all_candidates), OUTPUT_FIELDS)
    atomic_write_csv(replacement_path, output_rows(eligible_variants), OUTPUT_FIELDS)
    atomic_write_csv(family_path, output_rows(family_representatives), OUTPUT_FIELDS)
    atomic_write_csv(final_path, output_rows(final_rows), OUTPUT_FIELDS)
    atomic_write_text(
        final_fasta_path,
        "".join(
            f">{row['candidate_id']} rank={row['final_rank']} source={row['source_type']} "
            f"parent={row['parent_candidate_id']} score={float(row['final_composite_score']):.6f}\n"
            f"{row['sequence']}\n"
            for row in final_rows
        ),
    )

    summary = {
        "schema_version": 1,
        "stage": "10_select_final_candidates",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "counts": {
            "input_parents": len(parents),
            "input_optimized_variants": len(variants),
            "optimized_absolute_structure_pass": sum(
                int(row["absolute_structure_pass"]) for row in variants
            ),
            "optimized_replacement_eligible": len(eligible_variants),
            "families_with_eligible_replacement": len(
                {row["parent_candidate_id"] for row in eligible_variants}
            ),
            "optimized_family_winners": sum(
                row["source_type"] == "afdesign" for row in family_representatives
            ),
            "parent_family_winners": sum(
                row["source_type"] == "parent" for row in family_representatives
            ),
            "family_representatives": len(family_representatives),
            "diversity_clusters": len(clusters),
            "final_candidates": len(final_rows),
        },
        "final_by_source_type": distribution(final_rows, "source_type"),
        "final_by_mutation_count": distribution(final_rows, "mutation_count"),
        "final_by_length": distribution(final_rows, "length"),
        "optimized_absolute_pass_by_mutation_count": dict(sorted(Counter(
            str(row["mutation_count"])
            for row in variants if row["absolute_structure_pass"]
        ).items())),
        "replacement_eligible_by_mutation_count": dict(sorted(Counter(
            str(row["mutation_count"])
            for row in eligible_variants
        ).items())),
        "policy": {
            "parent_is_family_fallback": True,
            "afdesign_loss_used_for_selection": False,
            "absolute_structure_thresholds": absolute_thresholds,
            "paired_improvement": paired,
            "ranking_metrics": config["ranking_metrics"],
            "diversity": config["diversity"],
        },
    }
    summary_path = output_dir / "final_selection_summary.json"
    atomic_write_json(summary_path, summary)

    manifest = {
        "schema_version": 1,
        "stage": "10_select_final_candidates",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": str(config_path.resolve()),
        "config_fingerprint": canonical_json_hash(config),
        "inputs": {
            "parent_metrics_csv": str(parent_path),
            "parent_metrics_csv_sha256": sha256_file(parent_path),
            "optimized_metrics_csv": str(optimized_path),
            "optimized_metrics_csv_sha256": sha256_file(optimized_path),
            "parent_variant_comparison_csv": str(comparison_path),
            "parent_variant_comparison_csv_sha256": sha256_file(comparison_path),
            "stage09_manifest": str(stage09_manifest_path),
            "stage09_manifest_sha256": sha256_file(stage09_manifest_path),
        },
        "counts": summary["counts"],
        "outputs": {
            "candidate_evaluation_all": str(evaluation_path),
            "optimized_replacement_candidates": str(replacement_path),
            "family_representatives": str(family_path),
            "final_candidates_csv": str(final_path),
            "final_candidates_fasta": str(final_fasta_path),
            "final_selection_summary": str(summary_path),
        },
    }
    atomic_write_json(output_dir / "final_selection_manifest.json", manifest)
    return output_dir


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    try:
        output_dir = run(args.config.resolve())
    except Exception as exc:
        print(f"10_select_final_candidates FAILED: {exc}", file=sys.stderr)
        return 1
    print(f"10_select_final_candidates PASS: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
