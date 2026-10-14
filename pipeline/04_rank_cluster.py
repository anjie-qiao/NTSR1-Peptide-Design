#!/usr/bin/env python3
"""Build a rule-defined, non-redundant PepMLM core candidate library."""

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
    CANONICAL_AAS,
    atomic_write_csv,
    atomic_write_json,
    atomic_write_text,
    canonical_json_hash,
    sha256_file,
)


ANNOTATION_FIELDS = [
    "ppl_rank_within_length",
    "ppl_percentile_within_length",
    "ppl_pool_pass",
    "cluster_id",
    "cluster_identity_threshold",
    "cluster_size",
    "cluster_representative_id",
    "is_cluster_representative",
    "selected_for_core",
    "selection_reason",
]


def resolve_path(project_root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def load_config(config_path: Path) -> tuple[dict, Path]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    required = {
        "input_csv", "input_manifest", "output_dir", "output_basename",
        "ppl_keep_fraction", "cluster_identity_threshold", "cluster_sensitivity_thresholds",
    }
    missing = required - set(config)
    if missing:
        raise ValueError(f"Missing config keys: {sorted(missing)}")
    if not 0 < float(config["ppl_keep_fraction"]) <= 1:
        raise ValueError("ppl_keep_fraction must be in (0, 1]")
    thresholds = [float(value) for value in config["cluster_sensitivity_thresholds"]]
    thresholds.append(float(config["cluster_identity_threshold"]))
    if any(not 0 < value <= 1 for value in thresholds):
        raise ValueError("All identity thresholds must be in (0, 1]")
    return config, config_path.resolve().parent.parent


def read_candidates(path: Path) -> tuple[list[str], list[dict]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {path}")
        required = {
            "candidate_id", "sequence", "length", "sampling_profile",
            "mean_log_probability", "pseudo_perplexity",
        }
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Input CSV is missing fields: {sorted(missing)}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise ValueError("Input candidate CSV is empty")
    seen_ids: set[str] = set()
    seen_sequences: set[str] = set()
    for row in rows:
        candidate_id = row["candidate_id"].strip()
        sequence = row["sequence"].strip().upper()
        if not candidate_id or candidate_id in seen_ids:
            raise ValueError(f"Missing or duplicate candidate_id: {candidate_id!r}")
        if not sequence or set(sequence) - set(CANONICAL_AAS) or sequence in seen_sequences:
            raise ValueError(f"Invalid or duplicate sequence for {candidate_id}")
        try:
            length = int(row["length"])
            ppl = float(row["pseudo_perplexity"])
            mean_log_probability = float(row["mean_log_probability"])
        except ValueError as exc:
            raise ValueError(f"Invalid numeric value for {candidate_id}") from exc
        if length != len(sequence) or not math.isfinite(ppl) or ppl <= 0:
            raise ValueError(f"Invalid length or pseudo-perplexity for {candidate_id}")
        if not math.isfinite(mean_log_probability):
            raise ValueError(f"Invalid mean_log_probability for {candidate_id}")
        row["candidate_id"] = candidate_id
        row["sequence"] = sequence
        row["_length"] = length
        row["_ppl"] = ppl
        row["_mean_log_probability"] = mean_log_probability
        seen_ids.add(candidate_id)
        seen_sequences.add(sequence)
    return list(reader.fieldnames), rows


def rank_within_length(rows: list[dict], keep_fraction: float) -> dict[int, int]:
    by_length: dict[int, list[dict]] = defaultdict(list)
    for row in rows:
        by_length[row["_length"]].append(row)
    keep_counts = {}
    for length, group in sorted(by_length.items()):
        group.sort(key=lambda row: (row["_ppl"], -row["_mean_log_probability"], row["candidate_id"]))
        keep_count = min(len(group), max(1, math.ceil(len(group) * keep_fraction)))
        keep_counts[length] = keep_count
        denominator = max(len(group) - 1, 1)
        for rank, row in enumerate(group, start=1):
            row["ppl_rank_within_length"] = rank
            row["ppl_percentile_within_length"] = (rank - 1) / denominator
            row["ppl_pool_pass"] = int(rank <= keep_count)
    return keep_counts


def required_matches(length: int, identity_threshold: float) -> int:
    return math.ceil(length * identity_threshold - 1e-12)


def same_length_identity(sequence_a: str, sequence_b: str) -> float:
    if len(sequence_a) != len(sequence_b):
        raise ValueError("Identity is only defined here for equal-length peptides")
    return sum(a == b for a, b in zip(sequence_a, sequence_b)) / len(sequence_a)


def greedy_clusters(rows: list[dict], identity_threshold: float) -> list[dict]:
    """Greedy equal-length Hamming clusters, ordered by best within-length PPL first."""
    by_length: dict[int, list[dict]] = defaultdict(list)
    for row in rows:
        by_length[row["_length"]].append(row)
    clusters = []
    for length, group in sorted(by_length.items()):
        group.sort(key=lambda row: (row["_ppl"], -row["_mean_log_probability"], row["candidate_id"]))
        length_clusters = []
        minimum_matches = required_matches(length, identity_threshold)
        for row in group:
            assignment = None
            for cluster in length_clusters:
                matches = sum(
                    aa == bb for aa, bb in zip(row["sequence"], cluster["representative"]["sequence"])
                )
                if matches >= minimum_matches:
                    assignment = cluster
                    break
            if assignment is None:
                assignment = {"representative": row, "members": []}
                length_clusters.append(assignment)
            assignment["members"].append(row)
        for index, cluster in enumerate(length_clusters, start=1):
            cluster["cluster_id"] = f"L{length:02d}-C{index:04d}"
            cluster["length"] = length
            cluster["identity_threshold"] = identity_threshold
            clusters.append(cluster)
    return clusters


def annotate_clusters(clusters: list[dict]) -> None:
    for cluster in clusters:
        representative = cluster["representative"]
        for row in cluster["members"]:
            row["cluster_id"] = cluster["cluster_id"]
            row["cluster_identity_threshold"] = cluster["identity_threshold"]
            row["cluster_size"] = len(cluster["members"])
            row["cluster_representative_id"] = representative["candidate_id"]
            row["is_cluster_representative"] = int(row is representative)


def cluster_summary_rows(clusters: list[dict], selected_ids: set[str]) -> list[dict]:
    output = []
    for cluster in clusters:
        representative = cluster["representative"]
        output.append({
            "cluster_id": cluster["cluster_id"],
            "length": cluster["length"],
            "identity_threshold": f"{cluster['identity_threshold']:.3f}",
            "cluster_size": len(cluster["members"]),
            "representative_id": representative["candidate_id"],
            "representative_sequence": representative["sequence"],
            "representative_ppl": f"{representative['_ppl']:.8f}",
            "representative_ppl_percentile": f"{representative['ppl_percentile_within_length']:.8f}",
            "representative_profile": representative["sampling_profile"],
            "eligible_for_core": representative["ppl_pool_pass"],
            "selected_for_core": int(representative["candidate_id"] in selected_ids),
        })
    return output


def sensitivity_rows(rows: list[dict], thresholds: list[float]) -> list[dict]:
    output = []
    for threshold in sorted(set(thresholds)):
        clusters = greedy_clusters(rows, threshold)
        sizes = [len(cluster["members"]) for cluster in clusters]
        per_length = Counter(cluster["length"] for cluster in clusters)
        for length in sorted({row["_length"] for row in rows}):
            length_sizes = [len(cluster["members"]) for cluster in clusters if cluster["length"] == length]
            output.append({
                "identity_threshold": f"{threshold:.3f}",
                "scope": f"length_{length}",
                "length": length,
                "cluster_count": per_length[length],
                "largest_cluster_size": max(length_sizes),
            })
        output.append({
            "identity_threshold": f"{threshold:.3f}",
            "scope": "all_lengths",
            "length": "",
            "cluster_count": len(clusters),
            "largest_cluster_size": max(sizes),
        })
    return output


def write_optional_figures(output_dir: Path, length_rows: list[dict], sensitivity: list[dict]) -> dict:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return {"generated": False, "reason": "matplotlib_not_installed", "files": []}

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 12,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 1.5,
        "legend.frameon": False,
        "svg.fonttype": "none",
    })
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    lengths = [int(row["length"]) for row in length_rows]
    ppl_scored = [int(row["ppl_scored"]) for row in length_rows]
    ppl_pool = [int(row["ppl_pool"]) for row in length_rows]
    selected = [int(row["core_representatives"]) for row in length_rows]
    fig, ax = plt.subplots(figsize=(13, 5.5))
    ax.plot(lengths, ppl_scored, marker="o", linewidth=2.4, color="#CFCECE", label="PPL scored")
    ax.plot(lengths, ppl_pool, marker="o", linewidth=2.4, color="#3775BA", label="Top PPL pool")
    ax.plot(lengths, selected, marker="o", linewidth=2.4, color="#0F4D92", label="Core representatives")
    ax.set(xlabel="Peptide length (aa)", ylabel="Candidate count", xticks=lengths)
    ax.legend(ncol=3, loc="upper right")
    fig.tight_layout(pad=2)
    files = []
    for suffix in ("png", "svg"):
        path = figure_dir / f"length_selection_funnel.{suffix}"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        files.append(str(path))
    plt.close(fig)

    totals = [row for row in sensitivity if row["scope"] == "all_lengths"]
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.bar(
        [float(row["identity_threshold"]) for row in totals],
        [int(row["cluster_count"]) for row in totals],
        width=0.055, color="#3775BA", edgecolor="black", linewidth=1.2,
    )
    ax.set(xlabel="Sequence identity threshold", ylabel="Number of clusters")
    fig.tight_layout(pad=2)
    for suffix in ("png", "svg"):
        path = figure_dir / f"cluster_threshold_sensitivity.{suffix}"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        files.append(str(path))
    plt.close(fig)
    return {"generated": True, "reason": "", "files": files}


def run(config_path: Path, skip_plots: bool = False) -> Path:
    config, project_root = load_config(config_path)
    input_path = resolve_path(project_root, config["input_csv"])
    manifest_path = resolve_path(project_root, config["input_manifest"])
    output_dir = resolve_path(project_root, config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    input_fields, rows = read_candidates(input_path)
    input_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected_count = input_manifest.get("counts", {}).get("ppl_computed")
    if expected_count is not None and int(expected_count) != len(rows):
        raise ValueError(f"Input manifest expects {expected_count} PPL rows; observed {len(rows)}")

    keep_counts = rank_within_length(rows, float(config["ppl_keep_fraction"]))
    threshold = float(config["cluster_identity_threshold"])
    clusters = greedy_clusters(rows, threshold)
    annotate_clusters(clusters)
    selected = [
        cluster["representative"] for cluster in clusters
        if cluster["representative"]["ppl_pool_pass"]
    ]
    selected_ids = {row["candidate_id"] for row in selected}
    if len(selected_ids) != len(selected):
        raise RuntimeError("Core selection contains duplicate candidate IDs")

    for row in rows:
        row["selected_for_core"] = int(row["candidate_id"] in selected_ids)
        if row["selected_for_core"]:
            row["selection_reason"] = "within_length_ppl_pool_cluster_representative"
        elif not row["ppl_pool_pass"]:
            row["selection_reason"] = "outside_within_length_ppl_pool"
        else:
            row["selection_reason"] = "sequence_cluster_redundant"

    output_fields = input_fields + [field for field in ANNOTATION_FIELDS if field not in input_fields]
    clean_rows = []
    for row in rows:
        clean = {field: row.get(field, "") for field in output_fields}
        clean["ppl_percentile_within_length"] = f"{row['ppl_percentile_within_length']:.8f}"
        clean["cluster_identity_threshold"] = f"{row['cluster_identity_threshold']:.3f}"
        clean_rows.append(clean)
    clean_rows.sort(key=lambda row: (int(row["length"]), int(row["ppl_rank_within_length"])))
    atomic_write_csv(output_dir / "ranked_all.csv", clean_rows, output_fields)

    selected.sort(key=lambda row: (row["_length"], row["_ppl"], row["candidate_id"]))
    selected_clean = []
    for row in selected:
        clean = {field: row.get(field, "") for field in output_fields}
        clean["ppl_percentile_within_length"] = f"{row['ppl_percentile_within_length']:.8f}"
        clean["cluster_identity_threshold"] = f"{row['cluster_identity_threshold']:.3f}"
        selected_clean.append(clean)
    output_basename = str(config["output_basename"])
    atomic_write_csv(output_dir / f"{output_basename}.csv", selected_clean, output_fields)
    atomic_write_text(
        output_dir / f"{output_basename}.fasta",
        "".join(
            f">{row['candidate_id']} length={row['_length']} profile={row['sampling_profile']} "
            f"ppl={row['_ppl']:.8f} cluster={row['cluster_id']}\n{row['sequence']}\n"
            for row in selected
        ),
    )

    cluster_fields = [
        "cluster_id", "length", "identity_threshold", "cluster_size", "representative_id",
        "representative_sequence", "representative_ppl", "representative_ppl_percentile",
        "representative_profile", "eligible_for_core", "selected_for_core",
    ]
    atomic_write_csv(
        output_dir / "cluster_summary.csv",
        cluster_summary_rows(clusters, selected_ids), cluster_fields,
    )

    selected_by_length = Counter(row["_length"] for row in selected)
    input_by_length = Counter(row["_length"] for row in rows)
    cluster_by_length = Counter(cluster["length"] for cluster in clusters)
    length_rows = []
    for length in sorted(input_by_length):
        profile_counts = Counter(
            row["sampling_profile"] for row in selected if row["_length"] == length
        )
        length_rows.append({
            "length": length,
            "ppl_scored": input_by_length[length],
            "ppl_pool": keep_counts[length],
            "clusters_all": cluster_by_length[length],
            "core_representatives": selected_by_length[length],
            "core_profile_counts": json.dumps(dict(sorted(profile_counts.items())), sort_keys=True),
        })
    length_fields = [
        "length", "ppl_scored", "ppl_pool", "clusters_all",
        "core_representatives", "core_profile_counts",
    ]
    atomic_write_csv(output_dir / "length_summary.csv", length_rows, length_fields)

    pool_rows = [row for row in rows if row["ppl_pool_pass"]]
    sensitivity = sensitivity_rows(
        pool_rows, [float(value) for value in config["cluster_sensitivity_thresholds"]]
    )
    sensitivity_fields = [
        "identity_threshold", "scope", "length", "cluster_count", "largest_cluster_size",
    ]
    atomic_write_csv(output_dir / "cluster_threshold_sensitivity.csv", sensitivity, sensitivity_fields)
    plot_status = (
        {"generated": False, "reason": "disabled_by_command", "files": []}
        if skip_plots or not bool(config.get("make_plots", True))
        else write_optional_figures(output_dir, length_rows, sensitivity)
    )

    profile_totals = Counter(row["sampling_profile"] for row in selected)
    selection_reason_totals = Counter(row["selection_reason"] for row in selected)
    selection_manifest = {
        "schema_version": 1,
        "stage": "04_rank_cluster",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": str(config_path.resolve()),
        "config_fingerprint": canonical_json_hash(config),
        "input": {
            "ppl_csv": str(input_path),
            "ppl_csv_sha256": sha256_file(input_path),
            "filter_manifest": str(manifest_path),
            "candidate_count": len(rows),
        },
        "method": {
            "ppl_ranking": "ascending within each peptide length",
            "ppl_keep_fraction": float(config["ppl_keep_fraction"]),
            "clustering": "greedy equal-length positional identity (Hamming identity)",
            "cluster_identity_threshold": threshold,
            "representative": "lowest PPL, then highest mean log probability, then candidate ID",
            "fixed_total_count": False,
            "length_quotas": False,
            "backfill_outside_ppl_pool": False,
        },
        "counts": {
            "input_ppl_scored": len(rows),
            "ppl_pool": sum(keep_counts.values()),
            "clusters_all": len(clusters),
            "core_cluster_representatives": len(selected),
        },
        "selected_by_length": {str(key): value for key, value in sorted(selected_by_length.items())},
        "selected_by_profile": dict(sorted(profile_totals.items())),
        "selected_by_reason": dict(sorted(selection_reason_totals.items())),
        "figures": plot_status,
        "outputs": {
            "ranked_all": str(output_dir / "ranked_all.csv"),
            "cluster_summary": str(output_dir / "cluster_summary.csv"),
            "length_summary": str(output_dir / "length_summary.csv"),
            "threshold_sensitivity": str(output_dir / "cluster_threshold_sensitivity.csv"),
            "candidates_csv": str(output_dir / f"{output_basename}.csv"),
            "candidates_fasta": str(output_dir / f"{output_basename}.fasta"),
        },
    }
    atomic_write_json(output_dir / "selection_manifest.json", selection_manifest)
    return output_dir


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--skip-plots", action="store_true")
    args = parser.parse_args()
    try:
        output_dir = run(args.config, skip_plots=args.skip_plots)
    except Exception as exc:
        print(f"04_rank_cluster FAILED: {exc}", file=sys.stderr)
        return 1
    print(f"04_rank_cluster PASS: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
