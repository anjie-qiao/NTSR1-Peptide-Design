#!/usr/bin/env python3
"""Filter PepMLM candidates and compute microbatched pseudo-perplexity.

Pseudo-perplexity follows the upstream generation.py definition exactly: each
binder residue is masked once and cross-entropy is averaged across binder
positions.  Rows from many candidates are streamed through a configurable GPU
microbatch to control memory use.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from pepmlm_model_runtime import load_local_model
from pepmlm_pipeline_common import (
    CANONICAL_AAS,
    allocate_integer,
    atomic_write_csv,
    atomic_write_json,
    atomic_write_text,
    canonical_json_hash,
    read_single_fasta,
    sequence_properties,
    sha256_file,
)


PROPERTY_FIELDS = [
    "net_charge_ph7_4",
    "hydrophobic_fraction",
    "max_hydrophobic_run",
    "max_homopolymer_run",
    "normalized_entropy",
    "cysteine_count",
]


def load_config(config_path: Path) -> tuple[dict, Path]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    required = {"input_csv", "target_fasta", "model_path", "output_dir", "hard_filters"}
    missing = required - set(config)
    if missing:
        raise ValueError(f"Missing config keys: {sorted(missing)}")
    project_root = config_path.resolve().parent.parent
    return config, project_root


def resolve_path(project_root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def read_rows(path: Path) -> tuple[list[str], list[dict]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {path}")
        rows = list(reader)
        return list(reader.fieldnames), rows


def filter_reason(sequence: str, row: dict, thresholds: dict, seen: set[str]) -> tuple[list[str], dict]:
    reasons = []
    if not sequence or set(sequence) - set(CANONICAL_AAS):
        reasons.append("noncanonical_or_empty")
        return reasons, {}
    try:
        recorded_length = int(row.get("length", len(sequence)))
    except ValueError:
        recorded_length = -1
    if recorded_length != len(sequence):
        reasons.append("recorded_length_mismatch")
    if sequence in seen:
        reasons.append("duplicate_sequence")
    else:
        seen.add(sequence)

    props = sequence_properties(sequence)
    if props["net_charge_ph7_4"] < float(thresholds["min_net_charge"]):
        reasons.append("net_charge_too_low")
    if props["net_charge_ph7_4"] > float(thresholds["max_net_charge"]):
        reasons.append("net_charge_too_high")
    if props["hydrophobic_fraction"] > float(thresholds["max_hydrophobic_fraction"]):
        reasons.append("hydrophobic_fraction_too_high")
    if props["max_hydrophobic_run"] > int(thresholds["max_hydrophobic_run"]):
        reasons.append("hydrophobic_run_too_long")
    if props["max_homopolymer_run"] > int(thresholds["max_homopolymer_run"]):
        reasons.append("homopolymer_run_too_long")
    if props["normalized_entropy"] < float(thresholds["min_normalized_entropy"]):
        reasons.append("sequence_complexity_too_low")
    if props["cysteine_count"] > int(thresholds["max_cysteine_count"]):
        reasons.append("too_many_cysteines")
    return reasons, props


def stratified_ppl_selection(rows: list[dict], max_per_length: int | None) -> list[dict]:
    if max_per_length is None or max_per_length <= 0:
        return list(rows)
    by_length: dict[int, list[dict]] = defaultdict(list)
    for row in rows:
        by_length[int(row["length"])].append(row)

    selected = []
    for length in sorted(by_length):
        length_rows = by_length[length]
        if len(length_rows) <= max_per_length:
            selected.extend(length_rows)
            continue
        by_profile: dict[str, list[dict]] = defaultdict(list)
        for row in length_rows:
            by_profile[row.get("sampling_profile", "unknown")].append(row)
        profile_names = sorted(by_profile)
        counts = allocate_integer(max_per_length, [len(by_profile[name]) for name in profile_names])
        for name, count in zip(profile_names, counts):
            group = sorted(
                by_profile[name],
                key=lambda row: float(row.get("mean_log_probability", "-inf")),
                reverse=True,
            )
            selected.extend(group[:count])
    return selected


def build_masked_rows(torch, tokenizer, target_sequence: str, candidate_rows: list[dict]):
    """Yield (input_ids, mask_position, target_token, candidate_index)."""
    for candidate_index, row in enumerate(candidate_rows):
        binder = row["sequence"]
        encoded = tokenizer.encode(target_sequence + binder, return_tensors="pt")[0].cpu()
        binder_positions = range(len(encoded) - len(binder) - 1, len(encoded) - 1)
        if len(list(binder_positions)) != len(binder):
            raise RuntimeError("Could not locate binder positions in tokenized sequence")
        for position in range(len(encoded) - len(binder) - 1, len(encoded) - 1):
            masked = encoded.clone()
            target_token = int(masked[position].item())
            masked[position] = tokenizer.mask_token_id
            yield masked, position, target_token, candidate_index


def compute_pseudo_perplexities(
    torch,
    tokenizer,
    model,
    device: str,
    target_sequence: str,
    candidate_rows: list[dict],
    row_batch_size: int,
) -> list[float]:
    import torch.nn.functional as F

    loss_sums = [0.0] * len(candidate_rows)
    loss_counts = [0] * len(candidate_rows)
    pending = []

    def flush(batch):
        if not batch:
            return
        input_ids = torch.stack([item[0] for item in batch]).to(device)
        attention_mask = torch.ones_like(input_ids, device=device)
        positions = torch.tensor([item[1] for item in batch], device=device)
        targets = torch.tensor([item[2] for item in batch], device=device)
        candidate_indices = [item[3] for item in batch]
        with torch.inference_mode():
            logits = model(input_ids=input_ids, attention_mask=attention_mask).logits
            row_indices = torch.arange(len(batch), device=device)
            masked_logits = logits[row_indices, positions].float()
            losses = F.cross_entropy(masked_logits, targets, reduction="none").cpu().tolist()
        for candidate_index, loss in zip(candidate_indices, losses):
            loss_sums[candidate_index] += float(loss)
            loss_counts[candidate_index] += 1

    for masked_row in build_masked_rows(torch, tokenizer, target_sequence, candidate_rows):
        if pending and len(masked_row[0]) != len(pending[0][0]):
            # Candidate lengths can differ; never stack token rows of different widths.
            flush(pending)
            pending = []
        pending.append(masked_row)
        if len(pending) >= row_batch_size:
            flush(pending)
            pending = []
    flush(pending)

    ppls = []
    for row, total, count in zip(candidate_rows, loss_sums, loss_counts):
        if count != len(row["sequence"]):
            raise RuntimeError(f"PPL mask count mismatch for {row['candidate_id']}: {count}")
        ppls.append(math.exp(total / count))
    return ppls


def validate_against_upstream(
    project_root: Path,
    torch,
    tokenizer,
    model,
    target_sequence: str,
    rows: list[dict],
    production_ppls: list[float],
    count: int,
) -> list[dict]:
    if count <= 0:
        return []
    upstream_path = Path(__file__).resolve().parents[1] / "pepmlm-main" / "scripts" / "generation.py"
    if not upstream_path.is_file():
        raise FileNotFoundError("Optional upstream comparison is disabled by default. Restore an authorized generation.py to " + str(upstream_path))
    spec = importlib.util.spec_from_file_location("pepmlm_upstream_generation", upstream_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load upstream generation.py: {upstream_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # The downloaded upstream file omits imports; inject only what its PPL function uses.
    module.torch = torch
    module.np = __import__("numpy")

    comparisons = []
    order = sorted(range(len(rows)), key=lambda idx: len(rows[idx]["sequence"]))[:count]
    for idx in order:
        reference = float(module.compute_pseudo_perplexity(
            model, tokenizer, target_sequence, rows[idx]["sequence"]
        ))
        production = float(production_ppls[idx])
        relative_error = abs(reference - production) / max(abs(reference), 1e-12)
        comparisons.append({
            "candidate_id": rows[idx]["candidate_id"],
            "upstream_ppl": reference,
            "production_ppl": production,
            "relative_error": relative_error,
            "passed": relative_error <= 5e-3,
        })
    if not all(item["passed"] for item in comparisons):
        raise RuntimeError("Microbatched PPL did not match the upstream implementation")
    return comparisons


def run(config_path: Path, model_override: str | None, skip_ppl: bool) -> Path:
    config, project_root = load_config(config_path)
    input_path = resolve_path(project_root, config["input_csv"])
    target_path = resolve_path(project_root, config["target_fasta"])
    model_path = resolve_path(project_root, model_override or config["model_path"])
    output_dir = resolve_path(project_root, config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    input_fields, input_rows = read_rows(input_path)
    _, target_sequence = read_single_fasta(target_path)

    seen: set[str] = set()
    annotated_rows = []
    passed_rows = []
    failure_counts: dict[str, int] = defaultdict(int)
    for row in input_rows:
        sequence = row.get("sequence", "").strip().upper()
        row = dict(row)
        row["sequence"] = sequence
        reasons, props = filter_reason(sequence, row, config["hard_filters"], seen)
        for field in PROPERTY_FIELDS:
            value = props.get(field, "")
            row[field] = f"{value:.6f}" if isinstance(value, float) else value
        row["filter_pass"] = int(not reasons)
        row["filter_reasons"] = ";".join(reasons)
        row["pseudo_perplexity"] = ""
        annotated_rows.append(row)
        if reasons:
            for reason in reasons:
                failure_counts[reason] += 1
        else:
            passed_rows.append(row)

    output_fields = input_fields + [field for field in PROPERTY_FIELDS if field not in input_fields]
    output_fields += ["filter_pass", "filter_reasons", "pseudo_perplexity"]
    atomic_write_csv(output_dir / "sequence_qc_all.csv", annotated_rows, output_fields)

    selected_rows = stratified_ppl_selection(
        passed_rows,
        config.get("ppl_max_candidates_per_length"),
    )
    ppls = []
    comparisons = []
    device = dtype = "not_loaded"
    if selected_rows and not skip_ppl and bool(config.get("compute_ppl", True)):
        torch, tokenizer, model, device, dtype = load_local_model(
            model_path,
            requested_device=config.get("device", "auto"),
            requested_dtype=config.get("dtype", "auto"),
        )
        ppls = compute_pseudo_perplexities(
            torch=torch,
            tokenizer=tokenizer,
            model=model,
            device=device,
            target_sequence=target_sequence,
            candidate_rows=selected_rows,
            row_batch_size=int(config.get("ppl_row_batch_size", 4)),
        )
        for row, ppl in zip(selected_rows, ppls):
            row["pseudo_perplexity"] = f"{ppl:.8f}"
        comparisons = validate_against_upstream(
            project_root=project_root,
            torch=torch,
            tokenizer=tokenizer,
            model=model,
            target_sequence=target_sequence,
            rows=selected_rows,
            production_ppls=ppls,
            count=int(config.get("validate_upstream_count", 0)),
        )

    scored_rows = sorted(
        selected_rows,
        key=lambda row: (
            int(row["length"]),
            float(row["pseudo_perplexity"]) if row["pseudo_perplexity"] else float("inf"),
            -float(row.get("mean_log_probability", "-inf")),
        ),
    )
    atomic_write_csv(output_dir / "ppl_selected.csv", scored_rows, output_fields)
    atomic_write_text(
        output_dir / "ppl_selected.fasta",
        "".join(
            f">{row['candidate_id']} length={row['length']} ppl={row['pseudo_perplexity'] or 'NA'}\n{row['sequence']}\n"
            for row in scored_rows
        ),
    )

    manifest = {
        "schema_version": 1,
        "stage": "03_sequence_filter",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": str(config_path.resolve()),
        "config_fingerprint": canonical_json_hash(config),
        "inputs": {
            "candidate_csv": str(input_path),
            "candidate_csv_sha256": sha256_file(input_path),
            "target_fasta": str(target_path),
            "target_length": len(target_sequence),
        },
        "counts": {
            "input": len(input_rows),
            "hard_filter_pass": len(passed_rows),
            "ppl_selected": len(selected_rows),
            "ppl_computed": len(ppls),
        },
        "hard_filter_failure_counts": dict(sorted(failure_counts.items())),
        "ppl": {
            "skipped": bool(skip_ppl),
            "row_batch_size": int(config.get("ppl_row_batch_size", 4)),
            "device": str(device),
            "dtype": str(dtype),
            "upstream_validation": comparisons,
        },
        "outputs": {
            "all_qc": str(output_dir / "sequence_qc_all.csv"),
            "ppl_selected_csv": str(output_dir / "ppl_selected.csv"),
            "ppl_selected_fasta": str(output_dir / "ppl_selected.fasta"),
        },
    }
    atomic_write_json(output_dir / "filter_manifest.json", manifest)
    return output_dir


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--model-path", help="Override local model directory from config")
    parser.add_argument("--skip-ppl", action="store_true", help="Run sequence QC without loading torch/model")
    args = parser.parse_args()
    try:
        output_dir = run(args.config, args.model_path, args.skip_ppl)
    except Exception as exc:
        print(f"03_sequence_filter FAILED: {exc}", file=sys.stderr)
        return 1
    print(f"03_sequence_filter PASS: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
