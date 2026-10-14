#!/usr/bin/env python3
"""Vectorized, restartable PepMLM candidate generation.

The upstream repository is retained unchanged.  This implementation reuses its
all-mask/top-k sampling definition while avoiding one identical model forward
pass per candidate: for a fixed target and peptide length, logits are computed
once and many sequences are sampled from those logits.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

from pepmlm_pipeline_common import (
    CANONICAL_AAS,
    allocate_integer,
    atomic_write_csv,
    atomic_write_json,
    atomic_write_text,
    canonical_json_hash,
    read_single_fasta,
    sha256_file,
    validate_sampling_profiles,
)
from pepmlm_model_runtime import canonical_token_ids, load_local_model


CSV_FIELDS = [
    "candidate_id",
    "sequence",
    "length",
    "sampling_profile",
    "top_k",
    "temperature",
    "sampling_seed",
    "raw_log_probability",
    "mean_log_probability",
    "generation_rank_within_length",
]


def load_config(config_path: Path) -> tuple[dict, Path]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    required = {
        "target_fasta", "model_path", "output_dir", "lengths",
        "sampling_profiles", "base_seed",
    }
    missing = required - set(config)
    if missing:
        raise ValueError(f"Missing config keys: {sorted(missing)}")
    lengths = [int(length) for length in config["lengths"]]
    if len(lengths) != len(set(lengths)) or any(length <= 0 for length in lengths):
        raise ValueError("lengths must be unique positive integers")
    count_modes = [key for key in ("candidates_per_length", "total_candidates") if key in config]
    if len(count_modes) != 1:
        raise ValueError("Specify exactly one of candidates_per_length or total_candidates")
    if int(config[count_modes[0]]) <= 0:
        raise ValueError(f"{count_modes[0]} must be positive")
    if "length_weights" in config:
        weights = config["length_weights"]
        if len(weights) != len(lengths) or any(float(weight) < 0 for weight in weights):
            raise ValueError("length_weights must match lengths and contain no negative values")
        if sum(float(weight) for weight in weights) <= 0:
            raise ValueError("length_weights must contain at least one positive value")
    validate_sampling_profiles(config["sampling_profiles"])
    project_root = config_path.resolve().parent.parent
    return config, project_root


def resolve_path(project_root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def generation_plan(config: dict, project_root: Path, model_override: str | None = None) -> dict:
    target_path = resolve_path(project_root, config["target_fasta"])
    model_path = resolve_path(project_root, model_override or config["model_path"])
    output_dir = resolve_path(project_root, config["output_dir"])
    header, sequence = read_single_fasta(target_path)
    lengths = [int(length) for length in config["lengths"]]
    if "total_candidates" in config:
        counts_by_length = allocate_integer(
            int(config["total_candidates"]),
            [float(weight) for weight in config.get("length_weights", [1.0] * len(lengths))],
        )
    else:
        counts_by_length = [int(config["candidates_per_length"])] * len(lengths)
    profile_weights = [float(profile["weight"]) for profile in config["sampling_profiles"]]
    length_plan = []
    for length, length_count in zip(lengths, counts_by_length):
        profile_counts = allocate_integer(length_count, profile_weights)
        length_plan.append({
            "length": length,
            "candidate_count": length_count,
            "profile_counts": {
                profile["name"]: count
                for profile, count in zip(config["sampling_profiles"], profile_counts)
            },
        })
    total = sum(counts_by_length)
    return {
        "target_fasta": str(target_path),
        "target_header": header,
        "target_length": len(sequence),
        "target_sequence_sha256": canonical_json_hash(sequence),
        "model_path": str(model_path),
        "output_dir": str(output_dir),
        "lengths": lengths,
        "candidate_count_mode": "total_candidates" if "total_candidates" in config else "candidates_per_length",
        "candidates_per_length": int(config["candidates_per_length"]) if "candidates_per_length" in config else None,
        "candidates_by_length": {
            str(length): count for length, count in zip(lengths, counts_by_length)
        },
        "length_plan": length_plan,
        "total_candidates": total,
        "profiles": config["sampling_profiles"],
        "base_seed": int(config["base_seed"]),
    }


def compute_mask_logits(torch, tokenizer, model, device: str, target_sequence: str, peptide_length: int):
    masked_peptide = tokenizer.mask_token * peptide_length
    inputs = tokenizer(target_sequence + masked_peptide, return_tensors="pt").to(device)
    mask_positions = (inputs["input_ids"] == tokenizer.mask_token_id).nonzero(as_tuple=True)[1]
    if int(mask_positions.numel()) != peptide_length:
        raise ValueError(
            f"Tokenizer produced {int(mask_positions.numel())} mask tokens; expected {peptide_length}"
        )
    with torch.inference_mode():
        logits = model(**inputs).logits[0, mask_positions].float().cpu()
    return logits


def sample_profile(
    torch,
    canonical_logits,
    canonical_aas: str,
    profile: dict,
    desired_count: int,
    sampling_seed: int,
    seen: set[str],
    max_rounds: int,
) -> list[dict]:
    if desired_count == 0:
        return []
    top_k = int(profile["top_k"])
    temperature = float(profile["temperature"])
    top_logits, top_indices = torch.topk(canonical_logits / temperature, k=top_k, dim=-1)
    probabilities = torch.softmax(top_logits, dim=-1)
    log_probabilities = torch.log(probabilities.clamp_min(1e-30))
    generator = torch.Generator(device="cpu")
    generator.manual_seed(sampling_seed)
    accepted: list[dict] = []

    for round_index in range(max_rounds):
        remaining = desired_count - len(accepted)
        if remaining <= 0:
            break
        draw_count = max(remaining * 2, 32)
        sampled_local = torch.multinomial(
            probabilities,
            num_samples=draw_count,
            replacement=True,
            generator=generator,
        )
        sampled_aa_indices = top_indices.gather(1, sampled_local)
        selected_log_probs = log_probabilities.gather(1, sampled_local).sum(dim=0)
        sampled_aa_indices = sampled_aa_indices.transpose(0, 1).tolist()
        selected_log_probs = selected_log_probs.tolist()

        for aa_indices, log_probability in zip(sampled_aa_indices, selected_log_probs):
            sequence = "".join(canonical_aas[index] for index in aa_indices)
            if sequence in seen:
                continue
            seen.add(sequence)
            accepted.append({
                "sequence": sequence,
                "sampling_profile": profile["name"],
                "top_k": top_k,
                "temperature": temperature,
                "sampling_seed": sampling_seed,
                "raw_log_probability": float(log_probability),
                "mean_log_probability": float(log_probability) / len(sequence),
            })
            if len(accepted) == desired_count:
                break
    if len(accepted) != desired_count:
        raise RuntimeError(
            f"Profile {profile['name']} produced only {len(accepted)}/{desired_count} unique sequences "
            f"after {max_rounds} rounds; increase top_k/temperature or max_sampling_rounds"
        )
    return accepted


def read_csv_rows(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def prepare_length_shard(
    torch,
    tokenizer,
    model,
    device: str,
    target_sequence: str,
    peptide_length: int,
    desired_count: int,
    config: dict,
    output_dir: Path,
    fingerprint: str,
    target_hash: str,
    force: bool,
) -> list[dict]:
    shard_path = output_dir / "shards" / f"length_{peptide_length:02d}.csv"
    meta_path = output_dir / "shards" / f"length_{peptide_length:02d}.meta.json"
    if shard_path.exists() and meta_path.exists() and not force:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if meta.get("config_fingerprint") != fingerprint or meta.get("target_sequence_sha256") != target_hash:
            raise RuntimeError(f"Existing shard metadata does not match current run: {shard_path}")
        rows = read_csv_rows(shard_path)
        if len(rows) != desired_count:
            raise RuntimeError(f"Existing shard has an unexpected row count: {shard_path}")
        print(f"RESUME length={peptide_length}: {len(rows)} candidates")
        return rows

    logits = compute_mask_logits(torch, tokenizer, model, device, target_sequence, peptide_length)
    token_ids = canonical_token_ids(tokenizer)
    canonical_logits = logits[:, token_ids]
    counts = allocate_integer(
        desired_count,
        [float(profile["weight"]) for profile in config["sampling_profiles"]],
    )
    seen: set[str] = set()
    rows: list[dict] = []
    for profile_index, (profile, count) in enumerate(zip(config["sampling_profiles"], counts)):
        seed = int(config["base_seed"]) + peptide_length * 10_000 + profile_index * 1_000_000
        rows.extend(sample_profile(
            torch=torch,
            canonical_logits=canonical_logits,
            canonical_aas=CANONICAL_AAS,
            profile=profile,
            desired_count=count,
            sampling_seed=seed,
            seen=seen,
            max_rounds=int(config.get("max_sampling_rounds", 25)),
        ))

    for index, row in enumerate(rows, start=1):
        row["candidate_id"] = f"PMLM-L{peptide_length:02d}-{index:06d}"
        row["length"] = peptide_length
        row["generation_rank_within_length"] = index
    rows.sort(key=lambda row: int(row["generation_rank_within_length"]))
    atomic_write_csv(shard_path, rows, CSV_FIELDS)
    atomic_write_json(meta_path, {
        "length": peptide_length,
        "candidate_count": len(rows),
        "config_fingerprint": fingerprint,
        "target_sequence_sha256": target_hash,
        "created_utc": datetime.now(timezone.utc).isoformat(),
    })
    print(f"GENERATED length={peptide_length}: {len(rows)} candidates")
    return rows


def run(config_path: Path, model_override: str | None, force: bool) -> Path:
    config, project_root = load_config(config_path)
    plan = generation_plan(config, project_root, model_override)
    target_path = Path(plan["target_fasta"])
    model_path = Path(plan["model_path"])
    output_dir = Path(plan["output_dir"])
    _, target_sequence = read_single_fasta(target_path)
    target_hash = plan["target_sequence_sha256"]
    fingerprint_payload = {
        "target_sequence_sha256": target_hash,
        "lengths": plan["lengths"],
        "candidates_by_length": plan["candidates_by_length"],
        "sampling_profiles": config["sampling_profiles"],
        "base_seed": plan["base_seed"],
        "canonical_aas": CANONICAL_AAS,
    }
    fingerprint = canonical_json_hash(fingerprint_payload)
    output_dir.mkdir(parents=True, exist_ok=True)

    missing_lengths = []
    for length in plan["lengths"]:
        shard = output_dir / "shards" / f"length_{length:02d}.csv"
        meta = output_dir / "shards" / f"length_{length:02d}.meta.json"
        if force or not (shard.exists() and meta.exists()):
            missing_lengths.append(length)

    torch = tokenizer = model = device = dtype = None
    if missing_lengths:
        torch, tokenizer, model, device, dtype = load_local_model(
            model_path,
            requested_device=config.get("device", "auto"),
            requested_dtype=config.get("dtype", "auto"),
        )

    all_rows: list[dict] = []
    for length in plan["lengths"]:
        all_rows.extend(prepare_length_shard(
            torch=torch,
            tokenizer=tokenizer,
            model=model,
            device=device,
            target_sequence=target_sequence,
            peptide_length=length,
            desired_count=int(plan["candidates_by_length"][str(length)]),
            config=config,
            output_dir=output_dir,
            fingerprint=fingerprint,
            target_hash=target_hash,
            force=force,
        ))

    if len(all_rows) != plan["total_candidates"]:
        raise RuntimeError(f"Expected {plan['total_candidates']} candidates; observed {len(all_rows)}")
    sequences = [row["sequence"] for row in all_rows]
    if len(sequences) != len(set(sequences)):
        raise RuntimeError("Duplicate sequences detected after shard collation")

    csv_path = output_dir / "raw_candidates.csv"
    fasta_path = output_dir / "raw_candidates.fasta"
    atomic_write_csv(csv_path, all_rows, CSV_FIELDS)
    atomic_write_text(
        fasta_path,
        "".join(f">{row['candidate_id']} length={row['length']} profile={row['sampling_profile']}\n{row['sequence']}\n" for row in all_rows),
    )

    upstream_path = Path(__file__).resolve().parents[1] / "pepmlm-main" / "scripts" / "generation.py"
    model_weight = model_path / "pytorch_model.bin"
    manifest = {
        "schema_version": 1,
        "stage": "02_pepmlm_generate",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": str(config_path.resolve()),
        "config_fingerprint": fingerprint,
        "target": {
            "fasta": str(target_path),
            "sequence_length": len(target_sequence),
            "sequence_sha256": target_hash,
        },
        "model": {
            "local_path": str(model_path),
            "weight_file_size_bytes": model_weight.stat().st_size if model_weight.exists() else None,
            "device": str(device) if device is not None else "resumed_without_model_load",
            "dtype": str(dtype) if dtype is not None else "resumed_without_model_load",
        },
        "upstream_reference": {
            "source_url": "https://github.com/programmablebio/pepmlm/blob/main/scripts/generation.py",
            "path": str(upstream_path),
            "sha256": sha256_file(upstream_path) if upstream_path.exists() else None,
            "reuse_note": "All-mask/top-k definition retained; repeated identical forwards replaced by vectorized sampling.",
        },
        "generation": plan,
        "outputs": {
            "csv": str(csv_path),
            "fasta": str(fasta_path),
            "candidate_count": len(all_rows),
            "unique_sequence_count": len(set(sequences)),
        },
    }
    atomic_write_json(output_dir / "generation_manifest.json", manifest)
    return output_dir


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--model-path", help="Override local model directory from config")
    parser.add_argument("--plan-only", action="store_true", help="Validate config and print the run plan without torch")
    parser.add_argument("--force", action="store_true", help="Regenerate known length shards for this stage")
    args = parser.parse_args()
    try:
        config, project_root = load_config(args.config)
        plan = generation_plan(config, project_root, args.model_path)
        if args.plan_only:
            print(json.dumps(plan, indent=2, ensure_ascii=False))
            return 0
        output_dir = run(args.config, args.model_path, args.force)
    except Exception as exc:
        print(f"02_pepmlm_generate FAILED: {exc}", file=sys.stderr)
        return 1
    print(f"02_pepmlm_generate PASS: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
