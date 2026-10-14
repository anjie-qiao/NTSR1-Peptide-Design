#!/usr/bin/env python3
"""Predict one peptide–NTSR1 complex with the project's Boltz settings."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import shutil
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
jobs = importlib.import_module("05_prepare_structure_jobs")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sequence", help="Peptide amino-acid sequence")
    parser.add_argument(
        "--prepare-only", action="store_true", help="Write the Boltz YAML without running prediction"
    )
    args = parser.parse_args()

    peptide = args.sequence.strip().upper()
    invalid = sorted(set(peptide) - set(jobs.CANONICAL_AAS))
    if not peptide or invalid:
        parser.error(f"Invalid peptide sequence; unsupported residues: {invalid}")

    config = json.loads((PROJECT_ROOT / "config/05_prepare_structure_jobs.json").read_text(encoding="utf-8"))
    _, receptor = jobs.read_single_fasta(PROJECT_ROOT / config["target_fasta"])
    msa = config["receptor_msa"]
    if not (PROJECT_ROOT / msa).is_file():
        parser.error(f"NTSR1 MSA not found: {PROJECT_ROOT / msa}")

    job_name = f"NTSR1_{hashlib.sha256(peptide.encode()).hexdigest()[:12]}"
    job_dir = PROJECT_ROOT / "work/single_peptide_boltz"
    input_path = job_dir / f"{job_name}.yaml"
    output_dir = job_dir / "predictions"
    jobs.atomic_write_text(
        input_path,
        jobs.boltz_yaml(
            receptor, peptide, config["receptor_chain_id"], config["peptide_chain_id"],
            msa, config["peptide_msa"],
        ),
    )
    print(f"Boltz input: {input_path}", flush=True)
    if args.prepare_only:
        return 0
    if shutil.which("boltz") is None:
        parser.error("boltz is unavailable; activate the project's Boltz environment first")
    cache_dir = PROJECT_ROOT / "models/boltz_cache"
    required_cache_files = ("mols.tar", "boltz2_conf.ckpt", "boltz2_aff.ckpt")
    missing = [
        name for name in required_cache_files
        if not (cache_dir / name).is_file() or (cache_dir / name).stat().st_size == 0
    ]
    if missing:
        parser.error(
            f"Boltz cache is incomplete at {cache_dir}: {', '.join(missing)}. "
            "Download these files on a machine with Hugging Face access, copy them "
            "to this cache directory, and see docs/INPUTS.md section 5."
        )

    settings = config["boltz_prediction"]
    command = [
        "boltz", "predict", str(input_path.relative_to(PROJECT_ROOT)),
        "--out_dir", str(output_dir.relative_to(PROJECT_ROOT)),
        "--cache", str(cache_dir),
        "--accelerator", "gpu",
        "--devices", str(settings.get("devices", 1)),
        "--recycling_steps", str(settings.get("recycling_steps", 3)),
        "--sampling_steps", str(settings.get("sampling_steps", 200)),
        "--diffusion_samples", str(settings.get("diffusion_samples", 1)),
        "--output_format", str(settings.get("output_format", "pdb")),
    ]
    for key, flag in (
        ("write_full_pae", "--write_full_pae"),
        ("use_potentials", "--use_potentials"),
        ("no_kernels", "--no_kernels"),
    ):
        if settings.get(key, key != "no_kernels"):
            command.append(flag)
    if settings.get("max_msa_seqs") is not None:
        command.extend(["--max_msa_seqs", str(settings["max_msa_seqs"])])

    subprocess.run(command, cwd=PROJECT_ROOT, check=True)
    print(f"Boltz results: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        print(f"Boltz prediction failed with exit code {exc.returncode}", file=sys.stderr)
        raise SystemExit(exc.returncode) from exc
