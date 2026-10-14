#!/usr/bin/env python3
"""Prepare traceable Boltz/ColabFold complex-prediction jobs from the core library."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shlex
import sys
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


TASK_FIELDS = [
    "job_id", "candidate_id", "sequence", "length", "sampling_profile",
    "pseudo_perplexity", "ppl_percentile_within_length", "cluster_id",
    "boltz_yaml", "receptor_chain_id", "peptide_chain_id",
]


def resolve_path(project_root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def build_boltz_flags(boltz_config: dict) -> str:
    """Build the Boltz CLI flags for full-library prediction."""
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


def load_config(config_path: Path) -> tuple[dict, Path]:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    required = {
        "core_candidates_csv", "core_manifest", "target_fasta", "target_manifest",
        "screening_masks", "output_dir",
        "receptor_chain_id", "peptide_chain_id", "receptor_msa", "peptide_msa",
        "boltz_prediction",
    }
    missing = required - set(config)
    if missing:
        raise ValueError(f"Missing config keys: {sorted(missing)}")
    if config["receptor_chain_id"] == config["peptide_chain_id"]:
        raise ValueError("Receptor and peptide chain IDs must differ")
    if not isinstance(config["boltz_prediction"], dict):
        raise ValueError("boltz_prediction must be an object")
    return config, config_path.resolve().parent.parent


def read_core_candidates(path: Path) -> tuple[list[str], list[dict]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {path}")
        required = {
            "candidate_id", "sequence", "length", "sampling_profile",
            "pseudo_perplexity", "ppl_percentile_within_length", "cluster_id",
            "selected_for_core",
        }
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Core CSV is missing fields: {sorted(missing)}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise ValueError("Core candidate CSV is empty")
    ids = set()
    sequences = set()
    for row in rows:
        candidate_id = row["candidate_id"].strip()
        sequence = row["sequence"].strip().upper()
        try:
            length = int(row["length"])
            ppl = float(row["pseudo_perplexity"])
            percentile = float(row["ppl_percentile_within_length"])
        except ValueError as exc:
            raise ValueError(f"Invalid numeric field for {candidate_id}") from exc
        if row["selected_for_core"] not in {"1", "true", "True"}:
            raise ValueError(f"Non-core row found in core library: {candidate_id}")
        if not candidate_id or candidate_id in ids:
            raise ValueError(f"Missing or duplicate candidate ID: {candidate_id!r}")
        if not sequence or sequence in sequences or set(sequence) - set(CANONICAL_AAS):
            raise ValueError(f"Invalid or duplicate peptide sequence for {candidate_id}")
        if length != len(sequence) or not math.isfinite(ppl) or not 0 <= percentile <= 1:
            raise ValueError(f"Invalid length/PPL for {candidate_id}")
        row["candidate_id"] = candidate_id
        row["sequence"] = sequence
        row["_length"] = length
        row["_ppl"] = ppl
        row["_percentile"] = percentile
        ids.add(candidate_id)
        sequences.add(sequence)
    return list(reader.fieldnames), rows


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


def run(config_path: Path) -> Path:
    config, project_root = load_config(config_path)
    core_path = resolve_path(project_root, config["core_candidates_csv"])
    core_manifest_path = resolve_path(project_root, config["core_manifest"])
    target_path = resolve_path(project_root, config["target_fasta"])
    target_manifest_path = resolve_path(project_root, config["target_manifest"])
    masks_path = resolve_path(project_root, config["screening_masks"])
    output_dir = resolve_path(project_root, config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    _, rows = read_core_candidates(core_path)
    core_manifest = json.loads(core_manifest_path.read_text(encoding="utf-8"))
    expected_core = core_manifest.get("counts", {}).get("core_cluster_representatives")
    if expected_core is not None and int(expected_core) != len(rows):
        raise ValueError(f"Core manifest expects {expected_core} rows; observed {len(rows)}")
    target_header, receptor_sequence = read_single_fasta(target_path)
    target_manifest = json.loads(target_manifest_path.read_text(encoding="utf-8"))
    masks = json.loads(masks_path.read_text(encoding="utf-8"))
    if len(receptor_sequence) != int(target_manifest["target"]["sequence_length"]):
        raise ValueError("Target FASTA length does not match target manifest")

    core_yaml_dir = output_dir / "boltz_inputs_core"
    core_yaml_dir.mkdir(parents=True, exist_ok=True)
    expected_core_yamls = {f"{row['candidate_id']}.yaml" for row in rows}
    validate_yaml_directory(core_yaml_dir, expected_core_yamls)

    receptor_chain = str(config["receptor_chain_id"])
    peptide_chain = str(config["peptide_chain_id"])
    receptor_msa = str(config["receptor_msa"])
    peptide_msa = str(config["peptide_msa"])
    task_rows = []
    for row in sorted(rows, key=lambda item: (item["_length"], item["_ppl"], item["candidate_id"])):
        job_id = f"NTSR1__{row['candidate_id']}"
        yaml_text = boltz_yaml(
            receptor_sequence, row["sequence"], receptor_chain, peptide_chain,
            receptor_msa, peptide_msa,
        )
        core_yaml = core_yaml_dir / f"{row['candidate_id']}.yaml"
        atomic_write_text(core_yaml, yaml_text)
        task_rows.append({
            "job_id": job_id,
            "candidate_id": row["candidate_id"],
            "sequence": row["sequence"],
            "length": row["_length"],
            "sampling_profile": row["sampling_profile"],
            "pseudo_perplexity": row["pseudo_perplexity"],
            "ppl_percentile_within_length": row["ppl_percentile_within_length"],
            "cluster_id": row["cluster_id"],
            "boltz_yaml": relative_posix(core_yaml, project_root),
            "receptor_chain_id": receptor_chain,
            "peptide_chain_id": peptide_chain,
    })
    atomic_write_csv(output_dir / "structure_tasks_all.csv", task_rows, TASK_FIELDS)

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
        "production_structure_screening_ready": target_manifest["readiness"][
            "production_structure_screening_ready"
        ],
        "structure_screening_blockers": target_manifest["readiness"]["structure_screening_blockers"],
    }
    atomic_write_json(output_dir / "receptor_context.json", receptor_context)

    boltz_config = config["boltz_prediction"]
    boltz_flags = build_boltz_flags(boltz_config)

    preflight_script = (
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        "command -v boltz >/dev/null || { echo 'ERROR: boltz is not installed in this environment' >&2; exit 1; }\n"
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
    atomic_write_text(output_dir / "run_boltz_preflight.sh", preflight_script)

    run_script = (
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        "bash work/05_prepare_structure_jobs/run_boltz_preflight.sh\n\n"
        "boltz predict work/05_prepare_structure_jobs/boltz_inputs_core \\\n"
        "  --out_dir work/05_prepare_structure_jobs/boltz_predictions_core \\\n"
        f"{boltz_flags}\n"
    )
    run_script = run_script.replace("work/05_prepare_structure_jobs/", shlex.quote(output_dir.as_posix()) + "/")
    atomic_write_text(output_dir / "run_boltz_core.sh", run_script)

    manifest = {
        "schema_version": 1,
        "stage": "05_prepare_structure_jobs",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config_path": str(config_path.resolve()),
        "config_fingerprint": canonical_json_hash(config),
        "inputs": {
            "core_candidates_csv": str(core_path),
            "core_candidates_sha256": sha256_file(core_path),
            "core_manifest": str(core_manifest_path),
            "target_fasta": str(target_path),
            "target_manifest": str(target_manifest_path),
            "screening_masks": str(masks_path),
        },
        "method": {
            "primary_predictor_input": "Boltz YAML v1",
            "receptor_input": "full-length O88319 sequence",
            "receptor_msa": receptor_msa,
            "peptide_msa": peptide_msa,
            "templates": "none",
            "existing_incomplete_pdb_used_as_prediction_template": False,
            "orthosteric_pocket_used_as_positive_constraint": False,
        },
        "counts": {
            "core_tasks": len(task_rows),
        },
        "readiness": {
            "full_library_prediction_jobs_ready": True,
            "production_ensemble_screening_ready": False,
            "blockers": receptor_context["structure_screening_blockers"],
        },
        "outputs": {
            "all_tasks": str(output_dir / "structure_tasks_all.csv"),
            "boltz_core_inputs": str(core_yaml_dir),
            "receptor_context": str(output_dir / "receptor_context.json"),
            "boltz_core_script": str(output_dir / "run_boltz_core.sh"),
            "boltz_preflight_script": str(output_dir / "run_boltz_preflight.sh"),
        },
    }
    atomic_write_json(output_dir / "structure_job_manifest.json", manifest)
    report = (
        "# Stage 05 full-library structure-job readiness\n\n"
        f"- Core library tasks: {len(task_rows)}\n"
        "- Primary input: Boltz YAML v1, receptor chain A and peptide chain B.\n"
        f"- Receptor MSA: {receptor_msa}.\n"
        f"- Peptide MSA: {peptide_msa}.\n"
        "- Existing incomplete 6Z8N/6ZIN-derived structures are not forced as templates.\n"
        "- Orthosteric pocket residues and missing-segment masks are exported for post-prediction checks only.\n\n"
        "Outputs remain computational predictions and require downstream interface, extracellular-accessibility, "
        "membrane-orientation and missing-region artefact checks.\n"
    )
    atomic_write_text(output_dir / "structure_prediction_readiness.md", report)
    return output_dir


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    try:
        output_dir = run(args.config)
    except Exception as exc:
        print(f"05_prepare_structure_jobs FAILED: {exc}", file=sys.stderr)
        return 1
    print(f"05_prepare_structure_jobs PASS: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
