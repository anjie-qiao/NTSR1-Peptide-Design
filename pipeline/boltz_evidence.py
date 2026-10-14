"""Retain original Boltz prediction artifacts alongside the final candidate export."""
import json
import shutil
from pathlib import Path


def retain_boltz_evidence(candidate_id, source_pdb, output, sha256):
    """Copy bytes unchanged; a standalone example PDB need not have Boltz artifacts."""
    output = Path(output)
    root = Path(__file__).resolve().parents[1]
    source_pdb = Path(source_pdb)
    prediction_dir = source_pdb.parent
    destination = output / "boltz_evidence" / candidate_id
    canonical_pdb = output / "structures" / f"{candidate_id}_model_0.pdb"
    previous_path = output / "boltz_evidence" / "manifest.json"
    # Re-export remains possible after bulk intermediate predictions are archived.
    if prediction_dir.name != candidate_id or prediction_dir.parent.name != "predictions":
        if previous_path.is_file():
            previous = json.loads(previous_path.read_text(encoding="utf-8"))
            entry = previous.get("candidates", {}).get(candidate_id)
            if entry and all((root / item["path"]).is_file()
                             and sha256(root / item["path"]) == item["sha256"]
                             for item in entry["files"]):
                return entry, previous.get("runs", {}).get(entry["run_key"], {})
        return None, None

    result_root = prediction_dir.parent.parent
    stage_root = result_root.parent.parent
    input_name = result_root.name.removeprefix("boltz_results_")
    input_yaml = stage_root / input_name / f"{candidate_id}.yaml"
    required = [prediction_dir / f"confidence_{candidate_id}_model_0.json",
                prediction_dir / f"pae_{candidate_id}_model_0.npz",
                prediction_dir / f"plddt_{candidate_id}_model_0.npz", input_yaml]
    missing = [path.name for path in required if not path.is_file()]
    destination.mkdir(parents=True, exist_ok=True)
    files = []
    for source in sorted(prediction_dir.iterdir()) + [input_yaml]:
        if not source.is_file() or source.suffix.lower() in (".pdb", ".cif"):
            continue
        target = destination / source.name
        if source.resolve() != target.resolve():
            shutil.copy2(source, target)
        files.append({"path": target.relative_to(root).as_posix(),
                      "sha256": sha256(target), "bytes": target.stat().st_size})
    files.append({"path": canonical_pdb.relative_to(root).as_posix(),
                  "sha256": sha256(canonical_pdb), "bytes": canonical_pdb.stat().st_size})

    run_key = stage_root.name
    run_destination = output / "boltz_evidence" / "runs" / run_key
    metadata_sources = list(stage_root.glob("run_boltz*.sh")) + list(stage_root.glob("boltz*.log"))
    metadata_sources += [stage_root / "structure_job_manifest.json",
                         stage_root / "receptor_context.json"]
    metadata_sources += list(stage_root.glob("structure_tasks*.csv"))
    config_root = stage_root.parent / "config"
    if not config_root.is_dir():
        config_root = root / "config"
    metadata_sources += list(config_root.glob(f"{run_key[:2]}_*.json"))
    metadata_sources += [stage_root.parent / "logs" / f"{run_key[:2]}_boltz.log",
                         stage_root.parent / "runtime_environment.json"]
    processed_manifest = result_root / "processed" / "manifest.json"
    if processed_manifest.is_file():
        metadata_sources.append(processed_manifest)
    metadata_sources += sorted((result_root / "lightning_logs").rglob("hparams.yaml"))
    run_files = []
    for source in metadata_sources:
        if not source.is_file():
            continue
        if source == processed_manifest:
            name = "processed_manifest.json"
        elif source.name == "hparams.yaml":
            name = f"{source.parent.name}_hparams.yaml"
        else:
            name = source.name
        run_destination.mkdir(parents=True, exist_ok=True)
        target = run_destination / name
        if not target.is_file() or sha256(target) != sha256(source):
            shutil.copy2(source, target)
        run_files.append({"path": target.relative_to(root).as_posix(), "sha256": sha256(target)})

    entry = {"directory": destination.relative_to(root).as_posix(), "run_key": run_key,
             "source_prediction_directory": prediction_dir.relative_to(root).as_posix(),
             "status": "complete" if not missing else "partial", "missing": missing, "files": files}
    return entry, {"files": run_files}
