"""Export final candidate tables, structures and original Boltz artifacts."""
import csv
import hashlib
import json
import re
import shutil
import zipfile
from datetime import datetime, timezone
from pathlib import Path

try:
    from .boltz_evidence import retain_boltz_evidence
except ImportError:
    from boltz_evidence import retain_boltz_evidence

ROOT = Path(__file__).resolve().parents[1]
AA3 = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCQEGHILKMFPSTWYV"))


def read_rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pdb_sequences(path):
    chains = {}
    residues = set()
    with Path(path).open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM  ") or line[16:17] not in (" ", "A"):
                continue
            key = (line[21:22], line[22:27])
            if key in residues:
                continue
            residues.add(key)
            chains.setdefault(key[0], []).append(AA3.get(line[17:20].strip(), "X"))
    return {chain: "".join(seq) for chain, seq in chains.items()}


def expected_receptor_sequence():
    with zipfile.ZipFile(ROOT / "models/mouse_NTSR1_6Z8N_6ZIN_pack.zip") as archive:
        fasta = archive.read("source_info/O88319_mouse_NTSR1.fasta").decode("utf-8")
    return "".join(line.strip() for line in fasta.splitlines() if not line.startswith(">"))


def export_results(final_csv, risk_csv, output, structure_roots=(),
                   mode="historical_export", runtime_file=None, track=None):
    output = Path(output).resolve()
    output.relative_to(ROOT)  # Tables use paths relative to the submitted project package.
    candidates, risks = read_rows(final_csv), read_rows(risk_csv)
    if not candidates:
        raise ValueError("Final candidate table is empty")
    ids = [r["candidate_id"] for r in candidates]
    if len(ids) != len(set(ids)) or any(not re.fullmatch(r"[A-Za-z0-9_-]+", cid) for cid in ids):
        raise ValueError("Invalid or duplicate candidate IDs")
    if [r["candidate_id"] for r in risks] != ids:
        raise ValueError("Stage 11 IDs/order must exactly match Stage 10")
    risk_index = {r["candidate_id"]: r for r in risks}
    if any(r["sequence"] != risk_index[r["candidate_id"]]["sequence"] for r in candidates):
        raise ValueError("Stage 10/11 sequences disagree")
    settings = json.loads((ROOT / "config/submission.json").read_text(encoding="utf-8"))
    risk_manifest_path = Path(risk_csv).parent / "oral_cleavage_manifest.json"
    risk_manifest = json.loads(risk_manifest_path.read_text(encoding="utf-8")) if risk_manifest_path.is_file() else {}
    receptor_sequence = expected_receptor_sequence()
    source_index = {}
    for directory in structure_roots:
        directory = Path(directory).resolve()
        if not directory.is_dir():
            raise FileNotFoundError(directory)
        for path in directory.rglob("*_model_0.pdb"):
            cid = path.name.removesuffix("_model_0.pdb")
            if cid in ids:
                source_index.setdefault(cid, []).append(path)
    planned, issues = {}, []
    for row in candidates:
        cid = row["candidate_id"]
        source = Path(row.get("structure_pdb", ""))
        if not source.is_absolute():
            source = ROOT / source
        paths = source_index.get(cid, []).copy()
        if source.is_file():
            paths.append(source)
        paths = list(dict.fromkeys(p.resolve() for p in paths))
        status = "missing"
        if paths:
            if len({sha256(p) for p in paths}) != 1:
                status = "conflicting_structures"
            else:
                chains = pdb_sequences(paths[0])
                if chains.get("B") != row["sequence"]:
                    status = "peptide_sequence_mismatch"
                elif chains.get("A") != receptor_sequence:
                    status = "receptor_chain_invalid"
                else:
                    planned[cid] = paths[0]
                    status = "verified"
        if status != "verified":
            issues.append({"candidate_id": cid, "sequence": row["sequence"], "status": status})
    output.mkdir(parents=True, exist_ok=True)
    write_rows(output / "structure_issues.csv", issues, ["candidate_id", "sequence", "status"])
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(), "mode": mode,
        "project_version": settings["project_version"], "candidate_count": len(ids),
        "structure_count": len(planned), "complete": not issues,
        "submission_ready": not issues and len(ids) == 100,
        "final_selection_stage": 10, "annotation_stage": 11, "annotation_changes_ranking": False,
        "path_base": "project/", "source_final_sha256": sha256(final_csv),
        "source_risk_sha256": sha256(risk_csv), "structures": {},
        "historical_runtime_versions": settings["historical_runtime_versions"],
        "annotation_method": risk_manifest.get("method", "unknown"),
    }
    if runtime_file and Path(runtime_file).is_file():
        manifest["execution_runtime"] = json.loads(Path(runtime_file).read_text(encoding="utf-8"))
    (output / "export_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    if issues:
        stale = output / "results.csv"
        if stale.exists():
            stale.replace(output / "results.csv.previous")
        raise ValueError(f"Incomplete structures: {len(issues)}/{len(ids)}. See {output / 'structure_issues.csv'}; no results.csv was generated.")
    rows = []
    evidence = {"path_base": "project/", "candidates": {}, "runs": {},
                "note": "Original Boltz files are retained unchanged. PDBs are stored once in results/structures; source directories may be archived."}
    for row in candidates:
        cid = row["candidate_id"]
        copy = dict(row)
        copy["structure_pdb"] = ""
        copy["structure_sha256"] = ""
        copy["structure_status"] = "missing_or_invalid"
        copy["boltz_evidence_dir"] = ""
        copy["boltz_evidence_status"] = "not_available"
        if cid in planned:
            target = output / "structures" / f"{cid}_model_0.pdb"
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.resolve() != planned[cid]:
                shutil.copy2(planned[cid], target)
            copy["structure_pdb"] = target.relative_to(ROOT).as_posix()
            copy["structure_sha256"] = sha256(target)
            copy["structure_status"] = "verified"
            manifest["structures"][cid] = {"path": copy["structure_pdb"], "sha256": copy["structure_sha256"]}
            entry, run = retain_boltz_evidence(cid, planned[cid], output, sha256)
            if entry:
                evidence["candidates"][cid] = entry
                evidence["runs"][entry["run_key"]] = run
                copy["boltz_evidence_dir"] = entry["directory"]
                copy["boltz_evidence_status"] = entry["status"]
        copy.update({k: v for k, v in risk_index[cid].items() if k not in ("candidate_id", "sequence")})
        copy.update({"track": track or settings["track"], "target": settings["target"],
                     "model_versions": json.dumps(settings["models"], ensure_ascii=False),
                     "run_version": settings["project_version"], "run_mode": mode,
                     "annotation_method": manifest["annotation_method"],
                     "runtime_versions": json.dumps(manifest.get("execution_runtime", {}).get("model_packages", settings["historical_runtime_versions"]), ensure_ascii=False),
                     "notes": "Computational candidate; Stage 11 provides relative oral proteolysis risk annotation."})
        rows.append(copy)
    table_name = "results.csv"
    write_rows(output / table_name, rows)
    if evidence["candidates"]:
        evidence_path = output / "boltz_evidence" / "manifest.json"
        evidence_path.parent.mkdir(parents=True, exist_ok=True)
        evidence_path.write_text(json.dumps(evidence, ensure_ascii=False, indent=2), encoding="utf-8")
        manifest["boltz_evidence_manifest"] = evidence_path.relative_to(ROOT).as_posix()
    manifest["boltz_evidence_count"] = len(evidence["candidates"])
    manifest["table"] = (output / table_name).relative_to(ROOT).as_posix()
    (output / "export_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return output / table_name
