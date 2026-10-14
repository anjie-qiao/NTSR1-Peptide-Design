#!/usr/bin/env python3
"""Prepare the canonical NTSR1 sequence and structural reference artifacts.

This stage is intentionally dependency-free.  It validates the supplied model
archive, extracts only the files needed downstream, creates receptor-only PDBs,
and writes machine-readable masks and a preparation report.  It does not model
missing loops or claim that the structures are production-docking ready.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path


AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

MODEL_SPECS = {
    "mouse_NTSR1_6Z8N": {
        "archive_pdb": "models/mouse_NTSR1_6Z8N.pdb",
        "archive_cif": "models/mouse_NTSR1_6Z8N.cif",
        "ligand_pdb": "ligands/6Z8N_SR9_original_coords.pdb",
        "pocket_csv": "pockets/mouse_NTSR1_6Z8N_pocket5A.csv",
        "missing_csv": "source_info/mouse_NTSR1_6Z8N_missing_residues.csv",
        "state": "agonist_bound_active_reference",
        "ligand": "SRI-9829/SR9",
    },
    "mouse_NTSR1_6ZIN": {
        "archive_pdb": "models/mouse_NTSR1_6ZIN.pdb",
        "archive_cif": "models/mouse_NTSR1_6ZIN.cif",
        "ligand_pdb": "ligands/6ZIN_Q6Q_original_coords.pdb",
        "pocket_csv": "pockets/mouse_NTSR1_6ZIN_pocket5A.csv",
        "missing_csv": "source_info/mouse_NTSR1_6ZIN_missing_residues.csv",
        "state": "inverse_agonist_bound_inactive_reference",
        "ligand": "SR48692/Q6Q",
    },
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def read_fasta(text: str) -> tuple[str, str]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise ValueError("Canonical target FASTA is missing a valid header")
    sequence = "".join(lines[1:]).upper()
    if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
        raise ValueError("Canonical target FASTA contains invalid residues")
    return lines[0][1:], sequence


def parse_checksum_manifest(text: str) -> dict[str, str]:
    checksums: dict[str, str] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        digest, name = line.split(maxsplit=1)
        checksums[name.strip()] = digest.lower()
    return checksums


def verify_zip_checksums(zf: zipfile.ZipFile) -> dict:
    expected = parse_checksum_manifest(zf.read("SHA256SUMS").decode("utf-8"))
    failures = []
    for name, digest in expected.items():
        actual = sha256_bytes(zf.read(name))
        if actual != digest:
            failures.append({"file": name, "expected": digest, "actual": actual})
    return {
        "manifest_entries": len(expected),
        "all_passed": not failures,
        "failures": failures,
    }


def pdb_receptor_only(pdb_text: str) -> str:
    kept = []
    for line in pdb_text.splitlines():
        record = line[:6].strip()
        if record == "HETATM":
            continue
        if record == "CONECT":
            continue
        if record == "END":
            continue
        kept.append(line)
    kept.insert(2, "REMARK 901 RECEPTOR-ONLY COPY; LIGAND REMOVED BY 01_PREPARE_TARGET.")
    kept.append("END")
    return "\n".join(kept) + "\n"


def parse_pdb_residues(pdb_text: str, chain: str = "A") -> dict[int, str]:
    residues: dict[int, str] = {}
    for line in pdb_text.splitlines():
        if line[:6].strip() != "ATOM" or len(line) < 27:
            continue
        if line[21].strip() != chain:
            continue
        resname = line[17:20].strip()
        if resname not in AA3_TO_AA1:
            continue
        resnum = int(line[22:26])
        residues.setdefault(resnum, AA3_TO_AA1[resname])
    return residues


def parse_pocket_csv(text: str) -> list[dict]:
    return list(csv.DictReader(io.StringIO(text)))


def parse_missing_ranges(text: str) -> list[dict]:
    rows = list(csv.DictReader(io.StringIO(text)))
    return [
        {
            "start": int(row["missing_mouse_residue_start"]),
            "end": int(row["missing_mouse_residue_end"]),
            "count": int(row["count"]),
        }
        for row in rows
    ]


def find_n_glyco_sequons(sequence: str) -> list[dict]:
    sites = []
    for match in re.finditer(r"N[^P][ST]", sequence):
        sites.append({"residue": match.start() + 1, "motif": match.group(0)})
    return sites


def residue_label(row: dict) -> str:
    return f"{row['one_letter']}{int(row['mouse_residue_number'])}"


def load_config(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    required = {"archive", "output_dir", "target_accession", "expected_length"}
    missing = sorted(required - set(data))
    if missing:
        raise ValueError(f"Missing config keys: {', '.join(missing)}")
    return data


def render_report(manifest: dict) -> str:
    lines = [
        "# 01_prepare_target report",
        "",
        f"- Target: **{manifest['target']['name']}**",
        f"- Accession: **{manifest['target']['accession']}**",
        f"- Canonical length: **{manifest['target']['sequence_length']} aa**",
        f"- Archive SHA256: `{manifest['source_archive']['sha256']}`",
        f"- Internal checksum verification: **{'PASS' if manifest['source_archive']['internal_checksums']['all_passed'] else 'FAIL'}**",
        f"- PepMLM sequence generation ready: **{manifest['readiness']['pepmlm_generation_ready']}**",
        f"- Production structure screening ready: **{manifest['readiness']['production_structure_screening_ready']}**",
        "",
        "## Prepared models",
        "",
        "| Model | State | Modeled residues | Sequence mismatches | Missing ranges |",
        "|---|---|---:|---:|---|",
    ]
    for model in manifest["models"]:
        ranges = "; ".join(f"{x['start']}-{x['end']}" for x in model["missing_ranges"])
        lines.append(
            f"| {model['model_id']} | {model['state']} | {model['modeled_residue_count']} | "
            f"{model['canonical_sequence_mismatch_count']} | {ranges} |"
        )
    lines.extend([
        "",
        "## Orthosteric pocket reference",
        "",
        f"- Union ({len(manifest['orthosteric_reference']['union'])}): "
        + " ".join(manifest["orthosteric_reference"]["union"]),
        f"- Intersection ({len(manifest['orthosteric_reference']['intersection'])}): "
        + " ".join(manifest["orthosteric_reference"]["intersection"]),
        "",
        "The supplied small-molecule grids are retained only as references. They are not peptide docking boxes.",
        "",
        "## Putative N-glycosylation sequons",
        "",
        ", ".join(f"{x['motif']}@{x['residue']}" for x in manifest["target"]["putative_n_glyco_sequons"]),
        "",
        "Extracellular steric-mask sites used downstream: N4, N38, N42 and N211. N292 is an intracellular sequon and is not treated as an extracellular glycan site.",
        "",
        "## Remaining blockers for production structural screening",
        "",
    ])
    lines.extend(f"- {item}" for item in manifest["readiness"]["structure_screening_blockers"])
    lines.extend([
        "",
        "The receptor-only PDBs are validated reference structures, not completed or membrane-relaxed production docking receptors.",
        "",
    ])
    return "\n".join(lines)


def prepare(config_path: Path) -> Path:
    cfg = load_config(config_path)
    project_root = config_path.resolve().parent.parent
    archive_path = (project_root / cfg["archive"]).resolve()
    outdir = (project_root / cfg["output_dir"]).resolve()
    if not archive_path.is_file():
        raise FileNotFoundError(f"Archive not found: {archive_path}")

    archive_sha = hashlib.sha256(archive_path.read_bytes()).hexdigest()
    with zipfile.ZipFile(archive_path) as zf:
        checksum_status = verify_zip_checksums(zf)
        if not checksum_status["all_passed"]:
            raise RuntimeError("Archive internal checksum verification failed")

        fasta_header, sequence = read_fasta(
            zf.read("source_info/O88319_mouse_NTSR1.fasta").decode("utf-8")
        )
        if len(sequence) != int(cfg["expected_length"]):
            raise ValueError(
                f"Expected target length {cfg['expected_length']}, observed {len(sequence)}"
            )
        if cfg["target_accession"] not in fasta_header:
            raise ValueError("Target accession does not match canonical FASTA header")

        target_fasta = f">{cfg['target_accession']}|{cfg['target_name']}|canonical_length={len(sequence)}\n"
        target_fasta += "\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n"
        write_text(outdir / "target" / "O88319_mouse_NTSR1.fasta", target_fasta)
        write_text(outdir / "target" / "pepmlm_target.fasta", target_fasta)

        models = []
        pocket_sets: dict[str, set[str]] = {}
        for model_id, spec in MODEL_SPECS.items():
            pdb_bytes = zf.read(spec["archive_pdb"])
            cif_bytes = zf.read(spec["archive_cif"])
            ligand_bytes = zf.read(spec["ligand_pdb"])
            pdb_text = pdb_bytes.decode("utf-8")

            write_bytes(outdir / "structures" / "reference_with_ligand" / f"{model_id}.pdb", pdb_bytes)
            write_bytes(outdir / "structures" / "reference_with_ligand" / f"{model_id}.cif", cif_bytes)
            write_text(
                outdir / "structures" / "receptor_only" / f"{model_id}_receptor_only.pdb",
                pdb_receptor_only(pdb_text),
            )
            write_bytes(outdir / "structures" / "ligand_reference" / Path(spec["ligand_pdb"]).name, ligand_bytes)

            residues = parse_pdb_residues(pdb_text)
            mismatches = [
                {"residue": resnum, "pdb": aa, "canonical": sequence[resnum - 1]}
                for resnum, aa in sorted(residues.items())
                if not (1 <= resnum <= len(sequence)) or sequence[resnum - 1] != aa
            ]
            pocket_rows = parse_pocket_csv(zf.read(spec["pocket_csv"]).decode("utf-8"))
            pocket_sets[model_id] = {residue_label(row) for row in pocket_rows}
            missing_ranges = parse_missing_ranges(zf.read(spec["missing_csv"]).decode("utf-8"))
            models.append({
                "model_id": model_id,
                "state": spec["state"],
                "ligand": spec["ligand"],
                "modeled_residue_count": len(residues),
                "modeled_residue_min": min(residues),
                "modeled_residue_max": max(residues),
                "canonical_sequence_mismatch_count": len(mismatches),
                "canonical_sequence_mismatches": mismatches,
                "missing_ranges": missing_ranges,
                "pocket_residues_5A": sorted(pocket_sets[model_id], key=lambda x: int(x[1:])),
                "files": {
                    "reference_pdb": f"structures/reference_with_ligand/{model_id}.pdb",
                    "reference_cif": f"structures/reference_with_ligand/{model_id}.cif",
                    "receptor_only_pdb": f"structures/receptor_only/{model_id}_receptor_only.pdb",
                    "ligand_reference": f"structures/ligand_reference/{Path(spec['ligand_pdb']).name}",
                },
            })

        union = set().union(*pocket_sets.values())
        intersection = set.intersection(*pocket_sets.values())
        residue_sort = lambda labels: sorted(labels, key=lambda x: int(x[1:]))

        manifest = {
            "schema_version": 1,
            "stage": "01_prepare_target",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "source_archive": {
                "path": str(archive_path),
                "sha256": archive_sha,
                "internal_checksums": checksum_status,
            },
            "target": {
                "name": cfg["target_name"],
                "accession": cfg["target_accession"],
                "species": cfg["species"],
                "sequence_length": len(sequence),
                "sequence_sha256": sha256_bytes(sequence.encode("ascii")),
                "canonical_fasta": "target/O88319_mouse_NTSR1.fasta",
                "pepmlm_input_fasta": "target/pepmlm_target.fasta",
                "putative_n_glyco_sequons": find_n_glyco_sequons(sequence),
            },
            "models": models,
            "orthosteric_reference": {
                "purpose": "negative_design_reference_not_positive_peptide_docking_box",
                "union": residue_sort(union),
                "intersection": residue_sort(intersection),
            },
            "readiness": {
                "pepmlm_generation_ready": True,
                "sequence_filtering_ready": True,
                "production_structure_screening_ready": False,
                "structure_screening_blockers": [
                    "Complete and ensemble-sample missing ECL2 residues in both receptor states.",
                    "Complete and ensemble-sample missing ECL3 residues 335-340 in the 6ZIN-derived state.",
                    "Represent the missing extracellular N-terminus and its glycan steric shield.",
                    "Assign a common membrane orientation and extracellular half-space.",
                    "Add protonation/hydrogens and perform membrane-aware restrained minimization.",
                    "Do not reuse the supplied small-molecule grids as peptide docking boxes.",
                ],
            },
            "reference_peptide_policy": {
                "excluded_from_generation": True,
                "excluded_from_generation_configuration": True,
                "evaluate_only_after_de_novo_candidate_freeze": True,
            },
        }

        manifest_path = outdir / "manifest" / "target_manifest.json"
        write_text(manifest_path, json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")

        masks = {
            "target_accession": cfg["target_accession"],
            "orthosteric_pocket_union": manifest["orthosteric_reference"]["union"],
            "orthosteric_pocket_intersection": manifest["orthosteric_reference"]["intersection"],
            "putative_extracellular_glycan_sites": [4, 38, 42, 211],
            "do_not_dock_through_missing_segments": {
                model["model_id"]: model["missing_ranges"] for model in models
            },
            "small_molecule_grids_are_peptide_docking_boxes": False,
        }
        write_text(
            outdir / "manifest" / "structure_screening_masks.json",
            json.dumps(masks, indent=2, ensure_ascii=False) + "\n",
        )

        pocket_csv_path = outdir / "manifest" / "orthosteric_pocket_union.csv"
        pocket_csv_path.parent.mkdir(parents=True, exist_ok=True)
        with pocket_csv_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["residue", "in_6Z8N", "in_6ZIN", "in_both"])
            for label in residue_sort(union):
                in_a = label in pocket_sets["mouse_NTSR1_6Z8N"]
                in_i = label in pocket_sets["mouse_NTSR1_6ZIN"]
                writer.writerow([label, int(in_a), int(in_i), int(in_a and in_i)])

        write_text(outdir / "reports" / "prepare_target_report.md", render_report(manifest))

    return outdir


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/01_prepare_target.json"),
        help="Path to the preparation JSON config (default: config/01_prepare_target.json)",
    )
    args = parser.parse_args()
    try:
        outdir = prepare(args.config)
    except Exception as exc:
        print(f"01_prepare_target FAILED: {exc}", file=sys.stderr)
        return 1
    print(f"01_prepare_target PASS: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
