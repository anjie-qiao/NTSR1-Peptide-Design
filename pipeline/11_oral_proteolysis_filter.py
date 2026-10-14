#!/usr/bin/env python3
"""Stage 11: candidate sequences to cohort-relative oral cleavage risk levels."""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

from pepmlm_pipeline_common import (
    CANONICAL_AAS, atomic_write_csv, atomic_write_json, atomic_write_text,
    canonical_json_hash, sha256_file,
)

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = ("gingipains", "DPPIV", "neutrophil_elastase", "cathepsins", "MMPs")
PROFILE_IDS = {"neutrophil_elastase": "S01.131", "cathepsin_B": "C01.060",
               "cathepsin_L": "C01.032", "MMP8": "M10.002", "MMP9": "M10.004"}
MOTIF_RULES = {"P.gingivalis_RgpA": "after_R", "P.gingivalis_RgpB": "after_R",
               "P.gingivalis_Kgp": "after_K", "human_DPP4": "preferred_Nterminal_XP_XA",
               "P.gingivalis_DPPIV": "preferred_Nterminal_XP_XA"}
TIER_CN = {"relative_low": "相对较低", "relative_medium": "相对中等", "relative_high": "相对较高"}
FAMILY_CN = {"gingipains": "牙龈蛋白酶（Rgp/Kgp）", "DPPIV": "DPPIV/DPP4",
             "neutrophil_elastase": "中性粒细胞弹性蛋白酶", "cathepsins": "组织蛋白酶B/L",
             "MMPs": "MMP-8/9"}


def resolve(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def validate_sequence(sequence: str) -> None:
    if not sequence or set(sequence) - set(CANONICAL_AAS):
        raise ValueError("Expected canonical uppercase amino acids; chemistry must not be silently normalized")


def potential_bonds(sequence: str, rule: str) -> list[int]:
    """Return 1-based P1 positions, excluding the nonexistent bond after the end.

    DPPIV checks only the intact peptide's initial N terminus and requires a
    free N terminus. A motif match is not evidence of actual bond cleavage.
    """
    validate_sequence(sequence)
    if rule not in set(MOTIF_RULES.values()):
        raise ValueError(f"Unsupported motif rule: {rule}")
    if rule.startswith("after_"):
        return [i + 1 for i, aa in enumerate(sequence[:-1]) if aa == rule[-1]]
    return [2] if len(sequence) >= 3 and sequence[1] in "PA" else []


def profile_windows(sequence: str, profile: dict, alpha: float) -> list[dict]:
    """Geometric mean of smoothed, column-max-normalized P4-P4' counts.

    Only full windows are scored. This compatibility is not a probability;
    missing terminal flanks and unobserved residues do not imply resistance.
    """
    validate_sequence(sequence)
    if not math.isfinite(alpha) or alpha <= 0:
        raise ValueError("pseudocount must be positive and finite")
    counts = profile["counts"]
    if set(counts) != set(CANONICAL_AAS):
        raise ValueError("Profile must contain all 20 amino acids")
    for values in counts.values():
        if len(values) != 8 or any(type(x) is not int or x < 0 for x in values):
            raise ValueError("Profile counts must be eight nonnegative integers")
    maxima = [max(counts[aa][j] for aa in CANONICAL_AAS) for j in range(8)]
    if any(x == 0 for x in maxima):
        raise ValueError("A profile column has no experimental observations")
    windows = []
    for start in range(len(sequence) - 7):
        window = sequence[start:start + 8]
        log_score = math.fsum(math.log((counts[aa][j] + alpha) / (maxima[j] + alpha))
                        for j, aa in enumerate(window)) / 8
        windows.append({"position_after_1based": start + 4,
                        "context": window[:4] + "|" + window[4:],
                        "compatibility": math.exp(log_score)})
    return windows


def percentiles(values: list[float]) -> list[float]:
    """Tie-preserving midrank; 12-decimal quantization suppresses platform noise."""
    if not values or any(not math.isfinite(x) for x in values):
        raise ValueError("Percentiles require finite values")
    if len(values) == 1:
        return [50.0]
    values = [round(value, 12) for value in values]
    return [100 * (sum(x < v for x in values) + (sum(x == v for x in values) - 1) / 2)
            / (len(values) - 1) for v in values]


def relative_tier(value: float, cutoffs: list[float]) -> str:
    return "relative_low" if value <= cutoffs[0] else "relative_medium" if value <= cutoffs[1] else "relative_high"


def run(config_path: Path) -> Path:
    config = json.loads(config_path.read_text(encoding="utf-8-sig"))
    input_path = resolve(config["input_csv"])
    snapshot_path = resolve(config["specificity_snapshot"])
    output = resolve(config["output_dir"])
    if output.resolve() == input_path.parent.resolve():
        raise ValueError("Use an independent output directory")
    candidates = read_csv(input_path)
    if not candidates or len(candidates) != config["expected_input_count"]:
        raise ValueError(f"Unexpected input candidate count: {len(candidates)}; expected {config['expected_input_count']}")
    if any(not {"candidate_id", "sequence"} <= row.keys() for row in candidates):
        raise ValueError("Input CSV requires candidate_id and sequence columns")
    sequences = {r["candidate_id"]: r["sequence"] for r in candidates}
    if len(sequences) != len(candidates) or any(not cid or not cid.strip() for cid in sequences):
        raise ValueError("Candidate IDs must be nonempty and unique")
    for seq in sequences.values():
        validate_sequence(seq)
        if len(seq) < 8:
            raise ValueError("At least 8 residues required for full-window matrix comparison")
    cutoffs = config["relative_tier_cutoffs"]
    if len(cutoffs) != 2 or not 0 < cutoffs[0] < cutoffs[1] < 100:
        raise ValueError("Relative tier cutoffs must increase within 0-100")
    if config["profile_mapping"] != PROFILE_IDS:
        raise ValueError("Profile mapping must match the five specified MEROPS enzyme IDs")
    panel = config["motif_panel"]
    if (len(panel) != len(MOTIF_RULES)
            or {e["enzyme"]: e["rule"] for e in panel} != MOTIF_RULES
            or any(not e.get("source") for e in panel)):
        raise ValueError("Motif panel must contain the five sourced, specified rules")
    snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
    if snapshot["positions"] != ["P4", "P3", "P2", "P1", "P1prime", "P2prime", "P3prime", "P4prime"]:
        raise ValueError("Snapshot must use the P4-P4prime position order")
    rows, site_rows, enzyme_rows, motif_rows = [], [], [], []
    for candidate in candidates:
        cid, seq = candidate["candidate_id"], candidate["sequence"]
        total_bonds = len(seq) - 1
        rgp = potential_bonds(seq, "after_R")
        kgp = potential_bonds(seq, "after_K")
        dpp_positions = potential_bonds(seq, "preferred_Nterminal_XP_XA")
        row = {
            "candidate_id": cid, "sequence": seq, "length": len(seq),
            "original_final_rank": candidate.get("final_rank", ""),
            "gingipains_burden": (len(rgp) + len(kgp)) / total_bonds,
            "DPPIV_burden": float(bool(dpp_positions)),
            "rgp_positions": ";".join(map(str, rgp)),
            "kgp_positions": ";".join(map(str, kgp)),
            "dppiv_preferred_initial_position": ";".join(map(str, dpp_positions)),
            "DPPIV_condition": "free_N_terminus_required_chemistry_unverified",
        }
        for entry in panel:
            positions = potential_bonds(seq, entry["rule"])
            condition = "free_N_terminus_required_chemistry_unverified" if "Nterminal" in entry["rule"] else "P1_rule_only_extended_specificity_unverified"
            for pos in positions:
                motif_rows.append({
                    "candidate_id": cid, "enzyme": entry["enzyme"], "rule": entry["rule"],
                    "position_after_1based": pos,
                    "bond": f"{seq[pos - 1]}{pos}|{seq[pos]}{pos + 1}",
                    "context": seq[max(0, pos - 4):pos] + "|" + seq[pos:pos + 4],
                    "condition": condition, "source": entry["source"],
                })
        for name, merops_id in config["profile_mapping"].items():
            profile = snapshot["entries"][merops_id]
            windows = profile_windows(seq, profile, float(config["pseudocount"]))
            peak = max(windows, key=lambda w: w["compatibility"])
            mean = math.fsum(w["compatibility"] for w in windows) / len(windows)
            row[name + "_mean_compatibility"] = mean
            row[name + "_peak_compatibility"] = peak["compatibility"]
            row[name + "_peak_position"] = peak["position_after_1based"]
            enzyme_rows.append({
                "candidate_id": cid, "enzyme_profile": name, "merops_id": merops_id,
                "mean_compatibility": mean, "peak_compatibility": peak["compatibility"],
                "peak_position_after_1based": peak["position_after_1based"],
                "peak_context": peak["context"], "full_windows_assessed": len(windows),
                "total_peptide_bonds": total_bonds, "terminal_bonds_unscored": total_bonds - len(windows),
                "profile_bond_coverage": len(windows) / total_bonds, "source": profile["source"],
            })
            site_rows.extend({"candidate_id": cid, "enzyme_profile": name,
                              "merops_id": merops_id, **window} for window in windows)
        row["profile_bond_coverage"] = (len(seq) - 7) / total_bonds
        row["unscored_terminal_bond_count_per_profile"] = 6
        row["unassessed_mechanisms"] = "aminopeptidase_trimming;cathepsin_B_exopeptidase;fragment_recleavage;nonpreferred_DPPIV;terminal_profile_windows"
        row["evidence_coverage"] = "partial_sequence_only"
        row["chemistry_assumption"] = config["chemistry_assumption"]
        rows.append(row)

    # Each matrix is ranked only against itself across the candidate cohort.
    for name in config["profile_mapping"]:
        for row, rank in zip(rows, percentiles([r[name + "_mean_compatibility"] for r in rows])):
            row[name + "_profile_percentile"] = rank
    for family in ("gingipains", "DPPIV"):
        for row, rank in zip(rows, percentiles([r[family + "_burden"] for r in rows])):
            row[family + "_family_percentile"] = rank
    for row in rows:
        row["neutrophil_elastase_family_percentile"] = row["neutrophil_elastase_profile_percentile"]
        row["cathepsins_family_percentile"] = max(row["cathepsin_B_profile_percentile"], row["cathepsin_L_profile_percentile"])
        row["MMPs_family_percentile"] = max(row["MMP8_profile_percentile"], row["MMP9_profile_percentile"])
    # Re-rank family maxima so two-profile families do not acquire higher weight.
    for family in ("cathepsins", "MMPs"):
        for row, rank in zip(rows, percentiles([r[family + "_family_percentile"] for r in rows])):
            row[family + "_family_percentile"] = rank
    for row in rows:
        ranks = {f: row[f + "_family_percentile"] for f in FAMILIES}
        row["sequence_panel_index"] = math.fsum(ranks.values()) / len(FAMILIES)
        # Zero motif burdens can have a nonzero tied percentile, but no motif warning.
        eligible = {f: v for f, v in ranks.items()
                    if f not in {"DPPIV", "gingipains"} or row[f + "_burden"] > 0}
        row["highest_family_percentile"] = max(eligible.values())
        row["highest_relative_family"] = ";".join(f for f, v in eligible.items() if v == max(eligible.values()))
    for row, rank in zip(rows, percentiles([r["sequence_panel_index"] for r in rows])):
        row["relative_sequence_risk_percentile"] = rank
        row["relative_sequence_risk_tier"] = relative_tier(rank, cutoffs)
        row["relative_sequence_risk_label_cn"] = TIER_CN[row["relative_sequence_risk_tier"]]
        row["risk_basis"] = "five_family_equal_weight_cohort_relative_v1_1_not_calibrated"

    # One row per input candidate, in input order; never auto-select experimental peptides.
    risk_rows = [{
        "candidate_id": r["candidate_id"], "sequence": r["sequence"],
        "oral_cleavage_risk_level": r["relative_sequence_risk_label_cn"],
        "oral_cleavage_risk_percentile": r["relative_sequence_risk_percentile"],
        "main_risk_source": ";".join(FAMILY_CN[f] for f in r["highest_relative_family"].split(";")),
        "risk_scope": "本批次100条相对等级" if len(rows) == 100 else "本批次相对等级",
        "evidence_coverage": "部分酶和位点；末端化学状态未核实",
    } for r in rows]
    counts = Counter(r["relative_sequence_risk_tier"] for r in rows)
    atomic_write_csv(output / "oral_cleavage_risk.csv", risk_rows, list(risk_rows[0]))
    atomic_write_csv(output / "evidence/oral_cleavage_metrics_all.csv", rows, list(rows[0]))
    atomic_write_csv(output / "evidence/enzyme_profile_summary.csv", enzyme_rows, list(enzyme_rows[0]))
    atomic_write_csv(output / "evidence/profile_window_scores.csv", site_rows, list(site_rows[0]))
    atomic_write_csv(output / "evidence/motif_cleavage_sites.csv", motif_rows,
                     ["candidate_id", "enzyme", "rule", "position_after_1based", "bond", "context", "condition", "source"])
    manifest = {
        "stage": "11_oral_proteolysis_filter", "method": "experimental_specificity_heuristic_v1.1",
        "numerical_policy": "math.fsum; percentile inputs rounded to 12 decimal places before tie comparisons",
        "candidate_count": len(rows), "trained_model_run": False, "calibrated_probability": False,
        "experimental_validation": False, "selection_automated": False,
        "input_csv": str(input_path), "input_sha256": sha256_file(input_path),
        "source_snapshot_sha256": sha256_file(snapshot_path),
        "config_sha256": canonical_json_hash(config), "families": list(FAMILIES),
        "family_weights": {f: 1 / len(FAMILIES) for f in FAMILIES},
        "relative_tier_cutoffs": cutoffs, "relative_tier_counts": dict(counts),
        "motif_panel": panel, "profile_mapping": config["profile_mapping"],
        "snapshot_scope": snapshot["scope"], "output_order": "same_as_input",
        "limitations": [
            "Equal family weights and tertiles are review policy, not biological cutoffs.",
            "Re-ranking another cohort can change all percentile-based strata.",
            "Only complete P4-P4prime windows scored; six terminal bonds per profile are unknown.",
            "Summarized entries may pool species, experimental conditions and redundant substrates.",
            "Mean compatibility can hide a susceptible bond; retain peak positions and worst family.",
            "No free-peptide MD, accessibility, enzyme activity or terminal chemistry in this index.",
            "Aminopeptidases, cathepsin B exopeptidase and fragment re-cleavage remain unassessed.",
            "Low relative tier is not a biological resistance claim.",
        ],
    }
    atomic_write_json(output / "oral_cleavage_manifest.json", manifest)
    lines = ["# 第11步：口腔酶切风险初筛结果", "",
             f"输入 {len(rows)} 条候选，输出 {len(risk_rows)} 条风险记录，顺序与输入一致。", "",
             "主结果：oral_cleavage_risk.csv。等级表示本批次的相对序列风险，百分位越高提示越多，不是切割概率。", "",
             "| 风险等级 | 数量 |", "|---|---:|"]
    for tier in ("relative_low", "relative_medium", "relative_high"):
        lines.append(f"| {TIER_CN[tier]} | {counts[tier]} |")
    lines.extend(["", "规则命中与 MEROPS 切点邻域统计在本步骤内完成；无需旧中间文件或 MD 表格。",
                  "每酶摘要、中央窗口和规则切点依据见 evidence/。算法及字段见 docs/11_口腔酶切风险初筛说明.md。",
                  "本结果可与独立分子动力学评估共同用于人工选出5条实验候选。相对较低仍需酶切实验验证。", ""])
    atomic_write_text(output / "口腔酶切风险初筛结果.md", "\n".join(lines))
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/11_oral_proteolysis_filter.json")
    args = parser.parse_args()
    print(run(args.config))
