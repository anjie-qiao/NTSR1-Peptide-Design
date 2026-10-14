"""Reuse stage results in work/ and run absent stages on the server."""
from __future__ import annotations

import argparse
import copy
import csv
import json
import os
import shutil
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TITLES = (
    "靶标准备", "多肽生成", "序列质量筛选", "聚类去冗余", "初始结构任务准备",
    "初始结构预测与筛选", "AfDesign 局部优化", "变体结构任务准备",
    "变体结构预测与复核", "最终候选优选", "酶切风险注释",
)
MANIFESTS = (
    "manifest/target_manifest.json", "generation_manifest.json", "filter_manifest.json",
    "selection_manifest.json", "structure_job_manifest.json", "structure_filter_manifest.json",
    "afdesign_manifest.json", "structure_job_manifest.json", "structure_audit_manifest.json",
    "final_selection_manifest.json", "oral_cleavage_manifest.json",
)
# Existing output files are assumed complete; archived raw artifacts are not checked.
OUTPUTS = (
    ("target/pepmlm_target.fasta",), ("raw_candidates.csv",), ("ppl_selected.csv",),
    ("candidates_core.csv",), ("structure_tasks_all.csv",),
    ("structure_metrics_all.csv", "structure_pass.csv"),
    ("afdesign_variants_pass.csv",), ("structure_tasks_afdesign.csv",),
    ("optimized_structure_metrics_all.csv", "parent_variant_comparison_all.csv"),
    ("final_candidates_100.csv",), ("oral_cleavage_risk.csv",),
)


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


@dataclass
class Record:
    complete: bool
    summary: str = ""


class Runner:
    def __init__(self, root=ROOT, work_dir=None, summary_only=False, config_dir=None,
                 pepmlm_env="pepmlm", afdesign_env="afdesign", boltz_env="boltz"):
        self.root = Path(root).resolve()
        self.work = Path(work_dir or self.root / "work").resolve()
        self.summary_only = summary_only
        self.envs = {"pepmlm": pepmlm_env, "afdesign": afdesign_env, "boltz": boltz_env}
        self.stages = {int(p.name[:2]): p for p in sorted((self.root / "pipeline").glob("[0-9][0-9]_*.py"))}
        config_dir = Path(config_dir or self.root / "config")
        self.config_names = {int(p.name[:2]): p.name for p in sorted(config_dir.glob("[0-9][0-9]_*.json"))}
        self.cfgs = {n: self.remap(read_json(config_dir / name)) for n, name in self.config_names.items()}
        self.executed = []
        self.reused = 0

    def remap(self, value):
        if isinstance(value, dict):
            return {key: self.remap(item) for key, item in value.items()}
        if isinstance(value, list):
            return [self.remap(item) for item in value]
        if isinstance(value, str):
            if value.startswith("work/"):
                return str(self.work / value[5:])
            if value.startswith(("models/", "data/")):
                return str(self.root / value)
        return value

    def result_paths(self, n):
        cfg = self.cfgs[n]
        names = OUTPUTS[n - 1] if n != 4 else (cfg["output_basename"] + ".csv",)
        return [Path(cfg["output_dir"]) / name for name in names]

    def has_results(self, n):
        return all(path.is_file() for path in self.result_paths(n))

    def plan(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        plan = dict(stages=[self.stages[n].name for n in range(1, 12)], configs=self.cfgs,
                    output=str(self.work), actions=[dict(stage=n, title=title,
                    action="reuse" if self.has_results(n) else "run") for n, title in enumerate(TITLES, 1)])
        path.write_text(json.dumps(plan, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(f"流程预览已生成：{Path(os.path.relpath(path, self.root)).as_posix()}", flush=True)
        print("Stages 01-11, two Boltz calls as needed; no computation started.", flush=True)
        return 0

    def inspect(self, n):
        cfg = self.cfgs[n]
        directory = Path(cfg["output_dir"])
        names = OUTPUTS[n - 1] if n != 4 else (cfg["output_basename"] + ".csv",)
        if not self.has_results(n):
            return Record(False)
        manifest = directory / MANIFESTS[n - 1]
        m = read_json(manifest) if manifest.is_file() else {}
        c = m.get("counts", {})
        if n == 1:
            target = m.get("target", {})
            sequence = "".join(line.strip() for line in (directory / names[0]).read_text().splitlines()
                               if not line.startswith(">"))
            text = f"靶标 {target.get('accession', cfg.get('target_name', 'NTSR1'))}，序列长度 {len(sequence)} aa"
        elif n == 2 and m:
            text = f"已生成 {m['outputs']['candidate_count']:,} 条多肽，去重后 {m['outputs']['unique_sequence_count']:,} 条"
        elif n == 3 and m:
            text = f"理化过滤通过 {c['hard_filter_pass']:,} 条，完成 PPL 计算 {c['ppl_computed']:,} 条"
        elif n == 4:
            size = c.get("core_cluster_representatives")
            size = len(rows(directory / names[0])) if size is None else size
            text = f"聚类去冗余后保留 {size:,} 条核心候选"
        elif n in (5, 8):
            size = c.get("core_tasks" if n == 5 else "boltz_tasks")
            size = len(rows(directory / names[0])) if size is None else size
            text = f"已准备 {size:,} 个{'初始复合物' if n == 5 else '变体'}结构任务"
        elif n == 6 and m:
            text = f"完成 {c['prediction_ok']:,} 个结构预测，筛选通过 {c['structure_pass']:,} 个亲本"
        elif n == 7 and m:
            text = f"优化 {c['input_parents']:,} 个亲本，生成 {c['generated_variants']:,} 个变体，序列筛选通过 {c['sequence_filter_pass_variants']:,} 个"
        elif n == 9:
            size = c.get("prediction_ok")
            size = len(rows(directory / names[0])) if size is None else size
            text = f"完成 {size:,} 个变体独立结构复核"
        elif n == 10:
            selected = rows(directory / names[0])
            kinds = Counter(row['source_type'] for row in selected)
            text = f"最终候选 {len(selected):,} 条：亲本 {kinds['parent']:,} 条，优化变体 {kinds['afdesign']:,} 条"
        elif n == 11 and m:
            tiers = m['relative_tier_counts']
            text = f"完成 {m['candidate_count']:,} 条注释，相对低 / 中 / 高风险：{tiers['relative_low']} / {tiers['relative_medium']} / {tiers['relative_high']}"
        else:
            text = f"结果表 {names[-1]}：{len(rows(directory / names[-1])):,} 条记录"
        return Record(True, text)

    def command(self, argv, label):
        log = self.work / "launcher/logs" / (label + ".log")
        log.parent.mkdir(parents=True, exist_ok=True)
        env = dict(os.environ, PYTHONUTF8="1", PYTHONUNBUFFERED="1", WANDB_MODE="disabled",
                   BOLTZ_CACHE=str(self.root / "models/boltz_cache"))
        print(f"    运行日志：{log}", flush=True)
        with log.open("w", encoding="utf-8") as handle:
            handle.write(subprocess.list2cmdline([str(a) for a in argv]) + "\n")
            handle.flush()
            result = subprocess.run(argv, cwd=self.root, env=env, stdout=handle, stderr=subprocess.STDOUT)
        if result.returncode:
            tail = "\n".join(log.read_text(encoding="utf-8", errors="replace").splitlines()[-12:])
            raise RuntimeError(f"{label} 运行失败（退出码 {result.returncode}），日志：{log}\n{tail}")

    def conda(self, kind):
        executable = shutil.which("conda") or os.environ.get("CONDA_EXE")
        if not executable:
            raise RuntimeError("未找到 conda，请先初始化服务器上的 Conda 环境")
        return [executable, "run", "--no-capture-output", "-n", self.envs[kind]]

    def effective_config(self, n):
        cfg = copy.deepcopy(self.cfgs[n])
        # Match the original full runner: derive downstream counts from actual inputs.
        keys = {
            7: {"expected_input_count": "input_csv"},
            8: {"expected_variant_count": "variants_csv"},
            9: {"expected_task_count": "tasks_csv", "expected_parent_metrics_count": "parent_metrics_csv"},
            10: {"expected_parent_count": "parent_metrics_csv", "expected_variant_count": "optimized_metrics_csv"},
            11: {"expected_input_count": "input_csv"},
        }
        for count_key, input_key in keys.get(n, {}).items():
            cfg[count_key] = len(rows(cfg[input_key]))
        cfg_path = self.work / "launcher/config" / self.config_names[n]
        cfg_path.parent.mkdir(parents=True, exist_ok=True)
        cfg_path.write_text(json.dumps(cfg, ensure_ascii=False, indent=2), encoding="utf-8")
        return cfg_path

    def execute_stage(self, n):
        config_path = self.effective_config(n)
        if n in (1, 11):
            command = [sys.executable]
        else:
            kind = "afdesign" if n == 7 else "boltz" if n in (6, 9) else "pepmlm"
            command = self.conda(kind) + ["python"]
        self.command(command + [str(self.stages[n]), "--config", str(config_path)], f"{n:02d}")

    def ensure_prediction(self, prep):
        phase = "core" if prep == 5 else "afdesign"
        output = Path(self.cfgs[prep + 1]["boltz_output_dir"])
        model = self.cfgs[prep + 1].get("model_index", 0)
        if next(output.rglob(f"*_model_{model}.pdb"), None) is not None:
            print("    已有 Boltz 预测结果，直接计算筛选指标", flush=True)
            return
        print("    Boltz 结构预测运行", flush=True)
        cfg = self.cfgs[prep]
        prediction = cfg["boltz_prediction"]
        command = self.conda("boltz") + ["boltz", "predict",
                   str(Path(cfg["output_dir"]) / ("boltz_inputs_" + phase)),
                   "--out_dir", str(output), "--cache", str(self.root / "models/boltz_cache"),
                   "--accelerator", "gpu", "--seed", str(prediction.get("seed", 42))]
        for flag in ("devices", "recycling_steps", "sampling_steps", "diffusion_samples", "max_msa_seqs", "output_format"):
            command += ["--" + flag, str(prediction[flag])]
        for flag in ("use_potentials", "write_full_pae", "no_kernels"):
            if prediction[flag]:
                command.append("--" + flag)
        self.command(command, f"{prep:02d}_boltz")

    def run(self):
        missing = []
        for n, title in enumerate(TITLES, 1):
            record = self.inspect(n)
            if record.complete:
                print(f"[{n}/11] {title}\n    {record.summary}", flush=True)
                self.reused += 1
                continue
            if self.summary_only:
                print(f"[{n}/11] {title}（待运行）", flush=True)
                missing.append(n)
                continue
            print(f"[{n}/11] {title}阶段运行", flush=True)
            if n in (6, 9):
                self.ensure_prediction(n - 1)
            self.execute_stage(n)
            record = self.inspect(n)
            if not record.complete:
                raise RuntimeError(f"阶段 {n:02d} 未生成对应结果，请查看 launcher/logs/{n:02d}.log")
            self.executed.append(n)
            print(f"    完成：{record.summary}", flush=True)
        if missing:
            print(f"待运行阶段：{', '.join(str(n) for n in missing)}", flush=True)
            return 2
        final = Path(os.path.relpath(Path(self.cfgs[10]['output_dir']) / 'final_candidates_100.csv', self.root)).as_posix()
        risk = Path(os.path.relpath(Path(self.cfgs[11]['output_dir']) / 'oral_cleavage_risk.csv', self.root)).as_posix()
        print(f"流程完成：运行了 {len(TITLES)} 个阶段。", flush=True)
        print(f"候选表：{final}", flush=True)
        print(f"酶切注释表：{risk}", flush=True)
        return 0


def main():
    parser = argparse.ArgumentParser(description="复用已有结果，运行缺失阶段，输出 11 阶段摘要。")
    parser.add_argument("--summary-only", action="store_true", help="只展示结果，不运行缺失阶段")
    parser.add_argument("--work-dir", default="work", help="结果目录，默认 work；新计算可指定 work/new_run")
    parser.add_argument("--config-dir", type=Path, help="读取指定配置目录，默认 config")
    parser.add_argument("--plan-only", action="store_true", help="生成执行计划，不运行模型")
    parser.add_argument("--plan-output", type=Path, help="执行计划文件路径")
    for name in ("pepmlm", "afdesign", "boltz"):
        parser.add_argument(f"--{name}-env", default=name)
    args = parser.parse_args()
    path = Path(args.work_dir)
    runner = Runner(work_dir=path if path.is_absolute() else ROOT / path, summary_only=args.summary_only,
                    config_dir=args.config_dir,
                    pepmlm_env=args.pepmlm_env, afdesign_env=args.afdesign_env, boltz_env=args.boltz_env)
    return runner.plan(args.plan_output or runner.work / "launcher/execution_plan.json") if args.plan_only else runner.run()


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (RuntimeError, ValueError, OSError, KeyError) as error:
        print(f"ERROR: {error}", file=sys.stderr, flush=True)
        sys.exit(1)
