"""Run the NTSR1 peptide design pipeline or export final results."""
import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

from pipeline.submission_export import export_results, read_rows, sha256

ROOT = Path(__file__).resolve().parent
STAGES = {int(p.name[:2]): p for p in sorted((ROOT / "pipeline").glob("[0-9][0-9]_*.py"))}
CONFIGS = {int(p.name[:2]): p for p in sorted((ROOT / "config").glob("[0-9][0-9]_*.json"))}


def save_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def local_path(value):
    path = Path(value)
    return (path if path.is_absolute() else ROOT / path).resolve()


def remap(value, base):
    if isinstance(value, dict):
        return {key: remap(item, base) for key, item in value.items()}
    if isinstance(value, list):
        return [remap(item, base) for item in value]
    if isinstance(value, str):
        if value.startswith("work/"):
            return str(base / value[5:])
        if value.startswith(("models/", "data/")):
            return str(ROOT / value)
    return value


def configs(base, args):
    result = {n: remap(json.loads(p.read_text(encoding="utf-8-sig")), base) for n, p in CONFIGS.items()}
    if args.input:
        result[1]["archive"] = str(local_path(args.input))
    for n in (5, 8):
        result[n]["boltz_prediction"].setdefault("seed", 42)
    return result


def run_command(command, log):
    log.parent.mkdir(parents=True, exist_ok=True)
    print("RUN:", subprocess.list2cmdline(command), flush=True)
    started = time.monotonic()
    env = os.environ.copy()
    env.update(PYTHONUTF8="1", PYTHONUNBUFFERED="1", BOLTZ_CACHE=str(ROOT / "models/boltz_cache"))
    env["WANDB_MODE"] = "disabled"
    with log.open("w", encoding="utf-8") as handle:
        with subprocess.Popen(command, cwd=ROOT, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, encoding="utf-8", errors="replace") as process:
            for line in process.stdout:
                print(line, end="", flush=True)
                handle.write(line)
                handle.flush()
            returncode = process.wait()
    if returncode:
        raise RuntimeError(f"Command failed ({returncode}): {command}; see {log}")
    return {"command": command, "seconds": time.monotonic() - started, "returncode": returncode}


def run_pipeline(base, args):
    cfgs = configs(base, args)
    plan = {"mode": args.mode, "stages": [STAGES[n].name for n in sorted(STAGES)],
            "output": str(base), "configs": cfgs,
            "note": "Two Boltz prediction calls occur after 05 and 08. Counts after 06 are derived from actual outputs; zero qualified candidates stops the run."}
    if args.plan_only:
        save_json(base / "execution_plan.json", plan)
        print(f"Execution plan: {base / 'execution_plan.json'}")
        print("Stages 01-11, two Boltz calls; no computation started.")
        return
    identity = {"mode": args.mode, "configs": {str(n): cfg for n, cfg in cfgs.items()}, "track": args.track}
    identity_path = base / "run_identity.json"
    if identity_path.exists() and json.loads(identity_path.read_text(encoding="utf-8")) != identity:
        raise ValueError("Output directory belongs to a different configuration; choose a new --output-dir")
    conda = shutil.which("conda")
    if platform.system() != "Linux":
        raise RuntimeError("GPU workflow requires Linux; use --plan-only to inspect the execution plan.")
    if not conda:
        raise RuntimeError("Conda not found; initialize Conda before running the pipeline")
    required = [ROOT / "models/PepMLM-650M/config.json", ROOT / "models/ntsr1_msa/ntsr1.csv",
                ROOT / "models/afdesign/params/params_model_1_multimer_v3.npz",
                ROOT / "models/boltz_cache/boltz2_conf.ckpt", ROOT / "models/boltz_cache/boltz2_aff.ckpt",
                ROOT / "models/boltz_cache/mols.tar"]
    missing = [str(p.relative_to(ROOT)) for p in required if not p.is_file() or p.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Prepare model assets first: " + ", ".join(missing))
    runtime = {"launcher_python": sys.version, "platform": platform.platform(), "model_packages": {}, "stage_runs": []}
    for name in (args.pepmlm_env, args.afdesign_env, args.boltz_env):
        run_command([conda, "run", "--no-capture-output", "-n", name, "python", "-m", "pip", "check"], base / f"logs/{name}_pip_check.txt")
        backend = "jax" if name == args.afdesign_env else "torch"
        result = subprocess.run([conda, "run", "--no-capture-output", "-n", name, "python",
                                 str(ROOT / "pipeline/runtime_probe.py"), "--backend", backend],
                                cwd=ROOT, text=True, encoding="utf-8", errors="replace", capture_output=True)
        probe_log = base / f"logs/{name}_probe.log"
        probe_log.parent.mkdir(parents=True, exist_ok=True)
        probe_log.write_text(result.stdout + result.stderr, encoding="utf-8")
        if result.returncode:
            raise RuntimeError(f"GPU/environment check failed for {name}; see {probe_log}: {result.stderr[-1500:]}")
        runtime["model_packages"][name] = json.loads(result.stdout.strip().splitlines()[-1])
        run_command([conda, "run", "--no-capture-output", "-n", name, "python", "-m", "pip", "freeze"], base / f"logs/{name}_pip_freeze.txt")
    if shutil.which("nvidia-smi"):
        result = subprocess.run(["nvidia-smi"], capture_output=True, text=True, errors="replace")
        runtime["nvidia_smi"] = result.stdout
    runtime_path = base / "runtime_environment.json"
    weight_files = list((ROOT / "models/PepMLM-650M").glob("*.safetensors")) + list((ROOT / "models/PepMLM-650M").glob("pytorch_model*.bin"))
    if not weight_files:
        raise FileNotFoundError("PepMLM weight files (.safetensors/.bin) are missing")
    runtime["asset_sha256"] = {str(p.relative_to(ROOT)): sha256(p) for p in dict.fromkeys(required + weight_files)}
    if runtime_path.exists():
        previous = json.loads(runtime_path.read_text(encoding="utf-8"))
        if any(previous.get(key) != runtime[key] for key in ("model_packages", "asset_sha256")):
            raise ValueError("Installed versions or model files changed within this run; choose a new output directory")
        runtime["stage_runs"] = previous.get("stage_runs", [])
    save_json(runtime_path, runtime)
    save_json(identity_path, identity)
    try:
        for n in range(args.resume_from, 12):
            cfg = cfgs[n]
            if n == 7:
                parents = read_rows(cfg["input_csv"])
                cfg["expected_input_count"] = len(parents)
                if not parents:
                    raise ValueError("No structural candidates passed")
            if n == 8:
                cfg["expected_variant_count"] = len(read_rows(cfg["variants_csv"]))
                if not cfg["expected_variant_count"]:
                    raise ValueError("No optimized variants passed; Stage 08 requires at least one variant")
            if n == 9:
                cfg["expected_task_count"] = len(read_rows(cfg["tasks_csv"]))
                cfg["parent_metrics_csv"] = cfgs[7]["input_csv"]
                cfg["expected_parent_metrics_count"] = len(read_rows(cfg["parent_metrics_csv"]))
            if n == 10:
                if args.resume_from <= 9:
                    cfg["parent_metrics_csv"] = cfgs[9]["parent_metrics_csv"]
                cfg["expected_parent_count"] = len(read_rows(cfg["parent_metrics_csv"]))
                cfg["expected_variant_count"] = len(read_rows(cfg["optimized_metrics_csv"]))
                if cfg["expected_parent_count"] < cfg["final_count"]:
                    raise ValueError("Insufficient parent families for the requested 100-candidate final selection")
            if n == 11:
                cfg["expected_input_count"] = len(read_rows(cfg["input_csv"]))
            config_path = base / "config" / CONFIGS[n].name
            save_json(config_path, cfg)
            env_name = args.afdesign_env if n == 7 else args.boltz_env if n in (6, 9) else args.pepmlm_env
            command = [conda, "run", "--no-capture-output", "-n", env_name, "python", str(STAGES[n]), "--config", str(config_path)]
            runtime["stage_runs"].append(run_command(command, base / f"logs/{n:02d}.log"))
            if n in (5, 8):
                phase = "core" if n == 5 else "afdesign"
                prediction = cfg["boltz_prediction"]
                destination = Path(cfg["output_dir"]) / f"boltz_predictions_{phase}"
                command = [conda, "run", "--no-capture-output", "-n", args.boltz_env, "boltz", "predict",
                           str(Path(cfg["output_dir"]) / f"boltz_inputs_{phase}"), "--out_dir", str(destination),
                           "--cache", str(ROOT / "models/boltz_cache"), "--accelerator", "gpu", "--seed", str(prediction["seed"])]
                for flag in ("devices", "recycling_steps", "sampling_steps", "diffusion_samples", "max_msa_seqs", "output_format"):
                    command.extend(["--" + flag, str(prediction[flag])])
                for flag in ("use_potentials", "write_full_pae", "no_kernels"):
                    if prediction[flag]:
                        command.append("--" + flag)
                runtime["stage_runs"].append(run_command(command, base / f"logs/{n:02d}_boltz.log"))
            save_json(runtime_path, runtime)
        path = export_results(Path(cfgs[10]["output_dir"]) / "final_candidates_100.csv",
                              Path(cfgs[11]["output_dir"]) / "oral_cleavage_risk.csv", base / "results",
                              [base / "05_prepare_structure_jobs", base / "08_prepare_optimized_structure_jobs"],
                              mode="full", runtime_file=runtime_path, track=args.track)
        runtime["status"] = "complete"
        print(f"Results: {path}")
    except Exception as error:
        runtime.update(status="failed", error=str(error))
        raise
    finally:
        save_json(runtime_path, runtime)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("full", "export"), default="full",
                        help="Run the complete pipeline or export existing final results")
    parser.add_argument("--input", help="NTSR1 target ZIP; defaults to the archive path in config/01_prepare_target.json")
    parser.add_argument("--output-dir", help="Path inside project; default work/<mode> or results for export")
    parser.add_argument("--track", default="AI大分子与多肽药物设计")
    parser.add_argument("--resume-from", type=int, choices=range(1, 12), default=1)
    parser.add_argument("--plan-only", action="store_true", help="Write the full execution plan without starting computation")
    parser.add_argument("--pepmlm-env", default="pepmlm")
    parser.add_argument("--afdesign-env", default="afdesign")
    parser.add_argument("--boltz-env", default="boltz")
    parser.add_argument("--final-csv", default="work/10_select_final_candidates/final_candidates_100.csv")
    parser.add_argument("--risk-csv", default="work/11_oral_proteolysis_filter/oral_cleavage_risk.csv")
    parser.add_argument("--structure-root", action="append", default=[])
    args = parser.parse_args()
    if args.mode == "export" and args.plan_only:
        parser.error("--plan-only is only supported with --mode full")
    base = local_path(args.output_dir or ("results" if args.mode == "export" else f"work/{args.mode}"))
    base.relative_to(ROOT)
    if base == ROOT or base == ROOT / "work":
        parser.error("Use a dedicated output subdirectory")
    if args.mode == "export":
        defaults = [ROOT / "work/05_prepare_structure_jobs", ROOT / "work/08_prepare_optimized_structure_jobs",
                    base / "structures"]
        result = export_results(local_path(args.final_csv), local_path(args.risk_csv), base,
                 [local_path(p) for p in args.structure_root] + [p for p in defaults if p.is_dir()],
                 track=args.track)
        print(f"Results: {result}")
    else:
        run_pipeline(base, args)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, FileNotFoundError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
