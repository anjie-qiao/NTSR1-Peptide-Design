"""Serve the NTSR1 workbench without changing the existing pipeline scripts."""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import hmac
import io
import json
import math
import os
import platform
import re
import secrets
import shutil
import subprocess
import sys
import threading
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from pipeline.cached_runner import Runner as CachedRunner

ROOT = Path(__file__).resolve().parent
JOB_ID = re.compile(r"^[0-9]{8}T[0-9]{6}-[a-f0-9]{8}$")
CONFIG_NAMES = {n: p.name for n, p in ((int(p.name[:2]), p) for p in sorted((ROOT / "config").glob("[0-9][0-9]_*.json")))}
BUSY = {"queued", "running"}


class WorkbenchError(Exception):
    def __init__(self, message, status=400):
        self.message, self.status = message, status


def now():
    return datetime.now(timezone.utc).isoformat()


def load_json(path, default=None):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8-sig"))
    except FileNotFoundError:
        return default


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + "." + secrets.token_hex(4) + ".tmp")
    try:
        temporary.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def config_revision(root):
    return hashlib.sha256(b"".join((root / "config" / name).read_bytes() for name in CONFIG_NAMES.values())).hexdigest()


def config_bundle(root):
    return {n: load_json(root / "config" / name) for n, name in CONFIG_NAMES.items()}


def nested(config, keys):
    for key in keys:
        config = config[key]
    return config


def assign(config, keys, value):
    for key in keys[:-1]:
        config = config[key]
    config[keys[-1]] = value


def fields():
    """Only expose numeric screening settings, never paths or executable commands."""
    result = []

    def add(key, label, stage, path, low, high, integer=False, group="基本参数", nullable=False, step=None):
        result.append(dict(key=key, label=label, stage=stage, path=path.split("."), min=low, max=high,
                           integer=integer, group=group, nullable=nullable,
                           step=step if step is not None else (1 if integer else .01)))

    add("total_candidates", "生成序列数量", 2, "total_candidates", 1, 1000000, True)
    add("length_min", "最短肽长 / aa", 2, "lengths", 12, 33, True)
    add("length_max", "最长肽长 / aa", 2, "lengths", 12, 33, True)
    add("final_count", "最终候选数量", 10, "final_count", 1, 100000, True)
    for key, label, low, high, integer in [
        ("min_net_charge", "最低净电荷", -100, 100, False),
        ("max_net_charge", "最高净电荷", -100, 100, False),
        ("max_hydrophobic_fraction", "最大疏水残基比例", 0, 1, False),
        ("max_hydrophobic_run", "最大连续疏水残基数", 1, 33, True),
        ("max_homopolymer_run", "最大同残基连续长度", 1, 33, True),
        ("min_normalized_entropy", "最低归一化序列熵", 0, 1, False),
        ("max_cysteine_count", "最多半胱氨酸数", 0, 33, True),
    ]:
        add("qc_" + key, label, 3, "hard_filters." + key, low, high, integer, "03 · 序列质量筛选")
    add("ppl_limit", "每个肽长最多进入 PPL 的数量", 3, "ppl_max_candidates_per_length", 1, 100000, True, "03 · PPL 计算")
    add("ppl_batch", "PPL 批大小", 3, "ppl_row_batch_size", 1, 1024, True, "03 · PPL 计算")
    add("ppl_keep", "各肽长 PPL 保留比例", 4, "ppl_keep_fraction", .01, 1, group="04 · 聚类去冗余")
    add("cluster_identity", "同长度聚类一致性阈值", 4, "cluster_identity_threshold", .01, 1, group="04 · 聚类去冗余")
    labels = [
        ("min_pair_iptm_mean", "最低平均 pair ipTM", 0, 1, False),
        ("min_pair_iptm_min", "最低最小 pair ipTM", 0, 1, False),
        ("min_peptide_ptm", "最低肽 pTM", 0, 1, False),
        ("min_peptide_mean_plddt", "最低平均肽 pLDDT", 0, 1, False),
        ("min_peptide_q10_plddt", "最低肽 pLDDT 第 10 百分位", 0, 1, False),
        ("max_interface_pae_mean", "最大平均界面 PAE / Å", 0, 1000, False),
        ("max_interface_pae_q90", "最大界面 PAE 第 90 百分位 / Å", 0, 1000, False),
        ("min_peptide_contact_fraction", "最低肽接触比例", 0, 1, False),
        ("min_interface_residue_pairs", "最少界面残基接触对数", 0, 100000, True),
        ("max_severe_clash_atom_pairs", "最多严重原子冲突对数", 0, 100000, True),
    ]
    for stage, prefix, group, config_key in [(6, "structure_", "06 · 初始结构筛选", "thresholds"),
                                              (10, "final_", "10 · 优化变体绝对结构门槛", "absolute_structure_thresholds")]:
        for key, label, low, high, integer in labels:
            add(prefix + key, label, stage, config_key + "." + key, low, high, integer, group,
                nullable=stage == 6 and key == "max_interface_pae_q90")
    add("final_identity", "最终集合聚类一致性阈值", 10, "diversity.same_length_identity_threshold", .01, 1,
        group="10 · 最终候选多样性")
    return result


FIELDS = fields()


def parameter_values(bundle):
    values = {f["key"]: nested(bundle[f["stage"]], f["path"]) for f in FIELDS if not f["key"].startswith("length_")}
    values.update(length_min=min(bundle[2]["lengths"]), length_max=max(bundle[2]["lengths"]))
    return values


def updated_configs(bundle, changes):
    if not isinstance(changes, dict) or not changes:
        raise WorkbenchError("请提交需要修改的参数。")
    allowed = {f["key"]: f for f in FIELDS}
    if set(changes) - set(allowed):
        raise WorkbenchError("包含不支持的参数。")
    values = parameter_values(bundle)
    for key, value in changes.items():
        f = allowed[key]
        if value is None and f["nullable"]:
            values[key] = None
            continue
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise WorkbenchError(f"{f['label']} 必须是有限数值。")
        if not f["min"] <= value <= f["max"] or (f["integer"] and int(value) != value):
            raise WorkbenchError(f"{f['label']} 的范围是 {f['min']}–{f['max']}，请检查输入。")
        values[key] = int(value) if f["integer"] else float(value)
    if values["length_min"] > values["length_max"]:
        raise WorkbenchError("最短肽长不能超过最长肽长。")
    if values["total_candidates"] < values["length_max"] - values["length_min"] + 1:
        raise WorkbenchError("生成数量至少应覆盖每个肽长一次。")
    if values["final_count"] > values["total_candidates"]:
        raise WorkbenchError("最终候选数不能超过生成数量；实际合格家族不足时原入口会停止。")
    if values["qc_min_net_charge"] > values["qc_max_net_charge"]:
        raise WorkbenchError("最低净电荷不能超过最高净电荷。")
    result = copy.deepcopy(bundle)
    for f in FIELDS:
        if f["key"] in changes and not f["key"].startswith("length_"):
            assign(result[f["stage"]], f["path"], values[f["key"]])
    if "length_min" in changes or "length_max" in changes:
        result[2]["lengths"] = list(range(values["length_min"], values["length_max"] + 1))
        if "length_weights" in result[2] and result[2]["lengths"] != bundle[2]["lengths"]:
            raise WorkbenchError("当前生成配置有自定义肽长权重，需先在配置文件中同步调整。")
    return result


@contextmanager
def state_lock(root, blocking=False):
    directory = root / "work" / "workbench"
    directory.mkdir(parents=True, exist_ok=True)
    handle = (directory / "state.lock").open("a+b")
    acquired = False
    try:
        if os.name == "nt":
            import msvcrt
            handle.seek(0, 2)
            if handle.tell() == 0:
                handle.write(b"0")
                handle.flush()
            handle.seek(0)
            msvcrt.locking(handle.fileno(), msvcrt.LK_LOCK if blocking else msvcrt.LK_NBLCK, 1)
        else:
            import fcntl
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB))
        acquired = True
        yield
    except (BlockingIOError, PermissionError) as error:
        if not acquired:
            raise WorkbenchError("筛选任务正在运行，当前不能保存参数或启动另一任务。", 409) from error
        raise
    finally:
        if acquired:
            if os.name == "nt":
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        handle.close()


def job_directory(root, job_id):
    if not JOB_ID.fullmatch(job_id):
        raise WorkbenchError("任务编号无效。", 404)
    directory = (root / "work" / "workbench" / job_id).resolve()
    directory.relative_to(root.resolve())
    if not (directory / "task.json").is_file():
        raise WorkbenchError("任务不存在。", 404)
    return directory


def pid_alive(pid):
    if not isinstance(pid, int) or pid <= 0:
        return False
    if os.name == "nt":
        import ctypes
        from ctypes import wintypes
        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
        kernel.OpenProcess.restype = wintypes.HANDLE
        kernel.GetExitCodeProcess.argtypes = [wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD)]
        kernel.CloseHandle.argtypes = [wintypes.HANDLE]
        handle = kernel.OpenProcess(0x1000, False, pid)
        if not handle:
            return False
        exit_code = wintypes.DWORD()
        try:
            return bool(kernel.GetExitCodeProcess(handle, ctypes.byref(exit_code))) and exit_code.value == 259
        finally:
            kernel.CloseHandle(handle)
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


def task(root, job_id):
    directory = job_directory(root, job_id)
    data = load_json(directory / "task.json")
    if data["status"] in BUSY and time.time() - data.get("launched_at", time.time()) > 15:
        if not pid_alive(data.get("worker_pid")) and not pid_alive(data.get("child_pid")):
            data.update(status="interrupted", error="运行进程已结束，未留下正常完成记录。请查看日志。")
    return data


def tasks(root):
    directory = root / "work" / "workbench"
    if not directory.exists():
        return []
    return [task(root, p.name) for p in sorted(directory.iterdir(), reverse=True)
            if JOB_ID.fullmatch(p.name) and (p / "task.json").is_file()][:50]


def assert_idle(root):
    if any(t["status"] in BUSY for t in tasks(root)):
        raise WorkbenchError("已有任务正在运行，请等待完成后再保存或启动。", 409)


def environment(root):
    required = ["models/PepMLM-650M/config.json", "models/ntsr1_msa/ntsr1.csv",
                "models/afdesign/params/params_model_1_multimer_v3.npz", "models/boltz_cache/boltz2_conf.ckpt",
                "models/boltz_cache/boltz2_aff.ckpt", "models/boltz_cache/mols.tar"]
    missing = [p for p in required if not (root / p).is_file() or (root / p).stat().st_size == 0]
    model_dir = root / "models" / "PepMLM-650M"
    if not list(model_dir.glob("*.safetensors")) + list(model_dir.glob("pytorch_model*.bin")):
        missing.append("models/PepMLM-650M/模型权重 (.safetensors/.bin)")
    reasons = ([] if platform.system() == "Linux" else ["全流程计算要求 Linux"])
    if not (shutil.which("conda") or os.environ.get("CONDA_EXE")):
        reasons.append("启动环境未找到 conda")
    if missing:
        reasons.append("部分模型资产尚未准备")
    launcher_ready = (root / "run_pipeline.sh").is_file() and (root / "pipeline/cached_runner.py").is_file()
    runner = CachedRunner(root=root)
    pending = [n for n in range(1, 12) if not runner.has_results(n)]
    needs_models = any(n not in (1, 11) for n in pending)
    if not launcher_ready:
        reasons.append("尚未安装运行脚本")
    return dict(platform=platform.system(), full_available=launcher_ready and (not needs_models or not reasons),
                reasons=reasons if needs_models or not launcher_ready else [], missing=missing if needs_models else [],
                note="运行脚本复用 work 中的已有结果，仅计算缺失阶段。")


def start_task(root, payload):
    mode = payload.get("mode")
    name = payload.get("name", "NTSR1 筛选任务")
    if mode not in {"full", "plan"} or not isinstance(name, str) or not 1 <= len(name.strip()) <= 80:
        raise WorkbenchError("请选择运行方式并填写 1–80 字的任务名称。")
    if mode == "full" and not environment(root)["full_available"]:
        raise WorkbenchError("请先准备缺失阶段所需的模型环境，或使用流程预览。")
    if not (root / "pipeline/cached_runner.py").is_file():
        raise WorkbenchError("请先上传 pipeline/cached_runner.py 和 run_pipeline.sh。")
    with state_lock(root):
        assert_idle(root)
        if payload.get("revision") != config_revision(root):
            raise WorkbenchError("参数已发生变化，请重新载入参数后再启动。", 409)
        job_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S") + "-" + secrets.token_hex(4)
        directory = root / "work" / "workbench" / job_id
        snapshot = directory / "config_snapshot"
        snapshot.mkdir(parents=True)
        for filename in CONFIG_NAMES.values():
            shutil.copy2(root / "config" / filename, snapshot / filename)
        data = dict(id=job_id, name=name.strip(), mode=mode, status="queued", created_at=now(),
                    launched_at=time.time(), revision=config_revision(root), active_stage=0, phase="preflight", runner="cached")
        atomic_json(directory / "task.json", data)
        command = [sys.executable, str(Path(__file__).resolve()), "--root", str(root), "--worker", job_id]
        env = os.environ.copy()
        env.update(PYTHONUTF8="1", PYTHONUNBUFFERED="1")
        try:
            with (directory / "worker.log").open("a", encoding="utf-8") as log:
                options = {"start_new_session": True} if os.name != "nt" else {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.CREATE_NO_WINDOW}
                process = subprocess.Popen(command, cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT,
                                           stdin=subprocess.DEVNULL, **options)
            data["worker_pid"] = process.pid
            atomic_json(directory / "task.json", data)
            threading.Thread(target=process.wait, daemon=True).start()
        except Exception:
            data.update(status="failed", error="无法创建运行进程。", finished_at=now())
            atomic_json(directory / "task.json", data)
            raise
        return data


def run_worker(root, job_id):
    directory = job_directory(root, job_id)
    with state_lock(root, blocking=True):
        data = load_json(directory / "task.json")
        log_path = directory / "launcher.log"
        try:
            if config_revision(root) != data["revision"]:
                raise RuntimeError("启动前项目配置发生变化，任务已停止；请重新创建任务。")
            command = ([sys.executable, "-u", str(root / "pipeline/cached_runner.py")] if os.name == "nt" else
                       ["bash", str(root / "run_pipeline.sh")])
            command += ["--config-dir", str(directory / "config_snapshot")]
            if data["mode"] == "plan":
                command += ["--plan-only", "--plan-output", str(directory / "execution_plan.json")]
            data.update(status="running", started_at=now(), command=command, worker_pid=os.getpid())
            atomic_json(directory / "task.json", data)
            env = os.environ.copy()
            env.update(PYTHONUTF8="1", PYTHONUNBUFFERED="1", PYTHON_BIN=sys.executable)
            with log_path.open("w", encoding="utf-8") as log:
                with subprocess.Popen(command, cwd=root, env=env, stdin=subprocess.DEVNULL,
                                      stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                      text=True, encoding="utf-8", errors="replace") as child:
                    data["child_pid"] = child.pid
                    atomic_json(directory / "task.json", data)
                    for line in child.stdout:
                        log.write(line)
                        log.flush()
                        match = re.match(r"\[(\d+)/11\]", line)
                        if match:
                            data.update(active_stage=int(match[1]), phase="stage")
                            atomic_json(directory / "task.json", data)
                        elif "Boltz 结构预测运行" in line:
                            data["phase"] = "boltz"
                            atomic_json(directory / "task.json", data)
                    code = child.wait()
            data.update(returncode=code, finished_at=now())
            runner = CachedRunner(root=root, config_dir=directory / "config_snapshot")
            expected = (directory / "execution_plan.json").is_file() if data["mode"] == "plan" else all(runner.has_results(n) for n in (10, 11))
            if code == 0 and expected:
                data.update(status="planned" if data["mode"] == "plan" else "complete", phase="finished")
            else:
                data.update(status="failed", error=f"入口退出码 {code}；请查看运行日志。" if code else "入口结束，但未发现预期结果文件。")
        except Exception as error:
            data.update(status="failed", error=str(error), finished_at=now())
            with log_path.open("a", encoding="utf-8") as log:
                log.write("\nWORKBENCH ERROR: " + str(error) + "\n")
        atomic_json(directory / "task.json", data)


def source_directory(root, source):
    return root if source == "historical" else job_directory(root, source)


def source_work(root, source):
    if source == "historical" or task(root, source).get("runner") == "cached":
        return root / "work"
    return job_directory(root, source)


def read_table(path):
    if not path.is_file():
        return []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def read_candidates(root, source):
    base = source_directory(root, source)
    work = source_work(root, source)
    final = work / "10_select_final_candidates/final_candidates_100.csv"
    if not final.is_file():
        return read_table(base / "results/results.csv")
    candidates = read_table(final)
    risks = {row["candidate_id"]: row for row in read_table(work / "11_oral_proteolysis_filter/oral_cleavage_risk.csv")}
    for row in candidates:
        row.update({key: value for key, value in risks.get(row["candidate_id"], {}).items()
                    if key not in ("candidate_id", "sequence")})
    return candidates


def candidate_structure(root, source, row):
    """Resolve the selected candidate's PDB within result/prediction directories."""
    base = source_directory(root, source)
    work = source_work(root, source)
    cid = row["candidate_id"]
    if not re.fullmatch(r"[A-Za-z0-9_-]+", cid):
        return None
    filename = cid + "_model_0.pdb"
    allowed = [base / "results/structures", root / "results/structures",
               work / "05_prepare_structure_jobs/boltz_predictions_core",
               work / "08_prepare_optimized_structure_jobs/boltz_predictions_afdesign"]
    value = row.get("structure_pdb", "").replace("\\", "/")
    if value.startswith("project/"):
        value = value[len("project/"):]
    candidates = [root / value] if value else []
    candidates += [allowed[0] / filename, allowed[1] / filename]
    for directory in allowed[2:]:
        candidates += list(directory.glob(f"boltz_results_*/predictions/{cid}/{filename}"))
    for candidate in candidates:
        path = candidate.resolve()
        if path.name != filename or not path.is_file():
            continue
        for directory in allowed:
            try:
                path.relative_to(directory.resolve())
                return path
            except ValueError:
                pass
    return None


def summary(root, source):
    base = source_directory(root, source)
    work = source_work(root, source)
    definitions = [
        (1, "靶标准备", "01 · 序列 / 结构 / 口袋掩码", "01_prepare_target/manifest/target_manifest.json", ["target", "sequence_length"]),
        (2, "多肽生成", "02 · PepMLM 条件生成", "02_pepmlm_generate_100k/generation_manifest.json", ["outputs", "unique_sequence_count"]),
        (3, "序列质量筛选", "03 · 理化过滤 / PPL 计算", "03_sequence_filter_100k/filter_manifest.json", ["counts", "hard_filter_pass"]),
        (4, "聚类去冗余", "04 · 肽长分层 / 排序 / 聚类", "04_rank_cluster_core/selection_manifest.json", ["counts", "core_cluster_representatives"]),
        (6, "结构预测与筛选", "05–06 · Boltz / 界面规则", "06_structure_filter/structure_metrics_summary.json", ["counts", "structure_pass"]),
        (9, "优化与独立复核", "07–09 · AfDesign / 第二轮 Boltz", "07_afdesign_optimize/afdesign_manifest.json", ["counts", "sequence_filter_pass_variants"]),
        (10, "最终候选优选", "10 · 改进 / 退化保护 / 多样性", "10_select_final_candidates/final_selection_summary.json", ["counts", "final_candidates"]),
        (11, "酶切风险注释", "11 · 风险注释 / 保留原始排序", "11_oral_proteolysis_filter/oral_cleavage_manifest.json", ["candidate_count"]),
    ]
    stages = []
    for n, title, subtitle, filename, keys in definitions:
        record = load_json(work / filename, {})
        try:
            count = nested(record, keys)
        except (KeyError, TypeError):
            count = None
        if n == 2 and count is None:
            count = record.get("counts", {}).get("total_unique")
        if n == 2 and count is None:
            count = record.get("counts", {}).get("generated_candidates")
        stages.append(dict(stage=n, title=title, subtitle=subtitle, count=count,
                           record=filename, counts=record.get("counts", {})))
    selection = load_json(work / "10_select_final_candidates/final_selection_summary.json", {})
    risks = load_json(work / "11_oral_proteolysis_filter/oral_cleavage_manifest.json", {})
    return dict(stages=stages, composition=selection.get("final_by_source_type", {}),
                risk_counts=risks.get("relative_tier_counts", {}),
                result_count=len(read_candidates(root, source)),
                parameters=parameter_values(config_bundle(root) if source == "historical" else
                                            {n: load_json(base / "config_snapshot" / name) for n, name in CONFIG_NAMES.items()}))


def create_app(root=ROOT, token=None):
    from flask import Flask, jsonify, render_template, request, send_file, session
    root = Path(root).resolve()
    app = Flask(__name__, template_folder=str(ROOT / "templates"), static_folder=str(ROOT / "workbench_static"))
    app.config.update(MAX_CONTENT_LENGTH=64 * 1024, SESSION_COOKIE_HTTPONLY=True, SESSION_COOKIE_SAMESITE="Strict")
    app.secret_key = secrets.token_hex(32)
    app.config["PROJECT_ROOT"] = root

    def csrf():
        if "csrf" not in session:
            session["csrf"] = secrets.token_hex(24)
        return session["csrf"]

    def body():
        value = request.get_json()
        if not isinstance(value, dict):
            raise WorkbenchError("请求内容必须是对象。")
        return value

    @app.before_request
    def protect():
        if request.path == "/" or request.path.startswith("/static/"):
            return None
        if request.method != "GET":
            if not hmac.compare_digest(session.get("csrf", "missing"), request.headers.get("X-Workbench-CSRF", "")):
                raise WorkbenchError("页面会话已失效，请刷新页面。", 403)
            if not request.is_json:
                raise WorkbenchError("请求格式必须为 JSON。", 415)
        if token and request.path != "/api/login" and not session.get("authenticated"):
            raise WorkbenchError("请先登录工作台。", 401)

    @app.after_request
    def response_headers(response):
        response.headers["Cache-Control"] = "no-store"
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["X-Frame-Options"] = "DENY"
        response.headers["Referrer-Policy"] = "same-origin"
        return response

    @app.errorhandler(WorkbenchError)
    def handle_error(error):
        return jsonify(error=error.message), error.status

    @app.errorhandler(400)
    def bad_request(error):
        return jsonify(error="请求内容无效。"), 400

    @app.get("/")
    def index():
        return render_template("workbench.html", csrf=csrf(), needs_login=bool(token and not session.get("authenticated")))

    @app.post("/api/login")
    def login():
        value = body().get("token", "")
        if token and (not isinstance(value, str) or not hmac.compare_digest(token, value)):
            raise WorkbenchError("访问口令不正确。", 401)
        session["authenticated"] = True
        return jsonify(ok=True)

    @app.get("/api/config")
    def get_config():
        revision = config_revision(root)
        values = parameter_values(config_bundle(root))
        if revision != config_revision(root):
            raise WorkbenchError("参数正在保存，请重新载入。", 409)
        return jsonify(fields=FIELDS, values=values, revision=revision)

    @app.post("/api/config")
    def save_config():
        payload = body()
        with state_lock(root):
            assert_idle(root)
            if payload.get("revision") != config_revision(root):
                raise WorkbenchError("配置已被其他操作修改，请重新载入参数。", 409)
            before = config_bundle(root)
            after = updated_configs(before, payload.get("values"))
            changed = [n for n in before if before[n] != after[n]]
            backup = root / "work" / "workbench" / "config_backups" / (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S") + "-" + secrets.token_hex(4))
            if changed:
                backup.mkdir(parents=True)
                for n in changed:
                    shutil.copy2(root / "config" / CONFIG_NAMES[n], backup / CONFIG_NAMES[n])
                try:
                    for n in changed:
                        atomic_json(root / "config" / CONFIG_NAMES[n], after[n])
                except Exception:
                    for n in changed:
                        shutil.copy2(backup / CONFIG_NAMES[n], root / "config" / CONFIG_NAMES[n])
                    raise
                atomic_json(backup / "changes.json", dict(saved_at=now(), changed_files=[CONFIG_NAMES[n] for n in changed]))
            return jsonify(values=parameter_values(after), revision=config_revision(root), changed_files=[CONFIG_NAMES[n] for n in changed],
                           backup=backup.relative_to(root).as_posix() if changed else None)

    @app.get("/api/state")
    def get_state():
        source = request.args.get("source", "historical")
        job = None if source == "historical" else task(root, source)
        return jsonify(source=source, task=job, tasks=tasks(root), environment=environment(root), **summary(root, source))

    @app.post("/api/tasks")
    def launch():
        return jsonify(start_task(root, body())), 202

    @app.get("/api/results")
    def results():
        return jsonify(rows=read_candidates(root, request.args.get("source", "historical")))

    @app.get("/api/log")
    def log():
        directory = job_directory(root, request.args.get("source", ""))
        path = directory / "launcher.log"
        if not path.exists():
            path = directory / "worker.log"
        if not path.exists():
            return jsonify(text="等待任务启动…", truncated=False)
        with path.open("rb") as handle:
            size = path.stat().st_size
            handle.seek(max(0, size - 65536))
            content = handle.read().decode("utf-8", errors="replace")
        return jsonify(text=content, truncated=size > 65536)

    @app.get("/api/download/<kind>")
    def download(kind):
        source = request.args.get("source", "historical")
        base = source_directory(root, source)
        mapping = {"results": base / "results/results.csv", "plan": base / "execution_plan.json", "log": base / "launcher.log"}
        if kind == "results":
            candidates = read_candidates(root, source)
            if not candidates:
                raise WorkbenchError("候选表尚未生成。", 404)
            text = io.StringIO(newline="")
            writer = csv.DictWriter(text, fieldnames=list(candidates[0]))
            writer.writeheader()
            writer.writerows(candidates)
            return send_file(io.BytesIO(text.getvalue().encode("utf-8-sig")), as_attachment=True,
                             download_name="final_candidates_100.csv", mimetype="text/csv")
        if kind == "structure":
            candidate_id = request.args.get("candidate", "")
            rows = read_candidates(root, source)
            row = next((r for r in rows if r["candidate_id"] == candidate_id), None)
            if not row:
                raise WorkbenchError("候选不存在。", 404)
            path = candidate_structure(root, source, row)
        else:
            path = mapping.get(kind)
        if path is None or not path.is_file():
            raise WorkbenchError("文件尚未生成。", 404)
        return send_file(path, as_attachment=True)

    return app


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8080)
    parser.add_argument("--root", type=Path, default=ROOT, help=argparse.SUPPRESS)
    parser.add_argument("--worker", help=argparse.SUPPRESS)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.worker:
        run_worker(root, args.worker)
        return
    token = os.environ.get("WORKBENCH_TOKEN")
    if args.host not in {"127.0.0.1", "localhost", "::1"} and not token:
        parser.error("对外监听时请设置 WORKBENCH_TOKEN 访问口令。")
    try:
        from waitress import serve
        app = create_app(root, token)
    except ImportError:
        parser.error("请先安装依赖：python -m pip install -r environment/requirements-workbench.txt")
    print(f"NTSR1 workbench: http://{args.host}:{args.port}", flush=True)
    serve(app, host=args.host, port=args.port, threads=8)


if __name__ == "__main__":
    main()
