"""Run with python -m unittest discover -s tests -v (Flask required)."""
import importlib.util
import json
import re
import shutil
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))
spec = importlib.util.spec_from_file_location("workbench", PROJECT / "workbench.py")
wb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wb)


class WorkbenchTests(unittest.TestCase):
    def setUp(self):
        (PROJECT / "work").mkdir(exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(prefix="workbench_test_", dir=PROJECT / "work")
        self.root = Path(self.temp.name)
        shutil.copytree(PROJECT / "config", self.root / "config")
        shutil.copytree(PROJECT / "pipeline", self.root / "pipeline", ignore=shutil.ignore_patterns("__pycache__"))
        shutil.copy2(PROJECT / "run_pipeline.sh", self.root / "run_pipeline.sh")
        self.app = wb.create_app(self.root)
        self.client = self.app.test_client()
        html = self.client.get("/").get_data(as_text=True)
        self.csrf = re.search(r'name="workbench-csrf" content="([^"]+)"', html)[1]
        self.headers = {"X-Workbench-CSRF": self.csrf}

    def tearDown(self):
        for t in wb.tasks(self.root):
            if t["status"] in wb.BUSY:
                self.wait_job(t["id"])
        self.root.resolve().relative_to((PROJECT / "work").resolve())
        self.temp.cleanup()

    def config(self):
        return self.client.get("/api/config").get_json()

    def post(self, url, payload):
        return self.client.post(url, json=payload, headers=self.headers)

    def wait_job(self, job_id):
        for _ in range(150):
            data = wb.task(self.root, job_id)
            if data["status"] not in wb.BUSY and not wb.pid_alive(data.get("worker_pid")):
                return data
            time.sleep(.1)
        self.fail("worker did not finish")

    def test_save_validated_parameters_backs_up_only_changed_configs(self):
        old = {n: (self.root / "config" / name).read_bytes() for n, name in wb.CONFIG_NAMES.items()}
        c = self.config()
        response = self.post("/api/config", dict(revision=c["revision"], values={"total_candidates": 2000, "length_min": 14, "length_max": 20, "structure_min_pair_iptm_mean": .9, "final_count": 25}))
        self.assertEqual(response.status_code, 200)
        data = response.get_json()
        self.assertEqual(set(data["changed_files"]), {wb.CONFIG_NAMES[n] for n in [2, 6, 10]})
        self.assertEqual(wb.load_json(self.root / "config" / wb.CONFIG_NAMES[2])["lengths"], list(range(14, 21)))
        for n, before in old.items():
            if n in [2, 6, 10]:
                self.assertEqual((self.root / data["backup"] / wb.CONFIG_NAMES[n]).read_bytes(), before)
            else:
                self.assertEqual((self.root / "config" / wb.CONFIG_NAMES[n]).read_bytes(), before)
        self.assertEqual(wb.load_json(self.root / "config" / wb.CONFIG_NAMES[10])["absolute_structure_thresholds"]["min_pair_iptm_mean"], .88)

    def test_invalid_and_stale_edits_do_not_change_files(self):
        c = self.config()
        for changes in [{"length_min": 32, "length_max": 12}, {"final_count": 1000001}, {"total_candidates": True}, {"structure_min_pair_iptm_mean": -1}, {"command": "rm"}, {"qc_min_net_charge": 10, "qc_max_net_charge": 5}, {"ppl_batch": 2.5}]:
            self.assertEqual(self.post("/api/config", dict(revision=c["revision"], values=changes)).status_code, 400)
            self.assertEqual(wb.config_revision(self.root), c["revision"])
        self.assertEqual(self.post("/api/config", dict(revision="stale", values={"final_count": 10})).status_code, 409)

    def test_plan_worker_uses_cached_runner_and_saved_settings_without_main(self):
        c = self.config()
        edited = self.post("/api/config", dict(revision=c["revision"], values={"total_candidates": 2048, "final_count": 25})).get_json()
        before = wb.config_revision(self.root)
        response = self.post("/api/tasks", dict(name="真实入口计划测试", mode="plan", revision=edited["revision"]))
        self.assertEqual(response.status_code, 202)
        job = self.wait_job(response.get_json()["id"])
        self.assertEqual(job["status"], "planned", job)
        directory = wb.job_directory(self.root, job["id"])
        plan = wb.load_json(directory / "execution_plan.json")
        self.assertEqual(len(plan["stages"]), 11)
        self.assertEqual(plan["configs"]["2"]["total_candidates"], 2048)
        self.assertEqual(plan["configs"]["10"]["final_count"], 25)
        self.assertEqual([a["action"] for a in plan["actions"]], ["run"] * 11)
        self.assertFalse((self.root / "main.py").exists())
        self.assertEqual(wb.config_revision(self.root), before)
        with self.client.get("/api/download/plan?source=" + job["id"]) as response:
            self.assertEqual(response.status_code, 200)
        self.assertIn("no computation started", self.client.get("/api/log?source=" + job["id"]).get_json()["text"])

    def test_concurrent_jobs_and_saving_are_blocked_but_reads_work(self):
        (self.root / "pipeline/cached_runner.py").write_text("import time\ntime.sleep(2)\nraise SystemExit(7)\n", encoding="utf-8")
        c = self.config()
        response = self.post("/api/tasks", dict(name="单任务锁测试", mode="plan", revision=c["revision"]))
        self.assertEqual(response.status_code, 202)
        job_id = response.get_json()["id"]
        self.assertEqual(self.post("/api/tasks", dict(name="第二个任务", mode="plan", revision=c["revision"])).status_code, 409)
        self.assertEqual(self.post("/api/config", dict(revision=c["revision"], values={"final_count": 10})).status_code, 409)
        self.assertEqual(self.client.get("/api/config").status_code, 200)
        job = self.wait_job(job_id)
        self.assertEqual(job["status"], "failed")
        self.assertEqual(job["returncode"], 7)
        self.assertEqual(wb.config_revision(self.root), c["revision"])

    def test_csrf_authentication_and_path_restrictions(self):
        c = self.config()
        self.assertEqual(self.client.post("/api/tasks", json=dict(name="x", mode="plan", revision=c["revision"])).status_code, 403)
        for url in ["/api/state?source=../../config", "/api/log?source=../../config", "/api/download/results?source=../../config"]:
            self.assertEqual(self.client.get(url).status_code, 404)
        app = wb.create_app(self.root, token="test-access")
        client = app.test_client()
        html = client.get("/").get_data(as_text=True)
        headers = {"X-Workbench-CSRF": re.search(r'name="workbench-csrf" content="([^"]+)"', html)[1]}
        self.assertEqual(client.get("/api/config").status_code, 401)
        self.assertEqual(client.post("/api/login", json={"token": "wrong"}, headers=headers).status_code, 401)
        self.assertEqual(client.post("/api/login", json={"token": "test-access"}, headers=headers).status_code, 200)
        self.assertEqual(client.get("/api/config").status_code, 200)
        self.assertEqual(self.post("/api/config", []).status_code, 400)


    def test_restart_can_read_persisted_task(self):
        c = self.config()
        response = self.post("/api/tasks", dict(name="重启后读取", mode="plan", revision=c["revision"]))
        job = self.wait_job(response.get_json()["id"])
        new_client = wb.create_app(self.root).test_client()
        data = new_client.get("/api/state?source=" + job["id"]).get_json()
        self.assertEqual(data["task"]["status"], "planned")
        self.assertEqual(data["result_count"], 0)

    def test_full_mode_checks_environment_then_records_success_and_progress(self):
        c = self.config()
        with patch.object(wb, "environment", return_value={"full_available": False}):
            self.assertEqual(self.post("/api/tasks", dict(name="环境检查", mode="full", revision=c["revision"])).status_code, 400)
        # Exercise orchestration using a fixture driver; no model inference occurs.
        (self.root / "pipeline/cached_runner.py").write_text(
            "from pathlib import Path\n"
            "base=Path('work')\n"
            "print('[5/11] 初始结构任务准备', flush=True)\n"
            "print('    Boltz 结构预测运行', flush=True)\n"
            "print('[11/11] 酶切风险注释', flush=True)\n"
            "(base/'10_select_final_candidates').mkdir(parents=True,exist_ok=True)\n"
            "(base/'11_oral_proteolysis_filter').mkdir(parents=True,exist_ok=True)\n"
            "(base/'10_select_final_candidates/final_candidates_100.csv').write_text('candidate_id,sequence,final_rank\\nTEST,AAAA,1\\n')\n"
            "(base/'11_oral_proteolysis_filter/oral_cleavage_risk.csv').write_text('candidate_id,oral_cleavage_risk_level\\nTEST,medium\\n')\n",
            encoding="utf-8")
        with patch.object(wb, "environment", return_value={"full_available": True}):
            response = self.post("/api/tasks", dict(name="完整任务状态测试", mode="full", revision=c["revision"]))
        self.assertEqual(response.status_code, 202)
        job = self.wait_job(response.get_json()["id"])
        self.assertEqual(job["status"], "complete", job)
        self.assertEqual(job["active_stage"], 11)
        self.assertEqual(job["returncode"], 0)
        self.assertNotIn("--plan-only", job["command"])
        self.assertEqual(self.client.get("/api/results?source=" + job["id"]).get_json()["rows"][0]["candidate_id"], "TEST")





if __name__ == "__main__":
    unittest.main()
