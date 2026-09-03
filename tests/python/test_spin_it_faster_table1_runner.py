from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "experiments/spin_it_reproduction/faster/run_paper_table1.py"


def load_runner():
    spec = importlib.util.spec_from_file_location("spin_it_faster_table1", RUNNER)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load Faster Table 1 runner")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SpinItFasterTable1RunnerTests(unittest.TestCase):
    def test_initial_points_are_deterministic_and_mode_constrained(self) -> None:
        runner = load_runner()
        first = runner.generate_initial_quadrics("YO_YO", seed=73, count=50)
        second = runner.generate_initial_quadrics("YO_YO", seed=73, count=50)
        self.assertEqual(first, second)
        self.assertEqual(len(first), 50)
        self.assertEqual(len(set(tuple(row) for row in first)), 50)
        for row in first:
            self.assertAlmostEqual(math.sqrt(sum(value * value for value in row)), 1.0)
            self.assertEqual(row[0], 0.0)
            self.assertLess(row[6], 0.0)
        for row in runner.generate_initial_quadrics("BUOYANT_STABILITY", 11, 4):
            self.assertEqual(row[4:], [0.0] * 6)

    def test_domain_aware_planes_intersect_the_design_bounds(self) -> None:
        runner = load_runner()
        corners = tuple(
            (x, y, z)
            for x in (-2.0, 3.0)
            for y in (-1.0, 4.0)
            for z in (-0.5, 0.75)
        )
        for mode in ("STATIC_BALANCE", "BUOYANT_STABILITY"):
            for row in runner.generate_initial_quadrics(mode, 91, 12, corners):
                values = [
                    row[0] + row[1] * x + row[2] * y + row[3] * z
                    for x, y, z in corners
                ]
                self.assertLess(min(values), 0.0)
                self.assertGreater(max(values), 0.0)

    def test_summary_uses_paper_statistics_and_hash_excludes_timings(self) -> None:
        runner = load_runner()
        attempts = [
            {"run": 0, "restart": 1, "status": "FAILED", "newton_iterations": 3,
             "function_evaluations": 8, "residual_seconds": 0.2, "derivative_seconds": 0.1},
            {"run": 0, "restart": 2, "status": "CONVERGED", "newton_iterations": 5,
             "function_evaluations": 9, "residual_seconds": 0.3, "derivative_seconds": 0.2},
            {"run": 1, "restart": 1, "status": "CONVERGED", "newton_iterations": 7,
             "function_evaluations": 11, "residual_seconds": 0.4, "derivative_seconds": 0.25},
        ]
        summary = runner.summarize_attempts(attempts, expected_runs=2)
        self.assertTrue(summary["all_runs_converged"])
        self.assertEqual(summary["average_restarts"], 1.5)
        self.assertEqual(summary["median_newton_iterations"], 6.0)
        self.assertEqual(summary["median_function_evaluations"], 10.0)
        first_hash = runner.scientific_hash({"status": "CONVERGED", "elapsed_seconds": 1.0})
        second_hash = runner.scientific_hash({"status": "CONVERGED", "elapsed_seconds": 9.0})
        self.assertEqual(first_hash, second_hash)

    def test_two_start_plan_records_all_reconstructed_cases(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "plan.json"
            subprocess.run(
                ["python3", str(RUNNER), "--plan-only", "--starts", "2", "--output", str(output)],
                check=True,
            )
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(payload["schema"], "spin-it-faster-table1-plan-v1")
            self.assertEqual(payload["starts_per_case"], 2)
            self.assertEqual(len(payload["cases"]), 8)
            self.assertTrue(all(len(case["initial_quadrics"]) == 2 for case in payload["cases"]))


class SpinItFasterTable1NativeRunnerTests(unittest.TestCase):
    def test_two_start_miniature_run_converges_within_paper_cutoff(self) -> None:
        runner = load_runner()
        manifest = json.loads(runner.DEFAULT_CASES.read_text(encoding="utf-8"))
        case = next(row for row in manifest["cases"] if row["id"] == "static-tilted-vase")
        payload = runner.run_case(case, starts=2, max_restarts=1)
        self.assertTrue(payload["summary"]["all_runs_converged"])
        self.assertEqual(payload["summary"]["successful_runs"], 2)
        self.assertEqual(payload["summary"]["average_restarts"], 1.0)
        self.assertEqual(len(payload["attempts"]), 2)
        for attempt in payload["attempts"]:
            self.assertLessEqual(attempt["function_evaluations"], 60)
            self.assertGreater(attempt["residual_seconds"], 0.0)
            self.assertGreater(attempt["derivative_seconds"], 0.0)
            self.assertTrue(attempt["trace"])


if __name__ == "__main__":
    unittest.main()
