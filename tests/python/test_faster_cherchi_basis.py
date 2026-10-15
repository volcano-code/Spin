"""Actual reference-Boolean finite differences, not synthetic Jacobian stubs."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT/"experiments/spin_it_reproduction/strict/run_faster_cherchi_differential_basis.py"
FIXTURE = ROOT/"experiments/spin_it_reproduction/faster/fixtures/cherchi-full-basis-method-v1/problem.json"
PROBE = ROOT/"build/cherchi-reference-singlethread-20260917/cherchi_probe"


class CherchiBasisTests(unittest.TestCase):
    def test_actual_paper_basis_rejects_missing_or_forged_source_provenance(self):
        runner = self.module()
        fixture = ROOT/"experiments/spin_it_reproduction/faster/fixtures/cherchi-paper-projective-basis-method-v1/problem.json"
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)/"output"
            result = runner.run(output, fixture=fixture, probe=PROBE)
            self.assertEqual(result["status"], "PASSED")
            self.assertEqual(result["evaluation_count"], 61)
            self.assertEqual(len(result["checks"]), 10)
            self.assertFalse(result["optimizer_executed"])
            records = [json.loads(line) for line in (output/"evaluations.jsonl").read_text().splitlines()]
            for defect in ("missing_features", "modified_features", "missing_map_identity"):
                bad = copy.deepcopy(records)
                for row in bad:
                    if defect == "missing_features":
                        row.pop("paper_source_features")
                    elif defect == "modified_features":
                        row["paper_source_features"]["normal_vertices"][0] += 1
                    else:
                        row["branch"].pop("paper_source_features_sha256")
                        row.pop("paper_source_features")
                with self.subTest(defect=defect), self.assertRaises(ValueError):
                    runner.check_records(result["protocol"], bad)

    def test_explicit_paper_mapping_replays_frozen_native_geometry(self):
        runner = self.module()
        fixture = ROOT/"experiments/spin_it_reproduction/faster/fixtures/cherchi-paper-projective-basis-method-v1/problem.json"
        problem = json.loads(fixture.read_text())
        native = runner.kernel.spin._native_module()
        source = native.evaluate_spin_it_faster_source_motion(
            problem["coefficients"], problem["design_vertices"], runner._settings(native, problem))
        self.assertTrue(source["defined"], source["diagnostic"])
        # Captured from the independent native dump before this adapter existed.
        self.assertEqual(source["triangulation_hash"],
                         "b038b65b251110b249512af5807417ab4fc5f987298b500079f908c147835b05")
        self.assertTrue(source["paper_projective_hyperboloid"])
        self.assertEqual(len(source["vertices"]), 70)
        self.assertEqual(len(source["faces"]), 132)
        self.assertEqual(source["truncation_radius"], 221.70250336881628)

    def test_invalid_paper_mapping_requests_are_not_silently_legacy(self):
        runner = self.module()
        native = runner.kernel.spin._native_module()
        problem = json.loads(FIXTURE.read_text())
        for bad in (1, "true", None):
            with self.subTest(flag=bad), self.assertRaises(ValueError):
                runner._settings(native, dict(problem, paper_projective_hyperboloid=bad))
        for radius in (True, -1.0, float("nan"), float("inf"), "2"):
            with self.subTest(radius=radius), self.assertRaises(ValueError):
                runner._settings(native, dict(problem, truncation_radius=radius))
        paraboloid = json.loads((FIXTURE.parent.parent/"cherchi-paraboloid-basis-method-v1/problem.json").read_text())
        with self.assertRaises(ValueError):
            runner._settings(native, dict(paraboloid, paper_projective_hyperboloid=True))

    def module(self):
        self.assertTrue(RUNNER.is_file(), "complete actual Cherchi basis checker is missing")
        sys.path.insert(0, str(RUNNER.parent))
        spec = importlib.util.spec_from_file_location("cherchi_basis_test_runner", RUNNER)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def test_timeout_keeps_partial_boolean_streams_and_receipt(self):
        runner = self.module()
        invoke = getattr(runner, "invoke_boolean", None)
        self.assertTrue(callable(invoke), "timeout-safe Boolean receipt adapter missing")
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            command = [sys.executable, "-u", "-c",
                       "import sys,time; print('partial-out'); print('partial-err',file=sys.stderr); time.sleep(5)"]
            with self.assertRaises(subprocess.TimeoutExpired):
                invoke(command, directory, "timeout", timeout=.5)
            self.assertEqual((directory/"boolean-timeout.stdout").read_text(), "partial-out\n")
            self.assertEqual((directory/"boolean-timeout.stderr").read_text(), "partial-err\n")
            receipt = json.loads((directory/"boolean-timeout.json").read_text())
            self.assertEqual(receipt["exception_type"], "TimeoutExpired")
            self.assertEqual(receipt["status"], "TIMEOUT")
            self.assertIsNone(receipt["returncode"])

    def test_comparison_failure_is_finalizable_and_not_a_pass(self):
        runner = self.module()
        capture = getattr(runner, "capture_comparison", None)
        self.assertTrue(callable(capture), "failed comparison finalization adapter missing")
        report = capture({"samples": []}, [])
        self.assertEqual(report["status"], "FAILED_COMPARISON")
        self.assertEqual(report["exception_type"], "ValueError")
        self.assertIn("missing or extra", report["diagnostic"])

    def test_actual_full_tangent_basis_and_radial_null_are_checked(self):
        runner = self.module()
        with tempfile.TemporaryDirectory() as directory:
            record = runner.run(Path(directory)/"output", fixture=FIXTURE, probe=PROBE)
            self.assertEqual(record["status"], "PASSED", json.dumps(record.get("checks")))
            self.assertFalse(record["optimizer_executed"])
            self.assertEqual(record["evaluation_count"], 61)
            self.assertEqual(len(record["checks"]), 10)
            for row in record["checks"]:
                self.assertEqual(len(row["steps"]), 3)
                self.assertTrue(all(step["same_branch"] for step in row["steps"]))
            basis = np.asarray(record["parameter_basis"])
            q = np.asarray(record["normalized_coefficients"])
            self.assertEqual(basis.shape, (10, 9))
            np.testing.assert_allclose(basis.T@basis, np.eye(9), atol=1e-14)
            np.testing.assert_allclose(basis@basis.T, np.eye(10)-np.outer(q,q), atol=1e-14)
            evaluations = [json.loads(line) for line in
                           (Path(directory)/"output/evaluations.jsonl").read_text().splitlines()]
            self.assertEqual(len(evaluations), 61)
            self.assertTrue(all(row["independent_moment_pass"] for row in evaluations))
            self.assertIn("LPI", evaluations[0]["provenance_types"])
            self.assertEqual(record["scope"], runner.SCOPE)
            # Relabelling every unchanged legacy source as the new point map
            # must be rejected even though its numerical differences agree.
            relabelled=copy.deepcopy(record["protocol"])
            relabelled["finite_point_map"]="ALGEBRAIC_V1"
            with self.assertRaisesRegex(ValueError,"finite.*map|algebraic"):
                runner.check_records(relabelled,evaluations)
            for defect in ("missing", "wrong_jacobian", "changed_branch", "nonfinite", "wrong_coefficients"):
                bad = copy.deepcopy(evaluations)
                if defect == "missing": bad.pop()
                elif defect == "wrong_jacobian": bad[0]["jacobian"][0][0] += .5
                elif defect == "changed_branch":
                    old = bad[1]["branch"]["operation"]
                    bad[1]["branch"]["operation"] = "intersection" if old == "subtraction" else "subtraction"
                elif defect == "nonfinite": bad[1]["raw_moments"][0] = float("nan")
                else: bad[1]["coefficients"][0] += .01
                with self.subTest(defect=defect):
                    try:
                        checked = runner.check_records(record["protocol"], bad)
                    except ValueError:
                        continue
                    self.assertNotEqual(checked["status"], "PASSED")


if __name__ == "__main__":
    unittest.main()
