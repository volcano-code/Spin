"""Fail-closed contracts for paper-simulation experiment entry points."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def _load_driver(name: str):
    path = REPOSITORY_ROOT / "experiments" / "paper_sim_v1" / f"{name}.py"
    specification = importlib.util.spec_from_file_location(
        f"paper_sim_test_{name}", path
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load experiment driver: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


class ExperimentDriverTests(unittest.TestCase):
    def test_official_case_flattens_the_surface_that_was_optimized(self):
        import numpy as np

        driver = _load_driver("run_official_case")
        source = np.zeros((9, 3), dtype=np.float64)
        replacement = np.ones((4, 3), dtype=np.float64)
        archive = SimpleNamespace(
            source_vertices=source,
            normalized_vertices=replacement,
        )
        self.assertIs(driver._flattening_reference_vertices(archive), replacement)

    def test_static_oracle_payload_normalizes_numpy_boolean(self):
        from dataclasses import dataclass
        import numpy as np

        driver = _load_driver("run_static_oracle")

        @dataclass
        class Result:
            inside_support: object
            margin: float

        payload = driver._serializable_result(Result(np.bool_(True), 0.5))
        self.assertIs(type(payload["inside_support"]), bool)
        json.dumps(payload)

    def test_static_oracle_uses_release_safe_region_tolerance(self):
        driver = _load_driver("run_static_oracle")
        results = {
            angle: [{"signed_margin": -5.0e-7}]
            for angle in ("0.0", "2.0", "5.0")
        }
        self.assertTrue(driver._accepted_required_angles(results, 1.0e-6))
        self.assertFalse(driver._accepted_required_angles(results, 1.0e-7))

    def test_lambda_candidate_requires_post_seal_safe_region(self):
        driver = _load_driver("run_lambda_search")
        run = {
            "completion": "stable",
            "sealed_artifact_status": "SEALED",
        }
        self.assertFalse(driver._candidate_is_eligible(0, run))
        run["sealed_safe_region_state"] = "OUTSIDE"
        self.assertFalse(driver._candidate_is_eligible(0, run))
        run["sealed_safe_region_state"] = "SAFE"
        self.assertTrue(driver._candidate_is_eligible(0, run))

    def test_lambda_search_writes_valid_boolean_when_no_candidate_succeeds(self):
        driver = _load_driver("run_lambda_search")
        failed = {
            "lambda": 5.0,
            "eligible": False,
            "process_returncode": 2,
            "completion": "not_converged",
            "sealed_artifact_status": "REJECTED",
            "normalized_rms_vertex_displacement": None,
            "final_projected_target_energy": None,
            "result_directory": "/unused",
        }
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "search"
            argv = [
                "run_lambda_search.py",
                "--beta-root", directory,
                "--case", "spheres",
                "--output", str(output),
            ]
            with mock.patch.object(sys, "argv", argv), mock.patch.object(
                driver, "_run_candidate", return_value=failed
            ):
                self.assertEqual(driver.main(), 2)
            report = json.loads(
                (output / "lambda-search.json").read_text(encoding="utf-8")
            )
            self.assertIs(report["full_resolution_only"], True)

    def test_toy_handle_reconstruction_is_deterministic_and_avoids_support(self):
        driver = _load_driver("run_toy_reconstruction")
        vertices = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ]
        first = driver.select_reconstructed_handles(
            vertices,
            gravity=(0.0, -1.0, 0.0),
            support_vertex_indices=(0, 1),
            count=3,
        )
        second = driver.select_reconstructed_handles(
            vertices,
            gravity=(0.0, -1.0, 0.0),
            support_vertex_indices=(0, 1),
            count=3,
        )
        self.assertEqual(first, second)
        self.assertEqual(len(first), 3)
        self.assertTrue(all(index not in {0, 1} for index in first))

    def test_table1_plan_runs_required_cases_in_diagnostic_order(self):
        driver = _load_driver("run_table1_matrix")
        plan = driver.build_table_plan(
            beta_root=Path("/archive"),
            output=Path("/results"),
            profile="paper",
            lambdas={
                "spheres": 5.0,
                "dinosaur": 10.0,
                "armadillo1": 20.0,
                "horse": 40.0,
                "toyrocking": 80.0,
            },
            repository_root=REPOSITORY_ROOT,
        )
        self.assertEqual(
            [item["case"] for item in plan],
            ["spheres", "dinosaur", "armadillo1", "horse", "toyrocking"],
        )
        horse = next(item for item in plan if item["case"] == "horse")
        self.assertIn("--replacement-normalized-surface", horse["command"])
        toy = next(item for item in plan if item["case"] == "toyrocking")
        self.assertTrue(toy["driver"].endswith("run_toy_reconstruction.py"))


if __name__ == "__main__":
    unittest.main()
