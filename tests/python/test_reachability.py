from __future__ import annotations

import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import make_it_stand.reachability as reachability
from make_it_stand.reachability import ReachabilityResult


class ReachabilitySchemaTests(unittest.TestCase):
    def test_native_result_has_strict_deterministic_schema(self) -> None:
        native = SimpleNamespace(
            input_hash="a" * 64,
            status="RELAXATION_REACHABLE",
            input_cell_count=2,
            variable_cell_count=2,
            fixed_mass=1.0,
            retained_fractions=[0.25, 0.75],
            final_mass=2.0,
            final_first_moment_2d=[0.5, -0.25],
            final_center_2d=[0.25, -0.125],
            maximum_primal_violation=0.0,
            primal_solver_status="Optimal",
            certificate_solver_status="",
            certificate_system_nonzero_count=6,
            certificate=SimpleNamespace(
                weights=[],
                minimum_weight=0.0,
                normalization_residual=0.0,
                stationarity_infinity_norm=0.0,
                contradiction=0.0,
                valid=False,
            ),
        )
        result = ReachabilityResult.from_native(native)
        payload = result.to_dict()
        self.assertEqual(
            payload["schema"],
            "make-it-stand-carving-reachability-v1",
        )
        self.assertIn(
            payload["status"],
            {
                "ALREADY_SAFE",
                "PROVEN_UNREACHABLE",
                "RELAXATION_REACHABLE",
                "INDETERMINATE",
            },
        )
        self.assertRegex(payload["input_hash"], r"^[0-9a-f]{64}$")
        self.assertEqual(payload["certificate_system_nonzero_count"], 6)
        self.assertEqual(payload, result.to_dict())

    def test_schema_rejects_nonfinite_native_number(self) -> None:
        native = SimpleNamespace(
            input_hash="b" * 64,
            status="INDETERMINATE",
            input_cell_count=0,
            variable_cell_count=0,
            fixed_mass=float("nan"),
            retained_fractions=[],
            final_mass=0.0,
            final_first_moment_2d=[0.0, 0.0],
            final_center_2d=[0.0, 0.0],
            maximum_primal_violation=0.0,
            primal_solver_status="API_ERROR",
            certificate_solver_status="",
            certificate_system_nonzero_count=0,
            certificate=SimpleNamespace(
                weights=[],
                minimum_weight=0.0,
                normalization_residual=0.0,
                stationarity_infinity_norm=0.0,
                contradiction=0.0,
                valid=False,
            ),
        )
        with self.assertRaises(ValueError):
            ReachabilityResult.from_native(native)

    def test_manifest_runner_is_deterministic_and_fail_closed(self) -> None:
        manifest = {
            "schema": "make-it-stand-reachability-cases-v1",
            "cases": [
                {"id": "zeta", "observed_optimizer_outcome": "FAILED"},
                {"id": "alpha", "observed_optimizer_outcome": "SUCCESS"},
            ],
        }

        def evaluator(case: dict[str, object]) -> dict[str, object]:
            status = (
                "RELAXATION_REACHABLE"
                if case["id"] == "zeta"
                else "ALREADY_SAFE"
            )
            return {
                "schema": reachability.SCHEMA,
                "schema_version": reachability.SCHEMA_VERSION,
                "input_hash": "c" * 64,
                "status": status,
                "input_cell_count": 2,
                "variable_cell_count": 2,
                "fixed_mass": 1.0,
                "retained_fractions": [0.5, 0.5],
                "final_mass": 2.0,
                "final_first_moment_2d": [0.0, 0.0],
                "final_center_2d": [0.0, 0.0],
                "maximum_primal_violation": 0.0,
                "primal_solver_status": "Optimal",
                "certificate_solver_status": "",
                "certificate_system_nonzero_count": 6,
                "certificate": {
                    "weights": [],
                    "minimum_weight": 0.0,
                    "normalization_residual": 0.0,
                    "stationarity_infinity_norm": 0.0,
                    "contradiction": 0.0,
                    "valid": False,
                },
            }

        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "first"
            second = Path(directory) / "second"
            payload_a = reachability.run_manifest(manifest, first, evaluator)
            payload_b = reachability.run_manifest(manifest, second, evaluator)
            self.assertEqual(payload_a, payload_b)
            self.assertEqual(payload_a["case_order"], ["alpha", "zeta"])
            self.assertEqual(
                payload_a["cases"][1]["classification"],
                "CURRENT_PLANE_CARVER_SOLVER_GAP",
            )
            self.assertFalse(payload_a["release_gate_passed"])
            compact_certificate = payload_a["cases"][0]["result"]["certificate"]
            self.assertNotIn("weights", compact_certificate)
            self.assertEqual(compact_certificate["weights_count"], 0)
            self.assertRegex(
                compact_certificate["weights_sha256"],
                r"^[0-9a-f]{64}$",
            )
            self.assertNotIn(
                "PERFECT_REPRODUCTION",
                json.dumps(payload_a, sort_keys=True),
            )
            self.assertEqual(
                (first / "summary.json").read_bytes(),
                (second / "summary.json").read_bytes(),
            )


if __name__ == "__main__":
    unittest.main()
