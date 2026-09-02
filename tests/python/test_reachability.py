from __future__ import annotations

from types import SimpleNamespace
import unittest

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


if __name__ == "__main__":
    unittest.main()
