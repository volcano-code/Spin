"""Independent value/derivative checks for the exposed Spin-It 2014 Givens form."""

from __future__ import annotations

import importlib
import math
import unittest

import numpy as np

from make_it_stand.verification.spin2014_givens_reference import evaluate_paper_givens


RAW = np.asarray((4.0, 0.4, -0.8, 2.0, 10.0, 14.0, 6.0, 0.3, -0.4, 0.2))
GAMMA_I, GAMMA_C, PHI = 2.5, 0.7, 0.37


class Spin2014GivensOracleTests(unittest.TestCase):
    def setUp(self) -> None:
        self.native = importlib.import_module("_make_it_stand")

    def test_every_objective_matches_independent_paper_formula(self) -> None:
        self.assertTrue(hasattr(self.native, "evaluate_spin_it_2014_givens"))
        for objective in ("YOYO", "TOP_INERTIA_ONLY", "TOP_LOWERING_ONLY",
                          "TOP_MASS_AND_LOWERING", "STATIC_LOWERING"):
            with self.subTest(objective=objective):
                phi = 0.0 if objective == "STATIC_LOWERING" else PHI
                actual = self.native.evaluate_spin_it_2014_givens(
                    getattr(self.native.Spin2014Objective, objective), RAW,
                    GAMMA_I, GAMMA_C, phi,
                )
                expected = evaluate_paper_givens(objective, RAW, GAMMA_I, GAMMA_C, phi)
                expected_value, expected_constraints = expected["objective"], expected["constraints"]
                self.assertEqual(actual["objective_token"], objective)
                np.testing.assert_allclose(actual["objective"], expected_value, rtol=0.0, atol=1e-13)
                np.testing.assert_allclose(actual["constraints"], expected_constraints, rtol=0.0, atol=1e-13)
                self.assertEqual(np.asarray(actual["objective_gradient"]).shape, (11,))
                self.assertEqual(np.asarray(actual["objective_hessian"]).shape, (11, 11))
                self.assertEqual(np.asarray(actual["constraint_jacobian"]).shape,
                                 (len(expected_constraints), 11))

    def test_gradient_and_hessian_are_consistent_at_three_steps(self) -> None:
        for objective in ("YOYO", "TOP_INERTIA_ONLY", "TOP_LOWERING_ONLY", "TOP_MASS_AND_LOWERING"):
            with self.subTest(objective=objective):
                state = np.concatenate((RAW, (PHI,)))
                evaluate = lambda x: self.native.evaluate_spin_it_2014_givens(
                    getattr(self.native.Spin2014Objective, objective), x[:10], GAMMA_I, GAMMA_C, x[10]
                )
                analytic = evaluate(state)
                gradient, hessian = np.asarray(analytic["objective_gradient"]), np.asarray(analytic["objective_hessian"])
                for h in (1e-4, 3e-5, 1e-5):
                    finite_gradient = np.zeros(11)
                    finite_hessian = np.zeros((11, 11))
                    for j in range(11):
                        plus, minus = state.copy(), state.copy()
                        plus[j] += h; minus[j] -= h
                        p, m = evaluate(plus), evaluate(minus)
                        finite_gradient[j] = (p["objective"] - m["objective"]) / (2 * h)
                        finite_hessian[:, j] = (np.asarray(p["objective_gradient"]) - np.asarray(m["objective_gradient"])) / (2 * h)
                    np.testing.assert_allclose(gradient, finite_gradient, rtol=2e-7, atol=2e-8)
                    np.testing.assert_allclose(hessian, finite_hessian, rtol=3e-6, atol=3e-7)

    def test_static_angle_and_invalid_physical_moments_are_rejected(self) -> None:
        static = self.native.Spin2014Objective.STATIC_LOWERING
        with self.assertRaisesRegex(ValueError, "phi=0"):
            self.native.evaluate_spin_it_2014_givens(static, RAW, GAMMA_I, GAMMA_C, PHI)
        bad_mass = RAW.copy(); bad_mass[0] = 0.0
        with self.assertRaises(ValueError):
            self.native.evaluate_spin_it_2014_givens(static, bad_mass, GAMMA_I, GAMMA_C, 0.0)
        bad_axial = RAW.copy(); bad_axial[4:6] = 0.0
        with self.assertRaises(ValueError):
            self.native.evaluate_spin_it_2014_givens(
                self.native.Spin2014Objective.TOP_INERTIA_ONLY,
                bad_axial, GAMMA_I, GAMMA_C, PHI,
            )


if __name__ == "__main__":
    unittest.main()
