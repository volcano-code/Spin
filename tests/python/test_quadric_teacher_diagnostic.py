"""Contract tests for the opt-in same-call legacy teacher diagnostic export."""

import unittest

import numpy as np

from test_exact_material_region import box
from make_it_stand import spin


COEFFICIENTS = (-.01, .01, .005, -.003, 1.1, .8, 1.3, .12, -.08, .17)


class QuadricTeacherDiagnosticTests(unittest.TestCase):
    def request(self):
        outer_vertices, outer_faces = box((-2.0,) * 3, (2.0,) * 3)
        design_vertices, design_faces = box((-1.0,) * 3, (1.0,) * 3)
        return spin.SpinItFasterRequest(
            outer_vertices.tolist(), outer_faces.tolist(),
            design_vertices.tolist(), design_faces.tolist(),
            mode="YO_YO",
            settings=spin.SpinSolverConfig(
                equation_profile="LEGACY_NATIVE_NORMALIZED_V1",
                maximum_iterations=1,
                maximum_restarts=1,
                maximum_function_evaluations_per_restart=1,
            ),
        )

    def test_diagnostic_exports_one_contract_and_does_not_change_default_api(self):
        request = self.request()
        native, native_request = spin._native_request(request)

        # The established API remains moments-and-hashes only.
        ordinary = native.evaluate_spin_it_faster_quadric(
            native_request, np.asarray(COEFFICIENTS, dtype=np.float64)
        )
        self.assertEqual(
            set(ordinary),
            {"variable_raw_moments", "total_raw_moments", "geometry_hash",
             "topology_hash", "interface_type", "interface_subtype"},
        )

        self.assertTrue(hasattr(spin, "evaluate_spin_it_faster_quadric_diagnostic"))
        diagnostic = spin.evaluate_spin_it_faster_quadric_diagnostic(
            request, COEFFICIENTS
        )
        self.assertEqual(diagnostic.variable_raw_moments, tuple(ordinary["variable_raw_moments"]))
        self.assertEqual(diagnostic.total_raw_moments, tuple(ordinary["total_raw_moments"]))
        self.assertEqual(diagnostic.geometry_hash, ordinary["geometry_hash"])
        self.assertEqual(diagnostic.topology_hash, ordinary["topology_hash"])
        self.assertEqual(diagnostic.equation_profile, "LEGACY_NATIVE_NORMALIZED_V1")
        self.assertEqual(diagnostic.derivative_method, "legacy-continuous-implicit-shape-v1")
        self.assertEqual(len(diagnostic.fixed_raw_moments), 10)
        np.testing.assert_allclose(
            diagnostic.fixed_raw_moments,
            np.asarray(diagnostic.total_raw_moments)
            - np.asarray(diagnostic.variable_raw_moments),
            rtol=0.0,
            atol=1e-14,
        )
        self.assertEqual(diagnostic.derivative_defined, diagnostic.moment_jacobian is not None)
        if diagnostic.moment_jacobian is not None:
            self.assertEqual(len(diagnostic.moment_jacobian), 10)
            self.assertTrue(all(len(row) == 10 for row in diagnostic.moment_jacobian))
            self.assertTrue(np.isfinite(diagnostic.moment_jacobian).all())
        self.assertGreaterEqual(diagnostic.residual_evaluation_seconds, 0.0)
        self.assertGreaterEqual(diagnostic.derivative_evaluation_seconds, 0.0)
        self.assertTrue(diagnostic.base_request_hash)
        self.assertTrue(diagnostic.evaluator_request_hash)

    def test_undefined_derivative_is_not_serialized_as_zero(self):
        request = self.request()
        diagnostic = spin._quadric_teacher_diagnostic_from_native(
            {
                "variable_raw_moments": np.zeros(10),
                "total_raw_moments": np.zeros(10),
                "moment_jacobian": None,
                "derivative_defined": False,
                "derivative_method": "legacy-continuous-implicit-shape-v1",
                "derivative_diagnostic": "legacy interface derivative undefined",
                "residual_evaluation_seconds": 0.0,
                "derivative_evaluation_seconds": 0.0,
                "geometry_hash": "geometry",
                "topology_hash": "topology",
                "base_request_hash": "base",
                "evaluator_request_hash": "evaluation",
                "derivative_scientific_hash": "scientific",
                "provenance_topology_hash": "provenance",
            },
            request,
            (.01,) + (0.0,) * 9,
        )
        self.assertFalse(diagnostic.derivative_defined)
        self.assertIsNone(diagnostic.moment_jacobian)
        self.assertTrue(diagnostic.derivative_diagnostic)

    def test_invalid_coefficients_fail_closed(self):
        request = self.request()
        for coefficients in ((1.0,) * 9, (1.0,) * 11, (float("nan"),) + (0.0,) * 9):
            with self.subTest(coefficients=coefficients):
                with self.assertRaises(ValueError):
                    spin.evaluate_spin_it_faster_quadric_diagnostic(request, coefficients)

    def test_native_claim_of_defined_but_invalid_jacobian_fails_closed(self):
        request = self.request()
        payload = {
            "variable_raw_moments": np.zeros(10),
            "total_raw_moments": np.zeros(10),
            "moment_jacobian": np.zeros((10, 9)),
            "derivative_defined": True,
            "derivative_method": "legacy-continuous-implicit-shape-v1",
            "derivative_diagnostic": "defined",
            "residual_evaluation_seconds": 0.0,
            "derivative_evaluation_seconds": 0.0,
            "geometry_hash": "geometry",
            "topology_hash": "topology",
            "base_request_hash": "base",
            "evaluator_request_hash": "evaluation",
            "derivative_scientific_hash": "scientific",
            "provenance_topology_hash": "provenance",
        }
        with self.assertRaisesRegex(RuntimeError, "non-10x10"):
            spin._quadric_teacher_diagnostic_from_native(
                payload, request, (.01,) + (0.0,) * 9
            )


if __name__ == "__main__":
    unittest.main()
