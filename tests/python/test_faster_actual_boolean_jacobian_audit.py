import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
AUDIT_PATH = ROOT / "experiments/spin_it_reproduction/strict/audit_faster_actual_boolean_jacobian.py"
sys.path.insert(0, str(AUDIT_PATH.parent))
SPEC = importlib.util.spec_from_file_location("audit_faster_actual_boolean_jacobian", AUDIT_PATH)
AUDIT = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = AUDIT
SPEC.loader.exec_module(AUDIT)


class FasterActualBooleanJacobianAuditTests(unittest.TestCase):
    def test_tangent_basis_spans_the_nine_dimensional_dual_sphere_tangent(self):
        q = np.asarray((.21, -.13, .07, .18, -.22, .11, .31, .12, -.19, .88), dtype=float)
        q /= np.linalg.norm(q)

        basis, sphere_point = AUDIT.dual_tangent_basis(q)

        self.assertEqual(basis.shape, (10, 9))
        np.testing.assert_allclose(basis.T @ basis, np.eye(9), atol=2e-14, rtol=0)
        np.testing.assert_allclose(sphere_point @ basis, np.zeros(9), atol=2e-14, rtol=0)
        np.testing.assert_allclose(np.linalg.norm(sphere_point), 1., atol=2e-15, rtol=0)

    def test_dual_geodesic_geometry_velocity_matches_independent_central_difference(self):
        q = np.asarray((-.0346, .0858, -.0500, .3437, -.0806, -.0004,
                        -.8900, .1476, -.2141, .0672), dtype=float)
        q /= np.linalg.norm(q)
        basis, _ = AUDIT.dual_tangent_basis(q)
        direction = basis[:, 3]
        step = 1e-6

        plus = AUDIT.geometry_from_dual_geodesic(q, direction, step)
        minus = AUDIT.geometry_from_dual_geodesic(q, direction, -step)
        analytic = AUDIT.geometry_direction_from_dual_tangent(q, direction)
        finite_difference = (plus - minus) / (2 * step)

        self.assertAlmostEqual(float(np.linalg.norm(plus)), 1., places=14)
        self.assertAlmostEqual(float(np.linalg.norm(minus)), 1., places=14)
        np.testing.assert_allclose(analytic, finite_difference, atol=2e-9, rtol=2e-8)

    def test_fixed_stratum_gate_rejects_cap_identity_changes(self):
        base_cap = {
            "surface_type": "HYPERBOLOID",
            "truncation_radius_hex": "0x1p+2",
            "initial_icosphere_subdivisions": 3,
            "effective_icosphere_subdivisions": 3,
            "maximum_radius_doublings": 4,
            "triangulation_hash": "triangulation-a",
        }
        moved_cap = {**base_cap, "triangulation_hash": "triangulation-b"}
        changed_cap_branch = {**base_cap, "truncation_radius_hex": "0x1.1p+2"}
        base_branch = AUDIT.cap_branch_identity(base_cap)
        self.assertEqual(base_branch, AUDIT.cap_branch_identity(moved_cap))
        self.assertNotEqual(base_branch, AUDIT.cap_branch_identity(changed_cap_branch))

        same = AUDIT.same_differential_stratum(
            "subtraction", "topology-a", base_branch,
            "subtraction", "topology-a", AUDIT.cap_branch_identity(moved_cap),
        )
        changed_cap = AUDIT.same_differential_stratum(
            "subtraction", "topology-a", base_branch,
            "subtraction", "topology-a", AUDIT.cap_branch_identity(changed_cap_branch),
        )
        changed_topology = AUDIT.same_differential_stratum(
            "subtraction", "topology-a", base_branch,
            "subtraction", "topology-b", base_branch,
        )
        self.assertTrue(same)
        self.assertFalse(changed_cap)
        self.assertFalse(changed_topology)


if __name__ == "__main__":
    unittest.main()
