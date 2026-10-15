"""Input-design invariants; these tests do not grant optimization success."""
import unittest
import numpy as np
import trimesh

from rounded_input import rounded_components, signed_moments, capacity_bounds


class RoundedInputTests(unittest.TestCase):
    def test_closed_nonprincipal_body_and_connected_pin_contact(self):
        body, inner, pin = rounded_components()
        for mesh in (body, inner, pin):
            self.assertTrue(mesh.is_watertight)
            self.assertTrue(mesh.is_winding_consistent)
            self.assertGreater(mesh.volume, 0)
        np.testing.assert_allclose(signed_moments(body)[1:4] / body.volume,
                                   [.1, 0, 1.12], atol=1e-13)
        self.assertGreater(abs(signed_moments(body)[8]), .01)
        self.assertGreater(np.ptp(body.vertices[:, 1]), 1.1)
        self.assertAlmostEqual(inner.volume / body.volume, .9**3, places=12)
        self.assertAlmostEqual(pin.bounds[0, 2], -.16, places=12)
        self.assertGreater(pin.bounds[1, 2], 2.25)
        self.assertGreater(pin.bounds[1, 0], .10)

    def test_capacity_is_a_necessary_bound_not_a_feasible_certificate(self):
        box = trimesh.creation.box(extents=[2, 2, 2])
        fixed = trimesh.creation.box(extents=[1, 1, 1])
        fixed.apply_translation([3, 0, 0])
        record = capacity_bounds(fixed.vertices, fixed.faces, box.vertices, box.faces)
        self.assertEqual(record['x']['design_min_exact'], '-2')
        self.assertEqual(record['x']['design_max_exact'], '2')
        self.assertEqual(record['status'], 'EXCLUDED_BY_FIRST_MOMENT_BOUND')
        self.assertFalse(record['feasible_witness'])
        fixed.apply_translation([-3, 0, 0])
        record = capacity_bounds(fixed.vertices, fixed.faces, box.vertices, box.faces)
        self.assertEqual(record['status'], 'UNKNOWN_NOT_EXCLUDED_BY_FIRST_MOMENT_BOUNDS')
        self.assertFalse(record['feasible_witness'])


if __name__ == '__main__':
    unittest.main()
