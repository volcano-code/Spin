import unittest
import numpy as np


class ClearanceDomainTests(unittest.TestCase):
    def test_analytic_ball_keeps_requested_clearance_and_excludes_pin(self):
        from clearance_domain import extract_distance_domain
        # Exact ball distance, with a fixed infinite cylinder along z.
        def distance(points):
            return np.minimum(1-np.linalg.norm(points, axis=1),
                              np.linalg.norm(points[:, :2], axis=1)-.15)
        mesh, report = extract_distance_domain(distance, [-1.1]*3, [1.1]*3,
                                               spacing=.04, clearance=.12)
        self.assertTrue(mesh.is_watertight)
        self.assertGreater(mesh.volume, 0)
        self.assertLessEqual(np.linalg.norm(mesh.vertices, axis=1).max(), .88)
        self.assertGreaterEqual(np.linalg.norm(mesh.vertices[:, :2], axis=1).min(), .27)
        self.assertFalse(report['independent_exact_admission'])

    def test_rejects_nonfinite_distance_before_constructing_geometry(self):
        from clearance_domain import extract_distance_domain
        with self.assertRaisesRegex(ValueError, 'NONFINITE_DISTANCE'):
            extract_distance_domain(lambda points: np.full(len(points), np.nan),
                                    [-1]*3, [1]*3, spacing=.2, clearance=.1)


if __name__ == '__main__':
    unittest.main()
