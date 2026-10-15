"""Catch dropping non-cancelling faces or admitting unbounded alias repairs."""
import json
from pathlib import Path
import unittest

from core_boolean_mesh import bounded_rational_core_mesh

HERE = Path(__file__).resolve().parent


class BoundedCoreMeshTests(unittest.TestCase):
    def test_frozen_73rd_boundary_closes_with_bounded_cancelled_alias_patch(self):
        packet = json.loads((HERE / 'noodles-dense-v5-core-failed-boolean-capture-v1'
            / 'capture-job-02/output/original-rational-output.json').read_text())
        mesh, report = bounded_rational_core_mesh(packet['vertices'], packet['faces'], length_scale=3.0)
        self.assertTrue(mesh.is_watertight)
        self.assertTrue(mesh.is_winding_consistent)
        self.assertGreater(mesh.volume, 0)
        self.assertEqual(report['removed_collapsed_faces'], 4)
        self.assertEqual(report['removed_opposite_face_pairs'], 2)
        self.assertLess(report['boundary_deviation_upper_bound'], 3e-6)
        self.assertLess(report['swept_volume_upper_bound'], 27e-6)
        self.assertFalse(report['independent_material_admission'])

    def test_same_orientation_duplicate_is_rejected(self):
        v = [['0','0','0'], ['1','0','0'], ['0','1','0'], ['0','0','1']]
        f = [[0,2,1], [0,1,3], [1,2,3], [2,0,3], [0,2,1]]
        with self.assertRaisesRegex(ValueError, 'NONCANCELLING_DUPLICATE'):
            bounded_rational_core_mesh(v, f, length_scale=2)

    def test_error_bound_is_checked_before_admission(self):
        # The exact tetrahedron has a vertex offset that rounds away at x=1.
        v = [['1','0','0'], ['100000000000000001/100000000000000000','0','0'],
             ['1','1','0'], ['1','0','1']]
        f = [[0,2,1], [0,1,3], [1,2,3], [2,0,3]]
        with self.assertRaises(ValueError):
            bounded_rational_core_mesh(v, f, length_scale=1e-15)


if __name__ == '__main__':
    unittest.main()
