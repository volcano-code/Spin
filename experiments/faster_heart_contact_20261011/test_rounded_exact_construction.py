"""Small exact-operation benchmark; not a reconstructed Heart solve."""
import unittest
import numpy as np
import trimesh
import _make_it_stand as native
from make_it_stand.verification.spin_problem import MeshSnapshot
from admit_rounded_input import exact_rows, union, construct
from rounded_input import exact_moments


class ExactConstructionTests(unittest.TestCase):
    def test_union_difference_and_reloaded_exact_partition(self):
        a = trimesh.creation.box(extents=[2, 2, 2])
        b = a.copy(); b.apply_translation([1, 0, 0])
        box = trimesh.creation.box(extents=[8, 8, 8])
        def operand(mesh):
            return exact_rows(mesh.vertices), mesh.faces.tolist()
        merged = union(native, operand(a), operand(b), operand(box))
        self.assertEqual(exact_moments(*merged)[0], 12)
        retained = construct(native, merged, operand(a), 'difference')
        self.assertEqual(exact_moments(*retained)[0], 4)
        snapshots = [MeshSnapshot.from_rationals(*x) for x in (merged, operand(a), retained)]
        report = native.check_exact_material_partition(snapshots[0].exact_region(),
                                                       [s.exact_region() for s in snapshots[1:]])
        self.assertTrue(report['completed'])
        self.assertTrue(report['complete_partition'])


if __name__ == '__main__':
    unittest.main()
