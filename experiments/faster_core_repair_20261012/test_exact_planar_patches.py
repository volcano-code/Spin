"""Catch treating a noncoplanar edge flip as harmless retriangulation."""
import unittest
from fractions import Fraction as F
from exact_planar_patches import patch_boundary_signature


class PlanarPatchTests(unittest.TestCase):
    def test_diagonal_flip_preserves_exact_patch(self):
        v=[[0,0,0],[1,0,0],[1,1,0],[0,1,0]]
        self.assertEqual(patch_boundary_signature(v,[[0,1,2],[0,2,3]],closed=False)['patches'],
                         patch_boundary_signature(v,[[0,1,3],[1,2,3]],closed=False)['patches'])

    def test_bent_diagonal_is_not_certified_as_coplanar(self):
        v=[[0,0,0],[1,0,0],[1,1,F(1,10**30)],[0,1,0]]
        self.assertNotEqual(patch_boundary_signature(v,[[0,1,2],[0,2,3]],closed=False)['patches'],
                            patch_boundary_signature(v,[[0,1,3],[1,2,3]],closed=False)['patches'])

    def test_true_zero_area_is_rejected_but_small_nonzero_is_not(self):
        with self.assertRaises(ValueError):
            patch_boundary_signature([[0,0,0],[1,0,0],[2,0,0]],[[0,1,2]],closed=False)
        result=patch_boundary_signature([[0,0,0],[1,0,0],[0,F(1,10**60),0]],[[0,1,2]],closed=False)
        self.assertEqual(len(result['patches']),1)

    def test_open_boundary_cannot_pass_closed_protocol(self):
        with self.assertRaises(ValueError):
            patch_boundary_signature([[0,0,0],[1,0,0],[0,1,0]],[[0,1,2]])


if __name__=='__main__': unittest.main()
