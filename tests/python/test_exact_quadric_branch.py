"""Independent analytic invariants; not Boolean or Newton success tests."""
import unittest
from make_it_stand.verification import quadric_branch


class ExactQuadricBranch(unittest.TestCase):
    def q(self,constant,diagonal,linear=(0,0,0),cross=(0,0,0)):
        return [f'{x}/1' for x in (constant,*linear,*diagonal,*cross)]

    def test_shifted_level_distinguishes_one_two_sheet_and_cone(self):
        for c,want,level in ((-1,'ONE_SHEET_HYPERBOLOID','-1/1'),
                             (1,'TWO_SHEET_HYPERBOLOID','1/1'),(0,'CONE','0/1')):
            value=quadric_branch.exact_full_rank_invariants(self.q(c,(1,1,-1)))
            self.assertEqual(value['subtype'],want)
            self.assertEqual(value['centered_level'],level)
            self.assertEqual(value['signature'],{'positive':2,'negative':1,'zero':0})
        translated=quadric_branch.exact_full_rank_invariants(self.q(3,(1,1,-1),(4,0,0)))
        self.assertEqual(translated['centered_level'],'-1/1')
        self.assertEqual(translated['subtype'],'ONE_SHEET_HYPERBOLOID')

    def test_offdiagonal_matrix_and_zero_diagonal_pivot_block(self):
        # q_xy is A_xy (the polynomial cross term is 2*q_xy), not A_xy/2.
        value=quadric_branch.exact_full_rank_invariants(self.q(-1,(0,0,1),cross=(0,0,1)))
        self.assertEqual(value['determinant'],'-1/1')
        self.assertEqual(value['signature'],{'positive':2,'negative':1,'zero':0})
        self.assertEqual(value['subtype'],'ONE_SHEET_HYPERBOLOID')

    def test_scale_sign_invariance_and_nonreal_ellipsoid(self):
        for sign in (1,-1):
            value=quadric_branch.exact_full_rank_invariants(self.q(-sign,(sign,sign,sign)))
            self.assertEqual(value['subtype'],'ELLIPSOID')
        self.assertEqual(quadric_branch.exact_full_rank_invariants(
            self.q(1,(1,1,1)))['subtype'],'NO_REAL_ZERO_SET')

    def test_rank_deficient_not_misclassified_as_infeasible(self):
        value=quadric_branch.exact_full_rank_invariants(self.q(-1,(1,1,0)))
        self.assertEqual(value['subtype'],'RANK_DEFICIENT_OUTSIDE_THIS_CHECK')
        self.assertIsNone(value['centered_level'])
        self.assertFalse(value['optimization_or_material_qualified'])
        with self.assertRaises(ValueError): quadric_branch.exact_full_rank_invariants(['NaN']*10)


if __name__=='__main__': unittest.main()
