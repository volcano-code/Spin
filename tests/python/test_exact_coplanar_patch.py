"""Reject geometric changes disguised as a triangulation-only branch event."""
from fractions import Fraction as F
import unittest

from make_it_stand.verification.coplanar_patch import certify_convex_quad_retriangulation


class ExactCoplanarPatch(unittest.TestCase):
    def setUp(self):
        self.points = [(F(0), F(0), F(0)), (F(2), F(0), F(0)),
                       (F(2), F(3), F(0)), (F(0), F(3), F(0))]
        self.first = [(0, 1, 2), (0, 2, 3)]
        self.second = [(0, 1, 3), (1, 2, 3)]

    def test_two_diagonals_cover_the_same_oriented_rectangle(self):
        proof = certify_convex_quad_retriangulation(self.points, self.first, self.second)
        self.assertEqual(proof['boundary_cycle'], [0, 1, 2, 3])
        self.assertEqual(proof['twice_area_vector'], ['0/1', '0/1', '12/1'])
        self.assertEqual(proof['plane'], ['0/1', '0/1', '1/1', '0/1'])
        reordered = [(2, 3, 1), (3, 0, 1)]
        self.assertEqual(certify_convex_quad_retriangulation(
            self.points, self.first, reordered), proof)

    def test_exact_nonplanarity_is_rejected_even_below_binary64_visibility(self):
        points = list(self.points)
        points[2] = (F(2), F(3), F(1, 2**512))
        with self.assertRaisesRegex(ValueError, 'coplanar'):
            certify_convex_quad_retriangulation(points, self.first, self.second)

    def test_wrong_winding_missing_face_overlap_and_degenerate_faces_reject(self):
        for bad in ([(0, 3, 1), (1, 2, 3)], [(0, 1, 3)],
                    [(0, 1, 3), (0, 1, 3)], [(0, 1, 1), (1, 2, 3)]):
            with self.subTest(faces=bad), self.assertRaises(ValueError):
                certify_convex_quad_retriangulation(self.points, self.first, bad)

    def test_concave_quad_cannot_inherit_a_convex_patch_certificate(self):
        points = list(self.points)
        points[2] = (F(1, 2), F(1, 2), F(0))
        with self.assertRaisesRegex(ValueError, 'convex'):
            certify_convex_quad_retriangulation(points, self.first, self.second)

    def test_float_impostor_wrong_index_and_extra_vertex_are_rejected(self):
        for points, faces in (([(float(x), y, z) for x, y, z in self.points], self.second),
                              (self.points, [(0, 1, 4), (1, 2, 3)]),
                              (self.points, [(True, 1, 3), (1, 2, 3)])):
            with self.assertRaises(ValueError):
                certify_convex_quad_retriangulation(points, self.first, faces)


if __name__ == '__main__':
    unittest.main()
