import tempfile
import unittest
from pathlib import Path
from fractions import Fraction as Q


class BoundedCandidateTests(unittest.TestCase):
    def module(self):
        import bounded_request_candidate as module
        return module

    def test_exact_triangle_cube_volume(self):
        module = self.module()
        self.assertEqual(module.triangle_tube_volume(
            [(Q(0),Q(0),Q(0)),(Q(1),Q(0),Q(0)),(Q(0),Q(1),Q(0))], Q(1,10)), Q(47,250))
        self.assertEqual(module.triangle_tube_volume(
            [(Q(0),)*3,(Q(1),)*3,(Q(0),)*3], Q(1,10)), Q(16,125))

    def test_large_exact_bound_serializes_without_rounding_or_global_limit_change(self):
        import sys
        before = sys.get_int_max_str_digits()
        value = self.module().text(Q(10**4500+1, 3))
        self.assertTrue(value.endswith('/3'))
        self.assertEqual(len(value.split('/')[0]), 4501)
        self.assertEqual(sys.get_int_max_str_digits(), before)

    def test_same_authority_vertex_is_reused_without_welding_indices(self):
        mapping = self.module().GlobalVertexMap()
        p = ['1/3','0/1','1/2']
        self.assertEqual(mapping.convert(p), mapping.convert(p))
        self.assertEqual(mapping.distinct_authority_vertices, 1)
        self.assertEqual(mapping.radius, abs(Q(1,3)-Q.from_float(float(Q(1,3)))))

    def test_distinct_vertices_colliding_in_float_are_rejected(self):
        mapping = self.module().GlobalVertexMap()
        mapping.convert(['1/1','0/1','0/1'])
        q = Q(1)+Q(1,2**60)
        with self.assertRaisesRegex(ValueError, 'VERTEX_COLLISION'):
            mapping.convert([f'{q.numerator}/{q.denominator}','0/1','0/1'])

    def test_rounding_degenerate_triangle_is_rejected(self):
        module = self.module()
        original = [(Q(0),)*3,(Q(1),Q(0),Q(0)),(Q(1,2),Q(1,2**1100),Q(0))]
        derived = [tuple(Q.from_float(float(x)) for x in row) for row in original]
        with self.assertRaisesRegex(ValueError, 'ZERO_AREA'):
            module.check_local_triangles(original, derived, [[0,1,2]])

    def test_local_success_never_grants_entity_or_solver_admission(self):
        module = self.module()
        v = [(Q(0),)*3,(Q(1),Q(0),Q(0)),(Q(0),Q(1),Q(0)),(Q(0),Q(0),Q(1))]
        faces = [[0,2,1],[0,1,3],[0,3,2],[1,2,3]]
        record = module.check_local_triangles(v, v, faces)
        self.assertTrue(record['closed_oriented_indexed_surface'])
        self.assertEqual(record['status'], 'PENDING_NATIVE_ENTITY_AND_PARTITION_CERTIFICATION')
        self.assertFalse(record['solver_admission_granted'])


if __name__ == '__main__': unittest.main()
