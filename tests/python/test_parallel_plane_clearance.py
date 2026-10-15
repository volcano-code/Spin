import math
import unittest

from experiments.teapot_spin2014.parallel_plane_clearance import audit_parallel_planes


class ParallelPlaneBoundTests(unittest.TestCase):
    def pair(self, gap=0.125):
        return [[0,0,0], [1,0,0], [0,1,0], [0,0,gap], [0,1,gap], [1,0,gap]], [[0,1,2], [3,4,5]]

    def test_endpoint_and_no_global_admission(self):
        v, f = self.pair()
        result = audit_parallel_planes(v, f, 0.125)
        self.assertEqual(result['status'], 'PARALLEL_PLANE_LOWER_BOUND_AT_LEAST_THRESHOLD')
        self.assertEqual(result['minimum_squared_plane_separation_exact'], '1/64')
        self.assertEqual(result['complete_wall_thickness'], 'UNKNOWN')
        self.assertEqual(result['geometry_admission'], 'NOT_GRANTED')
        self.assertEqual(audit_parallel_planes(v, f, math.nextafter(.125, math.inf))['status'], 'LOWER_BOUND_INCONCLUSIVE')

    def test_small_plane_gap_does_not_prove_wall(self):
        v, f = self.pair(.01)
        for row in v[3:]:
            row[0] += 100  # disjoint projected triangles, not a material wall
        result = audit_parallel_planes(v, f, .015)
        self.assertEqual(result['status'], 'LOWER_BOUND_INCONCLUSIVE')
        self.assertFalse(result['thin_material_witness'])

    def test_winding_and_exact_tilt(self):
        v, f = self.pair(.125)
        v = [[x,y,x+z] for x,y,z in v]
        result = audit_parallel_planes(v, [list(reversed(x)) for x in f], .01)
        self.assertEqual(result['minimum_squared_plane_separation_exact'], '1/128')

    def test_coplanar_count_and_minimum(self):
        v, f = self.pair(.5)
        v += [[0,0,.125], [1,0,.125], [0,1,.125]]
        result = audit_parallel_planes(v, f + [[0,2,1], [6,7,8]], .125)
        self.assertEqual(result['coplanar_triangle_pairs_excluded'], 1)
        self.assertEqual(result['distinct_plane_triangle_pairs'], 5)
        self.assertEqual(result['minimum_squared_plane_separation_exact'], '1/64')
        self.assertEqual(result['supporting_plane_face_indices'], [0,3])

    def test_nonparallel_does_not_pass(self):
        v, f = self.pair()
        v[5][2] += .01
        result = audit_parallel_planes(v, f, .015)
        self.assertEqual(result['status'], 'NO_DISTINCT_PARALLEL_PLANES')
        self.assertEqual(result['complete_wall_thickness'], 'UNKNOWN')

    def test_bad_geometry(self):
        v, f = self.pair()
        for bad_f in ([], [[0,0,1]], [[0,1,99]], [[False,1,2]], [[0.,1,2]]):
            with self.assertRaises(ValueError):
                audit_parallel_planes(v, bad_f, .015)
        for t in (0, -1, math.nan, math.inf, True):
            with self.assertRaises(ValueError):
                audit_parallel_planes(v, f, t)
        v[0][0] = math.nan
        with self.assertRaises(ValueError):
            audit_parallel_planes(v, f, .015)


if __name__ == '__main__':
    unittest.main()
