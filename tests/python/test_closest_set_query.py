import unittest
from fractions import Fraction as Q
from experiments.teapot_spin2014.closest_set_query import closest_set,cover
from experiments.teapot_spin2014.near_feature_witness import triangle


class ClosestSetTests(unittest.TestCase):
    def test_parallel_face_continuum(self):
        a=[[0,0,0],[1,0,0],[0,1,0]];b=[[0,0,.125],[1,0,.125],[0,1,.125]]
        result=closest_set(a,b,expected_squared_exact='1/64')
        self.assertEqual(result['dimension'],2)
        self.assertEqual(set(map(tuple,result['source_vertices_exact'])),{('0','0','0'),('1','0','0'),('0','1','0')})

    def test_skew_unique_and_parallel_edges(self):
        a=[[0,0,0],[1,0,0],[0,1,0]]
        b=[[.25,.25,.125],[.5,.25,.25],[.25,.5,.25]]
        self.assertEqual(closest_set(a,b,expected_squared_exact='1/64')['dimension'],0)
        b=[[0,0,.125],[1,0,.125],[0,-1,.125]]
        result=closest_set(a,b,expected_squared_exact='1/64')
        self.assertEqual(result['dimension'],1)
        self.assertEqual(result['source_vertices_exact'],[['0','0','0'],['1','0','0']])

    def test_nonparallel_edge_continuum(self):
        a=[[0,0,0],[1,0,0],[0,1,0]]
        b=[[0,0,.125],[1,0,.125],[0,0,1]]
        self.assertEqual(closest_set(a,b,expected_squared_exact='1/64')['dimension'],1)

    def test_segment_union_and_real_gap(self):
        poly=[['0','0','0'],['1','0','0']]
        left=triangle([[0,0,0],[.5,0,0],[0,1,0]])
        right=triangle([[.5,0,0],[1,0,0],[1,1,0]])
        self.assertTrue(cover(poly,[left,right])['all_global_minimizers_covered'])
        self.assertFalse(cover(poly,[left])['all_global_minimizers_covered'])

    def test_polygon_union_and_non_coplanar_does_not_cover_area(self):
        square=[['0','0','0'],['1','0','0'],['1','1','0'],['0','1','0']]
        first=triangle([[0,0,0],[1,0,0],[1,1,0]])
        second=triangle([[0,0,0],[1,1,0],[0,1,0]])
        self.assertTrue(cover(square,[first,second])['all_global_minimizers_covered'])
        self.assertFalse(cover(square,[first])['all_global_minimizers_covered'])
        line=triangle([[0,.5,-1],[1,.5,-1],[.5,.5,1]])
        result=cover(square,[line])
        self.assertFalse(result['all_global_minimizers_covered'])
        self.assertNotEqual(Q(result['uncovered_target_exact'][1]),Q(1,2))

    def test_excluded_corners_do_not_hide_plate_interior(self):
        top=[['0','0','1/8'],['1','0','1/8'],['0','1','1/8']]
        sides=[triangle([[0,0,0],[1,0,0],[0,0,.125]]),
               triangle([[0,0,0],[0,1,0],[0,0,.125]])]
        result=cover(top,sides)
        self.assertFalse(result['all_global_minimizers_covered'])
        self.assertFalse(result['other_local_extrema_excluded'])
        self.assertEqual(result['complete_wall_thickness'],'UNKNOWN')

    def test_malformed_polygon_rejected(self):
        with self.assertRaises(ValueError):cover([[0,0,0],[1,1,0],[0,1,0],[1,0,0]],[])
        with self.assertRaises(ValueError):cover([[0,0,0],[1,0,0],[1,1,0],[0,1,1]],[])


if __name__=='__main__':unittest.main()
