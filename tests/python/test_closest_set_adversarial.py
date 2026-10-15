from fractions import Fraction as Q
import unittest
from experiments.teapot_spin2014.closest_set_query import closest_set,cover
from experiments.teapot_spin2014.near_feature_witness import triangle


class ClosestSetAdversarialTests(unittest.TestCase):
    def test_hexagonal_parallel_overlap_and_reversed_winding(self):
        a=[[-2,-1,0],[2,-1,0],[0,2,0]]
        b=[[-2,1,.125],[2,1,.125],[0,-2,.125]]
        first=closest_set(a,b,expected_squared_exact='1/64')
        second=closest_set(a[::-1],b[::-1],expected_squared_exact='1/64')
        self.assertEqual(first['dimension'],2)
        self.assertEqual(len(first['source_vertices_exact']),6)
        self.assertEqual(set(map(tuple,first['source_vertices_exact'])),set(map(tuple,second['source_vertices_exact'])))
        p=[tuple(map(Q,v)) for v in first['source_vertices_exact']]
        area=abs(sum(x[0]*y[1]-x[1]*y[0] for x,y in zip(p,p[1:]+p[:1])))/2
        self.assertEqual(area,4)

    def test_small_hole_is_not_rounded_away(self):
        square=[[0,0,0],[1,0,0],[1,1,0],[0,1,0]]
        right=.5+2**-40
        covering=[triangle([[0,0,0],[.5,0,0],[.5,1,0]]),triangle([[0,0,0],[.5,1,0],[0,1,0]]),
                  triangle([[right,0,0],[1,0,0],[1,1,0]]),triangle([[right,0,0],[1,1,0],[right,1,0]])]
        result=cover(square,covering)
        self.assertFalse(result['all_global_minimizers_covered'])
        x=Q(result['uncovered_target_exact'][0])
        self.assertLess(Q(1,2),x);self.assertLess(x,Q(right))

    def test_coverage_independent_of_triangle_order(self):
        polygon=[[0,0,0],[1,0,0],[1,1,0],[0,1,0]]
        triangles=[triangle([[0,0,0],[1,0,0],[1,1,0]]),triangle([[0,0,0],[1,1,0],[0,1,0]])]
        self.assertTrue(cover(polygon,triangles)['all_global_minimizers_covered'])
        self.assertTrue(cover(polygon,triangles[::-1])['all_global_minimizers_covered'])


if __name__=='__main__':unittest.main()
