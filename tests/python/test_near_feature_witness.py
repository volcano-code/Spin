import unittest
from experiments.teapot_spin2014.near_feature_witness import disjoint_witness


class NearFeatureWitnessTests(unittest.TestCase):
    def test_parallel_interior_segment(self):
        a=[[0,0,0],[1,0,0],[0,1,0]]
        b=[[0,0,.125],[0,1,.125],[1,0,.125]]
        result=disjoint_witness(a,b,expected_squared_exact='1/64')
        self.assertEqual(result['local_halfspaces'],'BOTH_POSITIVE')
        reversed_result=disjoint_witness(a[::-1],b[::-1],expected_squared_exact='1/64')
        self.assertEqual(reversed_result['local_halfspaces'],'BOTH_NEGATIVE')
        self.assertEqual(reversed_result['open_segment_material_membership'],'UNKNOWN')

    def test_nonparallel_vertex_face(self):
        a=[[0,0,0],[1,0,0],[0,1,0]]
        b=[[.25,.25,.125],[.5,.25,.25],[.25,.5,.25]]
        result=disjoint_witness(a,b,expected_squared_exact='1/64')
        self.assertEqual(result['point_a_exact'],['1/4','1/4','0'])
        self.assertEqual(result['point_b_exact'],['1/4','1/4','1/8'])

    def test_skew_edge_edge(self):
        a=[[-1,0,0],[1,0,0],[0,-1,0]]
        b=[[0,-.5,.25],[0,.5,.25],[0,0,1]]
        result=disjoint_witness(a,b,expected_squared_exact='1/16')
        self.assertEqual(result['squared_exact'],'1/16')

    def test_incorrect_distance_and_intersection_rejected(self):
        a=[[0,0,0],[1,0,0],[0,1,0]]
        b=[[0,0,1],[1,0,1],[0,1,1]]
        for expected in ('0','-1','1/2'):
            with self.assertRaises(ValueError): disjoint_witness(a,b,expected_squared_exact=expected)


if __name__=='__main__': unittest.main()
