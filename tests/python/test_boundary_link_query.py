import unittest
import numpy as np
from experiments.teapot_spin2014.segment_material_query import ClosedBoundaryQueries
from experiments.teapot_spin2014.boundary_link_query import endpoint_link


class BoundaryLinkTests(unittest.TestCase):
    def test_vertex_link_two_switches_not_an_event(self):
        v=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],dtype=float)
        f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
        q=ClosedBoundaryQueries(v,f)
        result=endpoint_link(q,[.5,0,0],[0,0,0],[0,1,2])
        self.assertEqual(result['details']['cyclic_switches'],2)
        self.assertFalse(result['local_section_5_1_event'])
        self.assertEqual(result['wall_thickness_admission'],'NOT_GRANTED')

    def test_vertex_all_internal_and_incomplete_link(self):
        v=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],dtype=float)
        f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
        q=ClosedBoundaryQueries(v,f)
        result=endpoint_link(q,[1,1,1],[0,0,0],[0,1,2])
        self.assertTrue(result['local_section_5_1_event'])
        with self.assertRaises(ValueError): endpoint_link(q,[1,1,1],[0,0,0],[0,1])

    def test_edge_internal_external_and_missing_incidence(self):
        v=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],dtype=float)
        f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
        q=ClosedBoundaryQueries(v,f)
        positive=endpoint_link(q,[.5,.1,.1],[.5,0,0],[0,1])
        negative=endpoint_link(q,[.5,.1,-.1],[.5,0,0],[0,1])
        self.assertTrue(positive['local_section_5_1_event'])
        self.assertFalse(negative['local_section_5_1_event'])
        with self.assertRaises(ValueError):endpoint_link(q,[.5,.1,.1],[.5,0,0],[0])


if __name__=='__main__': unittest.main()
