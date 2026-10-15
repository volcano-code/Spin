import unittest
import numpy as np
from experiments.teapot_spin2014.segment_material_query import ClosedBoundaryQueries
from experiments.teapot_spin2014.source_star_query import SourceStarQuery


def thin_box():
    v=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],
                [0,0,.01],[1,0,.01],[1,1,.01],[0,1,.01]],dtype=float)
    f=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
                [1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])
    return SourceStarQuery(ClosedBoundaryQueries(v,f))


class SourceStarTests(unittest.TestCase):
    def test_excluded_minimizer_does_not_eliminate_real_thin_plate(self):
        q=thin_box()
        corner=q.audit(0,[0,0,0],[0,0,.01])
        interior=q.audit(0,[.75,.25,0],[.75,.25,.01])
        self.assertTrue(corner['target_in_initial_star'])
        self.assertFalse(interior['target_in_initial_star'])
        self.assertTrue(interior['moving_event_algebraic_condition'])
        self.assertEqual(interior['orthogonal_target_faces'],[2])
        self.assertFalse(corner['whole_triangle_pair_excluded'])
        self.assertEqual(interior['wall_thickness_admission'],'NOT_GRANTED')

    def test_edge_source_and_invalid_center(self):
        q=thin_box()
        result=q.audit(0,[.5,0,0],[.5,0,.01],source_edge=[0,1])
        self.assertTrue(result['target_in_initial_star'])
        self.assertEqual(result['source_mode'],'EDGE')
        with self.assertRaises(ValueError):q.audit(0,[.75,.25,0],[.75,.25,.01],source_edge=[0,1])
        with self.assertRaises(ValueError):q.audit(0,[.5,0,0],[.5,0,.01],source_edge=[0,4])
        with self.assertRaises(ValueError):q.audit(-1,[0,0,0],[0,0,.01])
        with self.assertRaises(ValueError):q.audit(True,[0,0,0],[0,0,.01])

    def test_full_incidence_and_nonperpendicular(self):
        q=thin_box()
        self.assertEqual(q.incident_faces([.75,.25,.01]),(2,))
        result=q.audit(0,[.7,.2,0],[.75,.25,.01])
        self.assertFalse(result['moving_event_algebraic_condition'])
        with self.assertRaises(ValueError):q.audit(0,[.75,.25,0],[.75,.25,.02])


if __name__=='__main__':unittest.main()
