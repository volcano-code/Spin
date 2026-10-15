from fractions import Fraction as Q
import unittest
import numpy as np
from experiments.teapot_spin2014.segment_material_query import ClosedBoundaryQueries
from experiments.teapot_spin2014.boundary_cap_query import BoundaryCapQuery


def box(height=1.):
    v=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],
                [0,0,height],[1,0,height],[1,1,height],[0,1,height]],dtype=float)
    f=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
                [1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])
    return BoundaryCapQuery(ClosedBoundaryQueries(v,f))


class BoundaryCapTests(unittest.TestCase):
    def test_ordinary_face_edge_corner_are_disk_caps(self):
        query=box()
        for p in ([.5,.5,0],[0,.5,0],[0,0,0],[.005,0,.005]):
            row=query.audit(p,Q(1,10000))
            self.assertEqual(row['state'],'SINGLE_DISK_ALL_POSITIVE_RADII_TO_THRESHOLD')
            self.assertEqual(row['complete_wall_thickness'],'UNKNOWN')

    def test_thin_slab_detected_without_corner_false_positive(self):
        query=box(.125)
        at=query.audit([.5,.5,0],Q(1,64))
        above=query.audit([.5,.5,0],Q(1,64)+Q(1,2**60))
        self.assertTrue(at['threshold_cap_single_disk'])
        self.assertEqual(above['state'],'NON_DISK_INTERVAL_FOUND')
        self.assertEqual(above['first_non_disk_interval']['after_squared_radius_exact'],'1/64')
        self.assertEqual(len(above['first_non_disk_interval']['components']),2)

    def test_transient_failure_not_hidden_by_final_disk(self):
        row=box(.01).audit([.03,.5,0],Q(1,100))
        self.assertTrue(row['threshold_cap_single_disk'])
        self.assertEqual(row['state'],'NON_DISK_INTERVAL_FOUND')

    def test_whole_sphere_boundary_is_not_disk(self):
        row=box().audit([0,0,0],Q(4))
        last=row['critical_intervals'][-1]
        self.assertFalse(last['single_disk'])
        self.assertEqual(last['components'][0]['euler'],2)

    def test_center_and_radius_rejected(self):
        query=box()
        with self.assertRaises(ValueError):query.audit([.5,.5,.5],Q(1,100))
        with self.assertRaises(ValueError):query.audit([0,0,0],Q(0))


if __name__=='__main__':unittest.main()
