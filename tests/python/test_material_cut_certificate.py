from fractions import Fraction as Q
import unittest
import numpy as np
from experiments.teapot_spin2014.segment_material_query import ClosedBoundaryQueries
from experiments.teapot_spin2014.boundary_cap_query import BoundaryCapQuery
from experiments.teapot_spin2014.material_cut_certificate import material_cut_witness


def mesh(z0=0.,z1=1.):
    v=np.array([[0,0,z0],[1,0,z0],[1,1,z0],[0,1,z0],[0,0,z1],[1,0,z1],[1,1,z1],[0,1,z1]],dtype=float)
    f=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
                [1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])
    return v,f


class MaterialCutCertificateTests(unittest.TestCase):
    def test_plate_through_cut_changes_euler(self):
        query=BoundaryCapQuery(ClosedBoundaryQueries(*mesh(0,.125)))
        record=material_cut_witness(query,[.5,.5,0],Q(1,16))
        self.assertEqual(record['state'],'MATERIAL_CUT_EULER_CHANGE_WITNESS')
        self.assertEqual(record['cut_euler_change'],-1)
        self.assertFalse(record['wall_metric_adopted'])
        self.assertEqual(record['wall_thickness_admission'],'NOT_GRANTED')

    def test_air_gap_is_not_material_cut(self):
        a,af=mesh(-1,0);b,bf=mesh(.125,1.125)
        query=BoundaryCapQuery(ClosedBoundaryQueries(np.vstack([a,b]),np.vstack([af,bf+8])))
        record=material_cut_witness(query,[.5,.5,0],Q(1,16))
        self.assertEqual(record['state'],'NO_MATERIAL_CUT_WITNESS')
        self.assertIsNone(record['cut_euler_change'])

    def test_corner_no_event(self):
        query=BoundaryCapQuery(ClosedBoundaryQueries(*mesh()))
        record=material_cut_witness(query,[.005,0,.005],Q(1,10000))
        self.assertEqual(record['state'],'NO_FIXED_CENTER_CAP_EVENT')


if __name__=='__main__':unittest.main()
