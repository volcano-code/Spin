from fractions import Fraction as Q
import unittest
import numpy as np
from experiments.teapot_spin2014.segment_material_query import ClosedBoundaryQueries,triangle_interval


def cube(low=0.,high=1.):
    v=np.array([[low,low,low],[high,low,low],[high,high,low],[low,high,low],
                [low,low,high],[high,low,high],[high,high,high],[low,high,high]],dtype=np.float64)
    f=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
                [1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])
    return v,f


class SegmentMaterialQueryTests(unittest.TestCase):
    def test_inside_and_boundary(self):
        query=ClosedBoundaryQueries(*cube())
        self.assertEqual(query.parity([.5,.5,.5])[0],'PARITY_INSIDE')
        self.assertEqual(query.parity([2,.5,.5])[0],'PARITY_OUTSIDE')
        self.assertEqual(query.parity([0,.5,.5])[0],'BOUNDARY')
        result=query.audit_segment([0,.25,.375],[1,.25,.375])
        self.assertEqual(result['state'],'OPEN_SEGMENT_PARITY_INSIDE')
        self.assertEqual(result['complete_wall_thickness'],'UNKNOWN')

    def test_coplanar_overlap(self):
        query=ClosedBoundaryQueries(*cube())
        result=query.audit_segment([0,.2,0],[1,.2,0])
        self.assertEqual(result['state'],'BOUNDARY_OVERLAP')

    def test_cavity_and_crossing(self):
        v,f=cube(); inner,faces=cube(.25,.75)
        query=ClosedBoundaryQueries(np.vstack([v,inner]),np.vstack([f,faces[:,::-1]+8]))
        self.assertEqual(query.parity([.5,.5,.5])[0],'PARITY_OUTSIDE')
        self.assertEqual(query.parity([.1,.5,.5])[0],'PARITY_INSIDE')
        self.assertEqual(query.audit_segment([0,.5,.5],[1,.5,.5])['state'],'INTERIOR_BOUNDARY_CONTACT')
        self.assertEqual(query.audit_segment([.25,.5,.5],[.75,.5,.5])['state'],'OPEN_SEGMENT_PARITY_OUTSIDE')

    def test_exact_endpoint_and_triangle_edge(self):
        t=((Q(0),Q(0),Q(0)),(Q(1),Q(0),Q(0)),(Q(0),Q(1),Q(0)))
        n=(Q(0),Q(0),Q(1))
        result=triangle_interval((Q(1,2),Q(0),Q(-1)),(Q(1,2),Q(0),Q(1)),t,n)
        self.assertEqual(result,(Q(1,2),Q(1,2),False))
        result=triangle_interval((Q(-1),Q(0),Q(0)),(Q(2),Q(0),Q(0)),t,n)
        self.assertEqual(result[:2],(Q(1,3),Q(2,3)))

    def test_missing_endpoint_and_invalid_mesh(self):
        query=ClosedBoundaryQueries(*cube())
        self.assertEqual(query.audit_segment([.2,.2,.2],[.8,.2,.2])['state'],'ENDPOINT_NOT_ON_BOUNDARY')
        v,f=cube()
        with self.assertRaises(ValueError): ClosedBoundaryQueries(v,f[:-1])
        with self.assertRaises(ValueError): query.audit_segment([0,0,0],[0,0,0])
        with self.assertRaises(ValueError): query.parity([float('nan'),0,0])


if __name__=='__main__': unittest.main()
