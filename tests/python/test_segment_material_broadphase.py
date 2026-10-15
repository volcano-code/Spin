from fractions import Fraction as Q
import unittest
import numpy as np
from experiments.teapot_spin2014.segment_material_query import ClosedBoundaryQueries,triangle_interval


class SegmentBroadphaseTests(unittest.TestCase):
    def test_outward_bounds_equal_exhaustive_rational_intersections(self):
        v=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],dtype=float)
        f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
        query=ClosedBoundaryQueries(v,f)
        rng=np.random.default_rng(914)
        pairs=[(tuple(Q(int(x),32) for x in row[0]),tuple(Q(int(x),32) for x in row[1]))
               for row in rng.integers(-32,65,size=(60,2,3))]
        tiny=Q(1,2**60)
        pairs += [((Q(1)+tiny,Q(0),Q(0)),(Q(1)-tiny,Q(0),Q(0))),
                  ((-tiny,Q(1,4),Q(1,4)),(tiny,Q(1,4),Q(1,4)))]
        for p,q in pairs:
            expected=[]
            for i in range(len(f)):
                t,n=query._triangle(i)
                found=triangle_interval(p,q,t,n)
                if found is not None:expected.append((i,*found))
            self.assertEqual(query.intersections(p,q),expected)

    def test_rational_inside_outside_near_face(self):
        v=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],dtype=float)
        f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
        query=ClosedBoundaryQueries(v,f)
        tiny=Q(1,2**100)
        self.assertEqual(query.parity([tiny,Q(1,4),Q(1,4)])[0],'PARITY_INSIDE')
        self.assertEqual(query.parity([-tiny,Q(1,4),Q(1,4)])[0],'PARITY_OUTSIDE')


if __name__=='__main__':unittest.main()
