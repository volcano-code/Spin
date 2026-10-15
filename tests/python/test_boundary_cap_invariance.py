from fractions import Fraction as Q
import unittest
import numpy as np
from experiments.teapot_spin2014.segment_material_query import ClosedBoundaryQueries
from experiments.teapot_spin2014.boundary_cap_query import BoundaryCapQuery


def cube():
    v=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,.125],[1,0,.125],[1,1,.125],[0,1,.125]],dtype=float)
    f=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
                [1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])
    return v,f


def states(row):
    output=[]
    for event in row['critical_intervals']:
        state=tuple(sorted(c['euler'] for c in event['components']))
        if not output or output[-1][1]!=state:output.append((event['after_squared_radius_exact'],state))
    return output


class BoundaryCapInvarianceTests(unittest.TestCase):
    def test_face_subdivision_preserves_topology_filtration(self):
        v,f=cube();vertices=v.tolist();faces=[]
        for a,b,c in f:
            index=len(vertices);vertices.append(((v[a]+v[b]+v[c])/3).tolist())
            faces.extend([[a,b,index],[b,c,index],[c,a,index]])
        original=BoundaryCapQuery(ClosedBoundaryQueries(v,f)).audit([.5,.5,0],Q(1,16))
        refined=BoundaryCapQuery(ClosedBoundaryQueries(np.array(vertices),np.array(faces))).audit([.5,.5,0],Q(1,16))
        self.assertEqual(states(original),states(refined))

    def test_closed_genus_one_surface_is_not_disk(self):
        outer=[[-2,-2],[2,-2],[2,2],[-2,2]];inner=[[-1,-1],[1,-1],[1,1],[-1,1]]
        v=np.array([[x,y,z] for z in [0,.25] for x,y in outer+inner],dtype=float);f=[]
        def quad(a,b,c,d):f.extend([[a,b,c],[a,c,d]])
        for i in range(4):
            j=(i+1)%4
            quad(i,j,j+8,i+8);quad(j+4,i+4,i+12,j+12)
            quad(i+8,j+8,j+12,i+12);quad(i+4,j+4,j,i)
        row=BoundaryCapQuery(ClosedBoundaryQueries(v,np.array(f))).audit(v[0],Q(100))
        self.assertEqual(row['critical_intervals'][-1]['components'][0]['euler'],0)
        self.assertFalse(row['threshold_cap_single_disk'])

    def test_pinched_vertex_rejected(self):
        v=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[-1,0,0],[0,-1,0],[0,0,-1]],dtype=float)
        f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3], [0,4,5],[0,6,4],[0,5,6],[4,6,5]])
        with self.assertRaises(ValueError):BoundaryCapQuery(ClosedBoundaryQueries(v,f))


if __name__=='__main__':unittest.main()
