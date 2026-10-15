"""A thin collapsed face can be close to a surviving edge, not its endpoints."""
from fractions import Fraction as Q
import unittest
import numpy as np
import trimesh
from core_boolean_mesh import bounded_rational_core_mesh

class EdgeErasureWitness(unittest.TestCase):
    def fixture(self):
        box=trimesh.creation.box(extents=(1,1,1))
        box.apply_translation((1,1,1))
        points=[list(map(lambda x:Q(float(x)),row)) for row in box.vertices]
        a,b=map(int,box.faces[0][:2]);t=Q(1,10**20)
        points.append([x+t*(y-x) for x,y in zip(points[a],points[b])]);p=len(points)-1
        faces=[]
        for face in box.faces:
            if a in face and b in face:
                for i in range(3):
                    u,v,w=map(int,(face[i],face[(i+1)%3],face[(i+2)%3]))
                    if {u,v}=={a,b}:faces.extend([[u,p,w],[p,v,w]]);break
            else:faces.append(list(map(int,face)))
        return [[str(x.numerator)+'/'+str(x.denominator) for x in row] for row in points],np.asarray(faces)

    def test_exact_edge_bound_accepts_thin_patch_without_changing_legacy(self):
        points,faces=self.fixture()
        with self.assertRaisesRegex(ValueError,'ERROR_EXCEEDED'):
            bounded_rational_core_mesh(points,faces,length_scale=2.)
        try:
            result=bounded_rational_core_mesh(points,faces,length_scale=2.,
                boundary_witness_policy='SURVIVING_EDGE_V2')
        except TypeError:
            self.fail('explicit exact surviving-edge witness policy missing')
        mesh,report=result
        self.assertTrue(mesh.is_watertight)
        self.assertAlmostEqual(mesh.volume,1.,places=12)
        self.assertLess(report['boundary_deviation_upper_bound'],1e-18)
        self.assertTrue(any('surviving_edge' in row for row in report['removed_face_witnesses']))
        self.assertFalse(report['independent_material_admission'])

    def test_edge_witness_does_not_relax_coordinate_rounding_budget(self):
        points,faces=self.fixture()
        with self.assertRaisesRegex(ValueError,'ERROR_EXCEEDED'):
            bounded_rational_core_mesh(points,faces,length_scale=1e-15,
                boundary_witness_policy='SURVIVING_EDGE_V2')

if __name__=='__main__':unittest.main()
