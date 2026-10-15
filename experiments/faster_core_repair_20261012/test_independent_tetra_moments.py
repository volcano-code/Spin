import unittest
import trimesh
import numpy as np
from independent_tetra_moments import integrate


class ExactTetraTests(unittest.TestCase):
    def test_translated_multidensity_box_against_analytic_moments(self):
        mesh=trimesh.creation.box(extents=(2.,2.,2.));mesh.apply_translation((1.,2.,3.))
        moments=integrate(mesh.vertices,mesh.faces,2.)
        expected=(16.,16.,32.,48.,16.*(1.+1./3),16.*(4.+1./3),16.*(9.+1./3),96.,48.,32.)
        np.testing.assert_allclose(moments,expected,rtol=0,atol=1e-12)


if __name__=='__main__':unittest.main()
