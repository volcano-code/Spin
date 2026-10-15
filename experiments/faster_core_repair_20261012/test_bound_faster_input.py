import unittest
import numpy as np
import trimesh
from shared_input_boundary import assemble_fixed_boundary

class BoundFasterInput(unittest.TestCase):
    def function(self):
        try:
            from bound_faster_input import normalize_bound_geometry
        except ImportError:
            normalize_bound_geometry=None
        self.assertTrue(callable(normalize_bound_geometry),'bound normalized input adapter missing')
        return normalize_bound_geometry

    def test_normalization_retains_target_and_density_is_applied_once(self):
        outer=trimesh.creation.box(extents=(4,4,4));design=trimesh.creation.box(extents=(2,2,2))
        for mesh in (outer,design):mesh.apply_translation((1,2,3))
        fixed=assemble_fixed_boundary(outer,design)
        result=self.function()(outer,design,fixed,[1,2,1],1200.)
        np.testing.assert_array_equal(result['target_origin'],[0,0,-.5])
        self.assertEqual(result['normalization']['L'],4.)
        self.assertEqual(result['normalization']['mass_unit_kg'],76800.)
        self.assertAlmostEqual(result['fixed_shell_moments'][0],.875,places=14)
        np.testing.assert_allclose(result['fixed_shell_moments'][1:4],[0,0,0],atol=1e-14)
        np.testing.assert_allclose(result['fixed_shell_moments'][4:7],[31/384]*3,atol=1e-14)
        other=self.function()(outer,design,fixed,[1,2,1],800.)
        np.testing.assert_array_equal(other['fixed_shell_moments'],result['fixed_shell_moments'])
        self.assertFalse(result['solver_admission'])

if __name__=='__main__':unittest.main()
