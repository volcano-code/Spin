"""Catches independent interior rounding and forgotten cavity orientation."""
import unittest
import numpy as np
import trimesh
import _make_it_stand as native

class SharedInputBoundary(unittest.TestCase):
    def test_same_design_boundary_defines_fixed_cavity_and_exact_partition(self):
        try:
            from shared_input_boundary import assemble_fixed_boundary
        except ImportError:
            assemble_fixed_boundary=None
        self.assertTrue(callable(assemble_fixed_boundary),'shared boundary assembly missing')
        outer=trimesh.creation.box(extents=(4,4,4));design=trimesh.creation.box(extents=(2,2,2))
        fixed=assemble_fixed_boundary(outer,design)
        self.assertTrue(fixed.is_watertight)
        self.assertAlmostEqual(fixed.volume,56,places=12)
        np.testing.assert_array_equal(fixed.vertices[len(outer.vertices):],design.vertices)
        np.testing.assert_array_equal(fixed.faces[len(outer.faces):]-len(outer.vertices),design.faces[:,::-1])
        region=lambda mesh:native.make_exact_material_region(mesh.vertices,mesh.faces.astype(np.int32))
        parts=native.check_exact_material_partition(region(outer),[region(fixed),region(design)])
        self.assertTrue(parts['completed'])
        self.assertTrue(parts['complete_partition'])

if __name__=='__main__':unittest.main()
