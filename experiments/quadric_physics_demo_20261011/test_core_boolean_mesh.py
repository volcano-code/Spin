import importlib.util
import numpy as np
import trimesh
import unittest


class CoreMeshTests(unittest.TestCase):
    def module(self):
        self.assertIsNotNone(importlib.util.find_spec('core_boolean_mesh'),'binary64 Boolean repair is missing')
        import core_boolean_mesh
        return core_boolean_mesh

    def test_exactly_identical_vertex_alias_is_collapsed_without_moving_geometry(self):
        box=trimesh.creation.box(extents=(1,1,1))
        vertices=np.vstack([box.vertices,box.vertices[0]])
        faces=box.faces.copy();faces[faces==0]=len(box.vertices)
        faces=np.vstack([faces,[0,len(box.vertices),1]])
        mesh,report=self.module().canonical_core_mesh(vertices,faces)
        self.assertTrue(mesh.is_watertight)
        self.assertTrue(mesh.is_winding_consistent)
        self.assertEqual(mesh.volume,1.)
        self.assertEqual(report['maximum_vertex_displacement'],0.)
        self.assertEqual(report['removed_collapsed_faces'],1)
        self.assertFalse(report['independent_material_admission'])

    def test_nearby_vertices_are_not_snapped_to_make_an_open_mesh_pass(self):
        box=trimesh.creation.box(extents=(1,1,1))
        vertices=np.vstack([box.vertices,box.vertices[0]+[1e-12,0,0]])
        faces=box.faces.copy();faces[0,faces[0]==0]=len(box.vertices)
        with self.assertRaises(ValueError):self.module().canonical_core_mesh(vertices,faces)

    def test_collinear_distinct_vertices_are_not_silently_deleted(self):
        box=trimesh.creation.box(extents=(1,1,1))
        vertices=np.vstack([box.vertices,[[0,0,0],[.1,0,0],[.2,0,0]]])
        faces=np.vstack([box.faces,[[8,9,10]]])
        with self.assertRaisesRegex(ValueError,'COLLINEAR'):
            self.module().canonical_core_mesh(vertices,faces)


if __name__=='__main__':unittest.main()
