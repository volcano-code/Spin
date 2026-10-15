import unittest
import importlib.util
import numpy as np
import trimesh


class InprocessBooleanTests(unittest.TestCase):
    def test_actual_box_intersection_and_constructor_binding(self):
        self.assertIsNotNone(importlib.util.find_spec('_cherchi_core_boolean'),'in-process Cherchi module missing')
        import _cherchi_core_boolean as native
        from make_it_stand.verification.cherchi_motion_differential import normalize_closed_component_orientation
        a=trimesh.creation.box(extents=(2,2,2));b=a.copy();b.apply_translation((1,0,0))
        result=native.boolean(a.vertices.tolist(),a.faces.tolist(),b.vertices.tolist(),b.faces.tolist(),'intersection')
        vertices=np.asarray(result['vertices']);faces=np.asarray(result['faces'])
        faces,_=normalize_closed_component_orientation(vertices,faces)
        mesh=trimesh.Trimesh(vertices,faces,process=False)
        self.assertTrue(mesh.is_watertight)
        self.assertAlmostEqual(mesh.volume,4.,places=10)
        np.testing.assert_allclose(mesh.center_mass,(.5,0,0),atol=1e-12)
        self.assertEqual(len(result['construction_records']),len(vertices))
        self.assertEqual([r['vertex_index'] for r in result['construction_records']],list(range(len(vertices))))
        with self.assertRaises(ValueError):
            native.boolean(a.vertices.tolist(),a.faces.tolist(),b.vertices.tolist(),b.faces.tolist(),'unsupported')


if __name__=='__main__':unittest.main()
