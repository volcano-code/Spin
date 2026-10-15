import unittest
from unittest.mock import patch
import numpy as np


class DenseNoodlesTests(unittest.TestCase):
    def test_dense_horizontal_array_has_frozen_shape_and_central_connection(self):
        from noodles_dense_input import SPEC, tube_layout, dense_centerline, identity
        layout=tube_layout()
        self.assertGreaterEqual(len(layout),40)
        self.assertEqual(identity(),identity())
        original=identity()
        with patch.dict(SPEC,{'density_kg_m3':SPEC['density_kg_m3']+1}):
            self.assertNotEqual(identity(),original)
        for index in range(len(layout)):
            points=dense_centerline(index)
            middle=points[len(points)//2]
            self.assertLess(np.linalg.norm(middle[:2]),SPEC['hub_radius_model'])
            self.assertGreater(middle[2],SPEC['hub_interval_model'][0])
            self.assertLess(middle[2],SPEC['hub_interval_model'][1])
            self.assertGreater(np.ptp(points[:,0]),2*np.ptp(points[:,2]))
            self.assertTrue(np.all(np.diff(points[:,0])>0))

    def test_closed_parts_and_no_degenerate_triangles(self):
        from noodles_dense_input import components
        from noodles_input import mesh_audit
        for name,mesh in components().items():
            with self.subTest(name=name):
                facts=mesh_audit(mesh)
                self.assertTrue(facts['closed'])
                self.assertTrue(facts['winding_consistent'])
                self.assertEqual(facts['components'],1)
                self.assertGreater(facts['signed_volume_model3'],0)
                self.assertGreater(facts['minimum_triangle_area_model2'],1e-12)


if __name__=='__main__': unittest.main()
