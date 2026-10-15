import unittest
import numpy as np


class DenseV4Tests(unittest.TestCase):
    def test_admission_requires_inherited_geometry_identity(self):
        import hashlib
        from pathlib import Path
        from noodles_dense_v4_input import SPEC,identity,tube_layout
        from noodles_admission import validate_manifest
        here=Path(__file__).parent
        fixed=['hub-fixed','pin-fixed']+[f'tube-{i:03d}-fixed' for i in range(len(tube_layout()))]
        manifest=dict(spec=SPEC,spec_sha256=identity(),construction=dict(fixed_parts=fixed),
            parts={n:{} for n in ['outer-ellipsoid','inner-ellipsoid-candidate']+fixed},
            creator_sha256=hashlib.sha256((here/'noodles_dense_v4_input.py').read_bytes()).hexdigest(),
            tube_mesh_helper_sha256=hashlib.sha256((here/'noodles_input.py').read_bytes()).hexdigest(),
            inherited_geometry_generator_sha256=hashlib.sha256((here/'noodles_dense_v3_input.py').read_bytes()).hexdigest())
        self.assertEqual(validate_manifest(manifest,'dense-v4'),SPEC)
        manifest['inherited_geometry_generator_sha256']='0'*64
        with self.assertRaisesRegex(ValueError,'INHERITED_GEOMETRY'):
            validate_manifest(manifest,'dense-v4')

    def test_true_core_is_smaller_but_visible_tube_geometry_is_not_relabelled(self):
        from noodles_dense_v4_input import SPEC,components,tube_layout,identity
        from noodles_dense_v3_input import components as v3_components,identity as v3_identity
        frozen_v3=v3_identity();parts=components();old=v3_components()
        self.assertEqual(len(tube_layout()),79)
        for name in parts:
            if name.startswith('tube-') or name in ('hub-fixed','pin-fixed'):
                np.testing.assert_array_equal(parts[name].vertices,old[name].vertices)
                np.testing.assert_array_equal(parts[name].faces,old[name].faces)
        self.assertLess(parts['outer-ellipsoid'].volume,old['outer-ellipsoid'].volume*.3)
        center=np.asarray(SPEC['ellipsoid_center_model'])
        np.testing.assert_allclose(parts['outer-ellipsoid'].bounds,
            [center-np.asarray(SPEC['core_ellipsoid_radii_model']),center+np.asarray(SPEC['core_ellipsoid_radii_model'])])
        self.assertEqual(v3_identity(),frozen_v3)
        self.assertNotEqual(identity(),frozen_v3)

    def test_actual_initial_material_includes_core_and_closed_source_parts(self):
        from noodles_dense_v4_input import components
        from noodles_input import mesh_audit
        parts=components()
        self.assertIn('outer-ellipsoid',parts)
        self.assertIn('inner-ellipsoid-candidate',parts)
        for name,mesh in parts.items():
            with self.subTest(name=name):
                a=mesh_audit(mesh)
                self.assertTrue(a['closed']);self.assertTrue(a['winding_consistent'])
                self.assertEqual(a['components'],1);self.assertGreater(a['minimum_triangle_area_model2'],1e-12)


if __name__=='__main__':unittest.main()
