import unittest
import numpy as np


class DenseV3Tests(unittest.TestCase):
    def test_v3_admission_adapter_keeps_frozen_part_identity(self):
        import hashlib
        from pathlib import Path
        from noodles_dense_v3_input import SPEC,identity,tube_layout
        from noodles_admission import validate_manifest
        here=Path(__file__).parent
        fixed=['hub-fixed','pin-fixed']+[f'tube-{i:03d}-fixed' for i in range(len(tube_layout()))]
        manifest=dict(spec=SPEC,spec_sha256=identity(),construction=dict(fixed_parts=fixed),
            parts={n:{} for n in ['outer-ellipsoid','inner-ellipsoid-candidate']+fixed},
            creator_sha256=hashlib.sha256((here/'noodles_dense_v3_input.py').read_bytes()).hexdigest(),
            tube_mesh_helper_sha256=hashlib.sha256((here/'noodles_input.py').read_bytes()).hexdigest())
        self.assertEqual(validate_manifest(manifest,'dense-v3'),SPEC)
        manifest['construction']['fixed_parts']=fixed[:-1]
        with self.assertRaisesRegex(ValueError,'FIXED_PARTS'):
            validate_manifest(manifest,'dense-v3')

    def test_dense_rounded_construction_keeps_actual_hub_connections(self):
        from noodles_dense_v3_input import SPEC,tube_layout,dense_centerline
        self.assertEqual(len(tube_layout()),79)
        for i in range(79):
            p=dense_centerline(i); middle=p[len(p)//2]
            self.assertLess(np.linalg.norm(middle[:2]),SPEC['hub_radius_model'])
            self.assertTrue(SPEC['hub_interval_model'][0]<middle[2]<SPEC['hub_interval_model'][1])
            self.assertTrue(np.all(np.diff(p[:,0])>0))
        body_min=SPEC['ellipsoid_center_model'][2]-SPEC['ellipsoid_radii_model'][2]
        body_max=SPEC['ellipsoid_center_model'][2]+SPEC['ellipsoid_radii_model'][2]
        self.assertGreaterEqual(SPEC['hub_interval_model'][0],body_min)
        self.assertLessEqual(SPEC['hub_interval_model'][1],body_max)

    def test_frozen_identity_does_not_modify_v2(self):
        from noodles_dense_input import identity as old_identity
        from noodles_dense_v3_input import identity,SPEC
        before=old_identity(); frozen=identity()
        original=SPEC['tube_radius_model']
        try:
            SPEC['tube_radius_model']=original+0.001
            self.assertNotEqual(identity(),frozen)
        finally:SPEC['tube_radius_model']=original
        self.assertEqual(old_identity(),before)

    def test_closed_rounded_parts_are_not_just_overlapping_visual_meshes(self):
        from noodles_dense_v3_input import components
        from noodles_input import mesh_audit
        parts=components()
        for name,mesh in parts.items():
            with self.subTest(name=name):
                f=mesh_audit(mesh)
                self.assertTrue(f['closed']);self.assertTrue(f['winding_consistent'])
                self.assertEqual(f['components'],1)
                self.assertGreater(f['minimum_triangle_area_model2'],1e-12)
                self.assertGreater(f['signed_volume_model3'],0)


if __name__=='__main__':unittest.main()
