"""Geometry and identity tests; no successful whole-object admission implied."""
import hashlib
import json
import unittest
import numpy as np
from noodles_input import SPEC, centerline, components, mesh_audit, spec_identity


class NoodlesInputTests(unittest.TestCase):
    def test_each_part_closed_oriented_positive_and_nondegenerate(self):
        for name, mesh in components().items():
            with self.subTest(part=name):
                audit = mesh_audit(mesh)
                self.assertTrue(audit['closed'])
                self.assertTrue(audit['winding_consistent'])
                self.assertEqual(audit['components'], 1)
                self.assertGreater(audit['signed_volume_model3'], 0)
                self.assertGreater(audit['minimum_triangle_area_model2'], 1e-12)

    def test_every_tube_has_both_end_disks_inside_central_ellipsoid(self):
        center = np.asarray(SPEC['ellipsoid_center_model'])
        axes = np.asarray(SPEC['ellipsoid_radii_model'])
        # Bounding each endpoint disk by a sphere gives a conservative overlap witness.
        for i in range(5):
            points = centerline(i)[[0, -1]]
            bound = np.linalg.norm((points-center)/axes, axis=1)+SPEC['tube_radius_model']/axes.min()
            self.assertTrue(np.all(bound < 1), (i, bound))

    def test_curve_nonlocal_sampled_separation_and_local_bend(self):
        # A screening check only: exact self-intersection remains an admission gate.
        for i in range(5):
            p = centerline(i)
            distances = np.linalg.norm(p[:, None, :]-p[None, :, :], axis=-1)
            separation = np.abs(np.arange(len(p))[:, None]-np.arange(len(p))[None, :])
            self.assertGreater(distances[separation >= 8].min(), 2*SPEC['tube_radius_model'])
            edge = np.diff(p, axis=0)
            tangent = edge/np.linalg.norm(edge, axis=1)[:, None]
            turn = np.arccos(np.clip(np.sum(tangent[1:]*tangent[:-1], axis=1), -1, 1))
            radius = np.minimum(np.linalg.norm(edge[:-1], axis=1), np.linalg.norm(edge[1:], axis=1))/(2*np.sin(turn/2))
            self.assertGreater(radius.min(), SPEC['tube_radius_model'])

    def test_spec_identity_is_deterministic_and_tracks_changes(self):
        self.assertEqual(spec_identity(), spec_identity())
        changed = dict(SPEC, density_kg_m3=SPEC['density_kg_m3']+1)
        digest = hashlib.sha256(json.dumps(changed, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
        self.assertNotEqual(digest, spec_identity())
        a, b = components(), components()
        for name in a:
            np.testing.assert_array_equal(a[name].vertices, b[name].vertices)
            np.testing.assert_array_equal(a[name].faces, b[name].faces)


if __name__ == '__main__': unittest.main()
