"""Prevent the diagnosed tube-radius/curvature fold without changing layout."""
import unittest
from pathlib import Path
import json
import hashlib
import numpy as np
import noodles_dense_v3_input as base


class GentleSweepTests(unittest.TestCase):
    def generator(self):
        try:
            import noodles_dense_v5_input as candidate
        except ImportError:
            candidate = base
        return candidate

    def curvature_product(self, points):
        first = np.gradient(points, axis=0)
        second = np.gradient(first, axis=0)
        curvature = np.linalg.norm(np.cross(first, second), axis=1) / np.linalg.norm(first, axis=1)**3
        return float(curvature[2:-2].max() * base.SPEC['tube_radius_model'])

    def test_diagnosed_tube003_does_not_fold_around_local_curvature(self):
        self.assertLessEqual(self.curvature_product(self.generator().dense_centerline(3)), .45)

    def test_every_sweep_preserves_layout_and_uses_safe_curvature(self):
        candidate = self.generator()
        self.assertEqual(candidate.tube_layout(), base.tube_layout())
        self.assertEqual(candidate.SPEC['tube_radius_model'], base.SPEC['tube_radius_model'])
        for index in range(len(base.tube_layout())):
            points = candidate.dense_centerline(index)
            self.assertTrue(np.all(np.diff(points[:, 0]) > 0))
            self.assertLessEqual(self.curvature_product(points), .45, msg=str(index))

    def test_admission_authenticates_v5_generator_and_rejects_changed_spec(self):
        from copy import deepcopy
        from noodles_admission import validate_manifest
        import noodles_dense_v5_input as candidate
        here=Path(__file__).resolve().parent
        record=json.loads((here/'noodles-dense-v3/manifest.json').read_text())
        record.update(spec=candidate.SPEC,spec_sha256=candidate.identity(),
            creator_sha256=hashlib.sha256(Path(candidate.__file__).read_bytes()).hexdigest(),
            inherited_geometry_generator_sha256=hashlib.sha256(Path(base.__file__).read_bytes()).hexdigest())
        self.assertEqual(validate_manifest(record,'dense-v5'),candidate.SPEC)
        changed=deepcopy(record)
        changed['spec']['tube_radius_model']*=2
        with self.assertRaisesRegex(ValueError,'STALE_FROZEN_SPEC'):
            validate_manifest(changed,'dense-v5')


if __name__ == '__main__':
    unittest.main()
