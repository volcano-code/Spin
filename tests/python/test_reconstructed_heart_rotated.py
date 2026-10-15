"""Input-construction contract for the non-principal reconstructed Heart probe."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HEART_DIR = ROOT / "experiments/reconstructed_heart_top"
for path in (str(HEART_DIR),):
    if path not in sys.path:
        sys.path.insert(0, path)

import heart_top  # noqa: E402

spec = importlib.util.spec_from_file_location(
    "reconstructed_heart_rotated", HEART_DIR / "run_rotated.py"
)
rotated = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(rotated)


class ReconstructedHeartRotatedInputTests(unittest.TestCase):
    def test_centroid_frame_keeps_target_origin_and_breaks_z_principal_alignment(self):
        vertices, faces = heart_top.heart_mesh()
        moments = heart_top.integrate(vertices, faces)
        centroid = moments[1:4] / moments[0]
        rotated_vertices = rotated._rotate_about_center(vertices, centroid, 25.0) - centroid
        rotated_moments = heart_top.integrate(rotated_vertices, faces)
        self.assertTrue(np.all(np.isfinite(rotated_moments)))
        np.testing.assert_allclose(rotated_moments[1:4], 0.0, atol=1e-12)
        _, _, inertia = heart_top.properties(rotated_moments)
        self.assertGreater(abs(float(inertia[0, 2])), 1e-4)
        # The source Heart is mirror-symmetric in Y; rotation about Y should
        # preserve the expected yz=0 symmetry while making xz nonzero.
        self.assertLess(abs(float(inertia[1, 2])), 1e-12)
        self.assertGreater(np.linalg.eigvalsh(inertia).min(), 0.0)


if __name__ == "__main__":
    unittest.main()
