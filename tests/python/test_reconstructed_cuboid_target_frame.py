"""Contracts for the frozen target-frame construction in the Cuboid runner."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "experiments/reconstructed_cuboid_yoyo/run.py"
SPEC = importlib.util.spec_from_file_location("reconstructed_cuboid_yoyo_runner", RUNNER)
assert SPEC is not None and SPEC.loader is not None
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)


class ReconstructedCuboidTargetFrameTests(unittest.TestCase):
    def test_nonprincipal_axis_has_deterministic_right_handed_frame(self):
        axis = np.array([1.0, 0.0, 1.0]) / np.sqrt(2.0)
        frame = np.asarray(runner.target_frame_for_axis(axis), dtype=float)

        self.assertEqual(frame.shape, (3, 3))
        np.testing.assert_allclose(frame.T @ frame, np.eye(3), rtol=0.0, atol=1e-12)
        self.assertAlmostEqual(float(np.linalg.det(frame)), 1.0, places=12)
        np.testing.assert_allclose(frame[:, 2], axis, rtol=0.0, atol=1e-12)
        np.testing.assert_allclose(
            frame,
            np.asarray(runner.target_frame_for_axis(axis), dtype=float),
            rtol=0.0,
            atol=0.0,
        )


if __name__ == "__main__":
    unittest.main()
