from __future__ import annotations

import math
import unittest

import numpy as np

import _make_it_stand as native
from make_it_stand.reachability import ReachabilityResult


def cube() -> tuple[np.ndarray, np.ndarray]:
    vertices = np.asarray(
        [
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
        ],
        dtype=np.float64,
    )
    faces = np.asarray(
        [
            [0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7],
            [0, 1, 5], [0, 5, 4], [3, 7, 6], [3, 6, 2],
            [0, 4, 7], [0, 7, 3], [1, 2, 6], [1, 6, 5],
        ],
        dtype=np.int32,
    )
    return vertices, faces


class NativeReachabilityTests(unittest.TestCase):
    def test_binding_matches_strict_schema_and_is_deterministic(self) -> None:
        vertices, faces = cube()
        config = native.BalanceModelConfig()
        config.voxel_resolution = 8
        model = native.BalanceModel(vertices, faces, config)
        model.add_standing_constraint_polygon(
            [[0, 0, 0], [1, 0, 0], [0, 0, 1]],
            [0, -1, 0],
            math.radians(5.0),
        )
        first = ReachabilityResult.from_native(
            model.analyze_carving_reachability(0)
        )
        second = ReachabilityResult.from_native(
            model.analyze_carving_reachability(0)
        )
        self.assertEqual(first.status, "RELAXATION_REACHABLE")
        self.assertGreater(first.variable_cell_count, 0)
        self.assertEqual(first.to_dict(), second.to_dict())


if __name__ == "__main__":
    unittest.main()
