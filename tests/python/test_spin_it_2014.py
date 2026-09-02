import json
import tempfile
import unittest

from make_it_stand import SpinIt2014Request, solve_spin_it_2014


class SpinIt2014PythonTests(unittest.TestCase):
    def test_sealed_deterministic_static_artifact_and_viewer_bundle(self) -> None:
        dimensions = (4, 4, 4)
        material = [1] * 64
        carveable = [0] * 64
        for z in range(1, 3):
            for y in range(1, 3):
                for x in range(1, 3):
                    carveable[x + 4 * (y + 4 * z)] = 1
        request = SpinIt2014Request(
            dimensions=dimensions,
            minimum=(-1.0, -1.0, 0.0),
            maximum=(1.0, 1.0, 2.0),
            material_mask=tuple(material),
            carveable_mask=tuple(carveable),
            mode="STATIC_BALANCE",
            initial_level=1,
            maximum_level=2,
            maximum_adaptive_passes=3,
        )
        first = solve_spin_it_2014(request)
        second = solve_spin_it_2014(request)
        self.assertTrue(first.succeeded)
        self.assertTrue(first.geometry_sealed)
        self.assertEqual(first.diagnostics["output_hash"], second.diagnostics["output_hash"])
        self.assertFalse(first.outer_vertices.flags.writeable)
        with tempfile.TemporaryDirectory() as directory:
            paths = first.save_viewer_bundle(directory)
            self.assertEqual(len(paths), 5)
            with open(paths[-1], encoding="utf-8") as stream:
                self.assertEqual(json.load(stream)["schema"], "spin-it-2014-artifact-v1")


if __name__ == "__main__":
    unittest.main()
