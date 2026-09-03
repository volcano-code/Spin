import json
import tempfile
import unittest

from make_it_stand import (
    SpinIt2014Request,
    solve_spin_it_2014,
    solve_spin_it_2014_cage,
)


class SpinIt2014PythonTests(unittest.TestCase):
    def test_cage_api_rejects_non_artifact_before_native_execution(self) -> None:
        request = SpinIt2014Request(
            dimensions=(1, 1, 1), minimum=(0.0, 0.0, 0.0), maximum=(1.0, 1.0, 1.0),
            material_mask=(1,), carveable_mask=(1,),
        )
        with self.assertRaises(TypeError):
            solve_spin_it_2014_cage(object(), request)  # type: ignore[arg-type]

    def test_binary_constraint_tolerance_must_be_positive(self) -> None:
        with self.assertRaises(ValueError):
            SpinIt2014Request(
                dimensions=(1, 1, 1), minimum=(0.0, 0.0, 0.0), maximum=(1.0, 1.0, 1.0),
                material_mask=(1,), carveable_mask=(1,),
                binary_constraint_tolerance=0.0,
            )

    def test_slqp_iteration_limit_must_be_positive(self) -> None:
        with self.assertRaises(ValueError):
            SpinIt2014Request(
                dimensions=(1, 1, 1), minimum=(0.0, 0.0, 0.0), maximum=(1.0, 1.0, 1.0),
                material_mask=(1,), carveable_mask=(1,),
                slqp_maximum_iterations=0,
            )

    def test_slqp_gradient_tolerance_must_be_positive(self) -> None:
        with self.assertRaises(ValueError):
            SpinIt2014Request(
                dimensions=(1, 1, 1), minimum=(0.0, 0.0, 0.0), maximum=(1.0, 1.0, 1.0),
                material_mask=(1,), carveable_mask=(1,),
                slqp_projected_gradient_tolerance=0.0,
            )

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
