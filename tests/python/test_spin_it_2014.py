import json
import tempfile
import unittest

import make_it_stand.spin2014 as spin2014_module
from make_it_stand import (
    SpinIt2014Request,
    solve_spin_it_2014,
    solve_spin_it_2014_cage,
)


class SpinIt2014PythonTests(unittest.TestCase):
    def test_objective_selector_must_match_mode(self) -> None:
        self.assertEqual(
            getattr(spin2014_module, "SPIN_2014_OBJECTIVES", None),
            frozenset({
                "YOYO",
                "TOP_INERTIA_ONLY",
                "TOP_LOWERING_ONLY",
                "TOP_MASS_AND_LOWERING",
                "STATIC_LOWERING",
            }),
        )
        common = dict(
            dimensions=(1, 1, 1),
            minimum=(0.0, 0.0, 0.0),
            maximum=(1.0, 1.0, 1.0),
            material_mask=(1,),
            carveable_mask=(1,),
        )
        with self.assertRaises(ValueError):
            SpinIt2014Request(
                **common,
                mode="STATIC_BALANCE",
                objective="TOP_INERTIA_ONLY",
            )
        with self.assertRaises(ValueError):
            SpinIt2014Request(
                **common,
                mode="YO_YO",
                objective="TOP_LOWERING_ONLY",
            )

    def test_objective_weights_must_be_finite_and_nonnegative(self) -> None:
        common = dict(
            dimensions=(1, 1, 1),
            minimum=(0.0, 0.0, 0.0),
            maximum=(1.0, 1.0, 1.0),
            material_mask=(1,),
            carveable_mask=(1,),
        )
        for field, value in (
            ("inertia_weight", -1.0),
            ("inertia_weight", float("nan")),
            ("center_of_mass_weight", float("inf")),
        ):
            with self.subTest(field=field, value=value):
                with self.assertRaises(ValueError):
                    SpinIt2014Request(**common, **{field: value})

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
            objective="STATIC_LOWERING",
            inertia_weight=2.5,
            initial_level=1,
            maximum_level=2,
            maximum_adaptive_passes=3,
        )
        first = solve_spin_it_2014(request)
        second = solve_spin_it_2014(request)
        self.assertTrue(first.succeeded)
        self.assertTrue(first.geometry_sealed)
        self.assertEqual(first.diagnostics["objective"], "STATIC_LOWERING")
        self.assertEqual(first.diagnostics["objective_resolution"], "EXPLICIT")
        self.assertEqual(first.diagnostics["inertia_weight"], 2.5)
        self.assertEqual(first.diagnostics["center_of_mass_weight"], 1.0)
        self.assertEqual(
            first.diagnostics["objective_semantics"],
            "ANGLE_ELIMINATED_TRANSVERSE_INERTIA",
        )
        self.assertFalse(first.diagnostics["original_phi_trajectory"])
        self.assertEqual(first.diagnostics["output_hash"], second.diagnostics["output_hash"])
        self.assertFalse(first.outer_vertices.flags.writeable)
        with tempfile.TemporaryDirectory() as directory:
            paths = first.save_viewer_bundle(directory)
            self.assertEqual(len(paths), 5)
            with open(paths[-1], encoding="utf-8") as stream:
                self.assertEqual(json.load(stream)["schema"], "spin-it-2014-artifact-v1")


if __name__ == "__main__":
    unittest.main()
