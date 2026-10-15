from __future__ import annotations

import unittest

from make_it_stand.verification.digital_cut_plan import (
    rasterize_discrete_cut_plane,
    verify_discrete_cut_plan,
    verify_discrete_cut_plane,
)


def _shell_with_voids(*voids: tuple[int, int, int]) -> set[tuple[int, int, int]]:
    material = {
        (x, y, z)
        for x in range(5)
        for y in range(5)
        for z in range(5)
    }
    material.difference_update(voids)
    return material


class DigitalCutPlanTests(unittest.TestCase):
    def test_axis_aligned_plane_rasterizes_exact_intersected_voxel_cells(self) -> None:
        cells = rasterize_discrete_cut_plane(
            (4, 2, 2),
            grid_origin=(0.0, 0.0, 0.0),
            voxel_spacing=(1.0, 1.0, 1.0),
            plane_point=(2.0, 0.0, 0.0),
            plane_normal=(7.0, 0.0, 0.0),
            kerf_width=0.2,
        )
        self.assertEqual(cells, [
            (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1),
            (2, 0, 0), (2, 0, 1), (2, 1, 0), (2, 1, 1),
        ])

    def test_plane_cut_exposes_a_cavity_and_reports_the_explicit_plane(self) -> None:
        dimensions = (5, 5, 5)
        material = _shell_with_voids((2, 2, 2))
        record = verify_discrete_cut_plane(
            dimensions, material,
            grid_origin=(0.0, 0.0, 0.0),
            voxel_spacing=(1.0, 1.0, 1.0),
            plane_point=(2.5, 0.0, 0.0),
            plane_normal=(1.0, 0.0, 0.0),
            kerf_width=0.1,
        )
        self.assertEqual(record["status"], "VERIFIED_DISCRETE_CUT_PLANE")
        self.assertEqual(record["exposed_preexisting_cavity_count"], 1)
        self.assertEqual(record["cut_material_component_count"], 2)
        self.assertEqual(record["cut_plane"]["normalized_normal"], [1.0, 0.0, 0.0])
        self.assertEqual(record["cut_plane"]["kerf_width"], 0.1)
        self.assertEqual(record["continuous_geometry_status"], "UNKNOWN")

    def test_plane_cut_refuses_to_remove_fixed_shell_cells(self) -> None:
        dimensions = (5, 5, 5)
        material = _shell_with_voids((2, 2, 2))
        with self.assertRaisesRegex(ValueError, "fixed material"):
            verify_discrete_cut_plane(
                dimensions, material,
                grid_origin=(0.0, 0.0, 0.0),
                voxel_spacing=(1.0, 1.0, 1.0),
                plane_point=(2.5, 0.0, 0.0),
                plane_normal=(1.0, 0.0, 0.0),
                kerf_width=0.1,
                fixed_material_cells={(2, 1, 1)},
            )

    def test_plane_cut_rejects_degenerate_or_nonfinite_geometry(self) -> None:
        common = dict(
            grid_origin=(0.0, 0.0, 0.0),
            voxel_spacing=(1.0, 1.0, 1.0),
            plane_point=(0.0, 0.0, 0.0),
            kerf_width=0.1,
        )
        for normal in ((0.0, 0.0, 0.0), (float("nan"), 0.0, 1.0)):
            with self.subTest(normal=normal), self.assertRaisesRegex(ValueError, "plane_normal"):
                rasterize_discrete_cut_plane((2, 2, 2), plane_normal=normal, **common)
        with self.assertRaisesRegex(ValueError, "kerf_width"):
            rasterize_discrete_cut_plane(
                (2, 2, 2), plane_normal=(1.0, 0.0, 0.0),
                **{**common, "kerf_width": 0.0},
            )

    def test_boundary_to_cavity_cut_exposes_exactly_that_preexisting_void(self) -> None:
        material = _shell_with_voids((2, 2, 2), (2, 1, 2))
        record = verify_discrete_cut_plan(
            (5, 5, 5), material,
            cut_cells={(0, 2, 2), (1, 2, 2), (2, 2, 2)},
        )
        self.assertEqual(record["status"], "VERIFIED_DISCRETE_CUT_PLAN")
        self.assertEqual(record["preexisting_cavity_count"], 1)
        self.assertEqual(record["exposed_preexisting_cavity_count"], 1)
        self.assertEqual(record["unexposed_preexisting_cavity_count"], 0)
        self.assertEqual(record["cut_material_component_count"], 1)
        self.assertTrue(record["cavities"][0]["exposed"])
        self.assertEqual(record["cavities"][0]["cells"], [[2, 1, 2], [2, 2, 2]])

    def test_cut_that_does_not_reach_a_void_does_not_count_as_exposure(self) -> None:
        material = _shell_with_voids((2, 2, 2))
        record = verify_discrete_cut_plan(
            (5, 5, 5), material,
            cut_cells={(0, 2, 2)},
        )
        self.assertEqual(record["preexisting_cavity_count"], 1)
        self.assertEqual(record["exposed_preexisting_cavity_count"], 0)
        self.assertEqual(record["unexposed_preexisting_cavity_count"], 1)
        self.assertFalse(record["cavities"][0]["exposed"])

    def test_diagonal_empty_cells_do_not_create_a_false_access_path(self) -> None:
        material = _shell_with_voids((2, 2, 2))
        record = verify_discrete_cut_plan(
            (5, 5, 5), material,
            cut_cells={(0, 0, 0), (1, 1, 1)},
        )
        self.assertEqual(record["exposed_preexisting_cavity_count"], 0)

    def test_fixed_material_is_an_immutable_cut_exclusion(self) -> None:
        material = _shell_with_voids((2, 2, 2))
        with self.assertRaisesRegex(ValueError, "fixed material"):
            verify_discrete_cut_plan(
                (5, 5, 5), material,
                cut_cells={(0, 2, 2)}, fixed_material_cells={(0, 2, 2)},
            )

    def test_permutation_of_input_cells_has_the_same_canonical_result(self) -> None:
        material = _shell_with_voids((2, 2, 2))
        first = verify_discrete_cut_plan(
            (5, 5, 5), material,
            cut_cells=[(0, 2, 2), (1, 2, 2), (2, 2, 2)],
        )
        second = verify_discrete_cut_plan(
            (5, 5, 5), list(reversed(sorted(material))),
            cut_cells=[(2, 2, 2), (1, 2, 2), (0, 2, 2)],
        )
        self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
