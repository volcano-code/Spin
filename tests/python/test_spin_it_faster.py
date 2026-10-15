from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np
import make_it_stand.spin as spin_api

from make_it_stand.spin import (
    SpinEvaluatorConfig,
    SpinItFasterRequest,
    SpinSolverConfig,
    compute_hydrostatic_geometry,
    find_hydrostatic_equilibrium,
    evaluate_spin_it_faster_quadric,
    evaluate_spin_it_faster_fullrank_quadric,
    integrate_raw_mass_moments,
    solve_spin_it_faster,
    solve_spin_it_faster_planar_fidget,
    solve_spin_it_faster_fullrank_quadric,
    solve_spin_it_faster_with_external_moment_evaluator,
    with_measured_fixed_shell_geometry,
)


def cube(minimum: float, maximum: float):
    vertices = (
        (minimum, minimum, minimum),
        (maximum, minimum, minimum),
        (maximum, maximum, minimum),
        (minimum, maximum, minimum),
        (minimum, minimum, maximum),
        (maximum, minimum, maximum),
        (maximum, maximum, maximum),
        (minimum, maximum, maximum),
    )
    faces = (
        (0, 2, 1), (0, 3, 2),
        (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4),
        (3, 7, 6), (3, 6, 2),
        (0, 4, 7), (0, 7, 3),
        (1, 2, 6), (1, 6, 5),
    )
    return vertices, faces


def box(lower, upper):
    x0, y0, z0 = lower
    x1, y1, z1 = upper
    vertices = (
        (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
        (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
    )
    faces = (
        (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4), (3, 7, 6), (3, 6, 2),
        (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5),
    )
    return vertices, faces


def combine_meshes(*meshes):
    vertices = []
    faces = []
    offset = 0
    for mesh_vertices, mesh_faces in meshes:
        vertices.extend(mesh_vertices)
        faces.extend(tuple(index + offset for index in face) for face in mesh_faces)
        offset += len(mesh_vertices)
    return tuple(vertices), tuple(faces)


class SpinSchemaTests(unittest.TestCase):
    def test_fixed_shell_adapter_integrates_real_mesh_at_its_declared_density(self) -> None:
        outer_vertices, outer_faces = cube(-5.0, 5.0)
        design_vertices, design_faces = cube(-1.0, 1.0)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            material_density=1.2,
        )
        fixed_vertices, fixed_faces = cube(-0.5, 0.5)
        transformed = tuple(
            (2.0*x + 0.4*y + 0.3, 0.25*x + 1.5*y - 0.2, 0.5*x + 0.75*y + 1.1*z + 0.1)
            for x, y, z in fixed_vertices
        )

        measured = with_measured_fixed_shell_geometry(
            request, transformed, fixed_faces, density=0.7,
        )
        expected = integrate_raw_mass_moments(transformed, fixed_faces, density=0.7)
        np.testing.assert_allclose(measured.fixed_shell_moments, expected, rtol=0.0, atol=0.0)
        self.assertAlmostEqual(measured.fixed_shell_moments[0], 0.7 * (2.0 * 1.5 - 0.4 * 0.25) * 1.1)
        self.assertNotAlmostEqual(measured.fixed_shell_moments[0], 1.2 * (2.0 * 1.5 - 0.4 * 0.25) * 1.1)
        self.assertNotEqual(request.fixed_shell_moments, measured.fixed_shell_moments)
        same_material = with_measured_fixed_shell_geometry(request, transformed, fixed_faces)
        self.assertAlmostEqual(same_material.fixed_shell_moments[0], 1.2 * (2.0 * 1.5 - 0.4 * 0.25) * 1.1)
        with self.assertRaisesRegex(ValueError, "density"):
            with_measured_fixed_shell_geometry(request, transformed, fixed_faces, density=0.0)

    def test_hydrostatic_geometry_is_measured_from_closed_mesh_not_hand_filled(self) -> None:
        vertices, faces = cube(-1.0, 1.0)
        measured = compute_hydrostatic_geometry(
            vertices, faces, waterline_origin=(0.0, 0.0, 0.0),
            waterline_normal=(0.0, 0.0, 1.0),
        )
        self.assertAlmostEqual(measured.displaced_volume, 4.0, places=10)
        np.testing.assert_allclose(measured.center_of_buoyancy, (0.0, 0.0, -0.5), atol=1e-10)
        self.assertAlmostEqual(measured.waterline_area, 4.0, places=10)
        np.testing.assert_allclose(measured.waterline_centroid, (0.0, 0.0, 0.0), atol=1e-10)
        np.testing.assert_allclose(
            measured.waterline_central_second_moment,
            ((4.0 / 3.0, 0.0), (0.0, 4.0 / 3.0)), atol=1e-10,
        )
        with self.assertRaisesRegex(ValueError, "strictly intersect"):
            compute_hydrostatic_geometry(
                vertices, faces, waterline_origin=(0.0, 0.0, 2.0),
                waterline_normal=(0.0, 0.0, 1.0),
            )

    def test_hydrostatic_equilibrium_solves_only_closed_boundary_volume(self) -> None:
        vertices, faces = cube(-1.0, 1.0)
        equilibrium = find_hydrostatic_equilibrium(
            vertices, faces, reference_origin=(0.0, 0.0, 0.0),
            waterline_normal=(0.0, 0.0, 1.0), target_displaced_volume=6.0,
        )
        self.assertAlmostEqual(equilibrium.waterline_signed_offset, 0.5, places=9)
        self.assertAlmostEqual(equilibrium.geometry.displaced_volume, 6.0, places=8)
        self.assertLessEqual(abs(equilibrium.volume_residual), 8e-10)
        with self.assertRaisesRegex(ValueError, "strictly interior"):
            find_hydrostatic_equilibrium(
                vertices, faces, reference_origin=(0.0, 0.0, 0.0),
                waterline_normal=(0.0, 0.0, 1.0), target_displaced_volume=8.0,
            )

    def test_buoyant_request_is_built_from_measured_closed_geometries(self) -> None:
        self.assertTrue(hasattr(spin_api, "with_measured_buoyancy_geometry"))
        outer_vertices, outer_faces = cube(-5.0, 5.0)
        design_vertices, design_faces = cube(-1.0, 1.0)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            material_density=1.5,
        )
        fixed_vertices, fixed_faces = cube(-0.5, 0.5)
        displacement_vertices, displacement_faces = cube(-1.0, 1.0)

        measured = spin_api.with_measured_buoyancy_geometry(
            request,
            displacement_vertices,
            displacement_faces,
            fixed_vertices,
            fixed_faces,
            waterline_origin=(0.0, 0.0, 0.0),
            waterline_normal=(0.0, 0.0, 1.0),
            liquid_density=1.0,
        )

        self.assertEqual(measured.mode, "BUOYANT_STABILITY")
        self.assertAlmostEqual(measured.displaced_liquid_volume, 4.0, places=10)
        np.testing.assert_allclose(measured.center_of_buoyancy_target, (0.0, 0.0, -0.5), atol=1e-10)
        self.assertAlmostEqual(measured.waterline_area, 4.0, places=10)
        np.testing.assert_allclose(
            measured.waterline_central_second_moment,
            ((4.0 / 3.0, 0.0), (0.0, 4.0 / 3.0)), atol=1e-10,
        )
        self.assertAlmostEqual(measured.fixed_shell_volume, 1.0, places=10)
        self.assertAlmostEqual(measured.fixed_shell_moments[0], 1.5, places=10)
        self.assertEqual(request.mode, "STATIC_BALANCE")

    def test_measured_buoyant_fixed_material_uses_its_own_density(self) -> None:
        outer_vertices, outer_faces = cube(-5.0, 5.0)
        design_vertices, design_faces = cube(-1.0, 1.0)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            material_density=1.5,
        )
        fixed_vertices, fixed_faces = cube(-0.5, 0.5)
        displacement_vertices, displacement_faces = cube(-1.0, 1.0)

        measured = spin_api.with_measured_buoyancy_geometry(
            request,
            displacement_vertices,
            displacement_faces,
            fixed_vertices,
            fixed_faces,
            waterline_origin=(0.0, 0.0, 0.0),
            waterline_normal=(0.0, 0.0, 1.0),
            liquid_density=1.0,
            fixed_shell_density=0.5,
        )

        self.assertAlmostEqual(measured.fixed_shell_volume, 1.0, places=10)
        self.assertAlmostEqual(measured.fixed_shell_moments[0], 0.5, places=10)

    def test_buoyant_solver_accepts_measured_heterogeneous_shell_and_rejects_wrong_density(self) -> None:
        outer_vertices, outer_faces = cube(-2.0, 2.0)
        design_vertices, design_faces = cube(-1.9, 1.9)
        design_faces_array = np.asarray(design_faces, dtype=np.int32)
        fixed_shell = combine_meshes(
            (outer_vertices, outer_faces),
            (design_vertices, design_faces_array[:, ::-1]),
        )
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            material_density=1.0,
        )
        measured = spin_api.with_measured_buoyancy_geometry(
            request,
            outer_vertices,
            outer_faces,
            np.asarray(fixed_shell[0]).tolist(),
            np.asarray(fixed_shell[1]).tolist(),
            waterline_origin=(0.0, 0.0, 0.0),
            waterline_normal=(0.0, 0.0, 1.0),
            liquid_density=1.0,
            fixed_shell_density=0.1,
        )

        solved = spin_api.solve_spin_it_faster_buoyant_plane(
            measured,
            initial_parameters=(0.0, 0.0, 0.0),
            maximum_iterations=60,
            constraint_tolerance=1e-10,
        )
        self.assertEqual(solved["status"], "CONVERGED", solved["diagnostic"])
        self.assertTrue(solved["equilibrium_constraints_verified"])
        self.assertTrue(solved["stable_at_prescribed_waterline"])
        self.assertAlmostEqual(measured.fixed_shell_moments[0], 0.9128, places=10)

        wrong_density = replace(measured, fixed_shell_density=1.0)
        rejected = spin_api.solve_spin_it_faster_buoyant_plane(wrong_density)
        self.assertEqual(rejected["status"], "INVALID_REQUEST")
        self.assertIn("fixed-shell mass is inconsistent", rejected["diagnostic"])

    def test_buoyant_geometry_adapter_rejects_nonvertical_waterline(self) -> None:
        outer_vertices, outer_faces = cube(-5.0, 5.0)
        design_vertices, design_faces = cube(-1.0, 1.0)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            material_density=1.0,
        )
        vertices, faces = cube(-1.0, 1.0)
        self.assertTrue(hasattr(spin_api, "with_measured_buoyancy_geometry"))
        with self.assertRaisesRegex(ValueError, "target-frame vertical"):
            spin_api.with_measured_buoyancy_geometry(
                request,
                vertices,
                faces,
                vertices,
                faces,
                waterline_origin=(0.0, 0.0, 0.0),
                waterline_normal=(1.0, 0.0, 0.0),
                liquid_density=1.0,
            )

    def test_buoyant_geometry_adapter_transforms_hydro_but_keeps_fixed_model_moments(self) -> None:
        local_outer, outer_faces = cube(-5.0, 5.0)
        local_design, design_faces = cube(-1.0, 1.0)
        axes = np.asarray(((1.0, 0.0, 0.0), (0.0, 0.0, -1.0), (0.0, 1.0, 0.0)))
        origin = np.asarray((3.0, -2.0, 4.0))

        def to_world(vertices):
            return tuple(tuple(row) for row in (np.asarray(vertices) @ axes.T + origin))

        request = SpinItFasterRequest(
            outer_vertices=to_world(local_outer),
            outer_faces=outer_faces,
            design_vertices=to_world(local_design),
            design_faces=design_faces,
            target_origin=tuple(origin),
            target_axes=tuple(tuple(row) for row in axes),
            material_density=1.0,
        )
        fixed_vertices, fixed_faces = cube(-0.5, 0.5)
        world_fixed = to_world(fixed_vertices)
        displacement_vertices, displacement_faces = cube(-1.0, 1.0)
        measured = spin_api.with_measured_buoyancy_geometry(
            request,
            to_world(displacement_vertices),
            displacement_faces,
            world_fixed,
            fixed_faces,
            waterline_origin=tuple(origin),
            waterline_normal=tuple(axes[:, 2]),
            liquid_density=1.0,
        )

        np.testing.assert_allclose(measured.center_of_buoyancy_target, (0.0, 0.0, -0.5), atol=1e-10)
        np.testing.assert_allclose(
            measured.fixed_shell_moments[1:4], origin, atol=1e-10,
        )
        # Fixed request moments are MODEL-frame; the actual evaluator applies
        # model→target once. Hydrostatic fields above are already target-frame.
        solved = spin_api.solve_spin_it_faster_buoyant_plane(measured)
        self.assertEqual(solved['status'], 'CONVERGED', solved['diagnostic'])
        np.testing.assert_allclose(solved['total_raw_moments'][1:3], (0., 0.), atol=1e-10)

    def test_buoyant_geometry_adapter_rejects_boolean_liquid_density(self) -> None:
        outer_vertices, outer_faces = cube(-5.0, 5.0)
        design_vertices, design_faces = cube(-1.0, 1.0)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            material_density=1.0,
        )
        vertices, faces = cube(-1.0, 1.0)
        with self.assertRaisesRegex(ValueError, "liquid_density"):
            spin_api.with_measured_buoyancy_geometry(
                request,
                vertices,
                faces,
                vertices,
                faces,
                waterline_origin=(0.0, 0.0, 0.0),
                waterline_normal=(0.0, 0.0, 1.0),
                liquid_density=True,
            )

    def test_measured_buoyant_input_runs_the_prescribed_plane_solver(self) -> None:
        self.assertTrue(hasattr(spin_api, "solve_spin_it_faster_buoyant_plane"))
        outer_vertices, outer_faces = box((-3.0, -3.0, -2.0), (4.0, 3.0, 4.0))
        design_vertices, design_faces = box((-1.0, -1.0, -1.0), (2.0, 1.0, 1.0))
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            material_density=1.0,
        )
        fixed = combine_meshes(
            box((-1.0, 1.0, -1.0), (2.0, 1.125, 1.0)),
            box((-1.0, -1.125, -1.0), (2.0, -1.0, 1.0)),
        )
        displacement = box((-1.0, -2.0, -1.0), (2.0, 2.0, 1.0))
        measured = spin_api.with_measured_buoyancy_geometry(
            request,
            *displacement,
            *fixed,
            waterline_origin=(0.0, 0.0, 0.0),
            waterline_normal=(0.0, 0.0, 1.0),
            liquid_density=0.625,
        )

        solved = spin_api.solve_spin_it_faster_buoyant_plane(
            measured,
            initial_parameters=(0.1, 0.2, -0.1),
            maximum_iterations=40,
            constraint_tolerance=1e-11,
        )
        self.assertEqual(solved["status"], "CONVERGED")
        self.assertTrue(solved["equilibrium_constraints_verified"])
        self.assertTrue(solved["eq7_numerical_conditions_verified"])
        self.assertGreater(solved["accepted_steps"], 0)
        self.assertLessEqual(solved["normalized_constraints_norm"], 1e-11)
        np.testing.assert_allclose(solved["plane_coefficients"], (0.0, 0.0, 0.0, -1.0), atol=1e-9)
        self.assertEqual(np.asarray(solved["variable_material_faces"]).shape[1], 3)
        total = np.asarray(solved["total_raw_moments"])
        self.assertAlmostEqual(total[0], 7.5, places=9)
        self.assertAlmostEqual(total[1], 3.75, places=9)
        self.assertAlmostEqual(total[2], 0.0, places=9)

    def test_request_is_frozen_and_rejects_unknown_mode(self) -> None:
        outer_vertices, outer_faces = cube(-5.0, 5.0)
        design_vertices, design_faces = cube(-1.0, 1.0)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            minimum_wall_thickness=4.0,
        )
        with self.assertRaises(FrozenInstanceError):
            request.mode = "YO_YO"  # type: ignore[misc]
        with self.assertRaisesRegex(ValueError, "unsupported Spin-It Faster mode"):
            replace(request, mode="UNKNOWN")

    def test_buoyant_request_normalizes_public_paper_inputs(self) -> None:
        outer_vertices, outer_faces = cube(-5.0, 5.0)
        design_vertices, design_faces = cube(-1.0, 1.0)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            mode="BUOYANT_STABILITY",
            material_density=0.5,
            minimum_wall_thickness=4.0,
            liquid_density=1.0,
            displaced_liquid_volume=1.5,
            center_of_buoyancy_target=(0.1, -0.2, 0.0),
            waterline_area=2.0,
            waterline_central_second_moment=((0.4, 0.05), (0.05, 0.3)),
            fixed_shell_volume=0.5,
        )
        self.assertEqual(request.mode, "BUOYANT_STABILITY")
        self.assertEqual(request.center_of_buoyancy_target, (0.1, -0.2, 0.0))
        with self.assertRaisesRegex(ValueError, "liquid_density"):
            replace(request, liquid_density=0.0)
        with self.assertRaisesRegex(ValueError, "symmetric"):
            replace(
                request,
                waterline_central_second_moment=((0.4, 0.3), (0.05, 0.3)),
            )

    def test_static_pareto_weight_is_mode_checked(self) -> None:
        outer_vertices, outer_faces = cube(-5.0, 5.0)
        design_vertices, design_faces = cube(-1.0, 1.0)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            static_objective="PARETO_TRADEOFF",
            static_pareto_weight=0.25,
            minimum_wall_thickness=4.0,
        )
        self.assertEqual(request.static_pareto_weight, 0.25)
        with self.assertRaisesRegex(ValueError, "static_pareto_weight"):
            replace(request, static_pareto_weight=1.1)
        with self.assertRaisesRegex(ValueError, "only.*PARETO_TRADEOFF"):
            replace(request, static_objective="MAXIMUM_MASS")

    def test_static_plane_family_public_path_solves_all_objectives_from_measured_material(self) -> None:
        """Public static dispatch must solve each paper objective on real fixed geometry."""
        self.assertTrue(
            hasattr(spin_api, "solve_spin_it_faster_static_plane_family"),
            "the static paper plane-family solver is not exposed through the public facade",
        )
        outer_vertices, outer_faces = box((-2.0, -2.0, -2.0), (2.0, 2.0, 2.0))
        design_vertices, design_faces = box((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0))
        fixed_geometry = combine_meshes(
            box((-2.0, -2.0, -2.0), (-1.0, 2.0, 2.0)),
            box((1.0, -2.0, -2.0), (2.0, 2.0, 2.0)),
            box((-1.0, -2.0, -2.0), (1.0, -1.0, 2.0)),
            box((-1.0, 1.0, -2.0), (1.0, 2.0, 2.0)),
            box((-1.0, -1.0, -2.0), (1.0, 1.0, -1.0)),
            box((-1.0, -1.0, 1.0), (1.0, 1.0, 2.0)),
        )
        base = with_measured_fixed_shell_geometry(
            SpinItFasterRequest(
                outer_vertices=outer_vertices,
                outer_faces=outer_faces,
                design_vertices=design_vertices,
                design_faces=design_faces,
                minimum_wall_thickness=0.0,
            ),
            *fixed_geometry,
        )
        lower = spin_api.solve_spin_it_faster_static_plane_family(
            replace(base, static_objective="LOWER_CENTER_OF_MASS"),
            initial_coefficients=(0.0, 0.0, 0.0, -1.0),
            maximum_iterations=40,
            constraint_tolerance=2e-10,
        )
        self.assertEqual(lower["status"], "CONVERGED")
        self.assertTrue(lower["paper_optimality_conditions_verified"])
        self.assertGreater(lower["accepted_steps"], 0)
        self.assertLessEqual(lower["normalized_constraints_norm"], 2e-10)

        cases = (
            ("MAXIMUM_MASS", 0.0, (1.0, -5.0, 0.0, 0.0)),
            ("PARETO_TRADEOFF", 1.0, (1.0, -5.0, 0.0, 0.0)),
            ("PARETO_TRADEOFF", 0.0, tuple(lower["plane_coefficients"])),
        )
        for objective, weight, initial in cases:
            with self.subTest(objective=objective, weight=weight):
                request = replace(
                    base,
                    static_objective=objective,
                    static_pareto_weight=weight,
                )
                solved = spin_api.solve_spin_it_faster_static_plane_family(
                    request,
                    initial_coefficients=initial,
                    maximum_iterations=40,
                    constraint_tolerance=2e-10,
                )
                self.assertEqual(solved["status"], "CONVERGED", solved["diagnostic"])
                self.assertTrue(solved["paper_optimality_conditions_verified"])
                self.assertLessEqual(solved["normalized_constraints_norm"], 2e-10)
                vertices = np.asarray(solved["variable_material_vertices"], dtype=np.float64)
                faces = np.asarray(solved["variable_material_faces"], dtype=np.int64)
                self.assertGreater(len(faces), 0)
                variable = integrate_raw_mass_moments(
                    vertices, faces.tolist(), density=request.material_density,
                )
                independent_total = variable + np.asarray(request.fixed_shell_moments)
                np.testing.assert_allclose(
                    solved["total_raw_moments"], independent_total,
                    rtol=1e-10, atol=1e-10,
                )

        full_design = spin_api.solve_spin_it_faster_static_plane_family(
            replace(base, static_objective="MAXIMUM_MASS"),
            initial_coefficients=(1.0, 0.0, 0.0, 0.0),
            maximum_iterations=40,
            constraint_tolerance=2e-10,
        )
        self.assertEqual(full_design["status"], "CONVERGED", full_design["diagnostic"])
        self.assertTrue(full_design["paper_optimality_conditions_verified"])
        self.assertEqual(full_design["accepted_steps"], 0)
        self.assertAlmostEqual(full_design["variable_raw_moments"][0], 8.0, places=10)
        self.assertLessEqual(full_design["normalized_constraints_norm"], 2e-10)
        off_center_fixed = box((-0.5, 1.2, -0.5), (0.5, 1.7, 0.5))
        unbalanced_base = with_measured_fixed_shell_geometry(
            replace(base, fixed_shell_moments=(0.0,) * 10), *off_center_fixed,
        )
        unbalanced_full = spin_api.solve_spin_it_faster_static_plane_family(
            replace(unbalanced_base, static_objective="MAXIMUM_MASS"),
            initial_coefficients=(1.0, 0.0, 0.0, 0.0),
            maximum_iterations=40,
            constraint_tolerance=2e-10,
        )
        self.assertEqual(unbalanced_full["status"], "DERIVATIVE_UNAVAILABLE")
        self.assertFalse(unbalanced_full["paper_optimality_conditions_verified"])
        self.assertGreater(unbalanced_full["normalized_constraints_norm"], 2e-10)

        # This axis-aligned cube fixture has a reproducible intermediate
        # Pareto trial whose Boolean boundary develops zero-area contact faces.
        # The public adapter must preserve that fail-closed outcome rather
        # than promote it to a successful material result.
        interior = spin_api.solve_spin_it_faster_static_plane_family(
            replace(base, static_objective="PARETO_TRADEOFF", static_pareto_weight=0.01),
            initial_coefficients=tuple(lower["plane_coefficients"]),
            maximum_iterations=40,
            constraint_tolerance=2e-10,
        )
        self.assertEqual(interior["status"], "LINE_SEARCH_EXHAUSTED")
        self.assertIn("degenerate_faces", interior["diagnostic"])
        self.assertFalse(interior["paper_optimality_conditions_verified"])


class SpinNativeArtifactTests(unittest.TestCase):
    def test_plane_halfspaces_return_an_exact_complementary_material_partition(self) -> None:
        self.assertTrue(hasattr(spin_api, "evaluate_spin_it_faster_plane_halfspaces"))
        outer_vertices, outer_faces = cube(-2.0, 2.0)
        design_vertices, design_faces = cube(-1.0, 1.0)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            material_density=1.0,
        )

        partition = spin_api.evaluate_spin_it_faster_plane_halfspaces(
            request, coefficients=(0.0, 1.0, 0.0, 0.0)
        )

        self.assertEqual(partition["status"], "EXACT_COMPLEMENTARY_PARTITION")
        self.assertLessEqual(partition["normalized_moment_partition_error"], 1e-10)
        self.assertTrue(partition["exact_partition_audit"]["complete_partition"])
        self.assertTrue(partition["exact_partition_audit"]["no_positive_volume_overlap"])
        self.assertTrue(partition["interface_stitch"]["exact_coordinate_matching"])
        self.assertFalse(partition["interface_stitch"]["vertex_coordinates_moved"])
        self.assertLessEqual(partition["independent_mesh_moment_error"], 1e-10)
        self.assertLessEqual(partition["independent_partition_moment_error"], 1e-10)
        positive, negative = partition["positive"], partition["negative"]
        self.assertAlmostEqual(positive["raw_moments"][0], 4.0, places=10)
        self.assertAlmostEqual(negative["raw_moments"][0], 4.0, places=10)
        self.assertTrue(any(abs(vertex[0]) <= 1e-12 for vertex in positive["vertices"]))
        self.assertTrue(any(abs(vertex[0]) <= 1e-12 for vertex in negative["vertices"]))

    def test_planar_fidget_solver_builds_fixed_and_design_moments_from_disjoint_loops(self) -> None:
        """The dogleg must use geometry-derived moments for a hole-separated fixed polygon."""
        pi = math.pi
        outer = ((-3,-3),(3,-3),(3,3),(-3,3))
        patch = ((pi/4-0.1,-0.1),(pi/4+0.1,-0.1),(pi/4+0.1,0.1),(pi/4-0.1,0.1))
        hole = tuple(reversed(patch))
        result = solve_spin_it_faster_planar_fidget(
            design_region_loops=(outer,hole),
            fixed_material_loops=(patch,),
            fixed_areal_density=100.0,
            design_areal_density=1.0,
            initial_center=(0.75,0.2),
            extrusion_interval=(-0.7,0.8),
            extrusion_volume_density=2.0/3.0,
            fixed_extrusion_layers=(
                {"loops": (patch,), "z_interval": (-0.1,0.0), "volume_density": 500.0},
                {"loops": (patch,), "z_interval": (0.0,0.1), "volume_density": 500.0},
            ),
        )
        self.assertEqual(result["status"], "CONVERGED")
        np.testing.assert_allclose(result["center"], (1.0,0.0), rtol=0, atol=2e-7)
        self.assertAlmostEqual(result["radius"], 1.0, places=7)
        self.assertGreater(result["accepted_steps"], 0)
        self.assertGreater(result["initial_normalized_residual"], 1e-6)
        self.assertLessEqual(result["normalized_residual"], 1e-10)
        self.assertEqual(result["fixed_editable_region_overlap_status"],
                         "DISJOINT_BY_EXACT_REGULARIZED_INTERSECTION")
        remaining = np.asarray(result["remaining_material_moments"], dtype=np.float64)
        self.assertAlmostEqual(remaining[0], 40.0-pi, places=7)
        self.assertAlmostEqual(remaining[1], 0.0, places=7)
        self.assertAlmostEqual(remaining[2], 0.0, places=7)
        self.assertTrue(result["boundary_defined"], result["boundary_diagnostic"])
        self.assertEqual(len(result["boundary_cycles"]), 2)
        pieces = [piece for cycle in result["boundary_cycles"] for piece in cycle]
        self.assertTrue(any(piece["type"] == "CIRCULAR_ARC" for piece in pieces))
        self.assertEqual(result["geometry_scope"], "exact_2d_boundary_and_analytic_layered_extrusion")
        self.assertEqual(result["extruded_material_scope"], "design_plus_fixed_layers")
        spatial = np.asarray(result["extruded_design_material_raw_moments"], dtype=np.float64)
        self.assertEqual(spatial.shape, (10,))
        self.assertTrue(np.all(np.isfinite(spatial)))
        self.assertAlmostEqual(spatial[0], 36.0-pi, places=7)
        # Independent closed-form oracle for the actual remaining region:
        # square [-3,3]^2 minus the unit disk centered at (1,0), extruded
        # uniformly through [-0.7,0.8] with volume density 2/3.
        design_mass = 36.0 - pi
        design_x = -pi
        design_xx = 108.0 - 5.0*pi/4.0
        design_yy = 108.0 - pi/4.0
        design_z_mean = 0.05
        design_z_second = 0.19
        expected_design_3d = np.asarray((
            design_mass, design_x, 0.0, design_mass*design_z_mean,
            design_xx, design_yy, design_mass*design_z_second,
            0.0, design_x*design_z_mean, 0.0,
        ))
        np.testing.assert_allclose(spatial, expected_design_3d, rtol=1e-10, atol=1e-9)
        fixed_layers = result["extruded_fixed_layer_raw_moments"]
        self.assertEqual(len(fixed_layers), 2)
        self.assertAlmostEqual(sum(layer[0] for layer in fixed_layers), 4.0, places=12)
        self.assertAlmostEqual(result["extruded_complete_material_raw_moments"][0], 40.0-pi, places=7)
        patch_area = 0.2 * 0.2
        patch_center_x = pi/4.0
        patch_variance = 0.2**2 / 12.0
        expected_layers = []
        for z_mean in (-0.05, 0.05):
            layer_mass = 500.0 * 0.1 * patch_area
            expected_layers.append(np.asarray((
                layer_mass, layer_mass*patch_center_x, 0.0, layer_mass*z_mean,
                layer_mass*(patch_center_x**2+patch_variance),
                layer_mass*patch_variance, layer_mass/300.0,
                0.0, layer_mass*patch_center_x*z_mean, 0.0,
            )))
        for observed, expected in zip(fixed_layers, expected_layers):
            np.testing.assert_allclose(observed, expected, rtol=1e-12, atol=1e-10)
        expected_complete_3d = expected_design_3d + sum(expected_layers)
        np.testing.assert_allclose(
            result["extruded_complete_material_raw_moments"],
            expected_complete_3d, rtol=1e-10, atol=1e-9,
        )
        entities = result["extruded_material_entities"]
        self.assertEqual(tuple(entity["role"] for entity in entities),
                         ("editable_design", "fixed_material", "fixed_material"))
        self.assertEqual(entities[0]["boundary_cycles"], result["boundary_cycles"])
        np.testing.assert_allclose(entities[0]["z_interval"], (-0.7, 0.8), atol=0.0)
        self.assertEqual(entities[0]["volume_density"], 2.0/3.0)
        np.testing.assert_allclose(entities[0]["raw_mass_moments"], expected_design_3d,
                                   rtol=1e-10, atol=1e-9)
        for index, (entity, layer, expected) in enumerate(
                zip(entities[1:], result["fixed_extrusion_layers"], expected_layers)):
            self.assertEqual(entity["loops"], layer["loops"])
            np.testing.assert_array_equal(entity["z_interval"], layer["z_interval"])
            self.assertEqual(entity["volume_density"], layer["volume_density"])
            np.testing.assert_allclose(entity["raw_mass_moments"], expected,
                                       rtol=1e-12, atol=1e-10,
                                       err_msg=f"fixed material entity {index}")
        np.testing.assert_allclose(
            sum((entity["raw_mass_moments"] for entity in entities), np.zeros(10)),
            result["extruded_complete_material_raw_moments"], rtol=1e-12, atol=1e-9,
        )
        layer_driven = solve_spin_it_faster_planar_fidget(
            design_region_loops=(outer,hole),
            fixed_material_loops=(), fixed_areal_density=None,
            design_areal_density=1.0, initial_center=(0.75,0.2),
            extrusion_interval=(-0.7,0.8), extrusion_volume_density=2.0/3.0,
            fixed_extrusion_layers=(
                {"loops": (patch,), "z_interval": (-0.1,0.0), "volume_density": 500.0},
                {"loops": (patch,), "z_interval": (0.0,0.1), "volume_density": 500.0},
            ),
        )
        self.assertEqual(layer_driven["status"], "CONVERGED")
        np.testing.assert_allclose(layer_driven["center"], result["center"], rtol=0, atol=2e-10)
        np.testing.assert_allclose(
            layer_driven["fixed_material_moments"],
            result["fixed_material_moments"], rtol=1e-12, atol=1e-12,
        )
        np.testing.assert_allclose(
            layer_driven["extruded_complete_material_raw_moments"],
            result["extruded_complete_material_raw_moments"], rtol=1e-10, atol=1e-9,
        )
        left_patch = ((-1.2,-0.1),(-1.0,-0.1),(-1.0,0.1),(-1.2,0.1))
        right_patch = ((1.0,-0.1),(1.2,-0.1),(1.2,0.1),(1.0,0.1))
        same_height_disjoint = solve_spin_it_faster_planar_fidget(
            design_region_loops=(outer,hole),
            fixed_material_loops=(), fixed_areal_density=None,
            design_areal_density=1.0, initial_center=(0.75,0.2),
            extrusion_interval=(-0.7,0.8), extrusion_volume_density=2.0/3.0,
            fixed_extrusion_layers=(
                {"loops": (left_patch,), "z_interval": (1.0,1.1), "volume_density": 500.0},
                {"loops": (right_patch,), "z_interval": (1.0,1.1), "volume_density": 500.0},
            ),
        )
        np.testing.assert_allclose(
            same_height_disjoint["fixed_material_moments"],
            (4.0, 0.0, 0.0, 4.0*(1.1**2+0.2**2/12.0),
             4.0*0.2**2/12.0, 0.0), rtol=1e-12, atol=1e-12,
        )
        with self.assertRaisesRegex(ValueError, "projections do not match"):
            solve_spin_it_faster_planar_fidget(
                design_region_loops=(outer,hole), fixed_material_loops=(patch,),
                fixed_areal_density=100.0, design_areal_density=1.0,
                initial_center=(0.75,0.2), extrusion_interval=(-0.7,0.8),
                extrusion_volume_density=2.0/3.0,
                fixed_extrusion_layers=({
                    "loops": (patch,), "z_interval": (-0.1,0.0),
                    "volume_density": 250.0,
                },),
            )
        with self.assertRaisesRegex(ValueError, "positive-volume overlap"):
            solve_spin_it_faster_planar_fidget(
                design_region_loops=(outer,hole), fixed_material_loops=(patch,),
                fixed_areal_density=100.0, design_areal_density=1.0,
                initial_center=(0.75,0.2), extrusion_interval=(-0.7,0.8),
                extrusion_volume_density=2.0/3.0,
                fixed_extrusion_layers=(
                    {"loops": (patch,), "z_interval": (-0.1,0.05), "volume_density": 500.0},
                    {"loops": (patch,), "z_interval": (0.0,0.1), "volume_density": 500.0},
                ),
            )

    def test_planar_fidget_layer_geometry_changes_the_optimized_solution(self) -> None:
        """Real fixed-layer location, not a caller-supplied moment, drives dogleg."""
        outer = ((-3,-3),(3,-3),(3,3),(-3,3))
        # The large central hole is the common fixed-material domain; moving
        # the layer inside it leaves the editable design geometry unchanged.
        hole = ((-1.5,-0.3),(-1.5,0.3),(1.5,0.3),(1.5,-0.3))
        design = (outer, hole)
        layers = {
            "left": ((-1.2,-0.1),(-1.0,-0.1),(-1.0,0.1),(-1.2,0.1)),
            "right": ((1.0,-0.1),(1.2,-0.1),(1.2,0.1),(1.0,0.1)),
        }

        def solve_with_layer(loop):
            return solve_spin_it_faster_planar_fidget(
                design_region_loops=design,
                fixed_material_loops=(), fixed_areal_density=None,
                design_areal_density=1.0, initial_center=(0.75,0.2),
                extrusion_interval=(-0.7,0.8),
                extrusion_volume_density=2.0/3.0,
                fixed_extrusion_layers=(
                    {"loops": (loop,), "z_interval": (-0.1,0.0),
                     "volume_density": 400.0},
                    {"loops": (loop,), "z_interval": (0.0,0.1),
                     "volume_density": 600.0},
                ),
            )

        left = solve_with_layer(layers["left"])
        right = solve_with_layer(layers["right"])
        self.assertEqual(left["status"], "CONVERGED")
        self.assertEqual(right["status"], "CONVERGED")
        # Each layer is integrated from its actual polygon, height and density;
        # the two unequal-density slabs share only a legal interface.
        self.assertAlmostEqual(left["fixed_material_moments"][0], 4.0, places=12)
        self.assertAlmostEqual(right["fixed_material_moments"][0], 4.0, places=12)
        self.assertAlmostEqual(left["fixed_material_moments"][1], -4.4, places=12)
        self.assertAlmostEqual(right["fixed_material_moments"][1], 4.4, places=12)
        self.assertGreater(left["accepted_steps"], 0)
        self.assertGreater(right["accepted_steps"], 0)
        self.assertGreater(np.linalg.norm(left["center"] - right["center"]), 2.0)
        self.assertLessEqual(left["normalized_residual"], 1e-10)
        self.assertLessEqual(right["normalized_residual"], 1e-10)
        for result in (left, right):
            self.assertEqual(
                tuple(entity["role"] for entity in result["extruded_material_entities"]),
                ("editable_design", "fixed_material", "fixed_material"),
            )
            self.assertAlmostEqual(
                sum(layer[0] for layer in result["extruded_fixed_layer_raw_moments"]),
                4.0, places=12,
            )
            self.assertAlmostEqual(
                result["extruded_complete_material_raw_moments"][0],
                result["extruded_design_material_raw_moments"][0] + 4.0,
                places=9,
            )

    def test_planar_fidget_solver_rejects_boolean_coordinates(self) -> None:
        patch = ((0.6,-0.1),(0.8,-0.1),(0.8,0.1),(0.6,0.1))
        with self.assertRaisesRegex(ValueError, "initial_center"):
            solve_spin_it_faster_planar_fidget(
                design_region_loops=((((-2,-2),(2,-2),(2,2),(-2,2))),),
                fixed_material_loops=(patch,),
                fixed_areal_density=1.0,
                design_areal_density=1.0,
                initial_center=(True,0.0),
            )

    def test_planar_fidget_rejects_fixed_editable_area_overlap(self) -> None:
        pi = math.pi
        outer = ((-3,-3),(3,-3),(3,3),(-3,3))
        patch = ((pi/4-0.1,-0.1),(pi/4+0.1,-0.1),
                 (pi/4+0.1,0.1),(pi/4-0.1,0.1))
        with self.assertRaisesRegex(ValueError, "fixed and editable planar regions"):
            solve_spin_it_faster_planar_fidget(
                design_region_loops=(outer,), fixed_material_loops=(patch,),
                fixed_areal_density=100.0, design_areal_density=1.0,
                initial_center=(0.75,0.2),
            )

    def test_external_moment_solver_fails_closed_without_defined_geometry_jacobian(self) -> None:
        outer_vertices, outer_faces = cube(-3.0, 3.0)
        design_vertices, design_faces = cube(-1.0, 1.0)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            settings=SpinSolverConfig(
                equation_profile="DUAL_MONOMIAL_RAW_V1",
                maximum_restarts=1,
                maximum_function_evaluations_per_restart=1,
            ),
        )
        calls: list[tuple[float, ...]] = []

        def evaluator(coefficients: tuple[float, ...]) -> dict[str, object]:
            calls.append(coefficients)
            return {
                "raw_moments": (1.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.0, 0.0, 0.0),
                "moment_jacobian": np.zeros((10, 10)),
                "derivative_defined": False,
                "derivative_diagnostic": "fixture topology provenance unavailable",
                "geometry_hash": "fixture-geometry-v1",
                "topology_hash": "fixture-topology-v1",
                "evaluation_map_identity": "a" * 64,
                "boolean_backend_identity": "fixture-boolean-backend-v1",
                "derivative_identity": "fixture-derivative-v1",
            }

        result = solve_spin_it_faster_with_external_moment_evaluator(
            request, evaluator, initial_quadrics=((1.0, 0.0, 0.0, 0.0, -1.0,
                                                     -1.0, -1.0, 0.0, 0.0, 0.0),)
        )
        self.assertEqual(len(calls), 1)
        self.assertEqual(result["function_evaluations"], 1)
        self.assertNotEqual(result["status"], "CONVERGED")
        self.assertFalse(result["final_geometry_verified"])
        self.assertEqual(result["scope"], "numerical_newton_only_external_geometry_not_sealed")
        self.assertTrue(result["event_hash"])
        self.assertTrue(result["trace_hash"])
        evaluation_events = [item for item in result["events"] if item["type"] == "evaluation"]
        self.assertEqual(len(evaluation_events), 1)
        self.assertEqual(evaluation_events[0]["material_geometry_hash"], "fixture-geometry-v1")
        self.assertEqual(evaluation_events[0]["material_topology_hash"], "fixture-topology-v1")
        self.assertEqual(len(evaluation_events[0]["quadric_coefficients"]), 10)
        self.assertEqual(evaluation_events[0]["equation_profile"], "DUAL_MONOMIAL_RAW_V1")

    def test_external_moment_solver_binds_backend_derivative_and_map_identity_to_event(self) -> None:
        from make_it_stand import spin

        outer_vertices, outer_faces = cube(-3.0, 3.0)
        design_vertices, design_faces = cube(-1.0, 1.0)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            settings=SpinSolverConfig(
                equation_profile="DUAL_MONOMIAL_RAW_V1",
                maximum_restarts=1,
                maximum_function_evaluations_per_restart=1,
            ),
        )
        identities = {
            "evaluation_map_identity": "d" * 64,
            "boolean_backend_identity": "cherchi-reference:" + "c" * 64,
            "derivative_identity": "cherchi-provenance-reynolds-v1",
        }

        def evaluator(coefficients: tuple[float, ...]) -> dict[str, object]:
            return {
                "raw_moments": (1.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.0, 0.0, 0.0),
                "moment_jacobian": np.eye(10),
                "derivative_defined": True,
                "geometry_hash": "fixture-geometry-v1",
                "topology_hash": "fixture-topology-v1",
                **identities,
            }

        result = solve_spin_it_faster_with_external_moment_evaluator(
            request, evaluator,
            initial_quadrics=((1.0, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0),),
        )
        evaluation_event = next(item for item in result["events"] if item["type"] == "evaluation")
        native = spin._native_module()
        self.assertEqual(
            native.spin_it_faster_quadric_hash(evaluation_event["quadric_coefficients"]),
            evaluation_event["quadric_hash"],
        )
        for token, identity in (
            ("evaluation_map=", identities["evaluation_map_identity"]),
            ("backend=", identities["boolean_backend_identity"]),
            ("derivative=", identities["derivative_identity"]),
        ):
            self.assertIn(token + identity, evaluation_event["detail"])

    def test_cylinder_stratum_grid_settings_are_validated_and_part_of_identity(self) -> None:
        outer_vertices, outer_faces = cube(-2.0, 2.0)
        design_vertices, design_faces = cube(-1.0, 1.0)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
        )
        coefficients = (0.25, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0)
        baseline = evaluate_spin_it_faster_quadric(request, coefficients)
        for field in ("hyperbolic_cylinder_grid_resolution", "parabolic_cylinder_grid_resolution"):
            changed = evaluate_spin_it_faster_quadric(
                request, coefficients,
                evaluator_settings=SpinEvaluatorConfig(**{field: 25}),
            )
            self.assertNotEqual(baseline.request_hash, changed.request_hash, field)
        for field, invalid in (("hyperbolic_cylinder_grid_resolution", 1),
                               ("parabolic_cylinder_grid_resolution", 0)):
            with self.assertRaisesRegex(ValueError, field):
                SpinEvaluatorConfig(**{field: invalid})

    def test_unverified_fullrank_failure_does_not_expose_material_mesh(self) -> None:
        outer_vertices, outer_faces = cube(-3.0, 3.0)
        design_vertices, design_faces = cube(-1.0, 1.0)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            settings=SpinSolverConfig(
                equation_profile="DUAL_MONOMIAL_RAW_V1",
                maximum_restarts=1,
                maximum_function_evaluations_per_restart=2,
            ),
        )
        result = solve_spin_it_faster_fullrank_quadric(
            request,
            initial_quadrics=((0.371, 0.8, 0.123, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),),
        )
        self.assertNotEqual(result.status, "CONVERGED")
        self.assertFalse(result.final_geometry_verified)
        self.assertIsNone(result.material_vertices)
        self.assertIsNone(result.material_faces)
        self.assertTrue(result.events)
        self.assertEqual(result.events[-1].type, "attempt_terminal")
        self.assertEqual(result.events[-1].reason, "initial_evaluation_failed")
        self.assertEqual(len(result.evaluation_timings), result.function_evaluations)
        self.assertEqual(sum(result.restart_function_evaluations), result.function_evaluations)

    def test_fullrank_solver_consumes_discrete_jacobian_from_perturbed_top_seed(self) -> None:
        outer_vertices, outer_faces = cube(-1.0, 1.0)
        design_vertices, design_faces = cube(-0.8, 0.8)
        outer_vertices = tuple((x, y, z + 2.0) for x, y, z in outer_vertices)
        design_vertices = tuple((x, y, z + 2.0) for x, y, z in design_vertices)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            mode="SPINNING_TOP",
            center_of_mass_weight=1.0,
            minimum_wall_thickness=0.15,
            settings=SpinSolverConfig(
                equation_profile="DUAL_MONOMIAL_RAW_V1",
                maximum_iterations=60,
                maximum_restarts=1,
                maximum_function_evaluations_per_restart=60,
            ),
        )
        seed_coefficients = np.asarray(
            (-0.31, 0.07, -0.09, 0.11, 0.83, 1.17, -0.72, 0.06, -0.04, 0.03)
        )
        linear = seed_coefficients[1:4].copy()
        quadratic = np.asarray((
            (seed_coefficients[4], seed_coefficients[9], seed_coefficients[8]),
            (seed_coefficients[9], seed_coefficients[5], seed_coefficients[7]),
            (seed_coefficients[8], seed_coefficients[7], seed_coefficients[6]),
        ))
        center = np.asarray((0.0, 0.0, 2.0))
        target = seed_coefficients.copy()
        target[0] = seed_coefficients[0] - linear @ center + center @ quadratic @ center
        target[1:4] = linear - 2.0 * quadratic @ center
        target_moments = evaluate_spin_it_faster_fullrank_quadric(request, target)
        self.assertTrue(target_moments.derivative_defined, target_moments.derivative_diagnostic)
        outer_moments = np.asarray(integrate_raw_mass_moments(outer_vertices, outer_faces))
        design_moments = np.asarray(integrate_raw_mass_moments(design_vertices, design_faces))
        request = replace(
            request,
            fixed_shell_moments=tuple(outer_moments - design_moments),
        )
        seed = target.copy()
        seed[0] += 1e-4
        seed[6] -= 1e-4
        self.assertNotEqual(tuple(seed), tuple(target))

        result = solve_spin_it_faster_fullrank_quadric(request, initial_quadrics=(seed,))
        self.assertEqual(result.status, "CONVERGED", result.diagnostic)
        self.assertGreater(len(result.trace), 1, "the fixture must take a Newton step")
        self.assertTrue(result.final_geometry_verified, result.verification_diagnostic)
        self.assertTrue(result.final_derivative_defined, result.final_derivative_diagnostic)
        self.assertLessEqual(result.independent_raw_residual, request.settings.equation_tolerance)
        self.assertRegex(result.output_hash, r"^[0-9a-f]{64}$")
        self.assertTrue(result.events)
        self.assertEqual(result.events[0].type, "attempt_start")
        self.assertTrue(any(
            event.reason == "topology_transition_relinearized"
            for event in result.events
        ), "the perturbed fixture must cross a Boolean topology stratum and relinearize")
        self.assertEqual(result.events[-1].reason, "converged")
        self.assertEqual(len(result.evaluation_timings), result.function_evaluations)
        self.assertEqual(sum(result.restart_function_evaluations), result.function_evaluations)
        self.assertGreater(len(result.material_faces), 0)
        np.testing.assert_allclose(
            result.final_total_raw_moments,
            result.independent_total_raw_moments,
            rtol=1e-10,
            atol=1e-10,
        )
        repeated = solve_spin_it_faster_fullrank_quadric(request, initial_quadrics=(seed,))
        self.assertEqual(result.output_hash, repeated.output_hash)
        self.assertEqual(result.event_hash, repeated.event_hash)
    def test_public_fullrank_hyperboloid_derivative_matches_discrete_moment_map(self) -> None:
        outer_vertices, outer_faces = cube(-3.0, 3.0)
        design_vertices, design_faces = cube(-1.0, 1.0)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            minimum_wall_thickness=2.0,
        )
        coefficients = np.asarray(
            (-0.3, 0.0, 0.0, 0.0, 1.0, 1.0, -0.7, 0.0, 0.0, 0.0),
            dtype=np.float64,
        )
        base = evaluate_spin_it_faster_fullrank_quadric(request, coefficients)
        self.assertTrue(base.derivative_defined, base.derivative_diagnostic)
        self.assertRegex(base.boolean_backend_identity,
            r'^libigl::copyleft::cgal::mesh_boolean;libigl_source=revision=[0-9a-f]+,tree-sha256=[0-9a-f]{64}$')
        self.assertEqual(np.asarray(base.coefficient_jacobian).shape, (10, 10))
        self.assertRegex(base.geometry_hash, r"^[0-9a-f]{64}$")
        self.assertGreater(len(base.material_faces), 0)

        unit = coefficients / np.linalg.norm(coefficients)
        direction = np.eye(10)[4] - unit * unit[4]
        direction /= np.linalg.norm(direction)
        step = 1e-5
        plus_coefficients = np.cos(step) * unit + np.sin(step) * direction
        minus_coefficients = np.cos(step) * unit - np.sin(step) * direction
        plus = evaluate_spin_it_faster_quadric(request, plus_coefficients)
        minus = evaluate_spin_it_faster_quadric(request, minus_coefficients)
        finite_difference = (
            np.asarray(plus.variable_raw_moments)
            - np.asarray(minus.variable_raw_moments)
        ) / (2.0 * step)
        predicted = np.asarray(base.coefficient_jacobian) @ direction
        self.assertLess(
            np.linalg.norm(predicted - finite_difference)
            / max(1.0, np.linalg.norm(finite_difference)),
            2e-5,
        )

    def test_fullrank_adaptive_hyperboloid_boolean_moment_jacobian(self) -> None:
        outer_vertices, outer_faces = cube(-3.0, 3.0)
        design_vertices, design_faces = cube(-0.5, 0.5)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            minimum_wall_thickness=2.5,
        )
        # This two-sheet hyperboloid cuts the design box. The requested
        # truncation radius 0.3 doubles to 0.6 before the artificial caps lie
        # outside the box; this exercises the actual Boolean-to-moment path.
        coefficients = np.asarray(
            (0.15, 0.0, 0.0, 0.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0),
            dtype=np.float64,
        )
        evaluator_settings = SpinEvaluatorConfig(
            icosphere_subdivisions=2,
            truncation_radius=0.3,
        )
        base = evaluate_spin_it_faster_fullrank_quadric(
            request, coefficients, evaluator_settings=evaluator_settings
        )
        self.assertTrue(base.derivative_defined, base.derivative_diagnostic)
        self.assertGreater(len(base.material_faces or ()), 0)

        unit = coefficients / np.linalg.norm(coefficients)
        jacobian = np.asarray(base.coefficient_jacobian)
        previous = np.inf
        for step in (1e-4, 3e-5, 1e-5):
            error_squared = 0.0
            finite_difference_squared = 0.0
            for coordinate in range(10):
                direction = np.eye(10)[coordinate] - unit * unit[coordinate]
                direction_norm = np.linalg.norm(direction)
                if direction_norm <= 1e-12:
                    continue
                direction /= direction_norm
                plus_coefficients = np.cos(step) * unit + np.sin(step) * direction
                minus_coefficients = np.cos(step) * unit - np.sin(step) * direction
                plus = evaluate_spin_it_faster_quadric(
                    request, plus_coefficients, evaluator_settings=evaluator_settings
                )
                minus = evaluate_spin_it_faster_quadric(
                    request, minus_coefficients, evaluator_settings=evaluator_settings
                )
                finite_difference = (
                    np.asarray(plus.variable_raw_moments)
                    - np.asarray(minus.variable_raw_moments)
                ) / (2.0 * step)
                predicted = jacobian @ direction
                error_squared += float(np.dot(predicted - finite_difference,
                                              predicted - finite_difference))
                finite_difference_squared += float(np.dot(finite_difference,
                                                          finite_difference))
            relative_error = np.sqrt(error_squared) / max(
                1.0, np.sqrt(finite_difference_squared)
            )
            if step == 1e-5:
                self.assertLess(relative_error, 1e-7)
            if np.isfinite(previous):
                self.assertLess(relative_error, 0.2 * previous)
            previous = relative_error

    def test_fullrank_rejects_hyperboloid_cap_that_changes_design_material(self) -> None:
        outer_vertices, outer_faces = cube(-6.0, 6.0)
        design_vertices, design_faces = cube(-4.0, 4.0)
        outer_moments = integrate_raw_mass_moments(outer_vertices, outer_faces)
        design_moments = integrate_raw_mass_moments(design_vertices, design_faces)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            fixed_shell_moments=tuple(
                outer - design for outer, design in zip(outer_moments, design_moments)
            ),
            minimum_wall_thickness=2.0,
        )
        # q=4+x^2+y^2-z^2 is a two-sheet hyperboloid. At truncation radius
        # 3 its artificial end caps cross this design cube. With no allowed
        # radius refinement the evaluator must reject the uncertified source;
        # radius 8 places both caps outside the design region.
        coefficients = (4.0, 0.0, 0.0, 0.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0)
        with self.assertRaisesRegex(
            RuntimeError,
            "hyperboloid truncation failed after bounded radius refinement",
        ):
            evaluate_spin_it_faster_fullrank_quadric(
                request,
                coefficients,
                evaluator_settings=SpinEvaluatorConfig(
                    icosphere_subdivisions=2,
                    truncation_radius=3.0,
                    maximum_radius_doublings=0,
                ),
            )

        valid = evaluate_spin_it_faster_fullrank_quadric(
            request,
            coefficients,
            evaluator_settings=SpinEvaluatorConfig(
                icosphere_subdivisions=2,
                truncation_radius=8.0,
            ),
        )
        self.assertTrue(valid.derivative_defined, valid.derivative_diagnostic)
        self.assertGreater(len(valid.material_faces), 0)

        grown = evaluate_spin_it_faster_fullrank_quadric(
            request,
            coefficients,
            evaluator_settings=SpinEvaluatorConfig(
                icosphere_subdivisions=2,
                truncation_radius=3.0,
                cap_radius_policy="GROW_TO_CERTIFIED_V1",
            ),
        )
        self.assertTrue(grown.derivative_defined, grown.derivative_diagnostic)
        self.assertGreater(len(grown.material_faces), 0)

    def test_dual_profile_native_roundtrip_and_failed_serialization(self) -> None:
        from make_it_stand.spin import _native_request
        v, f = cube(-5.0, 5.0)
        d, g = cube(-1.0, 1.0)
        config = SpinSolverConfig(equation_profile="DUAL_MONOMIAL_RAW_V1",
                                  maximum_iterations=1, maximum_restarts=1,
                                  maximum_function_evaluations_per_restart=1)
        request = SpinItFasterRequest(v, f, d, g, minimum_wall_thickness=4.0, settings=config)
        native, converted = _native_request(request)
        self.assertEqual(converted.settings.equation_profile, native.SpinEquationProfile.DUAL_MONOMIAL_RAW_V1)
        with self.assertRaises(ValueError):
            replace(config, equation_profile="UNKNOWN")
        with self.assertRaises(ValueError):
            replace(config, maximum_function_evaluations_per_restart=61)
        initial = ((1, 0, 0, 0, 0, 0, 0, 0, 0, 0),)
        strict = solve_spin_it_faster(request, initial_quadrics=initial)
        repeated = solve_spin_it_faster(request, initial_quadrics=initial)
        legacy_request = replace(request, settings=replace(config,
            equation_profile="LEGACY_NATIVE_NORMALIZED_V1"))
        legacy = solve_spin_it_faster(legacy_request, initial_quadrics=initial)
        legacy_repeated = solve_spin_it_faster(legacy_request, initial_quadrics=initial)
        self.assertEqual(strict.equation_profile, "DUAL_MONOMIAL_RAW_V1")
        self.assertEqual(strict.solver_coordinate_identity, "literal_monomial_unit_sphere_v1")
        self.assertIsNone(strict.raw_equation_residual)
        self.assertEqual(strict.diagnostics["event_stream_version"], 2)
        self.assertEqual(strict.diagnostics, repeated.diagnostics)
        for event in strict.events:
            self.assertEqual(event.equation_profile, strict.equation_profile)
            self.assertEqual(event.solver_coordinate_identity, strict.solver_coordinate_identity)
            self.assertIsNone(event.raw_equation_residual)
        self.assertNotEqual(strict.diagnostics["input_hash"], legacy.diagnostics["input_hash"])
        self.assertNotEqual(strict.event_hash, legacy.event_hash)
        self.assertNotIn("equation_profile", legacy.diagnostics)
        self.assertEqual(legacy.equation_profile, "LEGACY_NATIVE_NORMALIZED_V1")
        self.assertEqual(legacy.events[1].equation_profile, legacy.equation_profile)
        # Canonical legacy constant-quadric failure. The preserved pre-endpoint
        # native module and new module emit identical bytes (see the scoped
        # compare_top_endpoint_legacy.py diagnostic). This baseline includes
        # the already-existing static-objective-measure metadata.
        # Independent evaluations must still emit identical bytes.
        self.assertEqual(legacy._native.diagnostics_json, legacy_repeated._native.diagnostics_json)
        self.assertEqual(hashlib.sha256((legacy._native.diagnostics_json + "\n").encode()).hexdigest(),
                         "0c7accec7aa9bb75cbde853e69c5a0a65a096bcedfe39dd701dcc4f53db843e0")
        json.dumps(strict.diagnostics, allow_nan=False)

        # A returned geometry evaluation has finite raw metrics and a native p
        # in each trace row, from which s=Dp/||Dp|| is recoverable.
        evaluated = solve_spin_it_faster(request, initial_quadrics=((.371,.8,.123,0,0,0,0,0,0,0),))
        self.assertTrue(evaluated.trace)
        self.assertIsNotNone(evaluated.raw_equation_residual)
        self.assertEqual(evaluated._native.raw_equation_residual, evaluated.raw_equation_residual)
        self.assertEqual(evaluated._native.equation_profile, evaluated.equation_profile)
        self.assertEqual(evaluated._native.solver_coordinate_identity, evaluated.solver_coordinate_identity)
        trace = evaluated.trace[0]
        self.assertIsInstance(trace.quadric_coefficients, tuple)
        self.assertEqual(trace.equation_profile, evaluated.equation_profile)
        self.assertEqual(trace.solver_coordinate_identity, evaluated.solver_coordinate_identity)
        self.assertIsNotNone(trace.raw_equation_residual)
        with self.assertRaises(TypeError):
            trace.quadric_coefficients[0] = 42  # type: ignore[index]

    def test_raw_moment_integration_matches_cube(self) -> None:
        vertices, faces = cube(-1.0, 1.0)
        moments = integrate_raw_mass_moments(vertices, faces)
        np.testing.assert_allclose(
            moments,
            (8.0, 0.0, 0.0, 0.0, 8.0 / 3.0, 8.0 / 3.0,
             8.0 / 3.0, 0.0, 0.0, 0.0),
            atol=1e-14,
        )

    def test_failed_evaluation_events_are_strict_json_safe(self) -> None:
        outer_vertices, outer_faces = cube(-5.0, 5.0)
        design_vertices, design_faces = cube(-1.0, 1.0)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            minimum_wall_thickness=4.0,
            settings=SpinSolverConfig(
                maximum_iterations=1,
                maximum_restarts=1,
                maximum_function_evaluations_per_restart=1,
            ),
        )
        failed = solve_spin_it_faster(
            request,
            initial_quadrics=((1.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0),),
        )
        self.assertEqual(failed.status, "GEOMETRY_EVALUATION_FAILED")
        self.assertEqual(failed.solve_function_evaluations, 1)
        self.assertEqual(len(failed.evaluation_timings), 1)
        self.assertFalse(failed.evaluation_timings[0]["evaluator_returned"])
        evaluation = failed.events[1]
        self.assertEqual(evaluation.outcome, "failed")
        self.assertEqual(evaluation.reason, "evaluation_exception")
        self.assertIsNotNone(evaluation.quadric_coefficients)
        self.assertIsInstance(evaluation.quadric_coefficients, tuple)
        self.assertEqual(len(evaluation.quadric_coefficients), 10)
        original_event_hash = failed.event_hash
        original_coefficient = evaluation.quadric_coefficients[0]
        with self.assertRaises(TypeError):
            evaluation.quadric_coefficients[0] = 42.0  # type: ignore[index]
        copied_diagnostics = failed.diagnostics
        copied_diagnostics["events"][1]["quadric_coefficients"][0] = 42.0
        self.assertEqual(
            failed.events[1].quadric_coefficients[0], original_coefficient
        )
        self.assertEqual(
            failed.diagnostics["events"][1]["quadric_coefficients"][0],
            original_coefficient,
        )
        self.assertEqual(failed.event_hash, original_event_hash)
        self.assertEqual(
            failed.events[1].quadric_hash,
            failed.diagnostics["events"][1]["quadric_hash"],
        )
        self.assertIsNone(evaluation.normalized_equation_residual)
        self.assertIsNone(evaluation.material_geometry_hash)
        self.assertIsNone(failed.events[-1].quadric_coefficients)
        self.assertEqual(
            failed.events[-1].reason,
            "initial_evaluation_failed",
        )
        json.dumps(failed.diagnostics, allow_nan=False)

    def test_planar_fixture_returns_immutable_sealed_artifact(self) -> None:
        outer_vertices, outer_faces = cube(-5.0, 5.0)
        design_vertices, design_faces = cube(-1.0, 1.0)
        quadric = (0.371, 0.8, 0.123, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        request = SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            minimum_wall_thickness=4.0,
            settings=SpinSolverConfig(maximum_restarts=1),
        )
        evaluated = evaluate_spin_it_faster_quadric(request, quadric)
        self.assertEqual(evaluated.interface_type, "PLANE")
        self.assertEqual(evaluated.interface_subtype, "PLANE")
        fixed = (
            10.0,
            -evaluated.variable_raw_moments[1],
            -evaluated.variable_raw_moments[2],
            -evaluated.variable_raw_moments[3],
            20.0,
            20.0,
            20.0,
            0.0,
            0.0,
            0.0,
        )
        request = replace(request, fixed_shell_moments=fixed)
        first = solve_spin_it_faster(request, initial_quadrics=(quadric,))
        second = solve_spin_it_faster(request, initial_quadrics=(quadric,))

        self.assertEqual(first.status, "CONVERGED")
        self.assertEqual(first.evidence, "PAPER_FORMULA_REPRODUCTION")
        self.assertTrue(first.succeeded)
        self.assertTrue(first.geometry_sealed)
        self.assertGreater(first.material_faces.shape[0], 0)
        self.assertGreater(first.cavity_faces.shape[0], 0)
        np.testing.assert_array_equal(first.design_vertices, design_vertices)
        self.assertEqual(first.diagnostics["mode"], "STATIC_BALANCE")
        self.assertEqual(len(first.diagnostics["center_of_mass_model"]), 3)
        self.assertEqual(len(first.diagnostics["principal_axes_model"]), 3)
        self.assertFalse(first.outer_vertices.flags.writeable)
        with self.assertRaises(ValueError):
            first.outer_vertices[0, 0] = 42.0
        self.assertEqual(
            first.diagnostics["output_hash"],
            second.diagnostics["output_hash"],
        )
        self.assertRegex(first.diagnostics["trace_hash"], r"^[0-9a-f]{64}$")
        self.assertEqual(
            first.diagnostics["trace_hash_semantics"],
            "legacy_accepted_states_v1",
        )
        self.assertRegex(first.event_hash, r"^[0-9a-f]{64}$")
        self.assertEqual(first.diagnostics["event_stream_version"], 1)
        self.assertEqual(first.diagnostics["event_hash"], first.event_hash)
        self.assertEqual(
            first.solve_function_evaluations,
            len(first.evaluation_timings),
        )
        self.assertEqual(
            first.solve_function_evaluations,
            first.diagnostics["solve_function_evaluations"],
        )
        self.assertEqual(first.post_solve_material_evaluation_count, 2)
        self.assertEqual(
            first.diagnostics["post_solve_material_evaluation_count"], 2
        )
        self.assertGreaterEqual(first.solve_evaluator_wall_seconds, 0.0)
        self.assertGreaterEqual(
            first.post_solve_material_evaluation_wall_seconds, 0.0
        )
        self.assertTrue(first.events)
        self.assertEqual(first.events[0].type, "attempt_start")
        self.assertEqual(first.events[-1].type, "attempt_terminal")
        self.assertEqual(first.events[-1].reason, "converged")
        json.dumps(first.diagnostics, allow_nan=False)
        self.assertGreaterEqual(len(first.trace), 1)
        self.assertGreaterEqual(first.residual_evaluation_seconds, 0.0)
        self.assertGreaterEqual(first.derivative_evaluation_seconds, 0.0)

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first.save_outer_mesh(root / "outer.obj")
            first.save_material_mesh(root / "material.obj")
            first.save_cavity_mesh(root / "cavity.obj")
            first.save_design_mesh(root / "design.obj")
            first.save_diagnostics(root / "diagnostics.json")
            self.assertTrue(all((root / name).stat().st_size > 0 for name in (
                "outer.obj", "design.obj", "material.obj", "cavity.obj",
                "diagnostics.json"
            )))
            payload = json.loads((root / "diagnostics.json").read_text())
            self.assertTrue(payload["geometry_sealed"])


if __name__ == "__main__":
    unittest.main()
