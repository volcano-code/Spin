from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from make_it_stand.spin import (
    SpinItFasterRequest,
    SpinSolverConfig,
    evaluate_spin_it_faster_quadric,
    integrate_raw_mass_moments,
    solve_spin_it_faster,
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


class SpinSchemaTests(unittest.TestCase):
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


class SpinNativeArtifactTests(unittest.TestCase):
    def test_raw_moment_integration_matches_cube(self) -> None:
        vertices, faces = cube(-1.0, 1.0)
        moments = integrate_raw_mass_moments(vertices, faces)
        np.testing.assert_allclose(
            moments,
            (8.0, 0.0, 0.0, 0.0, 8.0 / 3.0, 8.0 / 3.0,
             8.0 / 3.0, 0.0, 0.0, 0.0),
            atol=1e-14,
        )

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
