from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
import hashlib
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
        legacy = solve_spin_it_faster(replace(request, settings=replace(config,
            equation_profile="LEGACY_NATIVE_NORMALIZED_V1")), initial_quadrics=initial)
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
        # Captured from the actual f433a14 binary before Task6e rebuilds.
        self.assertEqual(hashlib.sha256((legacy._native.diagnostics_json + "\n").encode()).hexdigest(),
                         "c744e1daf10d8349b42801e787a9b54193464facfb62725dcc14b098fcf04dd8")
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
