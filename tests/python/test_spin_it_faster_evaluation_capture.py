from __future__ import annotations

import unittest

import numpy as np

from make_it_stand.spin import (
    SpinItFasterRequest,
    SpinSolverConfig,
    solve_spin_it_faster,
    solve_spin_it_faster_with_external_moment_evaluator,
)


def _cube(minimum: float, maximum: float):
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
        (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4), (3, 7, 6), (3, 6, 2),
        (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5),
    )
    return vertices, faces


class SpinItFasterEvaluationCaptureTests(unittest.TestCase):
    def _request(self, *, capture: bool,
                 equation_profile: str = "DUAL_MONOMIAL_RAW_V1") -> SpinItFasterRequest:
        outer_vertices, outer_faces = _cube(-3.0, 3.0)
        design_vertices, design_faces = _cube(-1.0, 1.0)
        return SpinItFasterRequest(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            design_vertices=design_vertices,
            design_faces=design_faces,
            settings=SpinSolverConfig(
                equation_profile=equation_profile,
                maximum_restarts=1,
                maximum_function_evaluations_per_restart=1,
                capture_evaluation_payload=capture,
            ),
        )

    def test_opt_in_payload_captures_the_effective_solver_jacobian(self) -> None:
        request = self._request(capture=True)
        native_jacobian = np.eye(10, dtype=np.float64)
        raw_moments = np.arange(1.0, 11.0, dtype=np.float64)
        coefficients = (
            0.31, -0.28, 0.17, 0.11, -0.43,
            -0.36, -0.52, 0.22, -0.13, 0.19,
        )

        def evaluator(_coefficients):
            return {
                "raw_moments": raw_moments,
                "variable_raw_moments": raw_moments * 0.5,
                "moment_jacobian": native_jacobian,
                "derivative_defined": True,
                "derivative_diagnostic": "",
                "geometry_hash": "capture-geometry-v1",
                "topology_hash": "capture-topology-v1",
                "evaluation_map_identity": "a" * 64,
                "boolean_backend_identity": "capture-backend-v1",
                "derivative_identity": "capture-operational-field-v1",
            }

        result = solve_spin_it_faster_with_external_moment_evaluator(
            request,
            evaluator,
            initial_quadrics=(coefficients,),
        )
        self.assertEqual(result["event_stream_version"], 6)
        event = next(row for row in result["events"] if row["type"] == "evaluation")
        payload = event["evaluation_payload"]
        self.assertEqual(payload["schema"], "spin-evaluation-payload-v1")
        np.testing.assert_allclose(payload["solver_raw_moments"], raw_moments)
        np.testing.assert_allclose(payload["variable_raw_moments"], raw_moments * 0.5)
        self.assertTrue(payload["derivative_defined"])
        self.assertEqual(
            payload["jacobian_coordinate_identity"],
            "literal_monomial_unit_sphere_v1",
        )
        self.assertEqual(payload["evaluation_map_identity"], "a" * 64)
        self.assertEqual(payload["boolean_backend_identity"], "capture-backend-v1")
        self.assertEqual(payload["derivative_identity"], "capture-operational-field-v1")

        diagonal = np.asarray((1.0,) * 7 + (2.0,) * 3)
        coordinate_scale = np.diag(diagonal)
        mapped = coordinate_scale @ np.asarray(coefficients, dtype=np.float64)
        solver_coordinates = mapped / np.linalg.norm(mapped)
        inverse_mapped = np.diag(1.0 / diagonal) @ solver_coordinates
        geometry_coordinates = inverse_mapped / np.linalg.norm(inverse_mapped)
        geometry_jacobian = (
            np.eye(10) - np.outer(geometry_coordinates, geometry_coordinates)
        ) @ np.diag(1.0 / diagonal) / np.linalg.norm(inverse_mapped)
        np.testing.assert_allclose(
            payload["solver_coordinates"],
            solver_coordinates,
            rtol=1e-12,
            atol=1e-12,
        )
        np.testing.assert_allclose(
            payload["solver_moment_jacobian"],
            native_jacobian @ geometry_jacobian,
            rtol=1e-12,
            atol=1e-12,
        )
        self.assertEqual(payload["geometry_hash"], "capture-geometry-v1")
        self.assertEqual(payload["topology_hash"], "capture-topology-v1")
        self.assertGreater(event["evaluation_index"], 0)

    def test_rank_deficient_terminal_event_captures_actual_newton_system(self) -> None:
        request = self._request(capture=True)
        moments = np.arange(1.0, 11.0, dtype=np.float64)
        jacobian = np.diag((1.0,) + (0.0,) * 9)
        initial = (0.31, -0.28, 0.17, 0.11, -0.43,
                   -0.36, -0.52, 0.22, -0.13, 0.19)

        def evaluator(_coefficients):
            return {
                "raw_moments": moments,
                "moment_jacobian": jacobian,
                "derivative_defined": True,
                "geometry_hash": "rank-system-geometry-v1",
                "topology_hash": "rank-system-topology-v1",
                "evaluation_map_identity": "d" * 64,
                "boolean_backend_identity": "rank-system-backend-v1",
                "derivative_identity": "rank-system-derivative-v1",
            }

        result = solve_spin_it_faster_with_external_moment_evaluator(
            request, evaluator, initial_quadrics=(initial,),
        )
        terminal = next(event for event in result["events"]
                        if event.get("reason") == "bordered_system_rank_deficient")
        payload = terminal["newton_system_payload"]
        self.assertEqual(payload["schema"], "spin-newton-system-payload-v2")
        matrix = np.asarray(payload["matrix"], dtype=np.float64)
        rhs = np.asarray(payload["right_hand_side"], dtype=np.float64)
        self.assertEqual(matrix.shape[0], matrix.shape[1])
        self.assertEqual(rhs.shape, (matrix.shape[0],))
        self.assertLess(payload["computed_rank"], matrix.shape[0])
        self.assertGreaterEqual(payload["equilibrated_rank"], 0)
        self.assertLessEqual(payload["equilibrated_rank"], matrix.shape[0])
        self.assertTrue(np.isfinite(matrix).all())
        self.assertTrue(np.isfinite(rhs).all())
        self.assertIsNotNone(payload["multiplier_mu"])
        self.assertEqual(np.asarray(payload["row_scales"]).shape, (matrix.shape[0],))
        self.assertIn(np.asarray(payload["column_scales"]).shape,
                      {(0,), (matrix.shape[0],)})

    def test_undefined_derivative_is_null_but_same_call_response_is_retained(self) -> None:
        request = self._request(capture=True)
        raw_moments = np.linspace(0.25, 2.5, 10, dtype=np.float64)

        def evaluator(_coefficients):
            return {
                "raw_moments": raw_moments,
                "variable_raw_moments": raw_moments,
                "moment_jacobian": np.zeros((10, 10), dtype=np.float64),
                "derivative_defined": False,
                "derivative_diagnostic": "topology provenance unavailable",
                "geometry_hash": "undefined-geometry-v1",
                "topology_hash": "undefined-topology-v1",
                "evaluation_map_identity": "b" * 64,
                "boolean_backend_identity": "capture-backend-v1",
                "derivative_identity": "undefined-operational-field-v1",
            }

        result = solve_spin_it_faster_with_external_moment_evaluator(
            request,
            evaluator,
            initial_quadrics=((1.0, 0.0, 0.0, 0.0, -1.0,
                               -1.0, -1.0, 0.0, 0.0, 0.0),),
        )
        event = next(row for row in result["events"] if row["type"] == "evaluation")
        payload = event["evaluation_payload"]
        self.assertFalse(payload["derivative_defined"])
        self.assertIsNone(payload["solver_moment_jacobian"])
        np.testing.assert_allclose(payload["solver_raw_moments"], raw_moments)
        np.testing.assert_allclose(payload["variable_raw_moments"], raw_moments)

    def test_external_solver_streams_finalized_events_before_return(self) -> None:
        request = self._request(capture=False)
        coefficients = (
            0.31, -0.28, 0.17, 0.11, -0.43,
            -0.36, -0.52, 0.22, -0.13, 0.19,
        )
        streamed = []
        solver_returned = False

        def evaluator(_coefficients):
            return {
                "raw_moments": np.arange(1.0, 11.0),
                "moment_jacobian": np.eye(10),
                "derivative_defined": True,
                "geometry_hash": "stream-geometry-v1",
                "topology_hash": "stream-topology-v1",
                "evaluation_map_identity": "c" * 64,
                "boolean_backend_identity": "stream-backend-v1",
                "derivative_identity": "stream-derivative-v1",
            }

        def on_event(event):
            self.assertFalse(solver_returned)
            streamed.append(event)

        result = solve_spin_it_faster_with_external_moment_evaluator(
            request, evaluator, initial_quadrics=(coefficients,), event_sink=on_event,
        )
        solver_returned = True

        self.assertTrue(streamed)
        self.assertEqual([row["type"] for row in streamed],
                         [row["type"] for row in result["events"]])
        self.assertEqual(len(streamed), len(result["events"]))
        for observed, returned in zip(streamed, result["events"]):
            for field in (
                "restart", "iteration", "evaluation_index",
                "restart_evaluation_index", "outcome", "reason", "detail",
                "quadric_hash", "normalized_equation_residual",
                "physical_relative_residual", "objective", "multiplier_mu",
                "material_geometry_hash", "material_topology_hash",
                "equation_profile", "solver_coordinate_identity",
                "raw_equation_residual",
            ):
                self.assertEqual(observed[field], returned[field], field)
            if returned["quadric_coefficients"] is not None:
                np.testing.assert_array_equal(
                    observed["quadric_coefficients"], returned["quadric_coefficients"],
                )
        self.assertEqual(streamed[-1]["type"], "attempt_terminal")
        self.assertTrue(any(row.get("outcome") == "accepted" for row in streamed))

    def test_native_capture_separates_fixed_shell_from_variable_response(self) -> None:
        request = SpinItFasterRequest(
            **{
                **self._request(capture=True).__dict__,
                "fixed_shell_moments": (7.0,) + (0.0,) * 9,
            }
        )
        result = solve_spin_it_faster(
            request,
            initial_quadrics=((1.0, 0.0, 0.0, 0.0, -1.0,
                               -1.0, -1.0, 0.0, 0.0, 0.0),),
        )
        event = next(row for row in result.events if row.type == "evaluation")
        payload = event.evaluation_payload
        self.assertIsNotNone(payload)
        assert payload is not None
        self.assertIsNotNone(payload.variable_raw_moments)
        assert payload.variable_raw_moments is not None
        self.assertAlmostEqual(
            payload.solver_raw_moments[0] - payload.variable_raw_moments[0],
            7.0,
            places=10,
        )
        np.testing.assert_allclose(
            np.asarray(payload.solver_raw_moments[1:]),
            np.asarray(payload.variable_raw_moments[1:]),
            rtol=0.0,
            atol=1e-10,
        )

    def test_capture_disabled_keeps_legacy_event_stream_and_has_no_payload(self) -> None:
        request = self._request(capture=False)
        identity = np.eye(10, dtype=np.float64)

        def evaluator(_coefficients):
            return {
                "raw_moments": np.ones(10, dtype=np.float64),
                "moment_jacobian": identity,
                "derivative_defined": True,
                "geometry_hash": "default-geometry-v1",
                "topology_hash": "default-topology-v1",
            }

        result = solve_spin_it_faster_with_external_moment_evaluator(
            request,
            evaluator,
            initial_quadrics=((1.0, 0.0, 0.0, 0.0, -1.0,
                               -1.0, -1.0, 0.0, 0.0, 0.0),),
        )
        self.assertEqual(result.get("event_stream_version", 2), 2)
        event = next(row for row in result["events"] if row["type"] == "evaluation")
        self.assertNotIn("evaluation_payload", event)


if __name__ == "__main__":
    unittest.main()
