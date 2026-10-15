"""Live authority events must survive interruption, not just native return."""
import json
from pathlib import Path
import tempfile
import unittest
import importlib.util
import numpy as np
import _make_it_stand as native


def inputs():
    vertices = np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
                         [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]], dtype=float)
    faces = np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
                      [3,7,6],[3,6,2],[0,4,7],[0,7,3],[1,2,6],[1,6,5]], dtype=np.int32)
    outer, design = native.SpinTriangleMesh(), native.SpinTriangleMesh()
    outer.vertices, outer.faces = vertices*3, faces
    design.vertices, design.faces = vertices*2, faces
    request = native.SpinDesignRequest()
    request.outer_surface, request.design_domain = outer, design
    request.material_density = 3
    request.mode = native.SpinMode.SPINNING_TOP
    settings = request.settings
    settings.maximum_restarts = 1
    request.settings = settings
    evaluation = native.SpinEvaluatorSettings()
    tri = evaluation.triangulation
    tri.icosphere_subdivisions = 0  # interface test, not functional level-2 evidence
    tri.rational_source_coordinates = True
    tri.compact_rational_source_coordinates = True
    tri.truncation_radius = 8
    tri.maximum_radius_doublings = 0
    evaluation.triangulation = tri
    evaluation.derivative_method = native.SpinDerivativeMethod.DISCRETE_FULL_RANK_QUADRIC
    evaluation.exact_binary64_source_boundary = True
    evaluation.exact_binary64_material_boundary = True
    return request, native.Spin2014QuadricObjectiveSettings(), evaluation


class DurableAuthorityEvents(unittest.TestCase):
    def test_public_policy_forwards_live_authority_observer(self):
        # Break caught: selecting authority native but forgetting its observer.
        from make_it_stand.spin import (SpinItFasterRequest, SpinSolverConfig,
            Spin2014QuadricObjectiveConfig, solve_spin2014_quadric_hollowing)
        request, _, _ = inputs()
        public = SpinItFasterRequest(
            tuple(tuple(float(x) for x in row) for row in request.outer_surface.vertices),
            tuple(tuple(int(x) for x in row) for row in request.outer_surface.faces),
            tuple(tuple(float(x) for x in row) for row in request.design_domain.vertices),
            tuple(tuple(int(x) for x in row) for row in request.design_domain.faces),
            mode="SPINNING_TOP", material_density=3.,
            settings=SpinSolverConfig(maximum_restarts=1, equation_profile="DUAL_MONOMIAL_RAW_V1"))
        seen = []
        def abort(event):
            seen.append(event)
            raise RuntimeError("public-authority-observer-stop")
        with self.assertRaisesRegex(RuntimeError, "public-authority-observer-stop"):
            solve_spin2014_quadric_hollowing(public, Spin2014QuadricObjectiveConfig(),
                fixed_material_density=2., coefficient_policy="COEFFICIENT_MPFR1024_V1",
                event_sink=abort)
        self.assertEqual(len(seen), 1)
        self.assertEqual(seen[0]["type"], "attempt_start")
        self.assertEqual(len(seen[0]["geometry_authority"]), 10)

    def test_observer_abort_preserves_start_before_initial_boolean(self):
        # Break caught: emitting only after native return or swallowing observer exceptions.
        path = Path(__file__).resolve().parents[2]/"experiments/spin2014_objectives_quadric_v1/run.py"
        spec = importlib.util.spec_from_file_location("authority_event_runner", path)
        runner = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(runner)
        request, objective, evaluation = inputs()
        with tempfile.TemporaryDirectory() as directory:
            journal = Path(directory)/"events.jsonl"
            persist = runner.durable_rotational_event_sink(journal)
            def abort(event):
                persist(event)
                saved = json.loads(journal.read_text())
                self.assertEqual(saved["event"]["type"], "attempt_start")
                self.assertEqual(saved["event"]["evaluation_index"], 0)
                self.assertEqual(len(saved["event"]["geometry_authority"]), 10)
                for key in ("context_identity", "request_identity", "objective_identity"):
                    self.assertEqual(len(saved["event"][key]), 64)
                self.assertEqual(saved["event"]["coefficient_policy"], "COEFFICIENT_MPFR1024_V1")
                raise RuntimeError("authority-observer-interrupted")
            with self.assertRaisesRegex(RuntimeError, "authority-observer-interrupted"):
                native.solve_spin2014_authority_quadric_hollowing(request, objective, 2.,
                    evaluation, evaluation_budget_per_attempt=1, event_sink=abort)
            rows = journal.read_text().splitlines()
            self.assertEqual(len(rows), 1)
            self.assertEqual(json.loads(rows[0])["sequence"], 0)

    def test_complete_failed_attempt_events_equal_live_stream(self):
        # Break caught: dropping failed attempts or altering science when observed.
        request, objective, evaluation = inputs()
        seen = []
        with tempfile.TemporaryDirectory() as directory:
            journal = Path(directory)/"completed.jsonl"
            with journal.open("x") as stream:
                def persist(event):
                    seen.append(event)
                    stream.write(json.dumps(dict(schema="spin2014-solver-event-journal-v1",
                        sequence=len(seen)-1,event=event))+"\n")
                    stream.flush()
                observed = native.solve_spin2014_authority_quadric_hollowing(request, objective, 2.,
                    evaluation, evaluation_budget_per_attempt=1, event_sink=persist)
            from make_it_stand.verification._spin2014_authority import check_authority_event_journal
            checked = check_authority_event_journal(journal,run=observed,
                expected_request_identity=observed["request_identity"],
                expected_context_identity=observed["context_identity"],
                expected_objective_identity=seen[0]["objective_identity"])
            self.assertEqual(checked["records"], len(seen))
            self.assertFalse(checked["sealed_artifact"])
        self.assertEqual(observed["events"], seen)
        self.assertEqual([event["sequence"] for event in seen], list(range(len(seen))))
        self.assertEqual(seen[0]["type"], "attempt_start")
        self.assertEqual(seen[-2]["type"], "attempt_terminal")
        self.assertEqual(seen[-1]["type"], "run_terminal")
        self.assertEqual(seen[-1]["evaluation_index"], observed["function_evaluations"])
        plain = native.solve_spin2014_authority_quadric_hollowing(request, objective, 2.,
            evaluation, evaluation_budget_per_attempt=1)
        self.assertEqual(plain["scientific_hash"], observed["scientific_hash"])
        self.assertEqual(plain["attempts"], observed["attempts"])
        self.assertEqual(plain["events"], observed["events"])


if __name__ == "__main__":
    unittest.main()
