"""Freeze guards and actual geometry, never author-case completion tests."""
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import subprocess
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]


class FrozenFunctionBatch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = ROOT / "experiments/spin2014_objectives_quadric_v1/run.py"
        spec = importlib.util.spec_from_file_location("fixed_shape_function_runner", path)
        cls.runner = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.runner)

    def test_heart_uses_frozen_source_resolution_and_original_units(self):
        r, objective, case = self.runner.load_frozen_rotation(
            "heart-original-axis-three-objectives", "TOP_LOWERING_ONLY", 0)
        self.assertEqual(len(r.outer_vertices), 530)
        self.assertEqual(len(r.outer_faces), 1056)
        self.assertEqual(r.target_origin, (0.,0.,0.))
        self.assertEqual(r.settings.deterministic_seed, 0)
        np.testing.assert_allclose(r.design_vertices,
            (np.asarray(r.outer_vertices)-[0.,0.,.95])*.72+[0.,0.,.95], atol=0, rtol=0)
        self.assertEqual(objective.center_weight, 1.)
        self.assertEqual(objective.solver_length_unit_in_declared_units, 1.)
        self.assertEqual(case["design_scale"], .72)

    def test_cuboid_has_real_initial_nonprincipal_inertia(self):
        from make_it_stand.spin import integrate_raw_mass_moments
        r, objective, _ = self.runner.load_frozen_rotation(
            "cuboid-offaxis-pairs-yoyo", "YOYO", 1)
        s = integrate_raw_mass_moments(r.outer_vertices, r.outer_faces)
        self.assertAlmostEqual(s[0], 3.84, places=12)
        expected_xz = (0.4608-1.28)*math.sin(math.radians(5))*math.cos(math.radians(5))
        self.assertAlmostEqual(s[8], expected_xz, places=12)
        self.assertGreater(abs(s[8]), .07)
        self.assertEqual(len(r.design_vertices), 16)
        self.assertEqual(len(r.design_faces), 24)
        self.assertEqual(r.settings.deterministic_seed, 1)
        self.assertEqual(objective.center_weight, 0.)

    def test_changed_protocol_or_unregistered_strategy_is_rejected(self):
        path = ROOT / "experiments/spin2014_objectives_quadric_v1/protocol.json"
        with tempfile.TemporaryDirectory() as tmp:
            altered = Path(tmp)/"protocol.json"
            value = json.loads(path.read_text())
            value["cases"][0]["design_scale"] = .71
            altered.write_text(json.dumps(value))
            with self.assertRaises(ValueError):
                self.runner.load_frozen_rotation("heart-original-axis-three-objectives",
                    "TOP_INERTIA_ONLY", 0, protocol_path=altered)
        for case, objective, seed in (("heart-original-axis-three-objectives", "YOYO", 0),
                                     ("cuboid-offaxis-pairs-yoyo", "YOYO", 2)):
            with self.subTest(case=case, objective=objective, seed=seed), self.assertRaises(ValueError):
                self.runner.load_frozen_rotation(case, objective, seed)

    def test_static_loads_authenticated_geometry_and_actual_bottom_hull(self):
        # Break caught: old maximum-mass dispatch, changed target or a bounding
        # box standing in for actual support-face vertices.
        request, objective, case, support = self.runner.load_frozen_static(
            "static-faceted-shell-lowering", "STATIC_LOWERING", 0)
        self.assertEqual(request.mode, "STATIC_BALANCE")
        self.assertEqual(request.target_origin, (.05,0.,-1.2))
        self.assertEqual(objective.objective, "STATIC_LOWERING")
        self.assertEqual(objective.inertia_weight, 0.)
        self.assertEqual(objective.center_weight, 1.)
        self.assertEqual(case["fixed_material_density"], 1.)
        self.assertEqual(len(support), 4)
        self.assertEqual(set(map(tuple,support)), {(-1.2,-1.2,-1.2),
            (1.2,-1.2,-1.2),(1.2,1.2,-1.2),(-1.2,1.2,-1.2)})
        self.assertTrue(all(tuple(p) in request.outer_vertices for p in support))
        with self.assertRaises(ValueError):
            self.runner.load_frozen_static("static-faceted-shell-lowering", "TOP_LOWERING_ONLY", 0)

    def test_static_cli_saves_fresh_squared_objective_not_rotational_dispatch(self):
        # Break caught: refusing the registered static strategy or routing it
        # through full-rank Newton and replacing failed outcomes by a pass.
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)/"static"
            result=subprocess.run([sys.executable,str(ROOT/
                "experiments/spin2014_objectives_quadric_v1/run.py"),
                "--case","static-faceted-shell-lowering","--objective","STATIC_LOWERING",
                "--seed","0","--output",str(output)],capture_output=True,text=True,timeout=120)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            record=json.loads((output/"native-run.json").read_text())
            self.assertEqual(record["objective_settings"]["objective"],"STATIC_LOWERING")
            self.assertEqual(record["initialization"],"INPUT_FULL_MATERIAL_HEIGHT_V1")
            self.assertIn("objective_hessian",record["squared_height"])
            self.assertIn("trace",record["candidate"])
            self.assertTrue((output/"begin-evidence.json").is_file())
            self.assertTrue((output/"run-evidence.json").is_file())
            summary=json.loads((output/"summary.json").read_text())
            self.assertEqual(summary["solver_family"],"PLANE_FAMILY_FRESH_2014_SQUARE_V1")
            self.assertFalse(summary["goal_complete"])
            self.assertTrue(summary["material_export_complete"])
            self.assertTrue(summary["material_geometry_attachment_verified"])
            self.assertFalse(summary["manufacturing_approved"])
            from make_it_stand.verification.spin_common import digest
            from make_it_stand.verification.spin_bundle import load_spin2014_static_material_output
            geometry = load_spin2014_static_material_output(output/"material-geometry",
                expected_source_run_digest=digest(record))
            self.assertTrue(geometry.audit["design_partition"]["complete_partition"])
            self.assertFalse(geometry.audit["sealed_artifact"])
            inputs=json.loads((output/"frozen-input.json").read_text())
            self.assertIsNone(inputs["evaluator"])
            self.assertEqual(inputs["initialization"],"DETERMINISTIC_FROM_INPUT")

    def test_saved_static_material_export_does_not_run_optimizer_again(self):
        # Break caught: repeating the optimization just to export existing
        # material, or overwriting the original run/evidence to grant a pass.
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            source=root/"source"
            runner=ROOT/"experiments/spin2014_objectives_quadric_v1/run.py"
            first=subprocess.run([sys.executable,str(runner),"--case","static-faceted-shell-lowering",
                "--objective","STATIC_LOWERING","--seed","0","--output",str(source)],
                capture_output=True,text=True,timeout=120)
            self.assertEqual(first.returncode,0,first.stdout+first.stderr)
            original=(source/"native-run.json").read_bytes()
            result=subprocess.run([sys.executable,str(runner),"--saved-static-run",str(source),
                "--output",str(root/"export")],capture_output=True,text=True,timeout=120)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            report=json.loads((root/"export"/"summary.json").read_text())
            self.assertEqual(report["optimizer_invocations"],0)
            self.assertTrue(report["material_geometry_attachment_verified"])
            self.assertFalse(report["goal_complete"])
            self.assertFalse(report["manufacturing_approved"])
            self.assertEqual((source/"native-run.json").read_bytes(),original)
            self.assertTrue((root/"export"/"run-evidence.json").is_file())

    def test_run_parameters_use_evidence_json_values_not_dataclass_tuples(self):
        from make_it_stand.spin import SpinEvaluatorConfig, _native_module
        from make_it_stand.verification.run_evidence import _nonempty_object
        request, objective, case = self.runner.load_frozen_rotation(
            "heart-original-axis-three-objectives", "TOP_INERTIA_ONLY", 0)
        inputs = self.runner.frozen_run_inputs(request, objective, case, 0,
            Path(_native_module().__file__).resolve(), SpinEvaluatorConfig(
                derivative_method="DISCRETE_FULL_RANK_QUADRIC",
                icosphere_subdivisions=2, maximum_icosphere_subdivisions=2))
        # The real evidence validator rejects tuples, even though json.dump
        # would silently serialize them. Validate what begin_run receives.
        _nonempty_object(inputs, "parameters")
        self.assertIsInstance(inputs["request"]["outer_vertices"], list)
        self.assertEqual(inputs["objective"]["objective"], "TOP_INERTIA_ONLY")
        self.assertEqual(inputs["fixed_material_density"], 1.)

    def test_source_audit_opt_in_is_bound_without_changing_frozen_science(self):
        # Break caught: applying the new source policy silently to historical
        # runs, failing to bind its true value, or changing the frozen request
        # while choosing the admission policy. No optimizer is called.
        from make_it_stand.spin import _native_module
        request, objective, case = self.runner.load_frozen_rotation(
            "heart-original-axis-three-objectives", "TOP_INERTIA_ONLY", 0)
        native_path = Path(_native_module().__file__).resolve()
        old = self.runner.frozen_run_inputs(request, objective, case, 0, native_path,
            self.runner.frozen_rotation_evaluator())
        exact = self.runner.frozen_run_inputs(request, objective, case, 0, native_path,
            self.runner.frozen_rotation_evaluator(exact_source_boundary=True))
        self.assertNotIn("exact_binary64_source_boundary", old["evaluator"])
        self.assertTrue(exact["evaluator"]["exact_binary64_source_boundary"])
        for key in ("request", "objective", "seed", "protocol_sha256", "fixed_material_density"):
            self.assertEqual(old[key], exact[key], key)
        self.assertEqual(exact["evaluator"]["icosphere_subdivisions"], 2)
        self.assertEqual(exact["request"]["settings"]["maximum_restarts"], 8)

    def test_rational_production_policy_reaches_frozen_request_without_science_changes(self):
        # Break caught: the batch runner silently leaves the real producer on
        # binary64 source coordinates, omits its policy identity, or changes
        # the frozen geometry/weights/initialization with the precision opt-in.
        from make_it_stand.spin import _native_module
        request, objective, case = self.runner.load_frozen_rotation(
            "cuboid-offaxis-pairs-yoyo", "YOYO", 1)
        native_path = Path(_native_module().__file__).resolve()
        old = self.runner.frozen_run_inputs(request, objective, case, 1, native_path,
            self.runner.frozen_rotation_evaluator())
        settings = self.runner.frozen_rotation_evaluator(
            exact_source_boundary=True, exact_material_boundary=True,
            rational_source_coordinates=True, compact_rational_source_coordinates=True)
        new = self.runner.frozen_run_inputs(request, objective, case, 1, native_path, settings)
        self.assertNotIn("rational_source_coordinates", old["evaluator"])
        self.assertNotIn("compact_rational_source_coordinates", old["evaluator"])
        self.assertTrue(new["evaluator"]["rational_source_coordinates"])
        self.assertTrue(new["evaluator"]["compact_rational_source_coordinates"])
        for key in ("request", "objective", "seed", "protocol_sha256", "fixed_material_density", "initialization"):
            self.assertEqual(old[key], new[key], key)
        self.assertEqual(settings.icosphere_subdivisions, 2)
        self.assertEqual(settings.maximum_icosphere_subdivisions, 2)
        self.assertEqual(request.settings.maximum_restarts, 8)
        self.assertEqual(request.settings.maximum_function_evaluations_per_restart, 60)

    def test_compact_batch_policy_requires_explicit_rational_source(self):
        with self.assertRaises(ValueError):
            self.runner.frozen_rotation_evaluator(compact_rational_source_coordinates=True)

    def test_static_cli_rejects_unused_rational_policy_before_writing_outputs(self):
        # Break caught: a plane result borrows a precision profile that the
        # plane path never executes, or writes output before rejecting it.
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)/"static"
            result = subprocess.run([sys.executable, str(ROOT/
                "experiments/spin2014_objectives_quadric_v1/run.py"),
                "--case", "static-faceted-shell-lowering", "--objective", "STATIC_LOWERING",
                "--seed", "0", "--rational-source-coordinates", "--output", str(output)],
                capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 2)
            self.assertIn("rotational", result.stderr)
            self.assertNotIn("unrecognized arguments", result.stderr)
            self.assertFalse(output.exists())

    def test_material_audit_opt_in_is_bound_without_changing_frozen_science(self):
        # Break caught: admitting changed geometry/weights/budgets while opting
        # into the complete exact signed-shell material check, or changing old
        # saved identities simply because a new false default was introduced.
        from make_it_stand.spin import _native_module
        request, objective, case = self.runner.load_frozen_rotation(
            "heart-original-axis-three-objectives", "TOP_INERTIA_ONLY", 0)
        loaded = Path(_native_module().__file__).resolve()
        old = self.runner.frozen_run_inputs(request, objective, case, 0, loaded,
            self.runner.frozen_rotation_evaluator())
        exact = self.runner.frozen_run_inputs(request, objective, case, 0, loaded,
            self.runner.frozen_rotation_evaluator(exact_material_boundary=True))
        self.assertNotIn("exact_binary64_material_boundary", old["evaluator"])
        self.assertTrue(exact["evaluator"]["exact_binary64_material_boundary"])
        for key in ("request", "objective", "seed", "protocol_sha256", "fixed_material_density"):
            self.assertEqual(old[key], exact[key], key)
        self.assertEqual(exact["evaluator"]["icosphere_subdivisions"], 2)
        self.assertEqual(exact["request"]["settings"]["maximum_restarts"], 8)

    def test_timing_hashes_do_not_change_scientific_events_or_hide_failures(self):
        one = {"solve": {"output_hash": 100, "event_hash": 10,
            "evaluation_timings": [{"evaluator_wall_seconds": .1}],
            "status": "RESTARTS_EXHAUSTED", "function_evaluations": 2,
            "events": [{"outcome": "failed", "raw_residual": .25,
                        "newton_system_payload": {"matrix": [[1., 2.]]}}]}}
        two = json.loads(json.dumps(one))
        two["solve"]["output_hash"] = 200
        two["solve"]["event_hash"] = 20
        two["solve"]["evaluation_timings"][0]["evaluator_wall_seconds"] = .2
        expected = {"solve": {"status": "RESTARTS_EXHAUSTED",
            "function_evaluations": 2,
            "events": [{"outcome": "failed", "raw_residual": .25,
                        "newton_system_payload": {"matrix": [[1., 2.]]}}]}}
        self.assertEqual(self.runner._scientific(one), expected)
        self.assertEqual(self.runner._scientific(two), expected)

    def test_explicit_science_policy_normalizes_only_oriented_cavity_cycles(self):
        # Break caught: raw CGAL face start choices change scientific identity,
        # or the opt-in silently changes legacy output or actual input arrays.
        one = {"cavity_faces": [[60, 58, 35], [63, 64, 66]],
               "cavity_vertices": [[1., 2., 3.]], "solve": {"status": "CONVERGED"}}
        two = {"cavity_faces": [[35, 60, 58], [64, 66, 63]],
               "cavity_vertices": [[1., 2., 3.]], "solve": {"status": "CONVERGED"}}
        expected = {"cavity_faces": [[35, 60, 58], [63, 64, 66]],
                    "cavity_vertices": [[1., 2., 3.]], "solve": {"status": "CONVERGED"}}
        self.assertEqual(self.runner._scientific(one, canonical_cavity_cycles=True), expected)
        self.assertEqual(self.runner._scientific(two, canonical_cavity_cycles=True), expected)
        self.assertEqual(self.runner._scientific(one), one)
        self.assertEqual(one["cavity_faces"], [[60, 58, 35], [63, 64, 66]])

    def test_science_cycle_policy_keeps_winding_missing_faces_and_coordinates_visible(self):
        # Break caught: unoriented sorting or dropped geometry hides a real
        # change instead of only removing an arbitrary triangle start index.
        one = {"cavity_faces": [[2, 0, 1]], "cavity_vertices": [[1., 2., 3.]]}
        baseline = self.runner._scientific(one, canonical_cavity_cycles=True)
        for changed in ({"cavity_faces": [[2, 1, 0]], "cavity_vertices": [[1., 2., 3.]]},
                        {"cavity_faces": [], "cavity_vertices": [[1., 2., 3.]]},
                        {"cavity_faces": [[2, 0, 1]], "cavity_vertices": [[1., 2., 3.1]]}):
            with self.subTest(changed=changed):
                self.assertNotEqual(self.runner._scientific(
                    changed, canonical_cavity_cycles=True), baseline)

    def test_science_cycle_policy_rejects_invalid_triangles_without_erasing_them(self):
        for faces in ([[0, 1]], [[0, 1, 1]], [[0, True, 2]], [[0, -1, 2]], [[0, 1., 2]]):
            with self.subTest(faces=faces), self.assertRaises(ValueError):
                self.runner._scientific({"cavity_faces": faces}, canonical_cavity_cycles=True)


if __name__ == "__main__":
    unittest.main()
