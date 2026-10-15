"""Real public interface checks, not three-function completion evidence."""
import unittest
from dataclasses import replace
import numpy as np
from make_it_stand.spin import _immutable_spin_solver_event

from make_it_stand.spin import (
    SpinItFasterRequest, SpinSolverConfig, SpinEvaluatorConfig,
    Spin2014QuadricObjectiveConfig, solve_spin2014_quadric_hollowing,
    solve_spin2014_static_quadric_hollowing,
)


def box(a, b):
    v = ((a,a,a),(b,a,a),(b,b,a),(a,b,a),
         (a,a,b),(b,a,b),(b,b,b),(a,b,b))
    f = ((0,2,1),(0,3,2),(4,5,6),(4,6,7),(0,1,5),(0,5,4),
         (3,7,6),(3,6,2),(0,4,7),(0,7,3),(1,2,6),(1,6,5))
    return v, f


class Quadric2014Interface(unittest.TestCase):
    def test_explicit_authority_policy_uses_real_seed_restart_entry(self):
        request = self.request()
        request = replace(request, settings=replace(request.settings,
            maximum_function_evaluations_per_restart=60))
        result = solve_spin2014_quadric_hollowing(request,
            Spin2014QuadricObjectiveConfig(objective="TOP_INERTIA_ONLY"),
            fixed_material_density=3.,
            coefficient_policy="COEFFICIENT_MPFR1024_V1",
            authority_evaluation_budget_per_attempt=1)
        self.assertEqual(result["schema"], "spin2014-authority-restarts-v1")
        self.assertEqual(result["scope"], "NUMERICAL_CANDIDATE_NOT_FINAL_SEAL")
        self.assertEqual(result["coefficient_policy"], "COEFFICIENT_MPFR1024_V1")
        self.assertEqual(result["precision_bits"], 1024)
        self.assertEqual(len(result["attempts"]), 1)
        self.assertLessEqual(result["function_evaluations"], 1)
        self.assertEqual(len(result["attempts"][0]["initial_geometry_authority"]), 10)
        self.assertIsNotNone(result["attempts"][0]["initial_total_moments"])

    def test_authority_policy_rejects_answer_initialization_and_unknown_policy(self):
        kwargs = dict(fixed_material_density=3., coefficient_policy="COEFFICIENT_MPFR1024_V1")
        config = Spin2014QuadricObjectiveConfig(objective="TOP_INERTIA_ONLY")
        with self.assertRaisesRegex(ValueError, "seed-only"):
            solve_spin2014_quadric_hollowing(self.request(), config,
                initial_quadrics=((1.,)*10,), **kwargs)
        with self.assertRaisesRegex(ValueError, "coefficient_policy"):
            solve_spin2014_quadric_hollowing(self.request(), config,
                fixed_material_density=3., coefficient_policy="silently-fallback")

    def test_durable_rotational_journal_survives_observer_abort(self):
        # Break caught: buffering a completed event until native return, or
        # overwriting an earlier failed run's journal on re-entry.
        import importlib.util
        import json
        from pathlib import Path
        import tempfile
        path = Path(__file__).resolve().parents[2]/"experiments/spin2014_objectives_quadric_v1/run.py"
        spec = importlib.util.spec_from_file_location("durable_rotation_runner", path)
        runner = importlib.util.module_from_spec(spec); spec.loader.exec_module(runner)
        with tempfile.TemporaryDirectory() as tmp:
            journal = Path(tmp)/"solver-events.jsonl"
            sink = runner.durable_rotational_event_sink(journal)
            def persist_then_abort(event):
                sink(event)
                row = json.loads(journal.read_text())
                self.assertEqual(row["sequence"], 0)
                self.assertEqual(row["event"]["type"], "attempt_start")
                raise RuntimeError("literal-durable-interruption")
            with self.assertRaisesRegex(RuntimeError, "literal-durable-interruption"):
                solve_spin2014_quadric_hollowing(self.request(),
                    Spin2014QuadricObjectiveConfig(objective="TOP_INERTIA_ONLY"),
                    fixed_material_density=3., event_sink=persist_then_abort)
            row = json.loads(journal.read_text())
            self.assertEqual(row["schema"], "spin2014-solver-event-journal-v1")
            self.assertEqual(len(row["event"]["quadric_coefficients"]), 10)
            original = journal.read_bytes()
            with self.assertRaises(FileExistsError):
                runner.durable_rotational_event_sink(journal)
            self.assertEqual(journal.read_bytes(), original)

    def test_live_observer_receives_start_before_a_solve_can_be_interrupted(self):
        # Break caught: keeping every event in memory until the whole solve
        # returns, so a resource cut loses even the actual attempt start.
        seen = []
        def stop_after_start(event):
            seen.append(event)
            raise RuntimeError("literal-observer-stop")
        with self.assertRaisesRegex(RuntimeError, "literal-observer-stop"):
            solve_spin2014_quadric_hollowing(self.request(),
                Spin2014QuadricObjectiveConfig(objective="TOP_INERTIA_ONLY"),
                fixed_material_density=3., event_sink=stop_after_start)
        self.assertEqual(len(seen), 1)
        self.assertEqual(seen[0]["type"], "attempt_start")
        self.assertEqual(seen[0]["evaluation_index"], 0)
        self.assertEqual(seen[0]["restart"], 0)
        self.assertEqual(len(seen[0]["quadric_coefficients"]), 10)

    def test_live_observer_preserves_actual_budget_exit_and_complete_event_stream(self):
        # Break caught: forwarding the callback only on success, dropping
        # failed/budget events, or changing science when recording is enabled.
        from dataclasses import asdict
        from make_it_stand.verification.faster_material import _json_value
        q = (-.31,.07,-.09,.11,.83,1.17,-.72,.06,-.04,.03)
        settings = SpinEvaluatorConfig(derivative_method="DISCRETE_FULL_RANK_QUADRIC",
            icosphere_subdivisions=2, truncation_radius=8.)
        config = Spin2014QuadricObjectiveConfig(objective="TOP_MASS_AND_LOWERING")
        kwargs = dict(fixed_material_density=3., evaluator_settings=settings,
                      initial_quadrics=(q,))
        original = solve_spin2014_quadric_hollowing(self.request(), config, **kwargs)
        events = []
        observed = solve_spin2014_quadric_hollowing(self.request(), config,
                                                   event_sink=events.append, **kwargs)
        self.assertEqual(observed.solve.function_evaluations, 1)
        self.assertEqual([_json_value(asdict(_immutable_spin_solver_event(e))) for e in events],
                         [_json_value(asdict(e)) for e in observed.solve.events])
        self.assertEqual(events[0]["type"], "attempt_start")
        self.assertEqual(events[-1]["type"], "attempt_terminal")
        self.assertEqual(observed.solve.status, original.solve.status)
        self.assertEqual(observed.solve.input_hash, original.solve.input_hash)
        self.assertEqual(observed.solve.trace_hash, original.solve.trace_hash)
        def timing_free(stream):
            rows = [_json_value(asdict(e)) for e in stream]
            for row in rows:
                if row["evaluation_payload"] is not None:
                    for name in ("evaluator_wall_seconds", "residual_evaluation_seconds",
                                 "derivative_evaluation_seconds"):
                        row["evaluation_payload"].pop(name)
            return rows
        self.assertEqual(timing_free(observed.solve.events), timing_free(original.solve.events))
        self.assertFalse(observed.material_export_complete)
        with self.assertRaisesRegex(TypeError, "event_sink"):
            solve_spin2014_quadric_hollowing(self.request(), config,
                                            event_sink=1, **kwargs)

    def test_saved_unit_coefficients_replay_same_endpoint_without_renormalizing(self):
        # Break caught: replaying an already-normalized solver endpoint through
        # raw-input normalization changes its represented coefficients/geometry.
        # This literal is a saved endpoint regression, NEVER an initial guess.
        import importlib.util
        from pathlib import Path
        from make_it_stand.spin import evaluate_spin_it_faster_quadric
        path = Path(__file__).resolve().parents[2]/"experiments/spin2014_objectives_quadric_v1/run.py"
        spec = importlib.util.spec_from_file_location("cuboid_saved_unit_replay", path)
        loader = importlib.util.module_from_spec(spec); spec.loader.exec_module(loader)
        request, _, _ = loader.load_frozen_rotation(
            "cuboid-offaxis-pairs-yoyo", "YOYO", 1)
        q = (-8.13645371385572e-25, 2.614050534547747e-15,
             -3.034437291903922e-17, 3.2530544979661192e-15,
             .12168606307734217, .26311511023847906, -.9276739278498335,
             3.2035652688754162e-18, -.23533810657135806, -1.1194757484761892e-18)
        settings = loader.frozen_rotation_evaluator(
            exact_source_boundary=True, exact_material_boundary=True)
        replay = evaluate_spin_it_faster_quadric(request, q,
            evaluator_settings=settings, coefficients_are_normalized=True)
        self.assertEqual(replay.geometry_hash,
            "ce7b18b6a63d5c1ba9be011f10f34a997bc26da198e43392a1ff2ed8b1004a8d")
        with self.assertRaises(ValueError):
            evaluate_spin_it_faster_quadric(request, [2*x for x in q],
                evaluator_settings=settings, coefficients_are_normalized=True)
        legacy = evaluate_spin_it_faster_quadric(request, [2*x for x in q],
            evaluator_settings=settings)
        self.assertGreater(legacy.variable_raw_moments[0], 0.)

    def test_exact_motion_tape_follows_actual_material_vertex_permutation(self):
        # Break caught: canonicalizing the material but leaving its SAME-call
        # rational birth coordinates in their pre-sort order. Saved seed-1
        # evaluation 67 is only a regression, never solver initialization.
        import importlib.util
        from pathlib import Path
        from make_it_stand.spin import evaluate_spin_it_faster_fullrank_quadric
        path = Path(__file__).resolve().parents[2]/"experiments/spin2014_objectives_quadric_v1/run.py"
        spec = importlib.util.spec_from_file_location("heart_assignment_regression", path)
        loader = importlib.util.module_from_spec(spec); spec.loader.exec_module(loader)
        request, _, _ = loader.load_frozen_rotation(
            "heart-original-axis-three-objectives", "TOP_INERTIA_ONLY", 1)
        q = (-.4040543998331367, -7.61118274051472e-15, 5.465270644114226e-16,
             .7721732076644878, .14774103340530337, .2873339753422354,
             -.368917813344626, -2.602763503339779e-16,
             3.3215920215046725e-15, -2.285645306079095e-18)
        settings = loader.frozen_rotation_evaluator(
            exact_source_boundary=True, exact_material_boundary=True)
        value = evaluate_spin_it_faster_fullrank_quadric(request, q, evaluator_settings=settings)
        self.assertGreater(value.variable_raw_moments[0], 0)
        self.assertTrue(value.derivative_defined, value.derivative_diagnostic)

    def test_exact_motion_preserves_actual_heart_birth_intersection(self):
        # Break caught: assembling facet motion at rounded Boolean coordinates
        # creates inconsistent equations although the actual rational birth
        # intersection has a full-rank, exactly consistent differential.
        # Saved trial coefficients are a regression, NEVER solver initialization.
        import importlib.util
        from pathlib import Path
        from make_it_stand.spin import evaluate_spin_it_faster_fullrank_quadric
        path = Path(__file__).resolve().parents[2]/"experiments/spin2014_objectives_quadric_v1/run.py"
        spec = importlib.util.spec_from_file_location("heart_motion_regression", path)
        loader = importlib.util.module_from_spec(spec); spec.loader.exec_module(loader)
        request, _, _ = loader.load_frozen_rotation(
            "heart-original-axis-three-objectives", "TOP_INERTIA_ONLY", 0)
        q = (-.4040543999518238, -5.297471429422767e-14,
             -8.395705814510158e-14, .7721732072931566,
             .14774103273330358, .2873339771134668, -.3689178128814416,
             3.467452526257984e-14, 2.1981390467341484e-14,
             5.552599095382224e-17)
        settings = loader.frozen_rotation_evaluator(
            exact_source_boundary=True, exact_material_boundary=True)
        value = evaluate_spin_it_faster_fullrank_quadric(request, q,
            evaluator_settings=settings)
        self.assertEqual(value.geometry_hash,
            "6c990e63d4416516d2d39678b1c0711f8d8d633b7a556115d6209dc47f4ac9a7")
        self.assertAlmostEqual(value.variable_raw_moments[0],
            .20331713298374085, places=12)
        self.assertTrue(value.derivative_defined, value.derivative_diagnostic)

    def test_exact_material_policy_recovers_actual_nonzero_heart_boundary(self):
        # Break caught: a global minimum-area heuristic rejects the actual
        # evaluation-33 boundary despite complete exact positive-boundary proof.
        # This is a saved coefficient diagnostic, never a solver initial guess.
        import importlib.util
        from pathlib import Path
        from make_it_stand.spin import (
            _native_request, _native_evaluator_settings,
            evaluate_spin_it_faster_fullrank_quadric)
        path = Path(__file__).resolve().parents[2]/"experiments/spin2014_objectives_quadric_v1/run.py"
        spec = importlib.util.spec_from_file_location("heart_material_audit_regression", path)
        loader = importlib.util.module_from_spec(spec); spec.loader.exec_module(loader)
        request, _, _ = loader.load_frozen_rotation(
            "heart-original-axis-three-objectives", "TOP_INERTIA_ONLY", 0)
        q = (-.4040544018473746, -7.223110007033643e-9,
             5.17806772037275e-10, .7721732004222513, .14774101952838603,
             .28733401035181794, -.3689178045869738,
             -1.6366680127308427e-10, 3.3806147526531734e-9,
             -5.533845411071279e-14)
        legacy = SpinEvaluatorConfig(derivative_method="DISCRETE_FULL_RANK_QUADRIC",
            icosphere_subdivisions=2, maximum_icosphere_subdivisions=2,
            exact_binary64_self_intersection=True, exact_binary64_source_boundary=True)
        native, nr = _native_request(request)
        with self.assertRaisesRegex(RuntimeError, "degenerate_faces=168"):
            native.evaluate_spin_it_faster_quadric(nr, q, _native_evaluator_settings(native, legacy))
        exact = replace(legacy, exact_binary64_material_boundary=True)
        packet = native.evaluate_spin_it_faster_quadric(nr, q, _native_evaluator_settings(native, exact))
        proof = native.audit_exact_triangle_boundary(
            packet["material_vertices"], packet["material_faces"], "OutwardPositive")
        heuristic = native.audit_triangle_boundary(packet["material_vertices"], packet["material_faces"])
        self.assertTrue(proof["boundary_valid"], proof)
        self.assertEqual(proof["degenerate_faces"], [])
        self.assertEqual(heuristic["degenerate_face_count"], 168)
        value = evaluate_spin_it_faster_fullrank_quadric(request, q, evaluator_settings=exact)
        self.assertTrue(value.derivative_defined, value.derivative_diagnostic)
        self.assertEqual(value.geometry_hash, packet["geometry_hash"])

    def test_exact_material_policy_is_explicit_and_preserves_legacy_payload(self):
        from dataclasses import asdict
        legacy = SpinEvaluatorConfig()
        payload = asdict(legacy)
        payload.pop("exact_binary64_source_boundary", None)
        payload.pop("exact_binary64_material_boundary", None)
        payload.pop("rational_source_coordinates", None)
        payload.pop("compact_rational_source_coordinates", None)
        self.assertEqual(legacy.identity_dict(), payload)
        exact = replace(legacy, exact_binary64_material_boundary=True)
        self.assertTrue(exact.identity_dict()["exact_binary64_material_boundary"])
        for value in (1, "true", None):
            with self.subTest(value=value), self.assertRaises(ValueError):
                replace(legacy, exact_binary64_material_boundary=value)

    def test_exact_source_policy_is_explicit_and_preserves_legacy_identity_payload(self):
        # Break caught: silently adding a false default to old saved identities,
        # or coercing a string/integer into permission to bypass a source audit.
        from dataclasses import asdict
        legacy = SpinEvaluatorConfig()
        payload = asdict(legacy)
        payload.pop("exact_binary64_source_boundary", None)
        payload.pop("exact_binary64_material_boundary", None)
        payload.pop("rational_source_coordinates", None)
        payload.pop("compact_rational_source_coordinates", None)
        self.assertEqual(legacy.identity_dict(), payload)
        exact = replace(legacy, exact_binary64_source_boundary=True)
        self.assertTrue(exact.identity_dict()["exact_binary64_source_boundary"])
        for value in (1, "true", None):
            with self.subTest(value=value), self.assertRaises(ValueError):
                replace(legacy, exact_binary64_source_boundary=value)

    def test_exact_source_audit_recovers_actual_nonzero_heart_operand(self):
        # Break caught: rejecting an EXACT-valid source only because a global
        # area heuristic fails. No mesh repair, goal optimization or known-solution
        # initialization is performed; this is the fixed evaluation-34 witness.
        import importlib.util
        from pathlib import Path
        from make_it_stand.spin import (
            _native_request, _native_evaluator_settings,
            evaluate_spin_it_faster_fullrank_quadric)
        path = Path(__file__).resolve().parents[2]/"experiments/spin2014_objectives_quadric_v1/run.py"
        spec = importlib.util.spec_from_file_location("heart_source_audit_regression", path)
        loader = importlib.util.module_from_spec(spec); spec.loader.exec_module(loader)
        request, _, _ = loader.load_frozen_rotation(
            "heart-original-axis-three-objectives", "TOP_INERTIA_ONLY", 0)
        q = (-.40405369737375213, -.0000570205537868547, -1.5357003516395616e-7,
             .7721734203024663, .1477409454475553, .28733334514912345,
             -.36891865822731423, 1.665328241522832e-7,
             .00002794693751502674, -1.6799271880358278e-8)
        legacy = SpinEvaluatorConfig(derivative_method="DISCRETE_FULL_RANK_QUADRIC",
            icosphere_subdivisions=2, maximum_icosphere_subdivisions=2,
            exact_binary64_self_intersection=True)
        native, nr = _native_request(request)
        packet = native.evaluate_spin_it_faster_quadric(
            nr, q, _native_evaluator_settings(native, legacy))
        proof = native.audit_exact_triangle_boundary(
            packet["interface_vertices"], packet["interface_faces"], "OutwardPositive")
        heuristic = native.audit_triangle_boundary(
            packet["interface_vertices"], packet["interface_faces"])
        self.assertTrue(proof["boundary_valid"])
        self.assertEqual(proof["degenerate_faces"], [])
        self.assertEqual(heuristic["degenerate_face_count"], 20)
        original = evaluate_spin_it_faster_fullrank_quadric(request, q, evaluator_settings=legacy)
        self.assertIsNone(original.coefficient_jacobian)
        exact = replace(legacy, exact_binary64_source_boundary=True)
        value = evaluate_spin_it_faster_fullrank_quadric(request, q, evaluator_settings=exact)
        self.assertTrue(value.derivative_defined, value.derivative_diagnostic)
        self.assertEqual(value.geometry_hash, original.geometry_hash)
        self.assertEqual(value.total_raw_moments, original.total_raw_moments)
        self.assertNotEqual(value.scientific_hash, original.scientific_hash)

    def test_fixed_shell_measurement_avoids_early_binary64_accumulation(self):
        # Break caught: summing hundreds of surface contributions in double
        # and rounding outer/design independently before their subtraction.
        # This is a one-evaluation input-measurement test, NOT a solved case.
        import importlib.util
        from pathlib import Path
        from fractions import Fraction
        from make_it_stand.verification.mass_moments import integrate_exact_rational_raw_moments
        path=Path(__file__).resolve().parents[2]/"experiments/spin2014_objectives_quadric_v1/run.py"
        spec=importlib.util.spec_from_file_location("fixed_input_moment_precision",path)
        loader=importlib.util.module_from_spec(spec);spec.loader.exec_module(loader)
        request,objective,_=loader.load_frozen_rotation(
            "heart-original-axis-three-objectives","TOP_MASS_AND_LOWERING",0)
        request=replace(request,settings=replace(request.settings,
            maximum_restarts=1,maximum_function_evaluations_per_restart=1))
        evaluator=SpinEvaluatorConfig(derivative_method="DISCRETE_FULL_RANK_QUADRIC",
            icosphere_subdivisions=2,maximum_icosphere_subdivisions=2,
            exact_binary64_self_intersection=True)
        run=solve_spin2014_quadric_hollowing(request,objective,
            fixed_material_density=1.,evaluator_settings=evaluator)
        self.assertEqual(run.solve.function_evaluations,1)
        outer=integrate_exact_rational_raw_moments(request.outer_vertices,
            [list(x) for x in request.outer_faces])
        design=integrate_exact_rational_raw_moments(request.design_vertices,
            [list(x) for x in request.design_faces])
        exact=[float((a-b)*Fraction(1)) for a,b in zip(outer,design)]
        actual=run.measured_request.fixed_shell_moments
        self.assertLessEqual(abs(actual[3]-exact[3]),2e-16)
        self.assertLessEqual(abs(actual[4]-exact[4]),5e-17)
        self.assertLessEqual(max(abs(actual[i]-exact[i]) for i in (1,2,7,8,9)),2e-19)

    def test_legacy_surface_moment_arithmetic_is_not_silently_replaced(self):
        # Break caught: routing legacy surface integration through the new
        # extended accumulator and changing historical binary64 results.
        import importlib.util
        from pathlib import Path
        from make_it_stand.spin import integrate_raw_mass_moments
        path=Path(__file__).resolve().parents[2]/"experiments/spin2014_objectives_quadric_v1/run.py"
        spec=importlib.util.spec_from_file_location("legacy_fixed_input_precision",path)
        loader=importlib.util.module_from_spec(spec);spec.loader.exec_module(loader)
        request,_,_=loader.load_frozen_rotation(
            "heart-original-axis-three-objectives","TOP_MASS_AND_LOWERING",0)
        outer=np.asarray(integrate_raw_mass_moments(request.outer_vertices,request.outer_faces))
        design=np.asarray(integrate_raw_mass_moments(request.design_vertices,request.design_faces))
        self.assertEqual(float((outer-design)[3]),.6325135187960063)

    def test_exact_intersection_policy_recovers_actual_heart_false_positive(self):
        # Break caught: using legacy shared-corner trimming for this explicit
        # strict map, or changing density/moments while fixing its audit. This
        # is a saved-trial regression, NOT answer-dependent optimization.
        import importlib.util
        from pathlib import Path
        from make_it_stand.spin import evaluate_spin_it_faster_fullrank_quadric, _native_module
        from make_it_stand.verification.mass_moments import (
            TriangleMesh, DensityRegion, aggregate_mass_properties)
        path=Path(__file__).resolve().parents[2]/"experiments/spin2014_objectives_quadric_v1/run.py"
        spec=importlib.util.spec_from_file_location("fixed_heart_regression",path)
        loader=importlib.util.module_from_spec(spec);spec.loader.exec_module(loader)
        request,_,_=loader.load_frozen_rotation("heart-original-axis-three-objectives","TOP_INERTIA_ONLY",0)
        q=(-.4020865727991198,-.09634794695033468,.03398247116777031,
           .764083344852022,.16085720423574648,.2885126086434484,
           -.3616562663362464,-.01934908526188745,.06147589016210791,.0013424654310866292)
        legacy=SpinEvaluatorConfig(derivative_method="DISCRETE_FULL_RANK_QUADRIC",
            icosphere_subdivisions=2,maximum_icosphere_subdivisions=2)
        with self.assertRaisesRegex(RuntimeError,"intersection_free=0"):
            evaluate_spin_it_faster_fullrank_quadric(request,q,evaluator_settings=legacy)
        exact=replace(legacy,exact_binary64_self_intersection=True)
        value=evaluate_spin_it_faster_fullrank_quadric(request,q,evaluator_settings=exact)
        self.assertTrue(value.derivative_defined,value.derivative_diagnostic)
        audit=_native_module().audit_exact_triangle_boundary(
            value.material_vertices,value.material_faces,"OutwardPositive")
        self.assertTrue(audit["boundary_valid"])
        mass=aggregate_mass_properties([DensityRegion(
            TriangleMesh(value.material_vertices,value.material_faces),1.)]).mass
        self.assertAlmostEqual(mass,.20533643678161223,places=12)
        self.assertAlmostEqual(mass,value.variable_raw_moments[0],places=12)
        # The audit policy must not detach J from the actual Boolean map.
        direction=np.array([.13,-.21,.17,.09,-.08,.11,-.16,.19,-.07,.23])
        coefficients=np.asarray(q)
        direction-=coefficients*np.dot(coefficients,direction)
        direction/=np.linalg.norm(direction)
        expected=np.asarray(value.coefficient_jacobian)@direction
        errors=[]
        # Near this narrow hyperboloid neck, coarse differences have large
        # O(h^2) truncation error. Retain those diagnostics and independently
        # test three finer steps; the 1e-6 accuracy gate is not changed.
        steps=(1e-4,5e-5,2.5e-5,1e-6,5e-7,2.5e-7)
        for h in steps:
            plus=evaluate_spin_it_faster_fullrank_quadric(request,coefficients+h*direction,
                evaluator_settings=exact)
            minus=evaluate_spin_it_faster_fullrank_quadric(request,coefficients-h*direction,
                evaluator_settings=exact)
            observed=(np.asarray(plus.variable_raw_moments)-np.asarray(minus.variable_raw_moments))/(2*h)
            errors.append(np.linalg.norm(observed-expected)/max(1.,np.linalg.norm(expected)))
        print("heart exact-policy derivative steps/errors:",list(zip(steps,map(float,errors))),flush=True)
        self.assertLess(errors[1],.3*errors[0],errors)
        self.assertLess(errors[2],.3*errors[1],errors)
        self.assertLess(max(errors[-3:]),1e-6,errors)
        self.assertLess(errors[-1],errors[-3],errors)

    def test_static_real_material_uses_fresh_square_and_request_derived_initial(self):
        # Break caught: borrowing the old h objective/pass flag, multiplying
        # density twice, or requiring an answer-dependent initial plane.
        outer, faces = box(0., 1.)
        outer = tuple((x,y,4*z) for x,y,z in outer)
        design, inner_faces = box(.1, .9)
        request = SpinItFasterRequest(outer, faces, design, inner_faces,
            mode="STATIC_BALANCE", material_density=3., target_origin=(.5,.5,0.))
        config = Spin2014QuadricObjectiveConfig(objective="STATIC_LOWERING",
            inertia_weight=0., center_weight=1.)
        run = solve_spin2014_static_quadric_hollowing(request, config,
            fixed_material_density=2.)
        self.assertEqual(run.objective_profile, "SPIN2014_OBJECTIVES_QUADRIC_V1")
        self.assertEqual(run.initialization, "INPUT_FULL_MATERIAL_HEIGHT_V1")
        self.assertEqual(run.measured_request.outer_vertices, request.outer_vertices)
        self.assertAlmostEqual(run.measured_request.fixed_shell_moments[0], 6.976, places=12)
        self.assertAlmostEqual(run.candidate["total_raw_moments"][0], 8.512, places=12)
        h = 16.256/8.512
        self.assertAlmostEqual(run.squared_height["objective"], h*h, places=12)
        self.assertAlmostEqual(run.squared_height["objective_gradient"][0], -2*h*h/8.512, places=12)
        self.assertTrue(run.squared_height["boundary_witness_verified"])
        self.assertTrue(run.squared_height["squared_height_conditions_verified"])
        self.assertEqual(run.candidate["accepted_steps"], 0)
        # A local numerical gate is NOT a material seal or nontrivial function.
        self.assertFalse(run.material_verified)
        scaled = solve_spin2014_static_quadric_hollowing(request, replace(config,
            solver_mass_unit_in_declared_units=2., solver_length_unit_in_declared_units=3.),
            fixed_material_density=2.)
        self.assertAlmostEqual(scaled.squared_height["objective"], 9*h*h, places=11)
        with self.assertRaises(ValueError):
            solve_spin2014_static_quadric_hollowing(replace(request,
                target_origin=(.5,.5,.01)), config, fixed_material_density=2.)
        with self.assertRaises(ValueError):
            solve_spin2014_static_quadric_hollowing(request, replace(config,
                objective="TOP_INERTIA_ONLY"), fixed_material_density=2.)

    def test_newton_payload_survives_python_event_conversion(self):
        # Break caught: a native matrix payload being rejected/dropped by the
        # immutable Python ledger. Literal serialization fixture, not a solve.
        entry = dict(type="newton_direction", restart=0, iteration=1,
            evaluation_index=2, restart_evaluation_index=2, step=None,
            outcome="rank_failure", reason="rank_deficient", detail="",
            quadric_coefficients=None, quadric_hash="fixture", objective=None,
            normalized_equation_residual=None, physical_relative_residual=None,
            multiplier_mu=None, material_geometry_hash=None, material_topology_hash=None,
            equation_profile="DUAL_MONOMIAL_RAW_V1",
            solver_coordinate_identity="literal_monomial_unit_sphere_v1",
            raw_equation_residual=None,
            newton_system_payload=dict(schema="spin-newton-system-payload-v2",
                matrix=[[2.,0.],[0.,0.]], right_hand_side=[1.,0.], computed_rank=1,
                equilibrated_rank=1, multiplier_mu=3., constraint_multipliers=[4.],
                row_scales=[.5,1.], column_scales=[1.,1.],
                linear_residual_infinity_norm=0., linear_backward_error=0.))
        event = _immutable_spin_solver_event(entry)
        self.assertEqual(event.newton_system_payload.matrix, ((2.,0.),(0.,0.)))
        self.assertEqual(event.newton_system_payload.right_hand_side, (1.,0.))
        self.assertEqual(event.newton_system_payload.computed_rank, 1)
        self.assertEqual(event.newton_system_payload.constraint_multipliers, (4.,))

    def request(self):
        outer, faces = box(-2., 2.)
        design, inner_faces = box(-1., 1.)
        return SpinItFasterRequest(
            outer, faces, design, inner_faces, mode="SPINNING_TOP",
            material_density=5., evidence="RECONSTRUCTED_ASSET",
            settings=SpinSolverConfig(maximum_restarts=1,
                maximum_function_evaluations_per_restart=1,
                equation_profile="DUAL_MONOMIAL_RAW_V1",
                capture_evaluation_payload=True))

    def test_config_rejects_missing_units_and_invalid_scales(self):
        for kwargs in ({"declared_length_unit":""},
                       {"solver_mass_unit_in_declared_units":0.},
                       {"solver_length_unit_in_declared_units":float("nan")},
                       {"objective":"not-an-objective"}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                Spin2014QuadricObjectiveConfig(**kwargs)

    def test_real_boolean_call_keeps_objective_and_complete_failure_ledger(self):
        # Break caught: trusting manual fixed moments, relabeling the old
        # objective, losing payloads, or treating a budget exit as material pass.
        request = self.request()
        q = (-.31,.07,-.09,.11,.83,1.17,-.72,.06,-.04,.03)
        config = Spin2014QuadricObjectiveConfig(objective="TOP_MASS_AND_LOWERING",
            declared_length_unit="test_length", declared_mass_unit="test_mass")
        evaluator = SpinEvaluatorConfig(derivative_method="DISCRETE_FULL_RANK_QUADRIC",
            icosphere_subdivisions=2, truncation_radius=8.)
        run = solve_spin2014_quadric_hollowing(request, config,
            fixed_material_density=3., evaluator_settings=evaluator,
            initial_quadrics=(q,))
        self.assertEqual(run.objective_profile, "SPIN2014_OBJECTIVES_QUADRIC_V1")
        self.assertEqual(run.objective_settings, config)
        self.assertAlmostEqual(run.measured_request.fixed_shell_moments[0], 168., places=11)
        self.assertEqual(request.fixed_shell_moments, (0.,)*10)
        self.assertEqual(run.measured_request.outer_vertices, request.outer_vertices)
        self.assertEqual(run.measured_request.outer_faces, request.outer_faces)
        self.assertEqual(run.solve.function_evaluations, 1)
        self.assertEqual(sum(run.solve.restart_function_evaluations), 1)
        self.assertFalse(run.solve.final_geometry_verified)
        self.assertFalse(run.material_export_complete)
        self.assertIsNone(run.cavity_vertices)
        events = [e for e in run.solve.events if e.evaluation_payload is not None]
        self.assertTrue(events, "actual evaluated moments are missing")
        payload = events[0].evaluation_payload
        s = np.asarray(payload.solver_raw_moments, dtype=float)
        mass, lx, ly, lz, qxx, qyy, qzz, qyz, qxz, qxy = s
        ix = qyy+qzz-(ly*ly+lz*lz)/mass
        iy = qxx+qzz-(lx*lx+lz*lz)/mass
        iz = qxx+qyy-(lx*lx+ly*ly)/mass
        ixy = -qxy+lx*ly/mass
        expected = (ix*ix+iy*iy+2*ixy*ixy)/(iz*iz)+lz*lz
        self.assertAlmostEqual(events[0].objective, expected, places=8)
        with self.assertRaises(ValueError):
            solve_spin2014_quadric_hollowing(replace(request,
                fixed_shell_moments=(336.,)+(0.,)*9), config,
                fixed_material_density=3., evaluator_settings=evaluator,
                initial_quadrics=(q,))

    def test_saved_material_rejects_relaxed_protocol_before_blessing_a_terminal_label(self):
        # Break caught: a saved CONVERGED label plus relaxed tolerance bypasses
        # the strict profile gate that the live native solver enforces. This
        # bounded negative run is NOT a solved material/functional success.
        from make_it_stand.verification.faster_material import prepare_spin2014_quadric_material_output
        request=self.request()
        evaluator=SpinEvaluatorConfig(derivative_method="DISCRETE_FULL_RANK_QUADRIC",
            icosphere_subdivisions=2,truncation_radius=8.)
        run=solve_spin2014_quadric_hollowing(request,Spin2014QuadricObjectiveConfig(),
            fixed_material_density=3.,evaluator_settings=evaluator,
            initial_quadrics=((- .31,.07,-.09,.11,.83,1.17,-.72,.06,-.04,.03),))
        self.assertNotEqual(run.solve.status,"CONVERGED")
        relaxed=replace(run.measured_request,settings=replace(run.measured_request.settings,
            equation_tolerance=1e-3))
        counterfeit=replace(run,measured_request=relaxed,
            solve=replace(run.solve,status="CONVERGED",final_geometry_verified=True))
        with self.assertRaisesRegex(ValueError,"strict|protocol|tolerance"):
            prepare_spin2014_quadric_material_output(counterfeit,evaluator_settings=evaluator)

    def test_saved_material_rejects_missing_ledger_before_geometry_replay(self):
        # Break caught: a terminal label is admitted even after failed attempts
        # or evaluations are deleted. This one-evaluation negative is not a
        # converged material witness and must fail at the evidence gate first.
        from make_it_stand.verification.faster_material import prepare_spin2014_quadric_material_output
        evaluator=SpinEvaluatorConfig(derivative_method="DISCRETE_FULL_RANK_QUADRIC",
            icosphere_subdivisions=2,truncation_radius=8.)
        run=solve_spin2014_quadric_hollowing(self.request(),Spin2014QuadricObjectiveConfig(),
            fixed_material_density=3.,evaluator_settings=evaluator,
            initial_quadrics=((- .31,.07,-.09,.11,.83,1.17,-.72,.06,-.04,.03),))
        counterfeit=replace(run,solve=replace(run.solve,status="CONVERGED",
            final_geometry_verified=True,events=()))
        with self.assertRaisesRegex(ValueError,"ledger"):
            prepare_spin2014_quadric_material_output(counterfeit,evaluator_settings=evaluator)

    def test_rotational_ledger_keeps_failed_attempt_and_rejected_evaluation(self):
        # Hand-counted serialized fixture: one failed attempt, then two
        # accepted evaluations plus a rejected trial, within the literal
        # 60-evaluation budget. It is not a method optimization success.
        from types import SimpleNamespace as NS
        from copy import deepcopy
        from make_it_stand.verification.faster_material import _check_spin2014_rotational_ledger
        settings=SpinSolverConfig(capture_evaluation_payload=False,
            equation_profile="DUAL_MONOMIAL_RAW_V1")
        q=(1.,)+(0.,)*9
        # Independently hashed literal C++ hexfloat serialization of (1,0,...,0).
        qhash="bad43883bc41c5b8fa364860e409cb62142634c183a18e59a802b61e16cd650f"
        def event(kind,restart,global_index,local_index,outcome,reason):
            return NS(type=kind,restart=restart,evaluation_index=global_index,
                restart_evaluation_index=local_index,outcome=outcome,reason=reason,
                quadric_coefficients=q,quadric_hash=qhash,evaluation_payload=None,
                equation_profile="DUAL_MONOMIAL_RAW_V1")
        events=(event("attempt_start",0,0,0,"started","restart_started"),
            event("evaluation",0,1,1,"failed","evaluation_exception"),
            event("attempt_terminal",0,1,1,"terminated","initial_evaluation_failed"),
            event("attempt_start",1,1,0,"started","restart_started"),
            event("evaluation",1,2,1,"accepted","initial_evaluation_succeeded"),
            event("evaluation",1,3,2,"rejected","no_residual_decrease"),
            event("evaluation",1,4,3,"accepted","residual_decreased"),
            event("attempt_terminal",1,4,3,"terminated","converged"))
        timings=tuple(NS(restart=r,evaluation_index=g,restart_evaluation_index=l)
            for r,g,l in ((0,1,1),(1,2,1),(1,3,2),(1,4,3)))
        solve=NS(status="CONVERGED",restart_count=2,restart_function_evaluations=(1,3),
            function_evaluations=4,events=events,evaluation_timings=timings)
        _check_spin2014_rotational_ledger(solve,settings)
        mutations=[]
        row=deepcopy(solve);row.events=row.events[3:];mutations.append(row)
        row=deepcopy(solve);row.events=tuple(e for e in row.events if e.outcome!="rejected");mutations.append(row)
        row=deepcopy(solve);row.events=row.events[:-1];mutations.append(row)
        row=deepcopy(solve);row.evaluation_timings=row.evaluation_timings[:-1];mutations.append(row)
        row=deepcopy(solve);row.restart_function_evaluations=(0,4);mutations.append(row)
        row=deepcopy(solve);row.events[4].restart_evaluation_index=2;mutations.append(row)
        row=deepcopy(solve);row.events[-1].reason="line_search_exhausted";mutations.append(row)
        row=deepcopy(solve);row.restart_function_evaluations=(1,61);mutations.append(row)
        for index,row in enumerate(mutations):
            with self.subTest(mutation=index),self.assertRaisesRegex(ValueError,"ledger"):
                _check_spin2014_rotational_ledger(row,settings)
        with self.assertRaisesRegex(ValueError,"ledger"):
            _check_spin2014_rotational_ledger(solve,replace(settings,capture_evaluation_payload=True))

    def test_rotational_ledger_rejects_coefficients_mismatched_to_identity(self):
        # Break caught: changing a failed/accepted trial's coefficients while
        # retaining a nonempty hash previously passed the saved-evidence gate.
        # This literal serialization fixture is not a successful optimization.
        from types import SimpleNamespace as NS
        from copy import deepcopy
        from make_it_stand.verification.faster_material import _check_spin2014_rotational_ledger
        settings=SpinSolverConfig(capture_evaluation_payload=False,
            equation_profile="DUAL_MONOMIAL_RAW_V1")
        q=(1.,)+(0.,)*9
        qhash="bad43883bc41c5b8fa364860e409cb62142634c183a18e59a802b61e16cd650f"
        def event(kind,index,outcome,reason):
            return NS(type=kind,restart=0,evaluation_index=index,
                restart_evaluation_index=index,outcome=outcome,reason=reason,
                quadric_coefficients=q,quadric_hash=qhash,evaluation_payload=None,
                equation_profile="DUAL_MONOMIAL_RAW_V1")
        events=(event("attempt_start",0,"started","restart_started"),
            event("evaluation",1,"accepted","initial_evaluation_succeeded"),
            event("evaluation",2,"failed","evaluation_exception"),
            event("evaluation",3,"accepted","residual_decreased"),
            event("attempt_terminal",3,"terminated","converged"))
        solve=NS(status="CONVERGED",restart_count=1,restart_function_evaluations=(3,),
            function_evaluations=3,events=events,
            evaluation_timings=tuple(NS(restart=0,evaluation_index=i,
                restart_evaluation_index=i) for i in (1,2,3)))
        _check_spin2014_rotational_ledger(solve,settings)
        for index in (0,1,2,3):
            for field,value in (("quadric_coefficients",(-1.,)+(0.,)*9),
                                ("quadric_hash","0"*64),
                                ("quadric_coefficients",(2.,)+(0.,)*9)):
                altered=deepcopy(solve)
                setattr(altered.events[index],field,value)
                with self.subTest(event=index,field=field,value=value), \
                        self.assertRaisesRegex(ValueError,"ledger.*identity"):
                    _check_spin2014_rotational_ledger(altered,settings)


if __name__ == "__main__":
    unittest.main()
