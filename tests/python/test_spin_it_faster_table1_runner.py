from __future__ import annotations

import importlib.util
import copy
import json
import math
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "experiments/spin_it_reproduction/faster/run_paper_table1.py"


def load_runner():
    spec = importlib.util.spec_from_file_location("spin_it_faster_table1", RUNNER)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load Faster Table 1 runner")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def planar_method_case(**controls):
    outer = [[-3., -3.], [3., -3.], [3., 3.], [-3., 3.]]
    patch = [[math.pi/4-.1, -.1], [math.pi/4+.1, -.1],
             [math.pi/4+.1, .1], [math.pi/4-.1, .1]]
    problem = {
        "design_region_loops": [outer, list(reversed(patch))],
        "design_areal_density": 1.,
        "extrusion_interval": [-.7, .8],
        "extrusion_volume_density": 2./3.,
        "fixed_extrusion_layers": [
            {"loops": [patch], "z_interval": [-.1, 0.], "volume_density": 500.},
            {"loops": [patch], "z_interval": [0., .1], "volume_density": 500.}],
        **controls,
    }
    return {"id": "planar-layered-method", "mode": "FIDGET_SPINNER",
            "asset_class": "SYNTHETIC_REGRESSION", "planar_problem": problem}


class SpinItFasterPlanarTableDispatchTests(unittest.TestCase):
    def test_actual_layered_dogleg_uses_preregistered_centers_and_keeps_entities(self):
        runner = load_runner()
        case = planar_method_case()
        plan = runner.plan_payload([case], starts=2, max_restarts=2)
        actual = runner.run_case(case, starts=2, max_restarts=2)
        self.assertEqual(actual["solver_family"], "PLANAR_DOGLEG_V1")
        self.assertEqual(actual["seed_sequence"], [0, 1])
        self.assertEqual(actual["initialization_plan"], plan["cases"][0]["runs"])
        self.assertEqual(plan["initialization_protocol"], "INPUT_BOUNDING_BOX_PLANAR_CENTERS_V1")
        self.assertEqual(plan["asset_class"], "SYNTHETIC_REGRESSION")
        self.assertEqual(actual["formal_paper_admission"], "NOT_GRANTED_BY_METHOD_RUNNER")
        for row in actual["attempts"]:
            self.assertEqual(row["solver_family"], "PLANAR_DOGLEG_V1")
            self.assertNotIn("raw_equation_residual", row)
            self.assertNotIn("initial_quadric", row)
            self.assertEqual(row["initial_center"],
                actual["initialization_plan"][row["run"]]["initial_centers"][row["restart"]-1])
            if row.get("solver_result") is not None:
                self.assertEqual(len(row["solver_result"]["extruded_material_entities"]), 3)
                self.assertEqual(row["solver_result"]["extruded_material_scope"],
                                 "design_plus_fixed_layers")
        runner.validate_attempt_ledger(actual["attempts"], expected_runs=2, maximum_restarts=2)
        json.dumps(actual, allow_nan=False)

    def test_failed_planar_budget_preserves_all_runs_restarts_and_auxiliary_cost(self):
        runner = load_runner()
        actual = runner.run_case(planar_method_case(maximum_function_evaluations=1),
                                 starts=2, max_restarts=2)
        self.assertEqual(len(actual["attempts"]), 4)
        self.assertEqual(actual["summary"]["failed_runs"], 2)
        self.assertEqual(actual["summary"]["resource_limited_runs"], 2)
        self.assertEqual(actual["summary"]["resource_limited_attempts"], 4)
        self.assertEqual(actual["summary"]["runs_with_resource_limited_attempts"], 2)
        self.assertEqual(actual["summary"]["total_function_evaluations_all_attempts"], 4)
        self.assertEqual(actual["total_auxiliary_residual_evaluations"], 8)
        self.assertEqual(actual["total_residual_evaluations"], 12)
        self.assertEqual(actual["summary"]["strict_terminal_numerics"]["numerically_qualified_runs"], 0)
        self.assertFalse(actual["passed"])
        self.assertTrue(all(row["status"] == "MAXIMUM_FUNCTION_EVALUATIONS"
                            for row in actual["attempts"]))

    def test_missing_planar_geometry_cannot_fall_back_to_legacy_3d(self):
        runner = load_runner()
        case = next(row for row in json.loads(runner.DEFAULT_CASES.read_text())["cases"]
                    if row["mode"] == "FIDGET_SPINNER")
        with self.assertRaisesRegex(ValueError, "planar_problem"):
            runner.run_case(case, starts=1, max_restarts=1)

    def test_supplied_answer_initialization_cannot_enter_seeded_planar_plan(self):
        runner = load_runner()
        with self.assertRaisesRegex(ValueError, "initial_center"):
            runner.plan_payload([planar_method_case(initial_center=[1., 0.])], starts=1)

    def test_planar_validator_rejects_dropped_calls_and_forged_convergence(self):
        runner = load_runner()
        actual = runner.run_case(planar_method_case(maximum_function_evaluations=1),
                                 starts=1, max_restarts=1)
        damaged = copy.deepcopy(actual["attempts"])
        damaged[0]["evaluation_ledger"]["solver_evaluations"] = []
        with self.assertRaises(ValueError):
            runner.validate_attempt_ledger(damaged, expected_runs=1, maximum_restarts=1)
        forged = copy.deepcopy(actual["attempts"][0])
        forged["status"] = "CONVERGED"
        forged["normalized_equation_residual"] = 0.
        forged["numerical_acceptance"] = {"passed": True}
        self.assertFalse(runner.numerical_attempt_acceptance(forged)["passed"])

    def test_planar_acceptance_rejects_invalid_thresholds_and_changed_terminal_value(self):
        runner = load_runner()
        rows = runner.run_case(planar_method_case(), starts=2, max_restarts=2)["attempts"]
        successful = next(row for row in rows if row["status"] == "CONVERGED")
        self.assertTrue(runner.numerical_attempt_acceptance(successful)["passed"])
        for value in (None, True, float("nan"), float("inf"), 1e-3):
            altered = copy.deepcopy(successful)
            altered["residual_tolerance"] = value
            with self.subTest(value=value):
                self.assertFalse(runner.numerical_attempt_acceptance(altered)["passed"])
        altered = copy.deepcopy(successful)
        altered["normalized_equation_residual"] = 0.
        altered["solver_result"]["normalized_residual"] = 0.
        self.assertFalse(runner.numerical_attempt_acceptance(altered)["passed"])


    def test_planar_acceptance_binds_terminal_center_and_successful_call_semantics(self):
        runner = load_runner()
        successful = next(row for row in runner.run_case(planar_method_case(), starts=2,
            max_restarts=2)["attempts"] if row["status"] == "CONVERGED")
        for mutation in ("center", "initial_failed", "missing_final_residual"):
            altered = copy.deepcopy(successful)
            if mutation == "center":
                altered["solver_result"]["center"] = [1e6, 1e6]
            elif mutation == "initial_failed":
                altered["evaluation_ledger"]["solver_evaluations"][0].update(
                    outcome="failed", has_residual=False, normalized_residual=None,
                    normalized_residual_vector=None)
            else:
                altered["evaluation_ledger"]["auxiliary_evaluations"][-1]["has_residual"] = False
            with self.subTest(mutation=mutation):
                self.assertFalse(runner.numerical_attempt_acceptance(altered)["passed"])


class SpinItFasterTable1RunnerTests(unittest.TestCase):
    def test_static_iteration_limit_is_not_geometry_failure_or_other_family_status(self):
        runner = load_runner()
        rows = [dict(run=i, restart=1, solver_family='STATIC_PLANE_FAMILY_V1',
                     status=status, function_evaluations=1)
                for i, status in enumerate(('ITERATION_LIMIT', 'GEOMETRY_EVALUATION_FAILED',
                    'LINE_SEARCH_EXHAUSTED', 'SINGULAR_REDUCED_SYSTEM'))]
        rows.append(dict(run=4, restart=1, solver_family='CHERCHI_RAW_EQ14_V1',
                         status='ITERATION_LIMIT', function_evaluations=1))
        summary = runner.summarize_attempts(rows, 5)
        self.assertEqual(summary['resource_limited_runs'], 1)
        self.assertEqual(summary['resource_limited_attempts'], 1)
        self.assertEqual(summary['runs_with_resource_limited_attempts'], 1)
        self.assertEqual(summary['failed_runs'], 5)

    def test_resource_statistics_keep_earlier_exhausted_attempts_after_restart(self):
        runner = load_runner()
        rows = [
            {'run':0, 'restart':1, 'status':'RESOURCE_LIMIT', 'function_evaluations':60},
            {'run':0, 'restart':2, 'status':'CONVERGED', 'function_evaluations':3},
            {'run':1, 'restart':1, 'status':'RESOURCE_LIMIT', 'function_evaluations':60},
            {'run':1, 'restart':2, 'status':'OPTIMIZATION_DID_NOT_CONVERGE', 'function_evaluations':1},
        ]
        summary = runner.summarize_attempts(rows, 2)
        # Terminal truncation and any-attempt truncation answer different
        # questions. An eventual native success must not erase earlier cost.
        self.assertEqual(summary['resource_limited_runs'], 0)
        self.assertEqual(summary.get('resource_limited_attempts'), 2)
        self.assertEqual(summary.get('runs_with_resource_limited_attempts'), 2)
        self.assertEqual(summary['total_function_evaluations_all_attempts'], 124)
        self.assertEqual(summary['successful_runs'], 1)
        self.assertEqual(summary['strict_terminal_numerics']['numerically_qualified_runs'], 0)

    def test_planar_budget_statuses_are_not_confused_with_geometric_failure(self):
        runner = load_runner()
        rows = [dict(run=i, restart=1, solver_family='PLANAR_DOGLEG_V1',
                     status=status, function_evaluations=1)
                for i, status in enumerate(('MAXIMUM_FUNCTION_EVALUATIONS',
                    'MAXIMUM_ITERATIONS', 'EVALUATION_FAILED', 'TRUST_REGION_COLLAPSED'))]
        summary = runner.summarize_attempts(rows, 4)
        self.assertEqual(summary['resource_limited_runs'], 2)
        self.assertEqual(summary['resource_limited_attempts'], 2)
        self.assertEqual(summary['runs_with_resource_limited_attempts'], 2)
        self.assertEqual(summary['failed_runs'], 4)

    def test_summary_separates_native_convergence_from_qualified_terminal_numerics(self):
        runner = load_runner()
        def attempt(run, restart, raw, mu):
            return {'run':run,'restart':restart,'status':'CONVERGED',
                'function_evaluations':3,'newton_iterations':2,
                'residual_seconds':0.,'derivative_seconds':0.,
                'equation_profile':'DUAL_MONOMIAL_RAW_V1',
                'raw_equation_residual':raw,'physical_relative_residual':1e-9,
                'multiplier_mu':mu,'numerical_acceptance':{'passed':True},
                'trace':[{'equation_profile':'DUAL_MONOMIAL_RAW_V1',
                          'raw_equation_residual':raw,'multiplier_mu':mu}]}
        rows=[attempt(0,1,1e-13,-1.),attempt(1,1,1e-13,2.)]
        summary=runner.summarize_attempts(rows,2)
        self.assertEqual(summary['successful_runs'],2)  # unchanged native status statistic
        qualified=summary.get('strict_terminal_numerics')
        self.assertIsInstance(qualified,dict,'qualified statistics are missing')
        self.assertEqual(qualified['numerically_qualified_runs'],1)
        self.assertEqual(qualified['unqualified_runs'],1)
        self.assertEqual(qualified['qualified_run_ids'],[1])
        self.assertFalse(qualified['all_runs_numerically_qualified'])
        self.assertEqual(qualified['native_converged_but_unqualified_attempts'],1)
        self.assertEqual(qualified['total_function_evaluations_all_attempts'],6)
        self.assertEqual(qualified['scope'],
            'TERMINAL_NUMERICS_STATISTICS_NOT_MATERIAL_OR_PAPER_ACCEPTANCE')
        self.assertIn('POSITIVE_MULTIPLIER_REQUIRED',qualified['rejected_attempts'][0]['reasons'])
        # Pre-filled passed=true must not bypass actual terminal checks.
        rows[0]['multiplier_mu']=2.;rows[0]['trace'][-1]['multiplier_mu']=2.
        rows[0]['raw_equation_residual']=1e-8;rows[0]['trace'][-1]['raw_equation_residual']=1e-8
        qualified=runner.summarize_attempts(rows,2)['strict_terminal_numerics']
        self.assertEqual(qualified['numerically_qualified_runs'],1)
        self.assertIn('RAW_EQ14_TOLERANCE',qualified['rejected_attempts'][0]['reasons'])

    def test_qualified_statistics_keep_missing_runs_and_reject_empty_population(self):
        runner=load_runner()
        summary=runner.summarize_attempts([],2)
        qualified=summary.get('strict_terminal_numerics')
        self.assertIsInstance(qualified,dict,'qualified statistics are missing')
        self.assertEqual(qualified['unqualified_runs'],2)
        self.assertEqual(qualified['numerically_qualified_runs'],0)
        self.assertFalse(qualified['all_runs_numerically_qualified'])
        with self.assertRaises(ValueError):runner.summarize_attempts([],0)

    def test_numerical_acceptance_requires_raw_terminal_residual_and_positive_multiplier(self):
        runner = load_runner()
        checker = getattr(runner, 'numerical_attempt_acceptance', None)
        self.assertTrue(callable(checker), 'raw terminal acceptance check is missing')
        attempt = {'status':'CONVERGED', 'equation_profile':'DUAL_MONOMIAL_RAW_V1',
            'raw_equation_residual':1e-13, 'normalized_equation_residual':1e-16,
            'physical_relative_residual':1e-8, 'multiplier_mu':2.,
            'trace':[{'equation_profile':'DUAL_MONOMIAL_RAW_V1',
                      'raw_equation_residual':1e-13, 'multiplier_mu':2.}]}
        approved = checker(attempt)
        self.assertTrue(approved['passed'])
        self.assertEqual(approved['scope'],
            'TERMINAL_NUMERICS_ONLY_NOT_MATERIAL_OR_PAPER_ACCEPTANCE')
        for change, reason in (
            ({'raw_equation_residual':1e-8}, 'RAW_EQ14_TOLERANCE'),
            ({'raw_equation_residual':float('nan')}, 'RAW_EQ14_TOLERANCE'),
            ({'multiplier_mu':0.}, 'POSITIVE_MULTIPLIER_REQUIRED'),
            ({'multiplier_mu':-1.}, 'POSITIVE_MULTIPLIER_REQUIRED'),
            ({'multiplier_mu':True}, 'POSITIVE_MULTIPLIER_REQUIRED'),
            ({'physical_relative_residual':float('inf')}, 'PHYSICAL_CONSTRAINT_TOLERANCE'),
            ({'equation_profile':'LEGACY_NATIVE_NORMALIZED_V1'}, 'RAW_EQ14_PROFILE_REQUIRED'),
            ({'trace':[]}, 'TERMINAL_TRACE_REQUIRED'),
            ({'status':'RESOURCE_LIMIT'}, 'NATIVE_CONVERGENCE_REQUIRED'),
        ):
            with self.subTest(change=change):
                result = checker({**attempt, **change})
                self.assertFalse(result['passed'])
                self.assertIn(reason, result['reasons'])

    def test_numerical_acceptance_rejects_mismatched_terminal_snapshot(self):
        runner = load_runner()
        checker = getattr(runner, 'numerical_attempt_acceptance', None)
        self.assertTrue(callable(checker), 'raw terminal acceptance check is missing')
        attempt = {'status':'CONVERGED', 'equation_profile':'DUAL_MONOMIAL_RAW_V1',
            'raw_equation_residual':1e-13, 'physical_relative_residual':1e-8,
            'multiplier_mu':2., 'trace':[{'equation_profile':'DUAL_MONOMIAL_RAW_V1',
                'raw_equation_residual':1e-8, 'multiplier_mu':-1.}]}
        result = checker(attempt)
        self.assertFalse(result['passed'])
        self.assertIn('TERMINAL_NUMERICS_MISMATCH', result['reasons'])

    def test_registered_seed_sequence_rejects_duplicate_and_coerced_runs(self):
        runner = load_runner()
        self.assertEqual(runner.registered_seeds(2),(0,1))
        for bad in [(0,0),(True,1),(0,-1),(0,),(.5,1)]:
            with self.subTest(seeds=bad):
                with self.assertRaises(ValueError):runner.registered_seeds(2,bad)
    def test_initial_points_are_deterministic_and_mode_constrained(self) -> None:
        runner = load_runner()
        first = runner.generate_initial_quadrics("YO_YO", seed=73, count=50)
        second = runner.generate_initial_quadrics("YO_YO", seed=73, count=50)
        self.assertEqual(first, second)
        self.assertEqual(len(first), 50)
        self.assertEqual(len(set(tuple(row) for row in first)), 50)
        for row in first:
            self.assertAlmostEqual(math.sqrt(sum(value * value for value in row)), 1.0)
            self.assertEqual(row[0], 0.0)
            self.assertLess(row[6], 0.0)
        for row in runner.generate_initial_quadrics("BUOYANT_STABILITY", 11, 4):
            self.assertEqual(row[4:], [0.0] * 6)
        for row in runner.generate_initial_quadrics("FIDGET_SPINNER", 17, 4):
            self.assertLess(row[0], 0.0)
            self.assertGreater(row[4], 0.0)
            self.assertEqual(row[4], row[5])

    def test_domain_aware_planes_intersect_the_design_bounds(self) -> None:
        runner = load_runner()
        corners = tuple(
            (x, y, z)
            for x in (-2.0, 3.0)
            for y in (-1.0, 4.0)
            for z in (-0.5, 0.75)
        )
        for mode in ("STATIC_BALANCE", "BUOYANT_STABILITY"):
            for row in runner.generate_initial_quadrics(mode, 91, 12, corners):
                values = [
                    row[0] + row[1] * x + row[2] * y + row[3] * z
                    for x, y, z in corners
                ]
                self.assertLess(min(values), 0.0)
                self.assertGreater(max(values), 0.0)

    def test_summary_uses_paper_statistics_and_hash_excludes_timings(self) -> None:
        runner = load_runner()
        attempts = [
            {"run": 0, "restart": 1, "status": "FAILED", "newton_iterations": 3,
             "function_evaluations": 8, "residual_seconds": 0.2, "derivative_seconds": 0.1},
            {"run": 0, "restart": 2, "status": "CONVERGED", "newton_iterations": 5,
             "function_evaluations": 9, "residual_seconds": 0.3, "derivative_seconds": 0.2},
            {"run": 1, "restart": 1, "status": "CONVERGED", "newton_iterations": 7,
             "function_evaluations": 11, "residual_seconds": 0.4, "derivative_seconds": 0.25},
        ]
        summary = runner.summarize_attempts(attempts, expected_runs=2)
        self.assertTrue(summary["all_runs_converged"])
        self.assertEqual(summary["average_restarts"], 1.5)
        self.assertEqual(summary["median_newton_iterations"], 6.0)
        self.assertEqual(summary["median_function_evaluations"], 10.0)
        first_hash = runner.scientific_hash({"status": "CONVERGED", "elapsed_seconds": 1.0})
        second_hash = runner.scientific_hash({"status": "CONVERGED", "elapsed_seconds": 9.0})
        self.assertEqual(first_hash, second_hash)

    def test_summary_keeps_failed_top_level_runs_in_its_denominator(self) -> None:
        runner = load_runner()
        attempts = [
            {"run": 0, "restart": 1, "status": "FAILED", "newton_iterations": 2,
             "function_evaluations": 5, "residual_seconds": 0.1, "derivative_seconds": 0.1},
            {"run": 0, "restart": 2, "status": "CONVERGED", "newton_iterations": 4,
             "function_evaluations": 7, "residual_seconds": 0.1, "derivative_seconds": 0.1},
            {"run": 1, "restart": 1, "status": "FAILED", "newton_iterations": 3,
             "function_evaluations": 8, "residual_seconds": 0.1, "derivative_seconds": 0.1},
            {"run": 1, "restart": 2, "status": "RESOURCE_LIMIT", "newton_iterations": 0,
             "function_evaluations": 0, "residual_seconds": 0.0, "derivative_seconds": 0.0},
        ]
        summary = runner.summarize_attempts(attempts, expected_runs=2)
        self.assertEqual(summary["successful_runs"], 1)
        self.assertEqual(summary["failed_runs"], 1)
        self.assertEqual(summary["resource_limited_runs"], 1)
        self.assertEqual(summary["average_attempts_all_runs"], 2.0)
        self.assertEqual(summary["average_restarts"], 2.0)

    def test_summary_counts_failed_attempt_evaluations_not_only_successful_terminal(self):
        runner = load_runner()
        attempts = [
            {"run": 0, "restart": 1, "status": "FAILED", "newton_iterations": 3,
             "function_evaluations": 8, "residual_seconds": .2, "derivative_seconds": .1},
            {"run": 0, "restart": 2, "status": "CONVERGED", "newton_iterations": 5,
             "function_evaluations": 9, "residual_seconds": .3, "derivative_seconds": .2},
            {"run": 1, "restart": 1, "status": "FAILED", "newton_iterations": 7,
             "function_evaluations": 11, "residual_seconds": .4, "derivative_seconds": .25},
            {"run": 1, "restart": 2, "status": "RESOURCE_LIMIT", "newton_iterations": 0,
             "function_evaluations": 0, "residual_seconds": 0., "derivative_seconds": 0.},
        ]
        summary = runner.summarize_attempts(attempts, expected_runs=2)
        self.assertEqual(summary.get("total_function_evaluations_all_attempts"), 28)
        self.assertEqual(summary.get("function_evaluations_by_run"), [17, 11])
        self.assertEqual(summary.get("average_function_evaluations_all_runs"), 14.)
        # Keep the legacy conditional median, but explicitly describe its
        # population: it must not be reported as complete top-level cost.
        self.assertEqual(summary["median_function_evaluations"], 9.)
        self.assertEqual(summary.get("conditional_statistics_population"),
                         "SUCCESSFUL_TERMINAL_ATTEMPTS_ONLY")

    def test_attempt_ledger_requires_complete_restart_events_and_budget_agreement(self) -> None:
        runner = load_runner()
        attempt = {
            "run": 0,
            "restart": 1,
            "status": "FAILED",
            "function_evaluations": 2,
            "solver_events": [
                {"type": "attempt_start", "restart": 0},
                {"type": "evaluation", "restart": 0,
                 "evaluation_index": 1, "restart_evaluation_index": 1},
                {"type": "evaluation", "restart": 0,
                 "evaluation_index": 2, "restart_evaluation_index": 2},
                {"type": "attempt_terminal", "restart": 0},
            ],
        }
        runner.validate_attempt_ledger([attempt], expected_runs=1, maximum_restarts=1)

        missing_terminal = {**attempt, "solver_events": attempt["solver_events"][:-1]}
        with self.assertRaisesRegex(ValueError, "one start and one terminal"):
            runner.validate_attempt_ledger(
                [missing_terminal], expected_runs=1, maximum_restarts=1
            )

        mismatched_evaluations = {**attempt, "function_evaluations": 1}
        with self.assertRaisesRegex(ValueError, "evaluation count"):
            runner.validate_attempt_ledger(
                [mismatched_evaluations], expected_runs=1, maximum_restarts=1
            )

        over_budget = {**attempt, "function_evaluations": 61}
        with self.assertRaisesRegex(ValueError, "per-attempt function budget"):
            runner.validate_attempt_ledger(
                [over_budget], expected_runs=1, maximum_restarts=1
            )

        with self.assertRaisesRegex(ValueError, "missing failed or unattempted"):
            runner.validate_attempt_ledger(
                [attempt], expected_runs=1, maximum_restarts=2
            )

        after_success = [
            {**attempt, "status": "CONVERGED",
             "solver_events":attempt["solver_events"][:-1]+[
                 {**attempt["solver_events"][-1],"outcome":"terminated","reason":"converged"}]},
            {**attempt, "restart": 2},
        ]
        with self.assertRaisesRegex(ValueError, "after its first convergence"):
            runner.validate_attempt_ledger(
                after_success, expected_runs=1, maximum_restarts=2
            )

    def test_solver_failure_events_contribute_to_scientific_hash(self) -> None:
        runner = load_runner()
        before = runner.scientific_hash({
            "status": "FAILED",
            "solver_events": [{"type": "evaluation", "reason": "boolean_failure"}],
        })
        after = runner.scientific_hash({
            "status": "FAILED",
            "solver_events": [{"type": "evaluation", "reason": "singular_jacobian"}],
        })
        self.assertNotEqual(before, after)

    def test_attempt_ledger_rejects_coerced_counters_and_foreign_event_order(self):
        runner = load_runner()
        attempt = {"run":0,"restart":1,"status":"FAILED","function_evaluations":1,
                   "solver_events":[{"type":"attempt_start","restart":0},
                                    {"type":"evaluation","restart":0,
                                     "evaluation_index":1,"restart_evaluation_index":1},
                                    {"type":"attempt_terminal","restart":0}]}
        for field, value in (("run",False),("run",0.9),("restart",1.8),
                             ("function_evaluations",1.9)):
            with self.subTest(field=field,value=value):
                corrupt = copy.deepcopy(attempt)
                corrupt[field] = value
                with self.assertRaisesRegex(ValueError,"integer"):
                    runner.validate_attempt_ledger([corrupt],expected_runs=1,maximum_restarts=1)
        foreign = copy.deepcopy(attempt)
        foreign["solver_events"][1]["restart"] = 1
        with self.subTest(event_failure="foreign_restart"):
            with self.assertRaisesRegex(ValueError,"restart identity"):
                runner.validate_attempt_ledger([foreign],expected_runs=1,maximum_restarts=1)
        reordered = copy.deepcopy(attempt)
        reordered["solver_events"] = list(reversed(reordered["solver_events"]))
        with self.subTest(event_failure="reversed_order"):
            with self.assertRaisesRegex(ValueError,"event order"):
                runner.validate_attempt_ledger([reordered],expected_runs=1,maximum_restarts=1)

    def test_aggregate_timing_summaries_do_not_change_scientific_hash(self) -> None:
        runner = load_runner()
        first = runner.scientific_hash({
            "status": "CONVERGED",
            "summary": {
                "median_residual_seconds": 0.5,
                "median_derivative_seconds": 0.2,
                "median_newton_iterations": 7,
            },
        })
        second = runner.scientific_hash({
            "status": "CONVERGED",
            "summary": {
                "median_residual_seconds": 0.8,
                "median_derivative_seconds": 0.3,
                "median_newton_iterations": 7,
            },
        })
        self.assertEqual(first, second)

    def test_two_start_plan_records_all_reconstructed_cases(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "plan.json"
            subprocess.run(
                ["python3", str(RUNNER), "--plan-only", "--starts", "2", "--output", str(output)],
                check=True,
            )
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(payload["schema"], "spin-it-faster-table1-plan-v2")
            self.assertEqual(payload["starts_per_case"], 2)
            self.assertEqual(len(payload["cases"]), 8)
            self.assertEqual(payload["initialization_protocol"], "PER_CASE_PROTOCOL_V1")
            for case in payload["cases"]:
                if case["mode"] in {"SPINNING_TOP","YO_YO"}:
                    self.assertEqual(case["status"],"BLOCKED_MISSING_CHERCHI_PROBLEM")
                    self.assertEqual(case["solver_family"],"CHERCHI_RAW_EQ14_V1")
                    self.assertEqual(case["runs"],[])
                    self.assertNotIn("initial_quadrics",case)
                    continue
                if case["mode"] == "FIDGET_SPINNER":
                    self.assertEqual(case["status"], "BLOCKED_MISSING_PLANAR_GEOMETRY")
                    self.assertEqual(case["runs"], [])
                    self.assertNotIn("initial_quadrics", case)
                    continue
                if case["mode"] == "BUOYANT_STABILITY":
                    self.assertEqual(case["status"], "BLOCKED_MISSING_BUOYANT_GEOMETRY")
                    self.assertEqual(case["solver_family"], "BUOYANT_PLANE_GEOMETRY_V1")
                    self.assertEqual(case["runs"], [])
                    self.assertNotIn("initial_quadrics", case)
                    continue
                if case["mode"] == "STATIC_BALANCE":
                    self.assertEqual(case["solver_family"], "STATIC_PLANE_FAMILY_V1")
                    self.assertEqual(len(case["runs"]), 2)
                    self.assertNotIn("initial_quadrics", case)
                    continue
                self.assertEqual(len(case["initial_quadrics"]), 2)
                self.assertEqual([row["seed"] for row in case["runs"]],[0,1])
                self.assertTrue(all(len(row["initial_quadrics"])==8 for row in case["runs"]))

    def test_plan_refuses_asset_bytes_that_disagree_with_manifest_digest(self) -> None:
        runner = load_runner()
        case = copy.deepcopy(json.loads(runner.DEFAULT_CASES.read_text(encoding="utf-8"))["cases"][0])
        case["design_domain"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "design_domain.*SHA-256"):
            runner.plan_payload([case], starts=1)

        case = copy.deepcopy(json.loads(runner.DEFAULT_CASES.read_text(encoding="utf-8"))["cases"][0])
        case["outer"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "outer.*SHA-256"):
            runner.plan_payload([case], starts=1)

    def test_plan_refuses_asset_paths_outside_manifest_directory(self) -> None:
        runner = load_runner()
        case = copy.deepcopy(json.loads(runner.DEFAULT_CASES.read_text(encoding="utf-8"))["cases"][0])
        case["design_domain"]["path"] = "../outside.obj"
        with self.assertRaisesRegex(ValueError, "design_domain.*inside"):
            runner.plan_payload([case], starts=1)


class SpinItFasterTable1NativeRunnerTests(unittest.TestCase):
    def test_duplicate_or_missing_evaluation_ordinals_cannot_replace_failed_trials(self):
        runner = load_runner()
        case = next(row for row in json.loads(runner.DEFAULT_CASES.read_text())["cases"]
                    if row["id"] == "static-tilted-vase")
        actual = runner.run_legacy_quadric_case(case, starts=1, max_restarts=1)
        runner.validate_attempt_ledger(actual["attempts"], expected_runs=1,
                                       maximum_restarts=1)
        evaluations = [event for event in actual["attempts"][0]["solver_events"]
                       if event["type"] == "evaluation"]
        self.assertGreaterEqual(len(evaluations), 2)
        for field in ("evaluation_index", "restart_evaluation_index"):
            for corruption in ("duplicate", "missing", "boolean"):
                with self.subTest(field=field, corruption=corruption):
                    forged = copy.deepcopy(actual["attempts"])
                    rows = [event for event in forged[0]["solver_events"]
                            if event["type"] == "evaluation"]
                    if corruption == "duplicate":
                        rows[1][field] = rows[0][field]
                    elif corruption == "missing":
                        rows[1].pop(field)
                    else:
                        rows[0][field] = True
                    with self.assertRaisesRegex(ValueError, "evaluation.*ordinal"):
                        runner.validate_attempt_ledger(forged, expected_runs=1,
                                                       maximum_restarts=1)

    def test_changed_success_label_cannot_override_actual_failed_native_terminal(self):
        runner = load_runner()
        case = next(row for row in json.loads(runner.DEFAULT_CASES.read_text())["cases"]
                    if row["id"] == "static-tilted-vase")
        actual = runner.run_legacy_quadric_case(case,starts=1,max_restarts=1)
        self.assertNotEqual(actual["attempts"][0]["status"],"CONVERGED")
        forged = copy.deepcopy(actual["attempts"])
        forged[0]["status"] = "CONVERGED"
        with self.assertRaisesRegex(ValueError,"convergence.*terminal"):
            runner.validate_attempt_ledger(forged,expected_runs=1,maximum_restarts=1)
        # Removing the contradictory terminal reason must not hide failure.
        forged[0]["solver_events"][-1].pop("reason")
        with self.assertRaisesRegex(ValueError,"convergence.*terminal"):
            runner.validate_attempt_ledger(forged,expected_runs=1,maximum_restarts=1)

    def test_plan_startpoints_match_actual_independent_top_level_runs(self):
        runner = load_runner()
        case = next(row for row in json.loads(runner.DEFAULT_CASES.read_text())["cases"]
                    if row["id"] == "static-tilted-vase")
        planned = runner.plan_payload([case],starts=2)["cases"][0]
        actual = runner.run_case(case,starts=2,max_restarts=1)
        first_attempts = [row["initial_plane"] for row in actual["attempts"]
                          if row["restart"] == 1]
        self.assertEqual([row["initial_planes"][0] for row in planned["runs"]],first_attempts)

    def test_failed_initial_evaluation_counts_even_without_accepted_trace(self):
        runner = load_runner()
        case = next(row for row in json.loads(runner.DEFAULT_CASES.read_text())["cases"]
                    if row["id"] == "static-tilted-vase")
        # Force a deliberately invalid constant initial shape, but execute the
        # real native solver: failed initial evaluations still spend budget.
        with patch.object(runner, "generate_initial_quadrics", return_value=[[-1.,0.,0.,0.,0.,0.,0.,0.,0.,0.]]):
            payload = runner.run_legacy_quadric_case(case,starts=1,max_restarts=1)
        attempt = payload["attempts"][0]
        self.assertNotEqual(attempt["status"], "CONVERGED")
        self.assertEqual(attempt["trace"], [])
        self.assertEqual(attempt["function_evaluations"], 1)
        self.assertEqual(sum(e["type"] == "evaluation" for e in attempt["solver_events"]), 1)

    def test_historical_two_start_miniature_retains_its_explicit_seed_regression(self) -> None:
        runner = load_runner()
        manifest = json.loads(runner.DEFAULT_CASES.read_text(encoding="utf-8"))
        case = next(row for row in manifest["cases"] if row["id"] == "static-tilted-vase")
        # Preserve the historical frozen regression, never substitute these
        # successful seeds for the new default 0/1 batch or its failures.
        payload = runner.run_legacy_quadric_case(case, starts=2, max_restarts=1,
                                  initial_seeds=(14375009,14376018))
        self.assertEqual(payload["seed_sequence"],[14375009,14376018])
        self.assertTrue(payload['legacy_numerical_passed'])
        self.assertFalse(payload['passed'])
        self.assertEqual(payload['material_validation'], 'NOT_VERIFIED_BY_THIS_LEGACY_DISPATCH')
        self.assertTrue(payload["summary"]["all_runs_converged"])
        self.assertEqual(payload["summary"]["successful_runs"], 2)
        self.assertEqual(payload["summary"]["average_restarts"], 1.0)
        self.assertEqual(len(payload["attempts"]), 2)
        for attempt in payload["attempts"]:
            self.assertFalse(attempt['numerical_acceptance']['passed'])
            self.assertIn('RAW_EQ14_PROFILE_REQUIRED', attempt['numerical_acceptance']['reasons'])
            self.assertLessEqual(attempt["function_evaluations"], 60)
            self.assertGreater(attempt["residual_seconds"], 0.0)
            self.assertGreater(attempt["derivative_seconds"], 0.0)
            self.assertTrue(attempt["trace"])
            self.assertIn("raw_equation_residual", attempt)
            for row in attempt["trace"]:
                self.assertIn("raw_equation_residual", row)
                self.assertIn("equation_profile", row)
                self.assertIn("solver_coordinate_identity", row)
            self.assertEqual(attempt["trace"][-1]["raw_equation_residual"],
                             attempt["raw_equation_residual"])
            self.assertEqual(
                sum(event["type"] == "evaluation" for event in attempt["solver_events"]),
                attempt["function_evaluations"],
            )
            self.assertEqual(
                sum(event["type"] == "attempt_start" for event in attempt["solver_events"]),
                1,
            )
            self.assertEqual(
                sum(event["type"] == "attempt_terminal" for event in attempt["solver_events"]),
                1,
            )


if __name__ == "__main__":
    unittest.main()
