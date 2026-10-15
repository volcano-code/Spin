import copy
import json
import math
import os
import subprocess
import unittest

import run_budget_campaign


HASH_A = "a" * 64
HASH_B = "b" * 64
HASH_C = "c" * 64
HASH_D = "d" * 64
HASH_E = "e" * 64


def valid_record(evaluation_budget=20):
    quadric = [0.1 * index for index in range(10)]
    moments = [1.0 + 0.1 * index for index in range(10)]
    accepted = {
        "type": "evaluation", "restart": 0, "iteration": 0,
        "evaluation_index": 1, "restart_evaluation_index": 1,
        "step": 0.0, "outcome": "accepted",
        "reason": "initial_evaluation_succeeded", "detail": "",
        "quadric_coefficients": quadric, "quadric_hash": HASH_A,
        "normalized_residual": 2.5e-13, "physical_residual": 5e-7,
        "objective": 1.25, "mu": 0.75, "geometry_hash": HASH_B,
        "topology_hash": HASH_C, "profile": "DUAL_MONOMIAL_RAW_V1",
        "raw_residual": 5e-13,
    }
    return {
        "schema": "m1-budget-canonical-v1", "z": 0, "weight": 0.1,
        "evaluation_budget": evaluation_budget,
        "numerical_status": "CONVERGED", "numerical_diagnostic": "",
        "final_verified": True, "final_diagnostic": "audit passed",
        "raw_residual": 5e-13, "independent_raw_residual": 4e-13,
        "physical_residual": 5e-7, "normalized_residual": 2.5e-13,
        "mu": 0.75, "multipliers": [0.1, 0.2, 0.3, 0.4],
        "initial_quadric": quadric, "quadric": quadric,
        "input_hash": HASH_A, "output_hash": HASH_B,
        "trace_hash": HASH_C, "event_hash": HASH_D,
        "final_geometry_hash": HASH_E, "final_moments": moments,
        "evaluations": 1, "restarts": 1, "restart_evaluations": [1],
        "event_stream_version": 2,
        "timing": {
            "solve_wrapper_seconds": 1.0,
            "coupled_evaluator_seconds": 0.7,
            "final_audit_seconds": 0.2,
            "solver_other_seconds": 0.1,
            "evaluator_reported_residual_seconds": 0.6,
            "evaluator_reported_derivative_seconds": 0.1,
            "entry_through_record_seconds": 1.01,
        },
        "events": [
            {
                "type": "attempt_start", "restart": 0, "iteration": 0,
                "evaluation_index": 0, "restart_evaluation_index": 0,
                "step": None, "outcome": "started", "reason": "restart_started",
                "detail": "", "quadric_coefficients": quadric,
                "quadric_hash": HASH_A, "normalized_residual": None,
                "physical_residual": None, "objective": None, "mu": None,
                "geometry_hash": None, "topology_hash": None,
                "profile": "DUAL_MONOMIAL_RAW_V1", "raw_residual": None,
            },
            accepted,
            {
                "type": "attempt_terminal", "restart": 0, "iteration": 0,
                "evaluation_index": 1, "restart_evaluation_index": 1,
                "step": None, "outcome": "terminated", "reason": "converged",
                "detail": "", "quadric_coefficients": None, "quadric_hash": "",
                "normalized_residual": None, "physical_residual": None,
                "objective": None, "mu": None, "geometry_hash": None,
                "topology_hash": None, "profile": "DUAL_MONOMIAL_RAW_V1",
                "raw_residual": None,
            },
        ],
        "trace": [{
            "restart": 0, "iteration": 0, "evaluations": 1,
            "normalized_residual": 2.5e-13, "physical_residual": 5e-7,
            "objective": 1.25, "mu": 0.75, "accepted_geodesic_step": 0.0,
            "quadric_hash": HASH_A, "geometry_hash": HASH_B,
            "topology_hash": HASH_C, "profile": "DUAL_MONOMIAL_RAW_V1",
            "raw_residual": 5e-13, "quadric_coefficients": quadric,
        }],
        "evaluation_timings": [{
            "evaluation_index": 1, "restart": 0,
            "restart_evaluation_index": 1, "returned": True,
            "wall_seconds": 0.1, "reported_residual_seconds": 0.08,
            "reported_derivative_seconds": 0.02,
        }],
    }


def make_unverified(record, status, terminal_reason):
    record["numerical_status"] = status
    record["final_verified"] = False
    record["final_diagnostic"] = "not audited"
    record["independent_raw_residual"] = None
    record["multipliers"] = []
    record["quadric"] = None
    record["final_geometry_hash"] = None
    record["final_moments"] = None
    record["timing"]["final_audit_seconds"] = 0.0
    record["timing"]["solver_other_seconds"] = 0.3
    record["events"][-1]["reason"] = terminal_reason
    return record


def early_native_failure(reason, evaluator_returned=None):
    record = make_unverified(
        valid_record(), "GEOMETRY_EVALUATION_FAILED",
        "initial_numerical_state_invalid"
        if reason == "invalid_numerical_state" else "initial_evaluation_failed")
    record["numerical_diagnostic"] = "native evaluation did not produce a finite Newton state"
    record["raw_residual"] = None
    record["physical_residual"] = "NONFINITE"
    record["normalized_residual"] = "NONFINITE"
    record["mu"] = 0.0
    event = record["events"][1]
    event["outcome"] = "failed"
    event["reason"] = reason
    event["detail"] = "invalid native state"
    for key in ("normalized_residual", "physical_residual", "raw_residual",
                "objective", "mu"):
        event[key] = None
    record["trace"] = []
    timing = record["evaluation_timings"][0]
    if evaluator_returned is None:
        evaluator_returned = reason == "invalid_numerical_state"
    timing["returned"] = evaluator_returned
    if not evaluator_returned:
        event["geometry_hash"] = None
        event["topology_hash"] = None
        timing["reported_residual_seconds"] = None
        timing["reported_derivative_seconds"] = None
    return record


class BudgetContract(unittest.TestCase):
    def test_accepts_each_registered_budget_and_binds_echo(self):
        for budget in (20, 40, 60):
            with self.subTest(budget=budget):
                run_budget_campaign.validate_record(
                    valid_record(budget), 0, 0.1, budget)

    def test_rejects_missing_wrong_and_boolean_budget(self):
        missing = valid_record()
        missing.pop("evaluation_budget")
        wrong = valid_record(40)
        boolean = valid_record()
        boolean["evaluation_budget"] = True
        cases = ((missing, 20), (wrong, 20), (boolean, 20),
                 (valid_record(), True), (valid_record(80), 80),
                 (valid_record(320), 320))
        for record, declared in cases:
            with self.subTest(record_budget=record.get("evaluation_budget"), declared=declared):
                with self.assertRaises(ValueError):
                    run_budget_campaign.validate_record(record, 0, 0.1, declared)

    def test_rejects_counter_above_declared_budget(self):
        record = valid_record(40)
        record["evaluations"] = 41
        with self.assertRaises(ValueError):
            run_budget_campaign.validate_record(record, 0, 0.1, 40)

    def test_rejects_resource_terminal_below_declared_budget(self):
        record = make_unverified(
            valid_record(40), "RESOURCE_LIMIT_REACHED",
            "maximum_evaluations_reached")
        with self.assertRaises(ValueError):
            run_budget_campaign.validate_record(record, 0, 0.1, 40)

    def test_budget20_comparison_strips_only_wrapper_and_timing(self):
        current = valid_record()
        legacy = copy.deepcopy(current)
        legacy["schema"] = "m1-six-case-canonical-v1"
        legacy.pop("evaluation_budget")
        legacy["timing"]["solve_wrapper_seconds"] = 2.0
        legacy["evaluation_timings"][0]["wall_seconds"] = 2.0
        self.assertEqual(
            run_budget_campaign.legacy_comparable_payload(current),
            run_budget_campaign.legacy_comparable_payload(legacy))
        legacy["events"][1]["objective"] = 9.0
        self.assertNotEqual(
            run_budget_campaign.legacy_comparable_payload(current),
            run_budget_campaign.legacy_comparable_payload(legacy))

    def test_cross_budget_prefix_uses_only_numeric_evaluations(self):
        left = valid_record(20)
        right = valid_record(40)
        right["input_hash"] = "f" * 64
        right["events"][-1]["reason"] = "maximum_iterations_reached"
        self.assertTrue(run_budget_campaign.compare_numeric_prefix(left, right)["match"])
        right["events"][1]["raw_residual"] = 6e-13
        mismatch = run_budget_campaign.compare_numeric_prefix(left, right)
        self.assertFalse(mismatch["match"])
        self.assertEqual(mismatch["first_mismatch_evaluation"], 1)
        self.assertEqual(mismatch["field"], "raw_residual")


class FrozenRegressionBehaviors(unittest.TestCase):
    def validate(self, record):
        run_budget_campaign.validate_record(record, 0, 0.1, 20)

    def assert_rejected(self, mutate):
        record = valid_record()
        mutate(record)
        with self.assertRaises(ValueError):
            self.validate(record)

    def test_accepts_complete_valid_record(self):
        self.validate(valid_record())

    def test_accepts_structured_failed_search_with_retained_trace(self):
        self.validate(make_unverified(
            valid_record(), "LINE_SEARCH_EXHAUSTED", "line_search_exhausted"))

    def test_accepts_native_initial_failures_without_an_accepted_trace(self):
        for reason, returned in (("invalid_numerical_state", True),
                                 ("evaluation_exception", False),
                                 ("evaluation_exception", True)):
            with self.subTest(reason=reason, returned=returned):
                self.validate(early_native_failure(reason, returned))

    def test_rejects_impossible_native_failure_timing_shapes(self):
        records = [
            early_native_failure("evaluation_exception", False),
            early_native_failure("invalid_numerical_state", True),
            early_native_failure("invalid_numerical_state", True),
        ]
        records[0]["events"][1]["geometry_hash"] = HASH_B
        records[1]["evaluation_timings"][0]["returned"] = False
        records[2]["evaluation_timings"][0]["reported_residual_seconds"] = None
        for record in records:
            with self.assertRaises(ValueError):
                self.validate(record)

    def test_rejects_nonfinite_tags_after_any_accepted_state(self):
        for key in ("physical_residual", "normalized_residual"):
            with self.subTest(key=key):
                self.assert_rejected(lambda record, key=key:
                                     record.__setitem__(key, "NONFINITE"))

    def test_rejects_status_terminal_reason_mismatches(self):
        def failed_with_converged(record):
            make_unverified(record, "LINE_SEARCH_EXHAUSTED", "converged")

        def resource_with_line_search(record):
            make_unverified(record, "RESOURCE_LIMIT_REACHED", "line_search_exhausted")

        for mutate in (
                lambda record: record["events"][-1].__setitem__(
                    "reason", "line_search_exhausted"),
                failed_with_converged, resource_with_line_search):
            self.assert_rejected(mutate)

    def test_rejects_reviewed_top_level_corruptions(self):
        corruptions = (
            lambda record: record.__setitem__("schema", "wrong-v9"),
            lambda record: record.__setitem__("restarts", 99),
            lambda record: record.__setitem__("raw_residual", "NONFINITE"),
            lambda record: record.__setitem__("raw_residual", 2e-12),
            lambda record: record.__setitem__("independent_raw_residual", 2e-12),
            lambda record: record.__setitem__("physical_residual", 1.0),
            lambda record: record.__setitem__("numerical_status", "MYSTERY"),
            lambda record: record.pop("final_diagnostic"),
            lambda record: record.__setitem__("surprise", 1),
            lambda record: record.__setitem__("evaluations", True),
        )
        for mutate in corruptions:
            self.assert_rejected(mutate)

    def test_rejects_deleted_or_malformed_event_stream(self):
        corruptions = (
            lambda record: record["events"].pop(1),
            lambda record: record["events"][1].pop("objective"),
            lambda record: record["events"][1].__setitem__("surprise", 1),
            lambda record: record["events"][1].__setitem__("outcome", "invented"),
            lambda record: record["events"][1].__setitem__("restart", 1),
            lambda record: record["events"][1].__setitem__("geometry_hash", "not-a-hash"),
        )
        for mutate in corruptions:
            self.assert_rejected(mutate)

    def test_rejects_deleted_or_inconsistent_accepted_trace(self):
        corruptions = (
            lambda record: record.__setitem__("trace", []),
            lambda record: record["trace"][0].pop("objective"),
            lambda record: record["trace"][0].__setitem__("objective", 9.0),
            lambda record: record["trace"][0].__setitem__("topology_hash", "bad"),
        )
        for mutate in corruptions:
            self.assert_rejected(mutate)

    def test_rejects_restart_evaluation_and_timing_accounting_corruption(self):
        corruptions = (
            lambda record: record.__setitem__("evaluations", 21),
            lambda record: record.__setitem__("restart_evaluations", [0]),
            lambda record: record["evaluation_timings"][0].__setitem__(
                "evaluation_index", 2),
            lambda record: record["evaluation_timings"][0].__setitem__(
                "returned", False),
            lambda record: record["timing"].pop("entry_through_record_seconds"),
            lambda record: record["timing"].__setitem__("solver_other_seconds", math.inf),
        )
        for mutate in corruptions:
            self.assert_rejected(mutate)

    def test_rejects_invalid_result_hashes_and_verified_material(self):
        corruptions = (
            lambda record: record.__setitem__("output_hash", "0" * 63),
            lambda record: record.__setitem__("final_geometry_hash", None),
            lambda record: record.__setitem__("final_moments", None),
            lambda record: record.__setitem__("multipliers", []),
            lambda record: record.__setitem__("quadric", [0.0] * 9),
        )
        for mutate in corruptions:
            self.assert_rejected(mutate)

    def test_collector_digest_binds_payload_and_excludes_timing(self):
        record = valid_record()
        original = run_budget_campaign.scientific_payload_sha256(record)
        timing_only = copy.deepcopy(record)
        timing_only["timing"]["solve_wrapper_seconds"] = 2.0
        timing_only["timing"]["solver_other_seconds"] = 1.1
        timing_only["timing"]["entry_through_record_seconds"] = 2.01
        timing_only["evaluation_timings"][0]["wall_seconds"] = 0.2
        self.assertEqual(original,
                         run_budget_campaign.scientific_payload_sha256(timing_only))
        changed_payload = copy.deepcopy(record)
        changed_payload["events"][1]["objective"] = 1.5
        self.assertEqual(record["event_hash"], changed_payload["event_hash"])
        self.assertNotEqual(
            original, run_budget_campaign.scientific_payload_sha256(changed_payload))

    def test_process_overhead_records_outer_partition_and_rejects_bad_values(self):
        record = valid_record()
        self.assertAlmostEqual(
            run_budget_campaign.process_overhead_seconds(1.25, record), 0.25)
        with self.assertRaises(ValueError):
            run_budget_campaign.process_overhead_seconds(0.5, record)
        with self.assertRaises(ValueError):
            run_budget_campaign.process_overhead_seconds(math.inf, record)


class ProbeArguments(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.binary = os.environ.get("M1_BUDGET_PROBE")
        if not cls.binary:
            raise unittest.SkipTest("M1_BUDGET_PROBE is not set")

    def invoke(self, *arguments):
        return subprocess.run(
            [self.binary, *arguments], capture_output=True, text=True,
            check=False, timeout=30)

    def test_accepts_budget20_and_echoes_it(self):
        completed = self.invoke("0", "0.1", "20")
        self.assertEqual(completed.returncode, 0, completed.stderr)
        lines = [line.removeprefix("M1_RESULT ")
                 for line in completed.stdout.splitlines()
                 if line.startswith("M1_RESULT ")]
        self.assertEqual(len(lines), 1)
        record = json.loads(lines[0])
        self.assertEqual(record["schema"], "m1-budget-canonical-v1")
        self.assertEqual(record["evaluation_budget"], 20)

    def test_recognizes_budget40_and60_before_case_validation_without_solving(self):
        for budget in ("40", "60"):
            with self.subTest(budget=budget):
                completed = self.invoke("3", "0.1", budget)
                self.assertEqual(completed.returncode, 2)
                self.assertNotIn("M1_RESULT ", completed.stdout)
                self.assertIn("case outside preregistered matrix", completed.stderr)
                self.assertNotIn("budget outside preregistered matrix", completed.stderr)

    def test_rejects_unknown_negative_bool_like_and_trailing_budget_without_result(self):
        for budget in ("21", "-20", "true", "20x", "+20", "020", "20.0",
                       "80", "320"):
            with self.subTest(budget=budget):
                completed = self.invoke("0", "0.1", budget)
                self.assertEqual(completed.returncode, 2)
                self.assertNotIn("M1_RESULT ", completed.stdout)
                self.assertIn("M1_HARNESS_ERROR", completed.stderr)

    def test_requires_exactly_three_arguments(self):
        for arguments in (("0", "0.1"), ("0", "0.1", "20", "extra")):
            with self.subTest(arguments=arguments):
                completed = self.invoke(*arguments)
                self.assertEqual(completed.returncode, 2)
                self.assertNotIn("M1_RESULT ", completed.stdout)


if __name__ == "__main__":
    unittest.main()
