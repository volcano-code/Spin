import copy
import math
import unittest

import run_campaign


HASH_A = "a" * 64
HASH_B = "b" * 64
HASH_C = "c" * 64
HASH_D = "d" * 64
HASH_E = "e" * 64


def valid_record():
    quadric = [0.1 * i for i in range(10)]
    moments = [1.0 + 0.1 * i for i in range(10)]
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
        "schema": "m1-six-case-canonical-v1", "z": 0, "weight": 0.1,
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


def early_native_failure(reason, evaluator_returned=None):
    record = valid_record()
    record["numerical_status"] = "GEOMETRY_EVALUATION_FAILED"
    record["numerical_diagnostic"] = "native evaluation did not produce a finite Newton state"
    record["final_verified"] = False
    record["final_diagnostic"] = "not audited: numerical search did not converge; no final material returned"
    record["raw_residual"] = None
    record["independent_raw_residual"] = None
    record["physical_residual"] = "NONFINITE"
    record["normalized_residual"] = "NONFINITE"
    record["mu"] = 0.0
    record["multipliers"] = []
    record["quadric"] = None
    record["final_geometry_hash"] = None
    record["final_moments"] = None
    record["timing"]["final_audit_seconds"] = 0.0
    record["timing"]["solver_other_seconds"] = 0.3
    event = record["events"][1]
    event["outcome"] = "failed"
    event["reason"] = reason
    event["detail"] = "invalid native state"
    event["normalized_residual"] = None
    event["physical_residual"] = None
    event["raw_residual"] = None
    event["objective"] = None
    event["mu"] = None
    record["trace"] = []
    timing = record["evaluation_timings"][0]
    if evaluator_returned is None:
        evaluator_returned = reason == "invalid_numerical_state"
    if evaluator_returned:
        timing["returned"] = True
    else:
        event["geometry_hash"] = None
        event["topology_hash"] = None
        timing["returned"] = False
        timing["reported_residual_seconds"] = None
        timing["reported_derivative_seconds"] = None
    if reason == "invalid_numerical_state":
        record["events"][-1]["reason"] = "initial_numerical_state_invalid"
    else:
        record["events"][-1]["reason"] = "initial_evaluation_failed"
    return record


class RecordValidation(unittest.TestCase):
    def assert_rejected(self, mutate):
        record = valid_record()
        mutate(record)
        with self.assertRaises(ValueError):
            run_campaign.validate_record(record, 0, 0.1)

    def test_accepts_complete_valid_record(self):
        run_campaign.validate_record(valid_record(), 0, 0.1)

    def test_accepts_structured_failed_search_with_retained_trace(self):
        record = valid_record()
        record["numerical_status"] = "LINE_SEARCH_EXHAUSTED"
        record["final_verified"] = False
        record["final_diagnostic"] = "not audited"
        record["independent_raw_residual"] = None
        record["multipliers"] = []
        record["quadric"] = None
        record["final_geometry_hash"] = None
        record["final_moments"] = None
        record["timing"]["final_audit_seconds"] = 0.0
        record["timing"]["solver_other_seconds"] = 0.3
        record["events"][-1]["reason"] = "line_search_exhausted"
        run_campaign.validate_record(record, 0, 0.1)

    def test_accepts_native_initial_failures_without_an_accepted_trace(self):
        cases = (
            ("invalid_numerical_state", True),
            ("evaluation_exception", False),
            ("evaluation_exception", True),
        )
        for reason, returned in cases:
            with self.subTest(reason=reason, returned=returned):
                try:
                    run_campaign.validate_record(
                        early_native_failure(reason, returned), 0, 0.1)
                except ValueError as error:
                    self.fail(f"legitimate native failure was rejected: {error}")

    def test_rejects_impossible_native_failure_timing_shapes(self):
        corruptions = {
            "pre-return exception cannot have a material hash": lambda r:
                r["events"][1].__setitem__("geometry_hash", HASH_B),
            "invalid state cannot precede evaluator return": lambda r:
                r["evaluation_timings"][0].__setitem__("returned", False),
            "returned invalid state must retain reported timings": lambda r:
                r["evaluation_timings"][0].__setitem__("reported_residual_seconds", None),
        }
        for name, mutate in corruptions.items():
            with self.subTest(name=name):
                record = early_native_failure("invalid_numerical_state")
                if name.startswith("pre-return"):
                    record = early_native_failure("evaluation_exception", False)
                mutate(record)
                with self.assertRaises(ValueError):
                    run_campaign.validate_record(record, 0, 0.1)

    def test_rejects_nonfinite_tags_after_any_accepted_state(self):
        for key in ("physical_residual", "normalized_residual"):
            with self.subTest(key=key):
                self.assert_rejected(lambda record, key=key:
                                     record.__setitem__(key, "NONFINITE"))

    def test_rejects_status_terminal_reason_mismatches(self):
        corruptions = {
            "converged status with failure terminal":
                lambda r: r["events"][-1].__setitem__("reason", "line_search_exhausted"),
            "failure status with converged terminal": lambda r: (
                r.__setitem__("numerical_status", "LINE_SEARCH_EXHAUSTED"),
                r.__setitem__("final_verified", False),
                r.__setitem__("independent_raw_residual", None),
                r.__setitem__("multipliers", []),
                r.__setitem__("quadric", None),
                r.__setitem__("final_geometry_hash", None),
                r.__setitem__("final_moments", None)),
            "resource status with line-search terminal": lambda r: (
                r.__setitem__("numerical_status", "RESOURCE_LIMIT_REACHED"),
                r.__setitem__("final_verified", False),
                r.__setitem__("independent_raw_residual", None),
                r.__setitem__("multipliers", []),
                r.__setitem__("quadric", None),
                r.__setitem__("final_geometry_hash", None),
                r.__setitem__("final_moments", None),
                r["events"][-1].__setitem__("reason", "line_search_exhausted")),
        }
        for name, mutate in corruptions.items():
            with self.subTest(name=name):
                self.assert_rejected(mutate)

    def test_rejects_reviewed_top_level_corruptions(self):
        corruptions = {
            "wrong schema": lambda r: r.__setitem__("schema", "wrong-v9"),
            "restart count 99": lambda r: r.__setitem__("restarts", 99),
            "nonfinite successful residual": lambda r: r.__setitem__("raw_residual", "NONFINITE"),
            "verified raw residual over tolerance": lambda r: r.__setitem__("raw_residual", 2e-12),
            "verified independent raw residual over tolerance": lambda r: r.__setitem__("independent_raw_residual", 2e-12),
            "verified physical residual over tolerance": lambda r: r.__setitem__("physical_residual", 1.0),
            "unknown status": lambda r: r.__setitem__("numerical_status", "MYSTERY"),
            "missing result field": lambda r: r.pop("final_diagnostic"),
            "extra result field": lambda r: r.__setitem__("surprise", 1),
            "boolean evaluation count": lambda r: r.__setitem__("evaluations", True),
        }
        for name, mutate in corruptions.items():
            with self.subTest(name=name):
                self.assert_rejected(mutate)

    def test_rejects_deleted_or_malformed_event_stream(self):
        corruptions = {
            "deleted evaluation event": lambda r: r["events"].pop(1),
            "missing event field": lambda r: r["events"][1].pop("objective"),
            "extra event field": lambda r: r["events"][1].__setitem__("surprise", 1),
            "wrong evaluation outcome": lambda r: r["events"][1].__setitem__("outcome", "invented"),
            "wrong event restart": lambda r: r["events"][1].__setitem__("restart", 1),
            "bad event hash": lambda r: r["events"][1].__setitem__("geometry_hash", "not-a-hash"),
        }
        for name, mutate in corruptions.items():
            with self.subTest(name=name):
                self.assert_rejected(mutate)

    def test_rejects_deleted_or_inconsistent_accepted_trace(self):
        corruptions = {
            "deleted trace": lambda r: r.__setitem__("trace", []),
            "missing trace field": lambda r: r["trace"][0].pop("objective"),
            "trace no longer matches accepted event": lambda r: r["trace"][0].__setitem__("objective", 9.0),
            "bad trace hash": lambda r: r["trace"][0].__setitem__("topology_hash", "bad"),
        }
        for name, mutate in corruptions.items():
            with self.subTest(name=name):
                self.assert_rejected(mutate)

    def test_rejects_restart_evaluation_and_timing_accounting_corruption(self):
        corruptions = {
            "evaluation budget exceeded": lambda r: r.__setitem__("evaluations", 21),
            "restart subtotal mismatch": lambda r: r.__setitem__("restart_evaluations", [0]),
            "wrong timing index": lambda r: r["evaluation_timings"][0].__setitem__("evaluation_index", 2),
            "returned flag disagrees with event": lambda r: r["evaluation_timings"][0].__setitem__("returned", False),
            "missing timing field": lambda r: r["timing"].pop("entry_through_record_seconds"),
            "nonfinite timing": lambda r: r["timing"].__setitem__("solver_other_seconds", math.inf),
        }
        for name, mutate in corruptions.items():
            with self.subTest(name=name):
                self.assert_rejected(mutate)

    def test_rejects_invalid_result_hashes_and_verified_material(self):
        corruptions = {
            "bad child hash": lambda r: r.__setitem__("output_hash", "0" * 63),
            "missing verified geometry hash": lambda r: r.__setitem__("final_geometry_hash", None),
            "missing verified moments": lambda r: r.__setitem__("final_moments", None),
            "wrong verified multiplier shape": lambda r: r.__setitem__("multipliers", []),
            "wrong quadric length": lambda r: r.__setitem__("quadric", [0.0] * 9),
        }
        for name, mutate in corruptions.items():
            with self.subTest(name=name):
                self.assert_rejected(mutate)

    def test_collector_digest_binds_payload_not_only_child_hashes_and_excludes_timing(self):
        self.assertTrue(hasattr(run_campaign, "scientific_payload_sha256"))
        record = valid_record()
        original = run_campaign.scientific_payload_sha256(record)
        timing_only = copy.deepcopy(record)
        timing_only["timing"]["solve_wrapper_seconds"] = 2.0
        timing_only["timing"]["solver_other_seconds"] = 1.1
        timing_only["timing"]["entry_through_record_seconds"] = 2.01
        timing_only["evaluation_timings"][0]["wall_seconds"] = 0.2
        self.assertEqual(original, run_campaign.scientific_payload_sha256(timing_only))
        changed_payload = copy.deepcopy(record)
        changed_payload["events"][1]["objective"] = 1.5
        self.assertEqual(record["event_hash"], changed_payload["event_hash"])
        self.assertNotEqual(original, run_campaign.scientific_payload_sha256(changed_payload))

    def test_process_overhead_records_outer_partition_and_rejects_bad_values(self):
        self.assertTrue(hasattr(run_campaign, "process_overhead_seconds"))
        record = valid_record()
        self.assertAlmostEqual(run_campaign.process_overhead_seconds(1.25, record), 0.25)
        with self.assertRaises(ValueError):
            run_campaign.process_overhead_seconds(0.5, record)
        with self.assertRaises(ValueError):
            run_campaign.process_overhead_seconds(math.inf, record)


if __name__ == "__main__":
    unittest.main()
