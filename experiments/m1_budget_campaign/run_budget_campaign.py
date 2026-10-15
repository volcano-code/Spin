"""Frozen bounded-budget M1 diagnostic; synthetic only, no B0/B1 claim."""
import argparse
import hashlib
import json
import math
import os
import re
import signal
import statistics
import subprocess
import time
from pathlib import Path


BUDGETS = (20, 40, 60)
CASES = tuple((z, weight) for z in (0, 2) for weight in (0.1, 1, 10))
TIMEOUT_SECONDS = 180
RESULT_KEYS = {
    "schema", "z", "weight", "evaluation_budget", "numerical_status",
    "numerical_diagnostic", "final_verified", "final_diagnostic",
    "raw_residual", "independent_raw_residual", "physical_residual",
    "normalized_residual", "mu", "multipliers", "initial_quadric",
    "quadric", "input_hash", "output_hash", "trace_hash", "event_hash",
    "final_geometry_hash", "final_moments", "evaluations", "restarts",
    "restart_evaluations", "event_stream_version", "timing", "events",
    "trace", "evaluation_timings",
}
TIMING_KEYS = {
    "solve_wrapper_seconds", "coupled_evaluator_seconds", "final_audit_seconds",
    "solver_other_seconds", "evaluator_reported_residual_seconds",
    "evaluator_reported_derivative_seconds", "entry_through_record_seconds",
}
EVENT_KEYS = {
    "type", "restart", "iteration", "evaluation_index",
    "restart_evaluation_index", "step", "outcome", "reason", "detail",
    "quadric_coefficients", "quadric_hash", "normalized_residual",
    "physical_residual", "objective", "mu", "geometry_hash",
    "topology_hash", "profile", "raw_residual",
}
TRACE_KEYS = {
    "restart", "iteration", "evaluations", "normalized_residual",
    "physical_residual", "objective", "mu", "accepted_geodesic_step",
    "quadric_hash", "geometry_hash", "topology_hash", "profile",
    "raw_residual", "quadric_coefficients",
}
EVALUATION_TIMING_KEYS = {
    "evaluation_index", "restart", "restart_evaluation_index", "returned",
    "wall_seconds", "reported_residual_seconds", "reported_derivative_seconds",
}
KNOWN_STATUSES = {
    "CONVERGED", "INFEASIBLE", "INVALID_REQUEST",
    "OPTIMIZATION_DID_NOT_CONVERGE", "LINE_SEARCH_EXHAUSTED",
    "GEOMETRY_EVALUATION_FAILED", "INVALID_MATERIAL_ARTIFACT",
    "RESOURCE_LIMIT_REACHED", "INDETERMINATE",
}
TERMINAL_REASONS = {
    "converged", "initial_evaluation_failed", "initial_numerical_state_invalid",
    "initial_sensitivity_vanished", "line_search_exhausted",
    "maximum_evaluations_reached", "maximum_iterations_reached",
    "newton_direction_invalid", "newton_direction_vanished",
    "newton_system_invalid", "bordered_system_rank_deficient",
}
EVALUATION_REASON_BY_OUTCOME = {
    "accepted": {"initial_evaluation_succeeded", "residual_decreased"},
    "rejected": {"no_residual_decrease"},
    "failed": {"evaluation_exception", "invalid_numerical_state"},
}
STATUS_TERMINAL_REASONS = {
    "CONVERGED": {"converged"},
    "GEOMETRY_EVALUATION_FAILED": {
        "initial_evaluation_failed", "initial_numerical_state_invalid",
    },
    "OPTIMIZATION_DID_NOT_CONVERGE": {
        "initial_sensitivity_vanished", "newton_system_invalid",
        "bordered_system_rank_deficient", "newton_direction_invalid",
        "newton_direction_vanished",
    },
    "LINE_SEARCH_EXHAUSTED": {"line_search_exhausted"},
    "RESOURCE_LIMIT_REACHED": {
        "maximum_evaluations_reached", "maximum_iterations_reached",
    },
}
PROFILE = "DUAL_MONOMIAL_RAW_V1"
HASH_PATTERN = re.compile(r"[0-9a-f]{64}")
SCIENTIFIC_TIMING_SIDECARS = {"timing", "evaluation_timings"}
NATIVE_HASH_KEYS = (
    "input_hash", "output_hash", "trace_hash", "event_hash",
    "final_geometry_hash",
)
PREFIX_KEYS = (
    "quadric_coefficients", "outcome", "reason", "step", "raw_residual",
    "normalized_residual", "physical_residual", "objective", "mu",
)
MATERIAL_NEGATIVE_OVERHEAD_SECONDS = -1e-6
PRIOR_SOURCE_IDENTITIES = {
    "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-domain-readiness-20260912/source/tests/canonical_campaign_probe.cpp":
        "53d4f91e1f8c7c0187504776a013825c2df2362bd2cae86d3007027fdc08d4c4",
    "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-domain-readiness-20260912/source/experiments/m1_canonical_campaign/CMakeLists.txt":
        "0ec4c51082759a9f05a93457dc8d6bec6798ef5d56e9fef939da4beb5f4a90bd",
    "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-domain-readiness-20260912/source/experiments/m1_canonical_campaign/run_campaign.py":
        "3950e4625b6572dffa614a321b873a9ff3a29406016dfdb65b06f8f2179aba94",
    "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-domain-readiness-20260912/source/experiments/m1_canonical_campaign/test_campaign_records.py":
        "ba084e47cffffcabf096c404e2c7c0bfe63e9c4721852db62de69a6c5ed60698",
}
REGISTRATION_IDENTITIES = {
    "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-budget-diagnostic-20260912/preregistration.json":
        "3c5b5768f7e0eada74076be56587add07f2109b766a8b6e6aca4ebd78e5b00df",
    "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-budget-diagnostic-20260912/preregistration-correction.json":
        "f4feb32c958bc69ee176ed048710e02ffba24bbf59e97096304447e01b90eb69",
}
FROZEN_NATIVE_IDENTITIES = {
    "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-replay-determinism-20260912/acceptance-build/libspin_canonical_replay.a":
        "cb5da5fd6d33c741a7b743d9ac26c5e49b07fe22e06f9854a5dcad46f5f48815",
    "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-general-differential-20260912/build/libspin_general_differential.a":
        "d9ed187b9603fef09bfa70765e05cbb8f2895af5de1dd3163feca57c473f06cc",
    "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-general-differential-20260912/build/native/libmake_it_stand_core.a":
        "7032f4ecad915edcdd2943a2f9d666a7ecb2841e62cc5830c14cd845468f3488",
    "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-plane-family-20260912/source/include/spin/SpinDesignRequest.h":
        "647f540a5d6fbf41bb00e8e5af1f569f9e3c80ab9f885a75de4e94d0ec6aed41",
    "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-plane-family-20260912/source/src/spin/SpinDesignRequest.cpp":
        "6c824059ec0cdd3cbdedc313f3cc8f7f95afb7076ca1f65e5314d5c99abaca47",
    "/data/xvxinyao/moment-guided-spin-v1/experiments/m1-plane-family-20260912/source/src/spin/QuadricHollowingSolver.cpp":
        "47d3f0f371d5766f8f1a2c6d0e1dc74ba39f94d5c90bdea018ec02aec5234490",
}
BASELINE_RECORDS_SHA256 = "af5072b19ea0d54a6ff4cfa0b8c691e30946a3f4b9c1f857c04de393fab32a32"


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _require_exact_keys(value, expected, label):
    _require(type(value) is dict, f"{label} must be an object")
    _require(set(value) == expected, f"{label} fields do not match schema")


def _is_int(value):
    return type(value) is int


def _is_finite_number(value):
    return type(value) in (int, float) and math.isfinite(value)


def _require_number(value, label, minimum=None):
    _require(_is_finite_number(value), f"{label} must be finite numeric")
    if minimum is not None:
        _require(value >= minimum, f"{label} is below its lower bound")


def _require_nullable_number(value, label, minimum=None):
    if value is not None:
        _require_number(value, label, minimum)


def _require_hash(value, label):
    _require(type(value) is str and HASH_PATTERN.fullmatch(value) is not None,
             f"{label} is not a lowercase SHA-256")


def _require_nullable_hash(value, label):
    if value is not None:
        _require_hash(value, label)


def _require_vector(value, length, label):
    _require(type(value) is list and len(value) == length,
             f"{label} must contain {length} numbers")
    for element in value:
        _require_number(element, label)


def _require_nullable_vector(value, length, label):
    if value is not None:
        _require_vector(value, length, label)


def _validate_event(event, count):
    _require_exact_keys(event, EVENT_KEYS, "event")
    for key in ("type", "outcome", "reason", "detail", "profile", "quadric_hash"):
        _require(type(event[key]) is str, f"event {key} must be a string")
    _require(event["profile"] == PROFILE, "event profile changed")
    for key in ("restart", "iteration", "evaluation_index", "restart_evaluation_index"):
        _require(_is_int(event[key]), f"event {key} must be an integer")
    _require(event["restart"] == 0, "event restart outside the single configured restart")
    _require(event["iteration"] >= 0, "event iteration is negative")
    _require(0 <= event["evaluation_index"] <= count, "event evaluation index out of range")
    _require(0 <= event["restart_evaluation_index"] <= count,
             "event restart evaluation index out of range")
    _require_nullable_number(event["step"], "event step", 0)
    _require_nullable_vector(event["quadric_coefficients"], 10, "event quadric")
    for key in ("normalized_residual", "physical_residual", "raw_residual"):
        _require_nullable_number(event[key], f"event {key}", 0)
    for key in ("objective", "mu"):
        _require_nullable_number(event[key], f"event {key}")
    _require_nullable_hash(event["geometry_hash"], "event geometry hash")
    _require_nullable_hash(event["topology_hash"], "event topology hash")
    event_type = event["type"]
    _require(event_type in {"attempt_start", "evaluation", "attempt_terminal"},
             "unknown event type")
    if event_type == "attempt_start":
        _require(event["outcome"] == "started" and event["reason"] == "restart_started",
                 "invalid attempt-start outcome")
        _require(event["evaluation_index"] == 0 and event["restart_evaluation_index"] == 0,
                 "attempt start has wrong evaluation accounting")
        _require(event["quadric_coefficients"] is not None, "attempt start lacks quadric")
        _require_hash(event["quadric_hash"], "attempt-start quadric hash")
        for key in ("step", "normalized_residual", "physical_residual", "raw_residual",
                    "objective", "mu", "geometry_hash", "topology_hash"):
            _require(event[key] is None, f"attempt start {key} must be null")
    elif event_type == "attempt_terminal":
        _require(event["outcome"] == "terminated" and event["reason"] in TERMINAL_REASONS,
                 "invalid attempt-terminal outcome")
        _require(event["evaluation_index"] == count
                 and event["restart_evaluation_index"] == count,
                 "attempt terminal has wrong evaluation accounting")
        _require(event["quadric_hash"] == "", "attempt terminal has unexpected quadric hash")
        for key in ("step", "quadric_coefficients", "normalized_residual",
                    "physical_residual", "raw_residual", "objective", "mu",
                    "geometry_hash", "topology_hash"):
            _require(event[key] is None, f"attempt terminal {key} must be null")
    else:
        outcome = event["outcome"]
        _require(outcome in EVALUATION_REASON_BY_OUTCOME, "unknown evaluation outcome")
        _require(event["reason"] in EVALUATION_REASON_BY_OUTCOME[outcome],
                 "evaluation outcome/reason mismatch")
        _require(event["evaluation_index"] >= 1
                 and event["restart_evaluation_index"] == event["evaluation_index"],
                 "evaluation has wrong accounting")
        _require(event["quadric_coefficients"] is not None, "evaluation lacks quadric")
        _require_hash(event["quadric_hash"], "evaluation quadric hash")
        if outcome == "failed":
            for key in ("normalized_residual", "physical_residual", "raw_residual",
                        "objective", "mu"):
                _require(event[key] is None, f"failed evaluation {key} must be null")
        else:
            for key in ("normalized_residual", "physical_residual", "raw_residual"):
                _require_number(event[key], f"evaluation {key}", 0)
            for key in ("objective", "mu", "step"):
                _require_number(event[key], f"evaluation {key}")
            _require_hash(event["geometry_hash"], "evaluation geometry hash")
            _require_hash(event["topology_hash"], "evaluation topology hash")


def _validate_trace(trace, accepted_events):
    _require(type(trace) is list, "trace must be an array")
    _require(bool(trace) == bool(accepted_events),
             "trace emptiness disagrees with accepted event count")
    _require(len(trace) == len(accepted_events), "accepted trace/event count mismatch")
    mapping = {
        "restart": "restart", "iteration": "iteration", "evaluations": "evaluation_index",
        "normalized_residual": "normalized_residual", "physical_residual": "physical_residual",
        "objective": "objective", "mu": "mu", "accepted_geodesic_step": "step",
        "quadric_hash": "quadric_hash", "geometry_hash": "geometry_hash",
        "topology_hash": "topology_hash", "profile": "profile",
        "raw_residual": "raw_residual", "quadric_coefficients": "quadric_coefficients",
    }
    for position, (item, event) in enumerate(zip(trace, accepted_events)):
        _require_exact_keys(item, TRACE_KEYS, "trace item")
        for key in ("restart", "iteration", "evaluations"):
            _require(_is_int(item[key]), f"trace {key} must be an integer")
        _require(item["restart"] == 0 and item["iteration"] >= 0
                 and item["evaluations"] >= 1, "trace index out of range")
        for key in ("normalized_residual", "physical_residual", "raw_residual"):
            _require_number(item[key], f"trace {key}", 0)
        for key in ("objective", "mu", "accepted_geodesic_step"):
            _require_number(item[key], f"trace {key}")
        for key in ("quadric_hash", "geometry_hash", "topology_hash"):
            _require_hash(item[key], f"trace {key}")
        _require(item["profile"] == PROFILE, "trace profile changed")
        _require_vector(item["quadric_coefficients"], 10, "trace quadric")
        for trace_key, event_key in mapping.items():
            _require(item[trace_key] == event[event_key],
                     f"trace item {position} does not match accepted event")


def validate_record(record, z, weight, evaluation_budget):
    _require(_is_int(evaluation_budget) and evaluation_budget in BUDGETS,
             "declared evaluation budget is not registered")
    _require_exact_keys(record, RESULT_KEYS, "result")
    _require(record["schema"] == "m1-budget-canonical-v1", "wrong result schema")
    _require(_is_int(record["evaluation_budget"])
             and record["evaluation_budget"] in BUDGETS,
             "record evaluation budget is not registered")
    _require(record["evaluation_budget"] == evaluation_budget,
             "record evaluation budget disagrees with declared budget")
    _require_number(record["z"], "z")
    _require_number(record["weight"], "weight")
    _require(record["z"] == z and record["weight"] == weight, "wrong case")
    for key in ("numerical_status", "numerical_diagnostic", "final_diagnostic"):
        _require(type(record[key]) is str, f"{key} must be a string")
    _require(record["numerical_status"] in KNOWN_STATUSES, "unknown numerical status")
    _require(record["numerical_status"] != "INFEASIBLE", "exhaustion is not infeasibility")
    _require(type(record["final_verified"]) is bool, "final_verified must be boolean")
    for key in ("raw_residual", "independent_raw_residual"):
        _require_nullable_number(record[key], key, 0)
    for key in ("physical_residual", "normalized_residual"):
        value = record[key]
        _require((_is_finite_number(value) and value >= 0) or value == "NONFINITE",
                 f"{key} must be finite numeric or the native NONFINITE failure tag")
    _require_number(record["mu"], "mu")
    _require(type(record["multipliers"]) is list and len(record["multipliers"]) in (0, 4),
             "multipliers must contain zero or four values")
    for value in record["multipliers"]:
        _require_number(value, "multiplier")
    _require_vector(record["initial_quadric"], 10, "initial quadric")
    _require_nullable_vector(record["quadric"], 10, "result quadric")
    _require_nullable_vector(record["final_moments"], 10, "final moments")
    for key in ("input_hash", "output_hash", "event_hash", "trace_hash"):
        _require_hash(record[key], key)
    _require_nullable_hash(record["final_geometry_hash"], "final geometry hash")
    count = record["evaluations"]
    _require(_is_int(count) and 1 <= count <= evaluation_budget,
             "evaluation count outside declared budget")
    _require(_is_int(record["restarts"]) and record["restarts"] == 1,
             "restart count differs from the single configured restart")
    restart_counts = record["restart_evaluations"]
    _require(type(restart_counts) is list and len(restart_counts) == 1
             and _is_int(restart_counts[0]) and restart_counts[0] == count,
             "restart evaluation accounting mismatch")
    _require(_is_int(record["event_stream_version"])
             and record["event_stream_version"] == 2, "unexpected event stream version")
    events = record["events"]
    _require(type(events) is list and len(events) == count + 2,
             "event stream length mismatch")
    for event in events:
        _validate_event(event, count)
    _require(events[0]["type"] == "attempt_start"
             and events[-1]["type"] == "attempt_terminal", "event framing is incomplete")
    evaluations = events[1:-1]
    _require(all(event["type"] == "evaluation" for event in evaluations),
             "non-evaluation event inside attempt")
    _require([event["evaluation_index"] for event in evaluations]
             == list(range(1, count + 1)), "evaluation sequence incomplete")
    accepted_events = [event for event in evaluations if event["outcome"] == "accepted"]
    _validate_trace(record["trace"], accepted_events)
    status = record["numerical_status"]
    terminal_reason = events[-1]["reason"]
    _require(status in STATUS_TERMINAL_REASONS
             and terminal_reason in STATUS_TERMINAL_REASONS[status],
             "numerical status and terminal reason disagree")
    if terminal_reason in {"initial_evaluation_failed", "initial_numerical_state_invalid"}:
        _require(count == 1 and not accepted_events,
                 "initial failure terminal has invalid event accounting")
        expected_reason = ("invalid_numerical_state"
                           if terminal_reason == "initial_numerical_state_invalid"
                           else "evaluation_exception")
        _require(evaluations[0]["outcome"] == "failed"
                 and evaluations[0]["reason"] == expected_reason,
                 "initial failure event and terminal reason disagree")
    if terminal_reason == "maximum_evaluations_reached":
        _require(count == evaluation_budget,
                 "maximum-evaluations terminal occurred below declared budget")
    if terminal_reason == "line_search_exhausted":
        _require(count < evaluation_budget,
                 "line-search terminal conflicts with exhausted evaluation budget")
    tagged_nonfinite = [record[key] == "NONFINITE"
                        for key in ("physical_residual", "normalized_residual")]
    if any(tagged_nonfinite):
        _require(all(tagged_nonfinite), "partial native NONFINITE top-level state")
        _require(status == "GEOMETRY_EVALUATION_FAILED" and not accepted_events
                 and record["raw_residual"] is None and record["mu"] == 0.0,
                 "NONFINITE tag outside the native pre-accepted-state failure shape")
    evaluation_timings = record["evaluation_timings"]
    _require(type(evaluation_timings) is list and len(evaluation_timings) == count,
             "evaluation timing count mismatch")
    for index, (item, event) in enumerate(zip(evaluation_timings, evaluations), start=1):
        _require_exact_keys(item, EVALUATION_TIMING_KEYS, "evaluation timing")
        for key in ("evaluation_index", "restart", "restart_evaluation_index"):
            _require(_is_int(item[key]), f"evaluation timing {key} must be an integer")
        _require(item["evaluation_index"] == index and item["restart"] == 0
                 and item["restart_evaluation_index"] == index,
                 "evaluation timing accounting mismatch")
        _require(type(item["returned"]) is bool, "evaluation returned flag must be boolean")
        if event["outcome"] != "failed":
            _require(item["returned"], "successful evaluation was not returned")
        if event["reason"] == "invalid_numerical_state":
            _require(item["returned"],
                     "invalid numerical state must follow a returned native evaluation")
        if not item["returned"]:
            _require(event["outcome"] == "failed"
                     and event["reason"] == "evaluation_exception"
                     and event["geometry_hash"] is None
                     and event["topology_hash"] is None,
                     "pre-return evaluator exception retained post-return material identity")
        _require_number(item["wall_seconds"], "evaluation wall time", 0)
        for key in ("reported_residual_seconds", "reported_derivative_seconds"):
            _require_nullable_number(item[key], f"evaluation {key}", 0)
        if item["returned"]:
            _require(item["reported_residual_seconds"] is not None
                     and item["reported_derivative_seconds"] is not None,
                     "returned evaluation lacks reported timing")
        else:
            _require(item["reported_residual_seconds"] is None
                     and item["reported_derivative_seconds"] is None,
                     "failed evaluation has reported timing")
    timing = record["timing"]
    _require_exact_keys(timing, TIMING_KEYS, "timing")
    for key in TIMING_KEYS:
        _require_number(timing[key], key, 0)
    parts = [timing[key] for key in
             ("coupled_evaluator_seconds", "final_audit_seconds", "solver_other_seconds")]
    _require(math.isclose(sum(parts), timing["solve_wrapper_seconds"],
                          rel_tol=1e-12, abs_tol=1e-12),
             "solve-wrapper timing double counted")
    _require(timing["entry_through_record_seconds"] >= timing["solve_wrapper_seconds"],
             "entry-through-record time is shorter than solve wrapper")
    if status != "CONVERGED":
        _require(not record["final_verified"], "failed search certified")
        _require(record["quadric"] is None and record["final_geometry_hash"] is None
                 and record["final_moments"] is None
                 and record["independent_raw_residual"] is None,
                 "failed search exposed final material")
    if record["final_verified"]:
        _require(status == "CONVERGED", "failed search certified")
        _require(record["raw_residual"] is not None and record["raw_residual"] <= 1e-12,
                 "verified raw residual exceeds 1e-12")
        _require(record["independent_raw_residual"] is not None
                 and record["independent_raw_residual"] <= 1e-12,
                 "independent verified raw residual exceeds 1e-12")
        _require(record["physical_residual"] <= 1e-6,
                 "verified physical residual exceeds 1e-6")
        _require(record["quadric"] is not None and record["final_moments"] is not None,
                 "verified final material is incomplete")
        _require(record["final_geometry_hash"] is not None,
                 "verified geometry identity missing")
        _require(len(record["multipliers"]) == 4,
                 "verified multiplier vector must contain four values")


def scientific_payload_sha256(record):
    payload = {key: value for key, value in record.items()
               if key not in SCIENTIFIC_TIMING_SIDECARS}
    serialized = json.dumps(payload, sort_keys=True, separators=(",", ":"),
                            ensure_ascii=False, allow_nan=False).encode("utf-8")
    return hashlib.sha256(serialized).hexdigest()


def process_overhead_seconds(process_wall_seconds, record):
    _require_number(process_wall_seconds, "process wall time", 0)
    solve_wrapper_seconds = record["timing"]["solve_wrapper_seconds"]
    _require_number(solve_wrapper_seconds, "solve wrapper time", 0)
    overhead = process_wall_seconds - solve_wrapper_seconds
    _require(math.isfinite(overhead), "process overhead is nonfinite")
    _require(overhead >= MATERIAL_NEGATIVE_OVERHEAD_SECONDS,
             "process overhead is materially negative")
    return overhead


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _identity(path, expected=None):
    path = Path(path).resolve(strict=True)
    digest = _sha256(path)
    if expected is not None:
        _require(digest == expected, f"frozen identity changed: {path}")
    return {"path": str(path), "sha256": digest, "bytes": path.stat().st_size}


def legacy_comparable_payload(record):
    return {key: value for key, value in record.items()
            if key not in SCIENTIFIC_TIMING_SIDECARS
            and key not in {"schema", "evaluation_budget"}}


def numeric_evaluation_prefix(record):
    return [{key: event[key] for key in PREFIX_KEYS}
            for event in record["events"] if event["type"] == "evaluation"]


def compare_numeric_prefix(left, right):
    left_prefix = numeric_evaluation_prefix(left)
    right_prefix = numeric_evaluation_prefix(right)
    common = min(len(left_prefix), len(right_prefix))
    for index in range(common):
        for key in PREFIX_KEYS:
            if left_prefix[index][key] != right_prefix[index][key]:
                return {
                    "match": False, "common_evaluations": common,
                    "first_mismatch_evaluation": index + 1, "field": key,
                    "left": left_prefix[index][key], "right": right_prefix[index][key],
                }
    return {"match": True, "common_evaluations": common,
            "left_evaluations": len(left_prefix),
            "right_evaluations": len(right_prefix)}


def _load_baseline(records_path):
    _require(_sha256(records_path) == BASELINE_RECORDS_SHA256,
             "baseline records identity changed")
    rows = [json.loads(line) for line in records_path.read_text().splitlines() if line]
    _require(len(rows) == 18, "baseline must contain exactly 18 records")
    groups = {}
    for row in rows:
        _require(type(row) is dict and type(row.get("result")) is dict,
                 "baseline row lacks a result")
        key = (row.get("z"), row.get("weight"))
        _require(key in CASES, "baseline row has an unknown case")
        groups.setdefault(key, []).append(row["result"])
    for case in CASES:
        records = groups.get(case, [])
        _require(len(records) == 3, "baseline case does not contain three runs")
        payloads = {_canonical_json(legacy_comparable_payload(record))
                    for record in records}
        _require(len(payloads) == 1,
                 "baseline case lacks one exact scientific payload")
    return groups


def _canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False)


def _source_identities(source_root):
    new_paths = (
        source_root / "tests/canonical_budget_probe.cpp",
        source_root / "experiments/m1_budget_campaign/CMakeLists.txt",
        source_root / "experiments/m1_budget_campaign/run_budget_campaign.py",
        source_root / "experiments/m1_budget_campaign/test_budget_records.py",
    )
    return {
        "budget_harness": [_identity(path) for path in new_paths],
        "prior_harness": [_identity(path, digest)
                          for path, digest in PRIOR_SOURCE_IDENTITIES.items()],
        "frozen_native": [_identity(path, digest)
                          for path, digest in FROZEN_NATIVE_IDENTITIES.items()],
        "registration": [_identity(path, digest)
                         for path, digest in REGISTRATION_IDENTITIES.items()],
    }


def _flatten_identities(identities):
    return {item["path"]: item["sha256"]
            for group in identities.values() for item in group}


def _schedule():
    schedule = []
    for budget in BUDGETS:
        for repeat in range(3):
            cases = tuple(reversed(CASES)) if repeat == 1 else CASES
            schedule.extend((budget, repeat, z, weight) for z, weight in cases)
    return schedule


def _baseline_comparison(rows, baseline_groups):
    comparisons = []
    for z, weight in CASES:
        reference = baseline_groups[(z, weight)][0]
        current = [row["result"] for row in rows
                   if row["evaluation_budget"] == 20 and row["z"] == z
                   and row["weight"] == weight and "result" in row]
        for repeat, record in enumerate(current):
            payload_match = (legacy_comparable_payload(record)
                             == legacy_comparable_payload(reference))
            native = {key: record[key] == reference[key] for key in NATIVE_HASH_KEYS}
            comparisons.append({
                "z": z, "weight": weight, "repeat_result_index": repeat,
                "raw_scientific_payload_match": payload_match,
                "native_hash_matches": native,
            })
    return comparisons


def _prefix_comparisons(rows):
    comparisons = []
    representatives = {}
    for row in rows:
        if "result" in row:
            representatives.setdefault(
                (row["evaluation_budget"], row["z"], row["weight"]),
                row["result"])
    for z, weight in CASES:
        for left_budget, right_budget in ((20, 40), (20, 60), (40, 60)):
            left = representatives.get((left_budget, z, weight))
            right = representatives.get((right_budget, z, weight))
            result = ({"match": False, "error": "missing valid representative"}
                      if left is None or right is None
                      else compare_numeric_prefix(left, right))
            comparisons.append({"z": z, "weight": weight,
                                "left_budget": left_budget,
                                "right_budget": right_budget, **result})
    return comparisons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--baseline-root", type=Path, required=True)
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    baseline_root = args.baseline_root.resolve(strict=True)
    baseline_records_path = baseline_root / "records.jsonl"
    baseline_records_path.resolve(strict=True)
    baseline_groups = _load_baseline(baseline_records_path)
    source_root = Path(__file__).resolve().parents[2]
    identities = _source_identities(source_root)
    initial_identities = _flatten_identities(identities)
    binary_identity = _identity(binary)
    baseline_identity = _identity(baseline_records_path, BASELINE_RECORDS_SHA256)
    schedule = _schedule()
    root = args.output_root.resolve()
    root.mkdir(parents=True, exist_ok=False)
    protocol = {
        "schema": "m1-budget-campaign-protocol-v1",
        "budgets": list(BUDGETS), "cases": CASES, "schedule": schedule,
        "exact_settings": {
            "maximum_function_evaluations_per_restart": list(BUDGETS),
            "maximum_iterations": 60, "maximum_restarts": 1,
            "equation_profile": PROFILE, "equation_tolerance": 1e-12,
            "physical_relative_tolerance": 1e-6,
            "line_search_reduction": 0.5, "maximum_geodesic_step": 0.5,
            "deterministic_seed": 0, "mode": "SPINNING_TOP",
            "icosphere_subdivisions": 1, "truncation_radius": 8,
            "density": 1, "minimum_wall_thickness": 0.15,
            "outer_box_radius": 1, "design_box_radius": 0.8,
            "z_values": [0, 2], "center_of_mass_weights": [0.1, 1, 10],
            "initial_quadric_seed":
                [-0.31, 0.07, -0.09, 0.11, 0.83, 1.17, -0.72, 0.06, -0.04, 0.03],
            "fixed_shell_moments":
                "outer surface moments minus design-domain surface moments",
            "z_translation": "same frozen coefficient translation as prior probe",
        },
        "binary": binary_identity, "source_identities": identities,
        "baseline_records": baseline_identity,
        "cache_policy": "fresh process, uncontrolled OS page cache",
        "timeout_seconds": TIMEOUT_SECONDS,
        "source_grade": "SYNTHETIC_REGRESSION", "formal_B0_B1": False,
        "collector_scientific_payload_sha256": {
            "serialization": "UTF-8 canonical JSON (sorted keys, compact separators)",
            "excluded_timing_sidecars": sorted(SCIENTIFIC_TIMING_SIDECARS),
            "purpose": "bind actual payload independently of probe-emitted hash claims",
        },
        "timing_partition": {
            "process_wall_seconds": "solve_wrapper_seconds + process_overhead_seconds",
            "solve_wrapper_seconds":
                "coupled_evaluator_seconds + final_audit_seconds + solver_other_seconds",
        },
    }
    (root / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    rows, errors = [], []
    for index, (budget, repeat, z, weight) in enumerate(schedule):
        cmd = [str(binary), str(z), str(weight), str(budget)]
        stdout_path = root / f"{index:02d}.stdout"
        stderr_path = root / f"{index:02d}.stderr"
        start = time.monotonic()
        with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
            proc = subprocess.Popen(cmd, stdout=stdout, stderr=stderr)
            timed_out = False
            while True:
                pid, status, usage = os.wait4(proc.pid, os.WNOHANG)
                if pid:
                    break
                if time.monotonic() - start > TIMEOUT_SECONDS:
                    os.kill(proc.pid, signal.SIGKILL)
                    pid, status, usage = os.wait4(proc.pid, 0)
                    timed_out = True
                    break
                time.sleep(0.02)
            proc.returncode = os.waitstatus_to_exitcode(status)
        row = {
            "index": index, "repeat": repeat, "z": z, "weight": weight,
            "evaluation_budget": budget, "command": cmd,
            "exit_code": proc.returncode, "timed_out": timed_out,
            "process_wall_seconds": time.monotonic() - start,
            "max_rss_kib": usage.ru_maxrss, "user_seconds": usage.ru_utime,
            "system_seconds": usage.ru_stime, "stdout": stdout_path.name,
            "stderr": stderr_path.name,
        }
        try:
            if timed_out:
                raise ValueError("resource timeout; raw diagnostics retained")
            if proc.returncode != 0:
                raise ValueError("probe exited unsuccessfully; raw diagnostics retained")
            lines = [line.removeprefix("M1_RESULT ")
                     for line in stdout_path.read_text().splitlines()
                     if line.startswith("M1_RESULT ")]
            if len(lines) != 1:
                raise ValueError("missing/duplicate result record")
            value = json.loads(lines[0])
            validate_record(value, z, weight, budget)
            row["process_overhead_seconds"] = process_overhead_seconds(
                row["process_wall_seconds"], value)
            row["collector_scientific_payload_sha256"] = scientific_payload_sha256(value)
            row["result"] = value
        except (ValueError, KeyError, TypeError, json.JSONDecodeError) as error:
            row["record_error"] = str(error)
            errors.append({"index": index, "error": str(error)})
        rows.append(row)
        with (root / "records.jsonl").open("a") as ledger:
            ledger.write(json.dumps(row, allow_nan=False) + "\n")
        print(f"{index + 1}/54 budget={budget} z={z} weight={weight} "
              f"exit={proc.returncode}", flush=True)
    if _sha256(binary) != binary_identity["sha256"]:
        errors.append({"error": "binary changed during campaign"})
    for path, digest in initial_identities.items():
        if _sha256(Path(path)) != digest:
            errors.append({"path": path, "error": "source changed during campaign"})
    summary = []
    for budget in BUDGETS:
        for z, weight in CASES:
            group = [row for row in rows if row["evaluation_budget"] == budget
                     and row["z"] == z and row["weight"] == weight
                     and "result" in row]
            payload_digests = {row["collector_scientific_payload_sha256"] for row in group}
            if len(group) != 3 or len(payload_digests) != 1:
                errors.append({
                    "evaluation_budget": budget, "z": z, "weight": weight,
                    "error": "incomplete or nondeterministic collector-bound scientific payload replay",
                })
            summary.append({
                "evaluation_budget": budget, "z": z, "weight": weight,
                "valid_runs": len(group),
                "unique_collector_scientific_payload_sha256": len(payload_digests),
                "collector_scientific_payload_sha256": sorted(payload_digests),
                "verified_runs": sum(row["result"]["final_verified"] for row in group),
                "statuses": sorted({row["result"]["numerical_status"] for row in group}),
                "median_process_seconds": statistics.median(
                    row["process_wall_seconds"] for row in group) if group else None,
                "peak_rss_kib": max((row["max_rss_kib"] for row in group), default=None),
            })
    baseline_comparisons = _baseline_comparison(rows, baseline_groups)
    for comparison in baseline_comparisons:
        if (not comparison["raw_scientific_payload_match"]
                or not all(comparison["native_hash_matches"].values())):
            errors.append({**comparison, "error": "budget20 differs from frozen baseline"})
    prefix_comparisons = _prefix_comparisons(rows)
    for comparison in prefix_comparisons:
        if not comparison["match"]:
            errors.append({**comparison, "error": "cross-budget numeric prefix mismatch"})
    report = {
        "status": "PASS" if not errors else "READINESS_FAILED", "errors": errors,
        "case_summary": summary, "budget20_baseline": baseline_comparisons,
        "cross_budget_numeric_prefixes": prefix_comparisons,
        "records": len(rows), "distinct_geometries": 2,
        "formal_B0_B1": False, "not_author_cases": True,
    }
    (root / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
