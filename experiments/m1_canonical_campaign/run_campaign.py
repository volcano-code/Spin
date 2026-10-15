"""Frozen six-case development diagnostic, no training and no B0/B1 claim."""
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


RESULT_KEYS = {
    "schema", "z", "weight", "numerical_status", "numerical_diagnostic",
    "final_verified", "final_diagnostic", "raw_residual",
    "independent_raw_residual", "physical_residual", "normalized_residual",
    "mu", "multipliers", "initial_quadric", "quadric", "input_hash",
    "output_hash", "trace_hash", "event_hash", "final_geometry_hash",
    "final_moments", "evaluations", "restarts", "restart_evaluations",
    "event_stream_version", "timing", "events", "trace",
    "evaluation_timings",
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
    "initial_sensitivity_vanished",
    "line_search_exhausted", "maximum_evaluations_reached",
    "maximum_iterations_reached", "newton_direction_invalid",
    "newton_direction_vanished", "newton_system_invalid",
    "bordered_system_rank_deficient",
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
MATERIAL_NEGATIVE_OVERHEAD_SECONDS = -1e-6


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
    event_key_by_trace_key = {
        "restart": "restart", "iteration": "iteration", "evaluations": "evaluation_index",
        "normalized_residual": "normalized_residual",
        "physical_residual": "physical_residual", "objective": "objective", "mu": "mu",
        "accepted_geodesic_step": "step", "quadric_hash": "quadric_hash",
        "geometry_hash": "geometry_hash", "topology_hash": "topology_hash",
        "profile": "profile", "raw_residual": "raw_residual",
        "quadric_coefficients": "quadric_coefficients",
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
        for trace_key, event_key in event_key_by_trace_key.items():
            _require(item[trace_key] == event[event_key],
                     f"trace item {position} does not match accepted event")


def validate_record(record, z, weight):
    _require_exact_keys(record, RESULT_KEYS, "result")
    _require(record["schema"] == "m1-six-case-canonical-v1", "wrong result schema")
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
    _require(_is_int(count) and 1 <= count <= 20, "evaluation count outside configured budget")
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
    _require([event["evaluation_index"] for event in evaluations] == list(range(1, count + 1)),
             "evaluation sequence incomplete")
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
        _require(count == 20, "maximum-evaluations terminal occurred below the configured budget")
    if terminal_reason == "line_search_exhausted":
        _require(count < 20, "line-search terminal conflicts with exhausted evaluation budget")

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
                          rel_tol=1e-12, abs_tol=1e-12), "solve-wrapper timing double counted")
    _require(timing["entry_through_record_seconds"] >= timing["solve_wrapper_seconds"],
             "entry-through-record time is shorter than solve wrapper")

    if record["numerical_status"] != "CONVERGED":
        _require(not record["final_verified"], "failed search certified")
        _require(record["quadric"] is None and record["final_geometry_hash"] is None
                 and record["final_moments"] is None
                 and record["independent_raw_residual"] is None,
                 "failed search exposed final material")
    if record["final_verified"]:
        _require(record["numerical_status"] == "CONVERGED", "failed search certified")
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
    """Collector digest of the complete result payload with timing sidecars removed."""
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    binary_sha = hashlib.sha256(binary.read_bytes()).hexdigest()
    root = args.output_root.resolve()
    root.mkdir(parents=True, exist_ok=False)  # Never overwrite a completed/partial campaign.
    cases = [(z, w) for z in (0, 2) for w in (.1, 1, 10)]
    schedule = [(repeat, z, w) for repeat in range(3)
                for z, w in (list(reversed(cases)) if repeat == 1 else cases)]
    protocol = {"cases": cases, "schedule": schedule, "binary": str(binary),
                "binary_sha256": binary_sha, "cache_policy": "fresh process, uncontrolled OS page cache",
                "timeout_seconds": 60, "source_grade": "SYNTHETIC_REGRESSION", "formal_B0_B1": False,
                "collector_scientific_payload_sha256": {
                    "serialization": "UTF-8 canonical JSON (sorted keys, compact separators)",
                    "excluded_timing_sidecars": sorted(SCIENTIFIC_TIMING_SIDECARS),
                    "purpose": "bind actual payload independently of probe-emitted hash claims",
                },
                "timing_partition": {
                    "process_wall_seconds": "solve_wrapper_seconds + process_overhead_seconds",
                    "solve_wrapper_seconds": "coupled_evaluator_seconds + final_audit_seconds + solver_other_seconds",
                }}
    (root / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    rows, errors = [], []
    for index, (repeat, z, weight) in enumerate(schedule):
        cmd = [str(binary), str(z), str(weight)]
        stdout_path, stderr_path = root / f"{index:02d}.stdout", root / f"{index:02d}.stderr"
        start = time.monotonic()
        with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
            proc = subprocess.Popen(cmd, stdout=stdout, stderr=stderr)
            timed_out = False
            while True:
                pid, status, usage = os.wait4(proc.pid, os.WNOHANG)
                if pid:
                    break
                if time.monotonic() - start > 60:
                    os.kill(proc.pid, signal.SIGKILL)
                    pid, status, usage = os.wait4(proc.pid, 0)
                    timed_out = True
                    break
                time.sleep(.02)
            proc.returncode = os.waitstatus_to_exitcode(status)
        row = {"index": index, "repeat": repeat, "z": z, "weight": weight, "command": cmd,
               "exit_code": proc.returncode, "timed_out": timed_out,
               "process_wall_seconds": time.monotonic() - start,
               "max_rss_kib": usage.ru_maxrss, "user_seconds": usage.ru_utime,
               "system_seconds": usage.ru_stime, "stdout": stdout_path.name, "stderr": stderr_path.name}
        try:
            if proc.returncode != 0:
                raise ValueError("probe exited unsuccessfully; raw diagnostics retained")
            lines = [line[len("M1_RESULT "):] for line in stdout_path.read_text().splitlines()
                     if line.startswith("M1_RESULT ")]
            if len(lines) != 1:
                raise ValueError("missing/duplicate result record")
            value = json.loads(lines[0])
            validate_record(value, z, weight)
            row["process_overhead_seconds"] = process_overhead_seconds(
                row["process_wall_seconds"], value)
            row["collector_scientific_payload_sha256"] = scientific_payload_sha256(value)
            row["result"] = value
        except (ValueError, KeyError, TypeError) as error:
            row["record_error"] = str(error)
            errors.append({"index": index, "error": str(error)})
        rows.append(row)
        with (root / "records.jsonl").open("a") as ledger:
            ledger.write(json.dumps(row, allow_nan=False) + "\n")
        print(f"{index+1}/18 z={z} weight={weight} exit={proc.returncode}", flush=True)
    if hashlib.sha256(binary.read_bytes()).hexdigest() != binary_sha:
        errors.append({"error": "binary changed during campaign"})
    summary = []
    for z, weight in cases:
        group = [r for r in rows if r["z"] == z and r["weight"] == weight and "result" in r]
        payload_digests = {r["collector_scientific_payload_sha256"] for r in group}
        if len(group) != 3 or len(payload_digests) != 1:
            errors.append({"z": z, "weight": weight,
                           "error": "incomplete or nondeterministic collector-bound scientific payload replay"})
        summary.append({"z": z, "weight": weight, "valid_runs": len(group),
                        "unique_collector_scientific_payload_sha256": len(payload_digests),
                        "collector_scientific_payload_sha256": sorted(payload_digests),
                        "verified_runs": sum(r["result"]["final_verified"] for r in group),
                        "statuses": sorted({r["result"]["numerical_status"] for r in group}),
                        "median_process_seconds": statistics.median(r["process_wall_seconds"] for r in group) if group else None,
                        "peak_rss_kib": max((r["max_rss_kib"] for r in group), default=None)})
    report = {"status": "PASS" if not errors else "READINESS_FAILED", "errors": errors,
              "case_summary": summary, "records": len(rows), "distinct_geometries": 2,
              "formal_B0_B1": False, "not_author_cases": True}
    (root / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
