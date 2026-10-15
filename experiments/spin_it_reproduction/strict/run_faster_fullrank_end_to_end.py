#!/usr/bin/env python3
"""Capture the native Faster full-rank integration chain as method evidence.

This runner exercises the current synthetic regression modes: adapter
construction, derivative/Boolean chain, solver guards, independent moment
verification, Boolean topology-transition relinearization, and process replay. It
does not claim an author asset, a CAD result, a global theorem, or a paper
Table-1 case.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import math
from pathlib import Path
import re
import subprocess
import time
import traceback

from make_it_stand.verification.run_evidence import (
    EXPLICIT_EFFECTIVE_SOURCE_POLICY,
    begin_run,
    check_freshness,
    finish_run,
)
from make_it_stand.verification.spin_common import canonical, parse


ROOT = Path(__file__).resolve().parents[3]
_BINARY = "spin-faster-canonical-fullrank-regression"
_SEED_SWEEP_BINARY = "spin-faster-unknown-seed-sweep"
_MAXIMUM_RESTARTS_PER_SEED = 2
_MAXIMUM_EVALUATIONS_PER_RESTART = 60
_MODES = (
    ("adapter", "[PASS] adapter"),
    ("chain", "[PASS] chain"),
    ("solver", "[PASS] solver"),
    ("integration", "[PASS] integration"),
    ("offtarget-yoyo", "[PASS] offtarget-yoyo"),
    ("audit", "[PASS] audit"),
    ("replay", "[PASS] replay"),
    ("negative", "[PASS] negative"),
)
_UNKNOWN_SEED_RECORD = re.compile(
    r"^seed=(?P<seed>\d+) generated_start=(?P<generated_start>[0-9a-f]{64}) "
    r"status=(?P<status>[A-Z_]+) evaluations=(?P<evaluations>\d+) "
    r"trace=(?P<trace>\d+) restarts=(?P<restarts>\d+) "
    r"restart_evaluations=(?P<restart_evaluations>\d+(?:,\d+)*) "
    r"residual=(?P<residual>[^ ]+) final_verified=(?P<final_verified>[01]) "
    r"input_hash=(?P<input_hash>[0-9a-f]{64}) "
    r"event_hash=(?P<event_hash>[0-9a-f]{64}) "
    r"trace_hash=(?P<trace_hash>[0-9a-f]{64}) "
    r"output_hash=(?P<output_hash>[0-9a-f]{64}) "
    r"initial_subtype=(?P<initial_subtype>[A-Z_]+) "
    r"initial_source_motion_defined=(?P<initial_source_motion_defined>[01]) "
    r"initial_source_diagnostic_hex=(?P<initial_source_diagnostic_hex>[0-9a-f]*) "
    r"first_derivative_failure_evaluation=(?P<first_derivative_failure_evaluation>-?\d+) "
    r"first_derivative_failure_iteration=(?P<first_derivative_failure_iteration>-?\d+) "
    r"first_derivative_failure_subtype=(?P<first_derivative_failure_subtype>[A-Z_]+) "
    r"first_derivative_failure_step=(?P<first_derivative_failure_step>[^ ]+) "
    r"first_derivative_failure_hex=(?P<first_derivative_failure_hex>[0-9a-f]*) "
    r"first_derivative_failure_q=(?P<first_derivative_failure_q>[0-9eE+.,-]+|none) "
    r"source_witness_vertex=(?P<source_witness_vertex>-?\d+) "
    r"source_witness_nearest_normal=(?P<source_witness_nearest_normal>-?\d+) "
    r"source_witness_nearest_distance=(?P<source_witness_nearest_distance>[^ ]+) "
    r"source_witness_radius_delta=(?P<source_witness_radius_delta>[^ ]+) "
    r"source_witness_q_residual=(?P<source_witness_q_residual>[^ ]+) "
    r"source_witness_point=(?P<source_witness_point>[^ ]+) "
    r"solver_diagnostic_hex=(?P<solver_diagnostic_hex>[0-9a-f]*) "
    r"verification_diagnostic_hex=(?P<verification_diagnostic_hex>[0-9a-f]*)$"
)
_SOLVER_STATUSES = {
    "CONVERGED",
    "INFEASIBLE",
    "INVALID_REQUEST",
    "OPTIMIZATION_DID_NOT_CONVERGE",
    "LINE_SEARCH_EXHAUSTED",
    "GEOMETRY_EVALUATION_FAILED",
    "INVALID_MATERIAL_ARTIFACT",
    "RESOURCE_LIMIT_REACHED",
    "INDETERMINATE",
}
_INITIAL_START = re.compile(
    r"^initial-start seed=(?P<seed>\d+) restart=(?P<restart>\d+) "
    r"hash=(?P<hash>[0-9a-f]{64}) coefficients=(?P<coefficients>[0-9eE+.,-]+)$"
)
_SOLVER_EVENT = re.compile(
    r"^solver-event seed=(?P<seed>\d+) restart=(?P<restart>\d+) "
    r"iteration=(?P<iteration>-?\d+) evaluation=(?P<evaluation>-?\d+) "
    r"type=(?P<type>\S+) outcome=(?P<outcome>\S+) reason=(?P<reason>\S+) "
    r"detail_hex=(?P<detail_hex>[0-9a-f]*) "
    r"coefficients=(?P<coefficients>none|[0-9eE+.,-]+)$"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write(path: Path, value: object) -> None:
    with path.open("x", encoding="ascii") as stream:
        stream.write(canonical(value) + "\n")


def parse_unknown_seed_records(
    transcript: str,
    *,
    expected_seeds: tuple[int, ...] = (0, 1, 2, 3, 4),
    require_completion: bool = True,
) -> list[dict[str, object]]:
    """Parse all five frozen-seed solver outcomes without dropping failures."""
    if require_completion and "[PASS] unknown-seed-sweep" not in transcript:
        raise ValueError("unknown-answer seed sweep lacks its completion marker")
    rows: list[dict[str, object]] = []
    for line in transcript.splitlines():
        match = _UNKNOWN_SEED_RECORD.fullmatch(line)
        if match is None:
            if line.startswith("seed="):
                raise ValueError("malformed unknown-answer seed record")
            continue
        values = match.groupdict()
        row: dict[str, object] = {
            "seed": int(values["seed"]),
            "generated_start": values["generated_start"],
            "status": values["status"],
            "evaluations": int(values["evaluations"]),
            "trace": int(values["trace"]),
            "restarts": int(values["restarts"]),
            "restart_evaluations": values["restart_evaluations"],
            "residual": values["residual"],
            "final_verified": values["final_verified"] == "1",
            "input_hash": values["input_hash"],
            "event_hash": values["event_hash"],
            "trace_hash": values["trace_hash"],
            "output_hash": values["output_hash"],
            "initial_subtype": values["initial_subtype"],
            "initial_source_motion_defined": values[
                "initial_source_motion_defined"
            ] == "1",
            "initial_source_diagnostic": bytes.fromhex(
                values["initial_source_diagnostic_hex"]
            ).decode("utf-8"),
            "first_derivative_failure_evaluation": int(
                values["first_derivative_failure_evaluation"]
            ),
            "first_derivative_failure_iteration": int(
                values["first_derivative_failure_iteration"]
            ),
            "first_derivative_failure_subtype": values[
                "first_derivative_failure_subtype"
            ],
            "first_derivative_failure_step": values[
                "first_derivative_failure_step"
            ],
            "first_derivative_failure": bytes.fromhex(
                values["first_derivative_failure_hex"]
            ).decode("utf-8"),
            "first_derivative_failure_q": values["first_derivative_failure_q"],
            "source_witness_vertex": int(values["source_witness_vertex"]),
            "source_witness_nearest_normal": int(
                values["source_witness_nearest_normal"]
            ),
            "source_witness_nearest_distance": values[
                "source_witness_nearest_distance"
            ],
            "source_witness_radius_delta": values[
                "source_witness_radius_delta"
            ],
            "source_witness_q_residual": values["source_witness_q_residual"],
            "source_witness_point": values["source_witness_point"],
            "solver_diagnostic": bytes.fromhex(
                values["solver_diagnostic_hex"]
            ).decode("utf-8"),
            "verification_diagnostic": bytes.fromhex(
                values["verification_diagnostic_hex"]
            ).decode("utf-8"),
        }
        if row["status"] not in _SOLVER_STATUSES:
            raise ValueError("unknown solver status in seed record")
        restarts = int(row["restarts"])
        evaluations = int(row["evaluations"])
        restart_evaluations = [
            int(value) for value in str(row["restart_evaluations"]).split(",")
        ]
        if restarts < 1 or restarts > _MAXIMUM_RESTARTS_PER_SEED:
            raise ValueError("restart count exceeds frozen limit")
        if (len(restart_evaluations) != restarts
                or any(value > _MAXIMUM_EVALUATIONS_PER_RESTART
                       for value in restart_evaluations)
                or sum(restart_evaluations) != evaluations):
            raise ValueError("function evaluation count exceeds frozen per-restart budget")
        row["restart_evaluations"] = restart_evaluations
        rows.append(row)
    seeds = [int(row["seed"]) for row in rows]
    if len(rows) != len(expected_seeds) or sorted(seeds) != list(expected_seeds):
        raise ValueError("seed records do not match the frozen seed selection")
    if len({row["generated_start"] for row in rows}) != len(expected_seeds):
        raise ValueError("registered seed protocol did not produce distinct starts")
    if len({row["input_hash"] for row in rows}) != len(expected_seeds):
        raise ValueError("seed-specific request identities are not distinct")
    return sorted(rows, key=lambda row: int(row["seed"]))


def parse_unknown_seed_diagnostics(
    transcript: str,
    *,
    expected_seeds: tuple[int, ...] = (0, 1, 2, 3, 4),
    require_completion: bool = True,
) -> dict[str, list[dict[str, object]]]:
    """Parse frozen starts and every emitted attempt/failure event."""
    starts: list[dict[str, object]] = []
    events: list[dict[str, object]] = []
    for line in transcript.splitlines():
        if line.startswith("initial-start=") or line.startswith("initial-start "):
            match = _INITIAL_START.fullmatch(line)
            if match is None:
                raise ValueError("malformed deterministic initial-start record")
            values = match.groupdict()
            coefficients = [float(value) for value in values["coefficients"].split(",")]
            if len(coefficients) != 10 or not all(map(math.isfinite, coefficients)):
                raise ValueError("initial start must contain ten finite coefficients")
            starts.append({
                "seed": int(values["seed"]),
                "restart": int(values["restart"]),
                "hash": values["hash"],
                "coefficients": values["coefficients"],
            })
        elif line.startswith("solver-event=") or line.startswith("solver-event "):
            match = _SOLVER_EVENT.fullmatch(line)
            if match is None:
                raise ValueError("malformed solver event record")
            values = match.groupdict()
            coefficients = values["coefficients"]
            if coefficients != "none":
                parsed = [float(value) for value in coefficients.split(",")]
                if len(parsed) != 10 or not all(map(math.isfinite, parsed)):
                    raise ValueError("solver event coefficients must be ten finite values")
            events.append({
                "seed": int(values["seed"]),
                "restart": int(values["restart"]),
                "iteration": int(values["iteration"]),
                "evaluation": int(values["evaluation"]),
                "type": values["type"],
                "outcome": values["outcome"],
                "reason": values["reason"],
                "detail_hex": values["detail_hex"],
                "coefficients": coefficients,
            })
    rows = parse_unknown_seed_records(
        transcript, expected_seeds=expected_seeds, require_completion=require_completion
    )
    expected_starts = {(seed, restart) for seed in expected_seeds
                       for restart in range(_MAXIMUM_RESTARTS_PER_SEED)}
    actual_starts = {(int(row["seed"]), int(row["restart"])) for row in starts}
    if len(starts) != len(expected_starts) or actual_starts != expected_starts:
        raise ValueError("each seed must have exactly two starts")
    rows_by_seed = {int(row["seed"]): row for row in rows}
    for event in events:
        seed = int(event["seed"])
        restart = int(event["restart"])
        if seed not in expected_seeds:
            raise ValueError("solver event references an unexpected seed")
        if restart < 0 or restart >= _MAXIMUM_RESTARTS_PER_SEED:
            raise ValueError("solver event restart exceeds frozen limit")
        if int(event["evaluation"]) < 0:
            raise ValueError("solver event evaluation index must be nonnegative")

    for seed in expected_seeds:
        seed_starts = [row for row in starts if row["seed"] == seed]
        if len({row["hash"] for row in seed_starts}) != 2:
            raise ValueError("restart starts must have distinct identities per seed")
        seed_events = [row for row in events if row["seed"] == seed]
        summary = rows_by_seed[seed]
        restart_count = int(summary["restarts"])
        restart_evaluations = [int(value) for value in summary["restart_evaluations"]]
        starts_by_restart: dict[int, list[dict[str, object]]] = {}
        terminals_by_restart: dict[int, list[dict[str, object]]] = {}
        for event in seed_events:
            restart = int(event["restart"])
            if event["type"] == "attempt_start":
                starts_by_restart.setdefault(restart, []).append(event)
            elif event["type"] == "attempt_terminal":
                terminals_by_restart.setdefault(restart, []).append(event)
        started_restarts = sorted(starts_by_restart)
        if started_restarts != list(range(restart_count)):
            raise ValueError("attempt-start events do not match reported restart count")
        if any(len(attempt_starts) != 1 for attempt_starts in starts_by_restart.values()):
            raise ValueError("each attempted restart must contain exactly one start event")
        for restart, start_events in starts_by_restart.items():
            frozen = next(row for row in seed_starts if row["restart"] == restart)
            event = start_events[0]
            if (event["coefficients"] != frozen["coefficients"]
                    or event["outcome"] != "started"):
                raise ValueError("attempt-start event does not match frozen start")
            expected_start_evaluation = sum(restart_evaluations[:restart])
            if int(event["evaluation"]) != expected_start_evaluation:
                raise ValueError("attempt-start evaluation index disagrees with prior budgets")

        if any(restart not in starts_by_restart for restart in terminals_by_restart):
            raise ValueError("attempt-terminal event has no matching start")
        if any(len(terminal_events) != 1
               for terminal_events in terminals_by_restart.values()):
            raise ValueError("each attempted restart must contain exactly one terminal event")
        if require_completion and set(terminals_by_restart) != set(range(restart_count)):
            raise ValueError("each complete attempt must contain exactly one terminal event")
        for restart, terminal_events in terminals_by_restart.items():
            event = terminal_events[0]
            expected_terminal_evaluation = sum(restart_evaluations[:restart + 1])
            if int(event["evaluation"]) != expected_terminal_evaluation:
                raise ValueError("attempt-terminal evaluation index disagrees with restart budget")
            if event["outcome"] != "terminated" or event["coefficients"] != "none":
                raise ValueError("attempt-terminal event has an invalid outcome or coefficient payload")

        # The text protocol emits each attempt's start, selected failures, then
        # its terminal record. Preserve that causal order and ensure no event
        # claims evaluations outside the attempt's reported budget interval.
        for restart in started_restarts:
            attempt_events = [event for event in seed_events
                              if int(event["restart"]) == restart]
            indexes = [int(event["evaluation"]) for event in attempt_events]
            if indexes != sorted(indexes):
                raise ValueError("solver events are not ordered by evaluation index")
            low = sum(restart_evaluations[:restart])
            high = sum(restart_evaluations[:restart + 1])
            if any(index < low or index > high for index in indexes):
                raise ValueError("solver event evaluation index exceeds its restart budget")
        if not seed_events:
            raise ValueError("solver attempt and failure events are missing")
    return {"initial_starts": starts, "solver_events": events}


def require_unknown_seed_replay(first: str, second: str) -> list[dict[str, object]]:
    """Require exact scientific record agreement between independent processes."""
    first_rows = parse_unknown_seed_records(first)
    second_rows = parse_unknown_seed_records(second)
    if first_rows != second_rows:
        raise ValueError("independent unknown-answer seed sweeps differ scientifically")
    if parse_unknown_seed_diagnostics(first) != parse_unknown_seed_diagnostics(second):
        raise ValueError("independent seed sweeps differ in starts or attempt events")
    return first_rows


def run(
    output_directory: Path,
    *,
    build_directory: Path,
    dependency_roots: dict[str, Path],
) -> dict[str, object]:
    output = output_directory.absolute()
    output.mkdir(parents=False, exist_ok=False)
    build = build_directory.resolve(strict=True)
    executable = (build / _BINARY).resolve(strict=True)
    if not executable.is_file() or not executable.stat().st_mode & 0o111:
        raise ValueError("required native full-rank regression executable is unavailable")
    seed_sweep_executable = (build / _SEED_SWEEP_BINARY).resolve(strict=True)
    if not seed_sweep_executable.is_file() or not seed_sweep_executable.stat().st_mode & 0o111:
        raise ValueError("required native unknown-answer seed-sweep executable is unavailable")

    started = datetime.now(timezone.utc).isoformat()
    tick = time.perf_counter()
    begin = None
    events: list[dict[str, object]] = []
    try:
        inputs = parse(canonical({
            "schema": "faster-fullrank-end-to-end-inputs-v2",
            "classification": "SYNTHETIC_REGRESSION",
            "binary": {"path": str(executable), "sha256": _sha256(executable)},
            "unknown_seed_sweep_binary": {
                "path": str(seed_sweep_executable),
                "sha256": _sha256(seed_sweep_executable),
            },
            "modes": [
                {"name": name, "required_marker": marker}
                for name, marker in _MODES
            ],
            "unknown_answer_initialization": {
                "mode": "unknown-seed-sweep",
                "seeds": [0, 1, 2, 3, 4],
                "starts": "deterministicInitialQuadrics from frozen request",
                "explicit_initial_guesses": False,
                "independent_process_repeats": 2,
                "maximum_restarts_per_seed": _MAXIMUM_RESTARTS_PER_SEED,
                "maximum_evaluations_per_restart": _MAXIMUM_EVALUATIONS_PER_RESTART,
                "maximum_total_evaluations_per_seed": (
                    _MAXIMUM_RESTARTS_PER_SEED * _MAXIMUM_EVALUATIONS_PER_RESTART
                ),
                "protocol": "two-start-development-method-check-v1",
                "evidence_mode": "development-method-check; not the five-seed repeated batch",
                "paper_case": False,
            },
            "strict_profile": {
                "effective_value_derivative": "discrete-fullrank-quadric-boolean-transport-v2",
                "value_and_newton_map": "same native Boolean material moments",
                "legacy_generic_solver": "kept separate for compatibility; not used by this runner",
            },
            "acceptance": {
                "positive_chain": "native adapter, Boolean/moment chain, solver guards, independent audit, and replay all return zero",
                "topology_transition": "off-target yo-yo accepts a residual-decreasing topology change only after fresh Boolean/moment evaluation; the native regression checks relinearization and final geometry verification",
            },
        }))
        input_path = output / "inputs.json"
        _write(input_path, inputs)
        begin = begin_run(
            repository_root=ROOT,
            runner_files=[Path(__file__)],
            checker_files=[
                ROOT / "tests/full_rank_hollowing_regression.cpp",
                ROOT / "src/spin/FullRankHollowingSolver.cpp",
                ROOT / "src/spin/FullRankMomentEvaluator.cpp",
                ROOT / "src/spin/PolyhedralMomentDifferential.cpp",
                ROOT / "src/spin/QuadricSourceMotion.cpp",
                ROOT / "src/spin/CanonicalFullRankHollowing.cpp",
                ROOT / "include/spin/FullRankHollowingSolver.h",
                ROOT / "include/spin/FullRankMomentEvaluator.h",
                ROOT / "include/spin/PolyhedralMomentDifferential.h",
                ROOT / "tests/full_rank_unknown_seed_sweep.cpp",
                ROOT / "tests/python/test_faster_unknown_seed_evidence.py",
                ROOT / "CMakeLists.txt",
            ],
            asset_files=[input_path, executable, seed_sweep_executable],
            build_configuration_files=[ROOT / "CMakeLists.txt", build / "CMakeCache.txt"],
            parameters=inputs,
            seed_protocol={"seed": 0, "protocol": "deterministic native full-rank mode sequence"},
            dependency_roots=dependency_roots,
            native_modules=[],
            pure_python=True,
            source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY,
        )
        _write(output / "begin.json", begin.as_dict())

        commands: list[dict[str, object]] = []
        for name, marker in _MODES:
            completed = subprocess.run(
                [str(executable), name],
                cwd=build,
                capture_output=True,
                text=True,
                check=False,
            )
            transcript = completed.stdout + completed.stderr
            transcript_path = output / f"fullrank-{name}.stdout"
            transcript_path.write_text(transcript, encoding="utf-8")
            if completed.returncode != 0 or marker not in transcript:
                raise RuntimeError(f"native full-rank mode failed: {name}")
            if name == "offtarget-yoyo" and (
                "status=CONVERGED" not in transcript
                or "final_verified=1" not in transcript
            ):
                raise RuntimeError(
                    "topology-transition regression did not converge with a verified final geometry"
                )
            commands.append({
                "name": name,
                "returncode": completed.returncode,
                "path": str(transcript_path),
                "sha256": _sha256(transcript_path),
                "required_marker": marker,
            })

        seed_sweep_transcripts: list[str] = []
        seed_sweep_records: list[dict[str, object]] | None = None
        seed_sweep_diagnostics: dict[str, list[dict[str, object]]] | None = None
        for repeat in range(2):
            name = f"unknown-seed-sweep-process-{repeat + 1}"
            completed = subprocess.run(
                [str(seed_sweep_executable)],
                cwd=build,
                capture_output=True,
                text=True,
                check=False,
            )
            transcript = completed.stdout + completed.stderr
            transcript_path = output / f"{name}.stdout"
            transcript_path.write_text(transcript, encoding="utf-8")
            if completed.returncode != 0:
                raise RuntimeError(f"native unknown-answer seed sweep failed: {name}")
            records = parse_unknown_seed_records(transcript)
            diagnostics = parse_unknown_seed_diagnostics(transcript)
            if seed_sweep_records is None:
                seed_sweep_records = records
                seed_sweep_diagnostics = diagnostics
            seed_sweep_transcripts.append(transcript)
            commands.append({
                "name": name,
                "returncode": completed.returncode,
                "path": str(transcript_path),
                "sha256": _sha256(transcript_path),
                "required_marker": "[PASS] unknown-seed-sweep",
            })
        assert seed_sweep_records is not None
        assert seed_sweep_diagnostics is not None
        require_unknown_seed_replay(
            seed_sweep_transcripts[0], seed_sweep_transcripts[1]
        )

        scientific_identity = {
            "commands": [
                {
                    "name": row["name"],
                    "returncode": row["returncode"],
                    "sha256": row["sha256"],
                    "required_marker": row["required_marker"],
                }
                for row in commands
            ],
            "topology_transition": {
                "mode": "offtarget-yoyo",
                "expected_status": "CONVERGED",
                "required_solver_event": "topology_transition_relinearized",
                "final_geometry_verified": True,
                "not_infeasibility_certificate": True,
            },
            "unknown_answer_seed_sweep": seed_sweep_records,
            "unknown_answer_seed_diagnostics": seed_sweep_diagnostics,
        }
        scientific_digest = hashlib.sha256(
            canonical(scientific_identity).encode("utf-8")
        ).hexdigest()
        scientific = {
            "schema": "faster-fullrank-end-to-end-scientific-v2",
            "classification": "SYNTHETIC_REGRESSION",
            "scope": "SYNTHETIC_FULLRANK_BOOLEAN_MOMENT_SOLVER_AND_REPLAY_METHODS_ONLY",
            "commands": commands,
            "scientific_digest": scientific_digest,
            "topology_transition": {
                "mode": "offtarget-yoyo",
                "expected_status": "CONVERGED",
                "required_solver_event": "topology_transition_relinearized",
                "final_geometry_verified": True,
                "not_infeasibility_certificate": True,
            },
            "unknown_answer_seed_sweep": seed_sweep_records,
            "unknown_answer_seed_diagnostics": seed_sweep_diagnostics,
            "seed_sweep_replay": "two independent process transcripts agree exactly on all seed records; failure statuses retained",
        }
        scientific_path = output / "scientific.json"
        _write(scientific_path, scientific)
        events.append({
            "iteration": 0,
            "status": "ok",
            "data": {
                "kind": "FASTER_FULLRANK_END_TO_END_NATIVE_MODES",
                "command_hashes": {row["name"]: row["sha256"] for row in commands},
            },
        })
        evidence = finish_run(
            begin,
            events=events,
            output_files=sorted(
                path for path in output.rglob("*")
                if path.is_file() and path.name != "begin.json"
            ),
            measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick},
        )
        _write(output / "evidence.json", evidence.as_dict())
        freshness = check_freshness(evidence)
        summary = {
            "status": "SYNTHETIC_FULLRANK_METHODS_AND_UNKNOWN_SEED_OUTCOMES_RECORDED",
            "classification": "SYNTHETIC_REGRESSION",
            "scientific_digest": scientific["scientific_digest"],
            "project_status": "PARTIAL_BLOCKED",
            "original_evidence_admission": "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS",
            "paper_case_admission": "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS",
            "unknown_seed_status_counts": {
                status: sum(row["status"] == status for row in seed_sweep_records)
                for status in sorted({str(row["status"]) for row in seed_sweep_records})
            },
            "freshness": {"status": freshness.status, "mismatches": list(freshness.mismatches)},
        }
        _write(output / "summary.json", summary)
        return summary
    except Exception as error:
        if begin is not None and not (output / "evidence.json").exists():
            events.append({"iteration": len(events), "status": "failed", "data": {"error": str(error)}})
            _write(output / "evidence.json", finish_run(
                begin,
                events=events,
                output_files=sorted(
                    path for path in output.rglob("*")
                    if path.is_file() and path.name != "begin.json"
                ),
                measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick},
            ).as_dict())
        _write(output / "failure.json", {"status": "FAILED", "error": str(error), "traceback": traceback.format_exc()})
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--build-directory", type=Path, required=True)
    parser.add_argument("--dependency-root", action="append", required=True, metavar="NAME=PATH")
    args = parser.parse_args()
    roots: dict[str, Path] = {}
    for item in args.dependency_root:
        name, separator, value = item.partition("=")
        if not name or not separator or name in roots:
            parser.error("invalid or duplicate dependency root")
        roots[name] = Path(value).resolve()
    result = run(args.output_directory, build_directory=args.build_directory, dependency_roots=roots)
    print(canonical(result))
    return 0 if result["freshness"]["status"] == "CURRENT" else 1


if __name__ == "__main__":
    raise SystemExit(main())
