#!/usr/bin/env python3
"""Capture the protected Faster Equation-14 Newton modes as strict evidence.

The native regression is deliberately composed of analytic/synthetic moment
maps.  That gives this runner a narrow and useful role: it binds execution of
the actual spherical Newton implementation, its raw-residual stopping rule,
the dual-monomial coordinate pairing, and restart accounting.  It cannot
establish that an author asset was available, that a Boolean geometry is
globally optimal, or that a published application was reproduced.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import subprocess
import sys
import time
import traceback

from make_it_stand.verification.run_evidence import (
    EXPLICIT_EFFECTIVE_SOURCE_POLICY, begin_run, check_freshness, finish_run,
)
from make_it_stand.verification.spin_common import canonical, parse


ROOT = Path(__file__).resolve().parents[3]
_BINARY = "quadric-hollowing-solver-regression"
_MODES = (
    ("default", (), ("[PASS] scenario=quadric-hollowing-solver iterations=", " residual=", " trace=")),
    ("raw-stopping", ("--raw-stopping",), ()),
    ("dual-pairing", ("--dual-pairing",), ()),
    ("dual-coordinates", ("--dual-coordinates",), ()),
)


def _write(path: Path, value: object) -> None:
    with path.open("x", encoding="ascii") as stream:
        stream.write(canonical(value) + "\n")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(output_directory: Path, *, build_directory: Path,
        dependency_roots: dict[str, Path]) -> dict[str, object]:
    output = output_directory.absolute()
    output.mkdir(parents=False, exist_ok=False)
    build = build_directory.resolve(strict=True)
    executable = (build / _BINARY).resolve(strict=True)
    if not executable.is_file() or not executable.stat().st_mode & 0o111:
        raise ValueError("required native Faster Equation-14 executable is unavailable")

    started, tick, begin, events = datetime.now(timezone.utc).isoformat(), time.perf_counter(), None, []
    try:
        inputs = parse(canonical({
            "schema": "faster-eq14-newton-inputs-v1",
            "classification": "SYNTHETIC_REGRESSION",
            "binary": {"path": str(executable), "sha256": _sha256(executable)},
            "modes": [{"name": name, "arguments": list(arguments), "markers": list(markers)}
                      for name, arguments, markers in _MODES],
            "paper_constants": {
                "line_search_reduction": 0.5,
                "maximum_geodesic_step": 0.5,
                "equation_tolerance": 1e-12,
                "maximum_function_evaluations_per_restart": 60,
            },
        }))
        input_path = output / "inputs.json"
        _write(input_path, inputs)
        begin = begin_run(
            repository_root=ROOT,
            runner_files=[Path(__file__)],
            checker_files=[
                ROOT / "tests/quadric_hollowing_solver_regression.cpp",
                ROOT / "tests/dual_monomial_newton_regression.cpp",
                ROOT / "src/spin/QuadricHollowingSolver.cpp",
                ROOT / "src/spin/DualMonomialCoordinates.cpp",
                ROOT / "include/spin/QuadricHollowingSolver.h",
                ROOT / "include/spin/DualMonomialCoordinates.h",
                ROOT / "include/spin/SpinDesignRequest.h",
            ],
            asset_files=[input_path, executable],
            build_configuration_files=[ROOT / "CMakeLists.txt", build / "CMakeCache.txt"],
            parameters=inputs,
            seed_protocol={"seed": 0, "protocol": "deterministic native Equation-14 regression modes"},
            dependency_roots=dependency_roots,
            native_modules=[],
            pure_python=True,
            source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY,
        )
        _write(output / "begin.json", begin.as_dict())
        commands = []
        for name, arguments, markers in _MODES:
            completed = subprocess.run([str(executable), *arguments], cwd=build,
                                       capture_output=True, text=True, check=False)
            transcript = completed.stdout + completed.stderr
            # Three protected native modes deliberately report success only
            # through their exit status.  Preserve that fact in a non-empty
            # evidence record; do not manufacture a solver-side marker.
            if not transcript:
                transcript = "[OBSERVED] native command returned 0 with no stdout or stderr\n"
            transcript_path = output / f"newton-{name}.stdout"
            with transcript_path.open("x", encoding="utf-8") as stream:
                stream.write(transcript)
            if completed.returncode != 0 or any(marker not in transcript for marker in markers):
                raise RuntimeError(f"native Faster Equation-14 mode failed: {name}")
            commands.append({
                "name": name, "arguments": list(arguments), "returncode": completed.returncode,
                "path": str(transcript_path), "sha256": _sha256(transcript_path),
                "required_markers": list(markers),
            })
        scientific = {
            "schema": "faster-eq14-newton-scientific-v1",
            "classification": "SYNTHETIC_REGRESSION",
            "scope": "SYNTHETIC_NATIVE_EQ14_NEWTON_GEODESIC_MULTIPLIER_AND_RESTART_METHODS_ONLY_NOT_AUTHOR_CASE_OR_GLOBAL_OPTIMALITY",
            "commands": commands,
        }
        scientific_path = output / "scientific.json"
        _write(scientific_path, scientific)
        events.append({"iteration": 0, "status": "ok", "data": {
            "kind": "FASTER_EQ14_NATIVE_REGRESSION_MODES",
            "command_hashes": {row["name"]: row["sha256"] for row in commands},
        }})
        evidence = finish_run(
            begin, events=events,
            output_files=sorted(path for path in output.rglob("*")
                                if path.is_file() and path.name != "begin.json"),
            measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick},
        )
        _write(output / "evidence.json", evidence.as_dict())
        freshness = check_freshness(evidence)
        summary = {
            "status": "SYNTHETIC_EQ14_NEWTON_METHOD_VALIDATED",
            "classification": "SYNTHETIC_REGRESSION",
            "project_status": "PARTIAL_BLOCKED",
            "original_evidence_admission": "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS",
            "paper_case_admission": "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS",
            "freshness": {"status": freshness.status, "mismatches": list(freshness.mismatches)},
        }
        _write(output / "summary.json", summary)
        return summary
    except Exception as error:
        if begin is not None and not (output / "evidence.json").exists():
            events.append({"iteration": len(events), "status": "failed", "data": {"error": str(error)}})
            _write(output / "evidence.json", finish_run(
                begin, events=events,
                output_files=sorted(path for path in output.rglob("*")
                                    if path.is_file() and path.name != "begin.json"),
                measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick},
            ).as_dict())
        _write(output / "failure.json", {"status": "FAILED", "error": str(error),
                                           "traceback": traceback.format_exc()})
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
    result = run(args.output_directory, build_directory=args.build_directory,
                 dependency_roots=roots)
    print(canonical(result))
    return 0 if result["freshness"]["status"] == "CURRENT" else 1


if __name__ == "__main__":
    raise SystemExit(main())
