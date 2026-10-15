#!/usr/bin/env python3
"""Capture narrow Faster raw-moment and literal-duality method evidence.

This runner intentionally does *not* certify Theorem 1, Theorem 2, global
optimality, a Boolean construction, or an author case.  It only binds the two
native regression programs that independently exercise the ten-moment
convention and the cross-term coefficient/raw-moment pairing.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import subprocess
import time
import traceback

from make_it_stand.verification.run_evidence import (
    EXPLICIT_EFFECTIVE_SOURCE_POLICY, begin_run, check_freshness, finish_run,
)
from make_it_stand.verification.spin_common import canonical, parse


ROOT = Path(__file__).resolve().parents[3]
_COMMANDS = (
    (
        "raw-moments",
        "quadric-moment-evaluator-regression",
        (),
        (
            "[PASS] raw ten-moment volume and first-moment convention",
            "[PASS] scenario=quadric-moment-evaluator derivative_error=",
        ),
    ),
    (
        "dual-monomial",
        "quadric-hollowing-solver-regression",
        ("--dual-coordinates",),
        (
            "[PASS] literal coefficient/raw-moment cross-term duality",
            "[PASS] literal coefficient-to-geometry differential",
            "[PASS] dual monomial coordinate Newton",
        ),
    ),
)
_SCOPE = (
    "SYNTHETIC_RAW_TEN_MOMENT_AND_LITERAL_CROSS_TERM_COEFFICIENT_DUALITY_"
    "METHODS_ONLY_NOT_FASTER_THEOREMS_OR_AUTHOR_CASE"
)


def _write(path: Path, value: object) -> None:
    with path.open("x", encoding="ascii") as stream:
        stream.write(canonical(value) + "\n")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(output_directory: Path, *, build_directory: Path,
        dependency_roots: dict[str, Path]) -> dict[str, object]:
    output = output_directory.absolute()
    output.mkdir(parents=False, exist_ok=False)
    build = build_directory.resolve(strict=True)
    executables = []
    for _, binary, _, _ in _COMMANDS:
        executable = (build / binary).resolve(strict=True)
        if not executable.is_file() or not executable.stat().st_mode & 0o111:
            raise ValueError("required Faster moment-foundation executable is unavailable: " + binary)
        executables.append(executable)

    started, tick, begin, events = datetime.now(timezone.utc).isoformat(), time.perf_counter(), None, []
    try:
        inputs = parse(canonical({
            "schema": "faster-moment-foundations-inputs-v1",
            "classification": "SYNTHETIC_REGRESSION",
            "commands": [
                {"name": name, "binary": {"path": str(executable), "sha256": _sha256(executable)},
                 "arguments": list(arguments), "required_markers": list(markers)}
                for (name, _, arguments, markers), executable in zip(_COMMANDS, executables)
            ],
        }))
        input_path = output / "inputs.json"
        _write(input_path, inputs)
        begin = begin_run(
            repository_root=ROOT, runner_files=[Path(__file__)],
            checker_files=[
                ROOT / "tests/quadric_moment_evaluator_regression.cpp",
                ROOT / "tests/dual_monomial_newton_regression.cpp",
                ROOT / "include/spin/DualMonomialCoordinates.h",
                ROOT / "src/spin/DualMonomialCoordinates.cpp",
                ROOT / "include/physics/RawMassMoments10.h",
                ROOT / "src/physics/RawMassMoments10.cpp",
            ], asset_files=[input_path, *executables],
            build_configuration_files=[ROOT / "CMakeLists.txt", build / "CMakeCache.txt"],
            parameters=inputs,
            seed_protocol={"seed": 0, "protocol": "deterministic native raw-moment and coefficient-duality regressions"},
            dependency_roots=dependency_roots, native_modules=[], pure_python=True,
            source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY,
        )
        _write(output / "begin.json", begin.as_dict())
        commands = []
        for (name, _, arguments, markers), executable in zip(_COMMANDS, executables):
            completed = subprocess.run([str(executable), *arguments], cwd=build,
                                       capture_output=True, text=True, check=False)
            transcript = output / f"{name}.stdout"
            transcript.write_text(completed.stdout + completed.stderr, encoding="utf-8")
            text = transcript.read_text(encoding="utf-8")
            if completed.returncode != 0 or any(marker not in text for marker in markers):
                raise RuntimeError("native Faster moment-foundation regression failed: " + name)
            commands.append({
                "name": name, "arguments": list(arguments), "returncode": completed.returncode,
                "required_markers": list(markers), "path": str(transcript), "sha256": _sha256(transcript),
            })
        scientific_path = output / "scientific.json"
        _write(scientific_path, {
            "schema": "faster-moment-foundations-scientific-v1",
            "classification": "SYNTHETIC_REGRESSION", "scope": _SCOPE,
            "commands": commands,
        })
        events.append({"iteration": 0, "status": "ok", "data": {
            "kind": "FASTER_RAW_MOMENT_AND_COEFFICIENT_DUALITY", "commands": commands,
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
            "status": "SYNTHETIC_MOMENT_FOUNDATIONS_VALIDATED",
            "classification": "SYNTHETIC_REGRESSION", "project_status": "PARTIAL_BLOCKED",
            "original_evidence_admission": "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS",
            "paper_case_admission": "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS",
            "theorem_coverage": "NOT_GRANTED_BY_THIS_RUNNER",
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
        _write(output / "failure.json", {
            "status": "FAILED", "error": str(error), "traceback": traceback.format_exc(),
        })
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
