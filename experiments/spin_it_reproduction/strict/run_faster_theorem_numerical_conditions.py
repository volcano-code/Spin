#!/usr/bin/env python3
"""Record narrow Faster Eq. (7) condition and canonical-replay evidence.

This runner deliberately records *supporting* evidence only.  A finite
numerical Eq. (7) witness and deterministic replay of one reconstructed
material fixture do not prove either Faster theorem, establish global
optimality, or reproduce an author-provided model.
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
    EXPLICIT_EFFECTIVE_SOURCE_POLICY,
    begin_run,
    check_freshness,
    finish_run,
)
from make_it_stand.verification.faster_theorem_exact import (
    reference_certificate,
    verify_certificate,
)
from make_it_stand.verification.spin_common import canonical, parse


ROOT = Path(__file__).resolve().parents[3]
_SCOPE = (
    "SYNTHETIC_EXACT_RATIONAL_THEOREM_INSTANCES_AND_CANONICAL_EQ7_"
    "MATERIAL_METHODS_ONLY_NOT_GENERAL_THEOREM_PROOF_OR_AUTHOR_CASE"
)
_COMMANDS = (
    (
        "theorem-conditions",
        "faster-theorem-conditions-regression",
        (),
        (
            "[PASS] Eq. 7 numerical hypotheses accept exact fixture",
            "[PASS] Eq. 7 numerical hypotheses reject infeasible and nonpositive-mu fixtures",
            "[PASS] finite signed variation certificate proves a strict discrete increase",
            "[PASS] finite beta-zero cells never stand in for a measure-zero proof",
            "[PASS] theorem audit scope remains numerical only",
        ),
    ),
    (
        "canonical-replay",
        "canonical_full_rank_regression",
        ("replay",),
        (
            "same-q replay geometry=1 moment-diff=0 J-diff=0",
            "[PASS] replay",
        ),
    ),
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


def _executable(build: Path, binary: str) -> Path:
    # The current CMake target has the project-qualified name below, while
    # older isolated replay builds used the shorter historical name.  Accept
    # either only for this canonical replay executable and record the actual
    # path in the evidence; never silently substitute another binary.
    candidates = (binary, "spin-faster-canonical-fullrank-regression") \
        if binary == "canonical_full_rank_regression" else (binary,)
    for candidate in candidates:
        path = (build / candidate).resolve()
        if path.is_file() and path.stat().st_mode & 0o111:
            return path
    raise ValueError("required Faster theorem-condition executable is unavailable: " + binary)


def run(output_directory: Path, *, build_directory: Path,
        canonical_build_directory: Path, dependency_roots: dict[str, Path]) -> dict[str, object]:
    output = output_directory.absolute()
    output.mkdir(parents=False, exist_ok=False)
    build = build_directory.resolve(strict=True)
    canonical_build = canonical_build_directory.resolve(strict=True)
    executables = (
        _executable(build, _COMMANDS[0][1]),
        _executable(canonical_build, _COMMANDS[1][1]),
    )
    started, tick, begin, events = datetime.now(timezone.utc).isoformat(), time.perf_counter(), None, []
    try:
        exact_certificate_path = output / "exact-certificate.json"
        exact_input = parse(canonical(reference_certificate()))
        _write(exact_certificate_path, exact_input)
        exact_result = verify_certificate(exact_input)
        exact_result_path = output / "exact-instances.json"
        _write(exact_result_path, exact_result)
        inputs = parse(canonical({
            "schema": "faster-theorem-numerical-conditions-inputs-v2",
            "classification": "SYNTHETIC_REGRESSION",
            "scope": _SCOPE,
            "exact_certificate_sha256": _sha256(exact_certificate_path),
            "commands": [
                {"name": name, "binary": {"path": str(executable), "sha256": _sha256(executable)},
                 "arguments": list(arguments), "required_markers": list(markers)}
                for (name, _, arguments, markers), executable in zip(_COMMANDS, executables)
            ],
        }))
        input_path = output / "inputs.json"
        _write(input_path, inputs)
        build_configuration_files = [ROOT / "CMakeLists.txt", build / "CMakeCache.txt"]
        canonical_cache = canonical_build / "CMakeCache.txt"
        if canonical_cache.resolve() != (build / "CMakeCache.txt").resolve():
            build_configuration_files.append(canonical_cache)
        begin = begin_run(
            repository_root=ROOT,
            runner_files=[
                Path(__file__),
                ROOT / "python/make_it_stand/verification/run_evidence.py",
                ROOT / "python/make_it_stand/verification/spin_common.py",
            ],
            checker_files=[
                ROOT / "tests/faster_theorem_conditions_regression.cpp",
                ROOT / "tests/full_rank_hollowing_regression.cpp",
                ROOT / "include/spin/FasterTheoremConditions.h",
                ROOT / "src/spin/FasterTheoremConditions.cpp",
                ROOT / "include/spin/FullRankHollowingSolver.h",
                ROOT / "src/spin/FullRankHollowingSolver.cpp",
                ROOT / "include/spin/CanonicalFullRankHollowing.h",
                ROOT / "src/spin/CanonicalFullRankHollowing.cpp",
                ROOT / "src/spin/CanonicalPolyhedralMomentDifferential.cpp",
                ROOT / "experiments/m1_canonical_replay/CMakeLists.txt",
                ROOT / "python/make_it_stand/verification/faster_theorem_exact.py",
                ROOT / "tests/python/test_faster_theorem_exact_instances.py",
            ],
            asset_files=[input_path, exact_certificate_path, *executables],
            build_configuration_files=build_configuration_files,
            parameters=inputs,
            seed_protocol={"seed": 0, "protocol": "deterministic native Eq. 7 condition and canonical material replay"},
            dependency_roots=dependency_roots,
            native_modules=[], pure_python=True,
            source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY,
        )
        _write(output / "begin.json", begin.as_dict())
        commands = []
        for (name, _, arguments, markers), executable in zip(_COMMANDS, executables):
            completed = subprocess.run([str(executable), *arguments], cwd=executable.parent,
                                       capture_output=True, text=True, check=False)
            transcript = output / f"{name}.stdout"
            transcript.write_text(completed.stdout + completed.stderr, encoding="utf-8")
            text = transcript.read_text(encoding="utf-8")
            if completed.returncode != 0 or any(marker not in text for marker in markers):
                raise RuntimeError("native Faster theorem-condition regression failed: " + name)
            commands.append({
                "name": name, "arguments": list(arguments), "returncode": completed.returncode,
                "required_markers": list(markers), "path": str(transcript), "sha256": _sha256(transcript),
            })
        scientific_path = output / "scientific.json"
        _write(scientific_path, {
            "schema": "faster-theorem-numerical-conditions-scientific-v2",
            "classification": "SYNTHETIC_REGRESSION", "scope": _SCOPE,
            "commands": commands,
            "exact_instance_result_sha256": _sha256(exact_result_path),
            "exact_instance_coverage": exact_result["status"],
            "theorem_coverage": "NOT_GRANTED_BY_THIS_RUNNER",
        })
        events.append({"iteration": 0, "status": "ok", "data": {
            "kind": "FASTER_EXACT_THEOREM_INSTANCES_EQ7_AND_CANONICAL_MATERIAL_REPLAY",
            "command_hashes": {row["name"]: row["sha256"] for row in commands},
            "exact_instance_result_sha256": _sha256(exact_result_path),
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
            "status": "SYNTHETIC_EXACT_THEOREM_INSTANCES_AND_NUMERICAL_CONDITIONS_VALIDATED",
            "classification": "SYNTHETIC_REGRESSION", "project_status": "PARTIAL_BLOCKED",
            "scope": _SCOPE, "theorem_coverage": "NOT_GRANTED_BY_THIS_RUNNER",
            "exact_instance_coverage": exact_result["status"],
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
        _write(output / "failure.json", {
            "status": "FAILED", "error": str(error), "traceback": traceback.format_exc(),
        })
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--build-directory", type=Path, required=True)
    parser.add_argument("--canonical-build-directory", type=Path, required=True)
    parser.add_argument("--dependency-root", action="append", required=True, metavar="NAME=PATH")
    args = parser.parse_args()
    roots: dict[str, Path] = {}
    for item in args.dependency_root:
        name, separator, value = item.partition("=")
        if not name or not separator or name in roots:
            parser.error("invalid or duplicate dependency root")
        roots[name] = Path(value).resolve()
    result = run(args.output_directory, build_directory=args.build_directory,
                 canonical_build_directory=args.canonical_build_directory, dependency_roots=roots)
    print(canonical(result))
    return 0 if result["freshness"]["status"] == "CURRENT" else 1


if __name__ == "__main__":
    raise SystemExit(main())
