#!/usr/bin/env python3
"""Capture the asset-independent Spin-It Faster Section 4.6 contract."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import subprocess
import time
import traceback

from make_it_stand.verification.run_evidence import EXPLICIT_EFFECTIVE_SOURCE_POLICY, begin_run, check_freshness, finish_run
from make_it_stand.verification.spin_common import canonical, parse

ROOT = Path(__file__).resolve().parents[3]
SCOPE = "SYNTHETIC_SPINNING_TOP_GRAVITY_AXIS_FREE_VERTICAL_CENTER_AND_EQ7_COEFFICIENT_CONDITIONS_ONLY_NOT_AUTHOR_CASE"
MARKERS = (
    "[PASS] spinning-top gravity-axis objective, free vertical center, and inertia-block constraints agree independently",
    "[PASS] spinning-top Eq. 7 coefficient condition A33<0 and p!=0 agrees independently",
    "[PASS] spinning-top application scope is synthetic method validation only",
)

def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def write(path: Path, value: object) -> None:
    path.write_text(canonical(value) + "\n", encoding="ascii")

def run(output_directory: Path, *, build_directory: Path, dependency_roots: dict[str, Path]) -> dict[str, object]:
    output = output_directory.absolute(); output.mkdir(parents=False, exist_ok=False)
    build = build_directory.resolve(strict=True); binary = (build / "spin-moment-problem-regression").resolve(strict=True)
    if not binary.is_file() or not binary.stat().st_mode & 0o111:
        raise ValueError("required spinning-top contract executable is unavailable")
    started, tick, begin, events = datetime.now(timezone.utc).isoformat(), time.perf_counter(), None, []
    try:
        inputs = parse(canonical({"schema": "faster-spinning-top-contract-inputs-v1", "classification": "SYNTHETIC_REGRESSION", "scope": SCOPE, "command": {"name": "spinning-top-contract", "binary": {"path": str(binary), "sha256": sha256(binary)}, "required_markers": list(MARKERS)}}))
        input_path = output / "inputs.json"; write(input_path, inputs)
        begin = begin_run(repository_root=ROOT, runner_files=[Path(__file__)], checker_files=[ROOT / "tests/spin_moment_problem_regression.cpp", ROOT / "include/spin/FasterTheoremConditions.h", ROOT / "src/spin/FasterTheoremConditions.cpp", ROOT / "include/spin/SpinMomentProblem.h", ROOT / "src/spin/SpinMomentProblem.cpp"], asset_files=[input_path, binary], build_configuration_files=[ROOT / "CMakeLists.txt", build / "CMakeCache.txt"], parameters=inputs, seed_protocol={"seed": 0, "protocol": "deterministic native spinning-top-application contract regression"}, dependency_roots=dependency_roots, native_modules=[], pure_python=True, source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY)
        write(output / "begin.json", begin.as_dict())
        completed = subprocess.run([str(binary)], cwd=build, capture_output=True, text=True, check=False)
        transcript = output / "spinning-top-contract.stdout"; transcript.write_text(completed.stdout + completed.stderr, encoding="utf-8")
        text = transcript.read_text(encoding="utf-8")
        if completed.returncode != 0 or any(marker not in text for marker in MARKERS):
            raise RuntimeError("native Faster spinning-top contract regression failed")
        command = {"name": "spinning-top-contract", "binary": str(binary), "binary_sha256": sha256(binary), "returncode": completed.returncode, "required_markers": list(MARKERS), "path": str(transcript), "sha256": sha256(transcript)}
        events.append({"iteration": 0, "status": "ok", "data": {"kind": "FASTER_SPINNING_TOP_APPLICATION_CONTRACT", "transcript_sha256": command["sha256"]}})
        write(output / "scientific.json", {"schema": "faster-spinning-top-contract-scientific-v1", "classification": "SYNTHETIC_REGRESSION", "scope": SCOPE, "command": command})
        evidence = finish_run(begin, events=events, output_files=sorted(p for p in output.rglob("*") if p.is_file() and p.name != "begin.json"), measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick})
        write(output / "evidence.json", evidence.as_dict()); freshness = check_freshness(evidence)
        summary = {"status": "SYNTHETIC_SPINNING_TOP_APPLICATION_CONTRACT_VALIDATED", "classification": "SYNTHETIC_REGRESSION", "scope": SCOPE, "project_status": "PARTIAL_BLOCKED", "original_evidence_admission": "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS", "paper_case_admission": "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS", "freshness": {"status": freshness.status, "mismatches": list(freshness.mismatches)}}
        write(output / "summary.json", summary); return summary
    except Exception as error:
        if begin is not None and not (output / "evidence.json").exists():
            events.append({"iteration": len(events), "status": "failed", "data": {"error": str(error)}})
            write(output / "evidence.json", finish_run(begin, events=events, output_files=sorted(p for p in output.rglob("*") if p.is_file() and p.name != "begin.json"), measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick}).as_dict())
        write(output / "failure.json", {"status": "FAILED", "error": str(error), "traceback": traceback.format_exc()}); raise

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--output-directory", type=Path, required=True); parser.add_argument("--build-directory", type=Path, required=True); parser.add_argument("--dependency-root", action="append", required=True, metavar="NAME=PATH")
    args = parser.parse_args(); roots: dict[str, Path] = {}
    for item in args.dependency_root:
        name, separator, value = item.partition("=")
        if not name or not separator or name in roots: parser.error("invalid or duplicate dependency root")
        roots[name] = Path(value).resolve()
    result = run(args.output_directory, build_directory=args.build_directory, dependency_roots=roots); print(canonical(result)); return 0 if result["freshness"]["status"] == "CURRENT" else 1

if __name__ == "__main__": raise SystemExit(main())
