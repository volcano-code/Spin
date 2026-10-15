#!/usr/bin/env python3
"""Capture the asset-independent Spin-It Faster two-dimensional method suite."""
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
SCOPE = "SYNTHETIC_PLANAR_DOGLEG_POLYGON_DISK_MOMENT_AND_EXTRUSION_METHODS_ONLY_NOT_AUTHOR_2D_CASES"
COMMANDS = (
    ("dogleg", "planar-fidget-dogleg-regression", (
        "FASTER_PLANAR_COMPONENT_JSON=",
    )),
    ("polygon_disk", "planar-disk-intersection-regression", (
        "FASTER_PLANAR_COMPONENT_JSON=",
    )),
    ("moment_problem", "planar-spin-moment-problem-regression", ()),
    ("mass_moments", "planar-mass-moments-regression", ()),
    ("extrusion", "planar-extrusion-regression", (
        "FASTER_PLANAR_COMPONENT_JSON=",
        "[PASS] planar ring moments equal independently integrated 3D extrusion",
        "[PASS] planar and 3D fidget objectives/constraints agree",
    )),
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
    binaries = [(name, (build / binary).resolve(strict=True), markers)
                for name, binary, markers in COMMANDS]
    if not all(binary.is_file() and binary.stat().st_mode & 0o111
               for _, binary, _ in binaries):
        raise ValueError("required Faster planar regression executable is unavailable")
    started, tick, begin, events = datetime.now(timezone.utc).isoformat(), time.perf_counter(), None, []
    try:
        inputs = parse(canonical({
            "schema": "faster-planar-methods-inputs-v1",
            "classification": "SYNTHETIC_REGRESSION", "scope": SCOPE,
            "commands": [{"name": name, "binary": {"path": str(binary), "sha256": _sha256(binary)},
                          "required_markers": list(markers)} for name, binary, markers in binaries],
        }))
        input_path = output / "inputs.json"
        _write(input_path, inputs)
        begin = begin_run(
            repository_root=ROOT, runner_files=[Path(__file__)],
            checker_files=[
                ROOT / "tests/planar_fidget_dogleg_regression.cpp",
                ROOT / "tests/planar_disk_intersection_regression.cpp",
                ROOT / "tests/planar_spin_moment_problem_regression.cpp",
                ROOT / "tests/planar_mass_moments_regression.cpp",
                ROOT / "tests/planar_extrusion_regression.cpp",
                ROOT / "src/spin/PlanarFidgetCircleProblem.cpp",
                ROOT / "src/spin/PlanarDiskIntersection.cpp",
                ROOT / "src/spin/PlanarSpinMomentProblem.cpp",
                ROOT / "src/spin/PlanarMassMoments.cpp",
            ],
            asset_files=[input_path, *(binary for _, binary, _ in binaries)],
            build_configuration_files=[ROOT / "CMakeLists.txt", build / "CMakeCache.txt"],
            parameters=inputs,
            seed_protocol={"seed": 0, "protocol": "deterministic native planar method regressions"},
            dependency_roots=dependency_roots, native_modules=[], pure_python=True,
            source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY,
        )
        _write(output / "begin.json", begin.as_dict())
        commands = []
        for name, binary, markers in binaries:
            completed = subprocess.run([str(binary)], cwd=build, capture_output=True, text=True, check=False)
            transcript = output / f"{name}.stdout"
            raw_output = completed.stdout + completed.stderr
            capture = "VERBATIM_STDOUT_STDERR_V1"
            if raw_output:
                transcript.write_text(raw_output, encoding="utf-8")
            else:
                capture = "EMPTY_STDOUT_STDERR_RECEIPT_V1"
                _write(transcript, {
                    "schema": "faster-planar-empty-output-receipt-v1",
                    "raw_stdout_stderr_bytes": 0,
                    "raw_stdout_stderr_sha256": hashlib.sha256(b"").hexdigest(),
                    "returncode": completed.returncode,
                })
            content = transcript.read_text(encoding="utf-8")
            if completed.returncode != 0 or any(marker not in raw_output for marker in markers):
                raise RuntimeError("native Faster planar regression failed: " + name)
            record = {"name": name, "binary": str(binary), "binary_sha256": _sha256(binary),
                      "path": str(transcript), "sha256": _sha256(transcript),
                      "returncode": completed.returncode, "required_markers": list(markers),
                      "capture": capture}
            commands.append(record)
            events.append({"iteration": len(events), "status": "ok", "data": {
                "kind": "FASTER_PLANAR_NATIVE_REGRESSION", "name": name,
                "transcript_sha256": record["sha256"],
            }})
        component_results = {}
        for name in ("dogleg", "polygon_disk", "extrusion"):
            transcript = output / f"{name}.stdout"
            lines = [line for line in transcript.read_text(encoding="utf-8").splitlines()
                     if line.startswith("FASTER_PLANAR_COMPONENT_JSON=")]
            if len(lines) != 1:
                raise RuntimeError("native planar component result missing or duplicated: " + name)
            result = parse(lines[0].split("=", 1)[1])
            if result.get("schema") != "faster-planar-component-result-v1" or \
                    result.get("classification") != "SYNTHETIC_REGRESSION" or \
                    result.get("component") != name:
                raise RuntimeError("native planar component result has invalid identity: " + name)
            component_results[name] = result
        _write(output / "scientific.json", {
            "schema": "faster-planar-methods-scientific-v1",
            "classification": "SYNTHETIC_REGRESSION", "scope": SCOPE,
            "commands": commands, "component_results": component_results,
        })
        evidence = finish_run(begin, events=events,
            output_files=sorted(path for path in output.rglob("*")
                                if path.is_file() and path.name != "begin.json"),
            measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick})
        _write(output / "evidence.json", evidence.as_dict())
        freshness = check_freshness(evidence)
        summary = {"status": "SYNTHETIC_PLANAR_METHODS_VALIDATED",
                   "classification": "SYNTHETIC_REGRESSION", "scope": SCOPE,
                   "project_status": "PARTIAL_BLOCKED",
                   "original_evidence_admission": "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS",
                   "paper_case_admission": "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS",
                   "freshness": {"status": freshness.status, "mismatches": list(freshness.mismatches)}}
        _write(output / "summary.json", summary)
        return summary
    except Exception as error:
        if begin is not None and not (output / "evidence.json").exists():
            events.append({"iteration": len(events), "status": "failed", "data": {"error": str(error)}})
            _write(output / "evidence.json", finish_run(begin, events=events,
                output_files=sorted(path for path in output.rglob("*")
                                    if path.is_file() and path.name != "begin.json"),
                measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick}).as_dict())
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
    result = run(args.output_directory, build_directory=args.build_directory, dependency_roots=roots)
    print(canonical(result))
    return 0 if result["freshness"]["status"] == "CURRENT" else 1


if __name__ == "__main__":
    raise SystemExit(main())
