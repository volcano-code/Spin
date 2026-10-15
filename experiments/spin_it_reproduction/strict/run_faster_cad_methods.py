#!/usr/bin/env python3
"""Capture a bounded, pinned-OCCT CAD transport regression for Spin-It Faster.

The fixture is deliberately a synthetic Box-minus-Sphere material object.  It
does not admit or substitute for a paper CAD asset or an original experiment.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import math
from pathlib import Path
import subprocess
import time
import traceback

from make_it_stand.verification.run_evidence import (
    EXPLICIT_EFFECTIVE_SOURCE_POLICY, begin_run, check_freshness, finish_run,
)
from make_it_stand.verification.spin_common import canonical, parse


ROOT = Path(__file__).resolve().parents[3]
SCOPE = "SYNTHETIC_PINNED_OCCT_7_9_3_BREP_STEP_ROUNDTRIP_AND_ANALYTIC_MASS_METHOD_ONLY_NOT_AUTHOR_CAD_CASE"
OCCT_ARCHIVE = ROOT / "external/reference-sources-20260917/archives/OCCT-a016080b.tar.gz"


def _write(path: Path, value: object) -> None:
    with path.open("x", encoding="ascii") as stream:
        stream.write(canonical(value) + "\n")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _validate(result: object) -> dict[str, object]:
    if not isinstance(result, dict):
        raise ValueError("CAD probe did not return a JSON object")
    expected_volume = 8.0 - math.pi / 6.0
    if result.get("classification") != "SYNTHETIC_REGRESSION" or \
            result.get("status") != "SYNTHETIC_OCCT_STEP_BREP_ROUNDTRIP_VALIDATED" or \
            not math.isclose(float(result.get("analytic_volume", math.nan)), expected_volume,
                             rel_tol=0.0, abs_tol=1e-14):
        raise ValueError("CAD probe has unexpected classification, status, or analytic volume")
    moments = result.get("moments")
    if not isinstance(moments, dict) or set(moments) != {"original", "brep", "step"}:
        raise ValueError("CAD probe omits an original, B-rep, or STEP moment")
    for name, row in moments.items():
        if not isinstance(row, dict) or not math.isclose(float(row.get("volume", math.nan)), expected_volume,
                                                          rel_tol=0.0, abs_tol=1e-10) or \
                row.get("center") != [1, 1, 1]:
            raise ValueError("CAD probe returns an invalid " + name + " mass moment")
    return parse(canonical(result))


def run(output_directory: Path, *, probe: Path,
        dependency_roots: dict[str, Path]) -> dict[str, object]:
    output = output_directory.absolute()
    output.mkdir(parents=False, exist_ok=False)
    probe = probe.resolve(strict=True)
    if not probe.is_file() or not probe.stat().st_mode & 0o111:
        raise ValueError("pinned OCCT CAD probe executable is unavailable")
    archive = OCCT_ARCHIVE.resolve(strict=True)
    started, tick, begin, events = datetime.now(timezone.utc).isoformat(), time.perf_counter(), None, []
    try:
        inputs = parse(canonical({
            "schema": "faster-cad-methods-inputs-v1", "classification": "SYNTHETIC_REGRESSION",
            "scope": SCOPE, "probe": {"path": str(probe), "sha256": _sha256(probe)},
            "geometry": "BOX_2x2x2_MINUS_CENTERED_SPHERE_RADIUS_0.5",
            "occt_version": "7.9.3", "occt_commit": "a016080bf6738d6aeae020badee4e888ad1540a5",
        }))
        input_path = output / "inputs.json"
        _write(input_path, inputs)
        begin = begin_run(
            repository_root=ROOT, runner_files=[Path(__file__)],
            checker_files=[ROOT / "experiments/occt_step_roundtrip/CMakeLists.txt",
                           ROOT / "experiments/occt_step_roundtrip/main.cpp",
                           ROOT / "tests/python/test_occt_step_roundtrip.py"],
            asset_files=[input_path, probe, archive],
            build_configuration_files=[probe.parent / "CMakeCache.txt"], parameters=inputs,
            seed_protocol={"seed": 0, "protocol": "deterministic single synthetic CAD transport fixture"},
            dependency_roots=dependency_roots, native_modules=[], pure_python=True,
            source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY,
        )
        _write(output / "begin.json", begin.as_dict())
        brep, step, transcript = output / "material.brep", output / "step.step", output / "probe.stdout"
        completed = subprocess.run([str(probe), str(brep), str(step)], cwd=probe.parent,
                                   capture_output=True, text=True, check=False)
        transcript.write_text(completed.stdout + completed.stderr, encoding="utf-8")
        if completed.returncode != 0:
            raise RuntimeError("pinned OCCT CAD probe failed")
        result = _validate(parse(completed.stdout))
        if not brep.is_file() or not step.is_file() or brep.stat().st_size == 0 or step.stat().st_size == 0:
            raise RuntimeError("CAD probe did not emit nonempty B-rep and STEP artifacts")
        _write(output / "result.json", result)
        scientific = {"schema": "faster-cad-methods-scientific-v1", "classification": "SYNTHETIC_REGRESSION",
                      "scope": SCOPE, "probe": {"path": str(probe), "sha256": _sha256(probe),
                      "returncode": completed.returncode, "transcript": {"path": str(transcript), "sha256": _sha256(transcript)}},
                      "artifacts": {"brep": {"path": str(brep), "sha256": _sha256(brep)},
                                    "step": {"path": str(step), "sha256": _sha256(step)}}, "result": result}
        _write(output / "scientific.json", scientific)
        events.append({"iteration": 0, "status": "ok", "data": {"kind": "FASTER_PINNED_OCCT_CAD_METHOD",
                       "probe_sha256": _sha256(probe), "result": result}})
        evidence = finish_run(begin, events=events,
            output_files=sorted(path for path in output.rglob("*") if path.is_file() and path.name != "begin.json"),
            measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick})
        _write(output / "evidence.json", evidence.as_dict())
        freshness = check_freshness(evidence)
        summary = {"status": "SYNTHETIC_PINNED_OCCT_CAD_METHOD_VALIDATED", "classification": "SYNTHETIC_REGRESSION",
                   "scope": SCOPE, "project_status": "PARTIAL_BLOCKED",
                   "original_evidence_admission": "BLOCKED_MISSING_AUTHOR_CAD_ASSETS_AND_PARAMETERS",
                   "paper_case_admission": "BLOCKED_MISSING_AUTHOR_CAD_ASSETS_AND_PARAMETERS",
                   "freshness": {"status": freshness.status, "mismatches": list(freshness.mismatches)}}
        _write(output / "summary.json", summary)
        return summary
    except Exception as error:
        if begin is not None and not (output / "evidence.json").exists():
            events.append({"iteration": len(events), "status": "failed", "data": {"error": str(error)}})
            _write(output / "evidence.json", finish_run(begin, events=events,
                output_files=sorted(path for path in output.rglob("*") if path.is_file() and path.name != "begin.json"),
                measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick}).as_dict())
        _write(output / "failure.json", {"status": "FAILED", "error": str(error), "traceback": traceback.format_exc()})
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--dependency-root", action="append", required=True, metavar="NAME=PATH")
    args = parser.parse_args()
    roots: dict[str, Path] = {}
    for item in args.dependency_root:
        name, separator, value = item.partition("=")
        if not name or not separator or name in roots:
            parser.error("invalid or duplicate dependency root")
        roots[name] = Path(value).resolve()
    result = run(args.output_directory, probe=args.probe, dependency_roots=roots)
    print(canonical(result))
    return 0 if result["freshness"]["status"] == "CURRENT" else 1


if __name__ == "__main__":
    raise SystemExit(main())
