#!/usr/bin/env python3
"""Run the complete source-cage and BBW native regressions as strict evidence.

The two regressions exercise all source-intersecting cells as cage handles and
the actual bounded-biharmonic QPs.  They are deliberately synthetic method
evidence, not author-model or paper-case evidence.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import subprocess
import time
import traceback

from make_it_stand.verification.run_evidence import (
    EXPLICIT_EFFECTIVE_SOURCE_POLICY, begin_run, check_freshness, finish_run,
)
from make_it_stand.verification.spin_common import canonical, native, parse

ROOT = Path(__file__).resolve().parents[3]


def _write(path: Path, value) -> None:
    with path.open("x", encoding="ascii") as stream:
        stream.write(canonical(value) + "\n")


def _run(executable: Path, output: Path, required: tuple[str, ...]) -> dict:
    completed = subprocess.run([str(executable)], cwd=executable.parent,
                               capture_output=True, text=True, check=False)
    text = completed.stdout + completed.stderr
    with output.open("x", encoding="utf-8") as stream:
        stream.write(text)
    if completed.returncode != 0 or any(token not in text for token in required):
        raise RuntimeError("native cage regression failed: " + executable.name)
    return {"path": str(output), "sha256": hashlib.sha256(text.encode()).hexdigest(),
            "returncode": completed.returncode, "required_markers": list(required)}


def run(output_directory: Path, *, build_directory: Path, dependency_roots: dict[str, Path]) -> dict:
    output = output_directory.absolute(); output.mkdir(parents=False, exist_ok=False)
    build = build_directory.resolve()
    surface = build / "spin-it-2014-source-surface-cage-regression"
    bbw = build / "spin-it-2014-source-cage-bbw-regression"
    if not all(path.is_file() and path.stat().st_mode & 0o111 for path in (surface, bbw)):
        raise ValueError("required native cage regression executable is unavailable")
    started, tick, begin, events = datetime.now(timezone.utc).isoformat(), time.perf_counter(), None, []
    try:
        inputs = parse(canonical({"schema": "source-spin2014-cage-inputs-v1",
            "classification": "SYNTHETIC_REGRESSION",
            "surface_binary": {"path": str(surface), "sha256": hashlib.sha256(surface.read_bytes()).hexdigest()},
            "bbw_binary": {"path": str(bbw), "sha256": hashlib.sha256(bbw.read_bytes()).hexdigest()}}))
        input_path = output / "inputs.json"; _write(input_path, inputs)
        extension = native()
        if Path(extension.__file__).resolve().parents[1] != build:
            raise ValueError("declared build directory does not contain the loaded native extension")
        begin = begin_run(repository_root=ROOT, runner_files=[Path(__file__)],
            checker_files=[ROOT / "tests/spin_it_2014_source_surface_cage_regression.cpp",
                           ROOT / "tests/spin_it_2014_source_cage_bbw_regression.cpp",
                           ROOT / "src/spin2014/SourceSurfaceCage.cpp", ROOT / "src/spin2014/SourceCageBbw.cpp"],
            asset_files=[input_path, surface, bbw],
            build_configuration_files=[ROOT / "CMakeLists.txt", build / "CMakeCache.txt"],
            parameters=inputs, seed_protocol={"seed": 0, "protocol": "deterministic native regressions"},
            dependency_roots=dependency_roots, native_modules=[extension], pure_python=False,
            source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY)
        _write(output / "begin.json", begin.as_dict())
        surface_result = _run(surface, output / "surface.stdout",
                              ("level9 GEOMETRY:", "complete geometry only"))
        bbw_result = _run(bbw, output / "bbw.stdout",
                          ("actual_702_variable_and_26_column_box_qps", "SYNTHETIC_REGRESSION: 0 failed groups"))
        scientific = {"schema": "source-spin2014-cage-scientific-v1", "classification": "SYNTHETIC_REGRESSION",
                      "surface": surface_result, "bbw": bbw_result}
        _write(output / "scientific.json", scientific)
        for kind, data in (("COMPLETE_SOURCE_CAGE", surface_result), ("ACTUAL_BBW_QPS", bbw_result)):
            events.append({"iteration": len(events), "status": "ok", "data": {"kind": kind, **data}})
        evidence = finish_run(begin, events=events,
            output_files=sorted(path for path in output.rglob("*") if path.is_file() and path.name != "begin.json"),
            measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick})
        _write(output / "evidence.json", evidence.as_dict())
        freshness = check_freshness(evidence)
        summary = {"status": "SYNTHETIC_METHOD_VALIDATED", "classification": "SYNTHETIC_REGRESSION",
                   "project_status": "PARTIAL_BLOCKED", "original_evidence_admission": "BLOCKED_UNKNOWN",
                   "freshness": {"status": freshness.status, "mismatches": list(freshness.mismatches)}}
        _write(output / "summary.json", summary); return summary
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
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--build-directory", type=Path, required=True)
    parser.add_argument("--dependency-root", action="append", required=True, metavar="NAME=PATH")
    args = parser.parse_args(); roots = {}
    for item in args.dependency_root:
        name, sep, value = item.partition("=")
        if not name or not sep or name in roots: parser.error("invalid or duplicate dependency root")
        roots[name] = Path(value).resolve()
    result = run(args.output_directory, build_directory=args.build_directory, dependency_roots=roots)
    print(canonical(result)); return 0 if result["freshness"]["status"] == "CURRENT" else 1


if __name__ == "__main__":
    raise SystemExit(main())
