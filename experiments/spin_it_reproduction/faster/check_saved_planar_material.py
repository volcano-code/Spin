#!/usr/bin/env python3
"""Reintegrate saved line/arc extrusions; never invoke an optimizer or native integrator."""
import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np

from make_it_stand.verification import planar_material
from make_it_stand.verification.run_evidence import (
    begin_run, finish_run, EXPLICIT_EFFECTIVE_SOURCE_POLICY)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--application", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[3]
    application = args.application.resolve(strict=True)
    # New output identity; never overwrite an earlier checker or failed run.
    args.output.mkdir(parents=True, exist_ok=False)
    checks = root/"python/make_it_stand/verification"
    start = time.perf_counter()
    begin = begin_run(repository_root=root, runner_files=[Path(__file__)],
        checker_files=[checks/name for name in ("planar_material.py", "spin_common.py", "run_evidence.py")],
        asset_files=[application], build_configuration_files=[root/"pyproject.toml"],
        parameters={"scope": planar_material.SCOPE, "operation": "READ_SAVED_PLANAR_GEOMETRY"},
        seed_protocol={"seed": 0, "protocol": "NO_RANDOMNESS_SAVED_OUTPUT_RECHECK"},
        dependency_roots={"numpy": Path(np.__file__).resolve().parent},
        native_modules=[], pure_python=True, source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY)
    (args.output/"begin.json").write_text(begin.to_json()+"\n")
    report = planar_material.audit_saved_planar_application(application)
    if "_make_it_stand" in sys.modules:
        raise RuntimeError("independent planar checker unexpectedly loaded the native solver")
    scientific = args.output/"scientific.json"
    scientific.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+"\n")
    events = [{"iteration": index, "status": "ok" if row["moments_match"] else "failed",
               "data": {"kind": "INDEPENDENT_SAVED_PLANAR_MOMENTS", "case": case["id"], **row}}
              for case in report["cases"] for index, row in enumerate(case["attempts"])]
    evidence = finish_run(begin, events=events, output_files=[scientific],
                          measurements={"elapsed_seconds": time.perf_counter()-start})
    (args.output/"evidence.json").write_text(evidence.to_json()+"\n")
    print(json.dumps({"scope": report["scope"], "cases": len(report["cases"]),
                      "native_solver_loaded": False, "evidence": str(args.output/"evidence.json")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
