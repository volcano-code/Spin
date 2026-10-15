#!/usr/bin/env python3
"""Read saved plane material; independent tetrahedron sums, no solver import."""
import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np

from make_it_stand.verification import static_material
from make_it_stand.verification.run_evidence import (
    begin_run, finish_run, EXPLICIT_EFFECTIVE_SOURCE_POLICY)


def main():
    parser = argparse.ArgumentParser()
    source=parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--application", type=Path)
    source.add_argument("--spin2014-function-directory",type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[3]
    directory=(args.spin2014_function_directory.resolve(strict=True)
               if args.spin2014_function_directory is not None else None)
    application = args.application.resolve(strict=True) if directory is None else None
    assets=[application] if directory is None else [directory/name for name in
        ("frozen-input.json","native-run.json","input-geometry.npz")]
    scope=static_material.SCOPE if directory is None else static_material.SPIN2014_STATIC_SCOPE
    args.output.mkdir(parents=True, exist_ok=False)
    checks = root/"python/make_it_stand/verification"
    start = time.perf_counter()
    begin = begin_run(repository_root=root, runner_files=[Path(__file__)],
        checker_files=[checks/name for name in ("static_material.py", "mass_moments.py",
            "rational_text.py", "spin_common.py", "run_evidence.py")]+(
            [root/"python/make_it_stand/stability_oracle.py"] if directory is not None else []),
        asset_files=assets, build_configuration_files=[root/"pyproject.toml"],
        parameters={"scope": scope, "operation": "READ_SAVED_STATIC_GEOMETRY"},
        seed_protocol={"seed": 0, "protocol": "NO_RANDOMNESS_SAVED_OUTPUT_RECHECK"},
        dependency_roots={"numpy": Path(np.__file__).resolve().parent},
        native_modules=[], pure_python=True, source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY)
    (args.output/"begin.json").write_text(begin.to_json()+"\n")
    report = (static_material.audit_saved_static_application(application) if directory is None
              else static_material.audit_spin2014_static_directory(directory))
    if "_make_it_stand" in sys.modules:
        raise RuntimeError("independent static checker unexpectedly loaded the native solver")
    scientific = args.output/"scientific.json"
    scientific.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+"\n")
    events = ([{"iteration": i, "status": "ok", "data": {
        "kind": "INDEPENDENT_SAVED_STATIC_MOMENTS", "case": case["id"], **row}}
        for case in report["cases"] for i,row in enumerate(case["attempts"])] if directory is None else
        [{"iteration":0,"status":"ok" if report["terminal_numerics_qualified"] else "failed",
          "data":{"kind":"INDEPENDENT_SAVED_2014_SQUARE_MOMENTS_AND_SUPPORT",**report}}])
    evidence = finish_run(begin, events=events, output_files=[scientific],
                          measurements={"elapsed_seconds": time.perf_counter()-start})
    (args.output/"evidence.json").write_text(evidence.to_json()+"\n")
    print(json.dumps({"scope": report["scope"], "cases": len(report["cases"]) if directory is None else 1,
                      "native_solver_loaded": False, "evidence": str(args.output/"evidence.json")}))


if __name__ == "__main__": main()
