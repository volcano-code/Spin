#!/usr/bin/env python3
"""Run deterministic reconstructed Spin-It 2014 mode fixtures twice."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

from make_it_stand import SpinIt2014Request, solve_spin_it_2014


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _fixture(mode: str) -> SpinIt2014Request:
    # Keep the cross-mode oracle deliberately small: it is a deterministic
    # formula/geometry contract, not a performance benchmark.  The separate
    # viewer-smoke fixture uses 16^3 cells for visual inspection.
    side = 4
    material = [1] * (side**3)
    carveable = [0] * (side**3)
    for z in range(1, side - 1):
        for y in range(1, side - 1):
            for x in range(1, side - 1):
                carveable[x + side * (y + side * z)] = 1
    minimum = (-1.0, -1.0, -1.0) if mode == "YO_YO" else (-1.0, -1.0, 0.0)
    maximum = (1.0, 1.0, 1.0) if mode == "YO_YO" else (1.0, 1.0, 2.0)
    return SpinIt2014Request(
        dimensions=(side, side, side),
        minimum=minimum,
        maximum=maximum,
        material_mask=tuple(material),
        carveable_mask=tuple(carveable),
        mode=mode,
        evidence="RECONSTRUCTED_ASSET",
        initial_level=1,
        maximum_level=2,
        maximum_adaptive_passes=3,
        initial_fill=0.5,
        fractional_epsilon=0.05,
        laplacian_weight=1e-3,
        center_of_mass_weight=1.0,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    source_commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], text=True
    ).strip()
    cases: list[dict[str, object]] = []
    for mode in ("STATIC_BALANCE", "SPINNING_TOP", "YO_YO"):
        request = _fixture(mode)
        start = time.monotonic()
        first = solve_spin_it_2014(request)
        repeated = solve_spin_it_2014(request)
        duration = time.monotonic() - start
        if not first.succeeded or not repeated.succeeded:
            raise RuntimeError(f"{mode} did not converge")
        first_hash = first.diagnostics["output_hash"]
        repeat_hash = repeated.diagnostics["output_hash"]
        if first_hash != repeat_hash:
            raise RuntimeError(f"{mode} repeat hash mismatch")
        case_root = args.output / mode.lower().replace("_", "-")
        saved = [Path(path) for path in first.save_viewer_bundle(case_root)]
        cases.append(
            {
                "mode": mode,
                "status": first.status,
                "evidence": first.evidence,
                "geometry_sealed": first.geometry_sealed,
                "output_hash": first_hash,
                "repeat_output_hash": repeat_hash,
                "repeat_identical": first_hash == repeat_hash,
                "binary_constraint_infinity_norm": first.diagnostics[
                    "binary_constraint_infinity_norm"
                ],
                "adaptive_pass_count": len(first.diagnostics["adaptive_passes"]),
                "center_of_mass_model": first.diagnostics["center_of_mass_model"],
                "files": {path.name: _sha256(path) for path in saved},
                "wall_seconds_excluded_from_hash": duration,
            }
        )
    summary = {
        "schema": "spin-it-2014-reconstructed-modes-v1",
        "schema_version": 1,
        "source_commit": source_commit,
        "claim_level": "RECONSTRUCTED_ASSET",
        "all_converged": all(case["status"] == "CONVERGED" for case in cases),
        "all_sealed": all(bool(case["geometry_sealed"]) for case in cases),
        "all_repeats_identical": all(bool(case["repeat_identical"]) for case in cases),
        "cases": cases,
    }
    summary_path = args.output / "summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
