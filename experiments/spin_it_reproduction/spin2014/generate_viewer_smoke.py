#!/usr/bin/env python3
"""Generate the 16^3 reconstructed static artifact used by the OpenGL smoke test."""

from __future__ import annotations

import argparse
from pathlib import Path

from make_it_stand import SpinIt2014Request, solve_spin_it_2014


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    side = 16
    material = [1] * (side**3)
    carveable = [0] * (side**3)
    for z in range(2, side - 2):
        for y in range(2, side - 2):
            for x in range(2, side - 2):
                carveable[x + side * (y + side * z)] = 1
    request = SpinIt2014Request(
        dimensions=(side, side, side),
        minimum=(-1.0, -1.0, 0.0),
        maximum=(1.0, 1.0, 2.0),
        material_mask=tuple(material),
        carveable_mask=tuple(carveable),
        mode="STATIC_BALANCE",
        evidence="RECONSTRUCTED_ASSET",
        initial_level=2,
        maximum_level=4,
        maximum_adaptive_passes=5,
    )
    first = solve_spin_it_2014(request)
    repeated = solve_spin_it_2014(request)
    if not first.succeeded or not first.geometry_sealed:
        raise RuntimeError(f"viewer smoke artifact failed: {first.status}")
    if first.diagnostics["output_hash"] != repeated.diagnostics["output_hash"]:
        raise RuntimeError("viewer smoke repeat output hash differs")
    first.save_viewer_bundle(args.output)
    print(first.diagnostics["output_hash"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
