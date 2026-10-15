#!/usr/bin/env python3
"""Diagnose stored partition algebra; this does not rerun geometry or a solver."""
import argparse
import hashlib
import json
import math
from pathlib import Path


def decompose(leaves, relaxed, binary, density):
    if not (len(leaves) == len(relaxed) == len(binary)) or not leaves:
        raise ValueError("partition and fill arrays differ or are empty")
    if not math.isfinite(density) or density <= 0:
        raise ValueError("invalid density")
    rows, groups = [], {}
    identities = set()
    for i, (leaf, a, b) in enumerate(zip(leaves, relaxed, binary)):
        if not math.isfinite(a) or not 0 <= a <= 1 or b != float(a >= .5):
            raise ValueError("invalid fill or violation of original half-threshold rounding")
        identity = (leaf["level"], *leaf["coordinates"])
        if identity in identities:
            raise ValueError("duplicate partition leaf")
        identities.add(identity)
        m = leaf["moments"]
        if len(m) != 10 or any(not math.isfinite(v) for v in m):
            raise ValueError("invalid leaf moments")
        delta = [-density * (b-a) * v for v in m]
        rows.append({"index": i, "level": leaf["level"], "coordinates": leaf["coordinates"],
                     "relaxed": a, "binary": b, "moment_delta": delta})
        groups.setdefault(str(leaf["level"]), []).append(rows[-1])
    return {"moment_delta": [math.fsum(r["moment_delta"][j] for r in rows) for j in range(10)],
            "by_level": {k: {"cells": len(v), "moment_delta": [math.fsum(r["moment_delta"][j] for r in v)
                          for j in range(10)]} for k,v in sorted(groups.items())}, "cells": rows}


def constraints(m, origin, axes, phi):
    # axes are target basis COLUMNS; target coordinates are axes^T*(p-origin).
    first = m[1:4]
    raw = [[m[4],m[9],m[8]],[m[9],m[5],m[7]],[m[8],m[7],m[6]]]
    centered = [[raw[i][j]-origin[i]*first[j]-first[i]*origin[j]+m[0]*origin[i]*origin[j]
                 for j in range(3)] for i in range(3)]
    l = [math.fsum(axes[k][i]*(first[k]-m[0]*origin[k]) for k in range(3)) for i in range(3)]
    q = [[math.fsum(axes[k][i]*centered[k][r]*axes[r][j] for k in range(3) for r in range(3))
          for j in range(3)] for i in range(3)]
    c,s = math.cos(phi),math.sin(phi)
    return [*l,q[0][2],q[1][2],c*s*(q[0][0]-q[1][1])+(c*c-s*s)*q[0][1]]


def analyze(path):
    path = Path(path)
    # These historical G4/G5 receipts are small; reject large trajectory files
    # rather than accidentally materialize a G7 full trace in memory.
    if path.stat().st_size > 64*1024*1024:
        raise ValueError("RESOURCE_LIMIT: this bounded receipt reader accepts at most 64 MiB")
    data = path.read_bytes()
    record = json.loads(data)
    native = record["native_run"]
    solve, settings = native["solve"], native["settings"]
    controls = settings["solver_settings"]
    if controls["objective"] != "YOYO":
        raise ValueError("diagnostic currently scoped to the frozen Cuboid yo-yo cases")
    density, side = settings["density"], settings["source_settings"]["root_side"]
    leaves, relaxed, binary = solve["final_partition"], solve["relaxed_fill"], solve["binary_fill"]
    result = decompose(leaves, relaxed, binary, density)
    base = [density*v for v in native["source"]["outer_surface_moments"]]
    def moments(fill):
        return [base[j]-density*math.fsum(v*cell["moments"][j] for cell,v in zip(leaves,fill)) for j in range(10)]
    rm,bm = moments(relaxed),moments(binary)
    degrees = [0,1,1,1,2,2,2,2,2,2]
    error = max(abs(a-b)/(density*side**(3+d)) for a,b,d in zip(bm,solve["binary_moments"],degrees))
    if error > 1e-10:
        raise ValueError("stored partition does not reproduce stored binary moments")
    rc = constraints(rm,controls["target_origin"],controls["target_axes"],solve["final_phi"])
    bc = constraints(bm,controls["target_origin"],controls["target_axes"],solve["final_phi"])
    direct = constraints(result["moment_delta"],controls["target_origin"],controls["target_axes"],solve["final_phi"])
    constraint_error = max(abs(a-b)/(density*side**(4 if i<3 else 5))
                           for i,(a,b) in enumerate(zip(bc,solve["binary_check"]["constraints"])))
    if constraint_error > 1e-10:
        raise ValueError("stored binary constraints are inconsistent with partition algebra")
    for group in result["by_level"].values():
        group["constraint_delta"] = constraints(group["moment_delta"],controls["target_origin"],controls["target_axes"],solve["final_phi"])
    for cell in result["cells"]:
        cell["constraint_delta"] = constraints(cell["moment_delta"],controls["target_origin"],controls["target_axes"],solve["final_phi"])
    return {"schema":"stored-rounding-breakdown-v1", "scope":"HISTORICAL_STORED_PARTITION_ALGEBRA_ONLY",
            "source_sha256":hashlib.sha256(data).hexdigest(), "source_status":solve["status"],
            "geometry_reintegrated":False, "new_solver_run":False, "infeasibility_certified":False,
            "reassembled_binary_moment_error":error, "reassembled_constraint_error":constraint_error,
            "relaxed_constraints":rc, "binary_constraints":bc, "constraint_delta":direct,
            "decomposition_identity_error":max(abs((b-a)-d) for a,b,d in zip(rc,bc,direct)),
            "normalized_binary_constraint_residual":max(abs(v)/(density*side**(4 if i<3 else 5)) for i,v in enumerate(bc)),
            **result}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input",type=Path)
    parser.add_argument("output",type=Path)
    args = parser.parse_args()
    report = analyze(args.input)
    with args.output.open("x") as stream:
        json.dump(report,stream,allow_nan=False,sort_keys=True,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k not in ("cells","by_level","moment_delta")},allow_nan=False))
