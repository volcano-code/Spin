#!/usr/bin/env python3
"""Capture independent mass-moment evidence on analytic and composite meshes.

The checker deliberately compares the independent signed-tetrahedron
integrator against closed-form fixtures written in this runner.  It does not
call the native solver's moment routine and it does not claim an original
Spin-It/Faster paper case.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import time
import traceback

from make_it_stand.verification import mass_moments
from make_it_stand.verification.mass_moments import DensityRegion, TriangleMesh
from make_it_stand.verification.run_evidence import (
    EXPLICIT_EFFECTIVE_SOURCE_POLICY, begin_run, check_freshness, finish_run,
)
from make_it_stand.verification.spin_common import canonical, parse


ROOT = Path(__file__).resolve().parents[3]
_SCOPE = (
    "SYNTHETIC_ANALYTIC_AND_COMPOSITE_INDEPENDENT_MASS_FIRST_SECOND_AND_"
    "CENTRAL_INERTIA_AUDIT_ONLY_NOT_AUTHOR_CASE"
)
_TOLERANCE = 5e-13


def _write(path: Path, value: object) -> None:
    with path.open("x", encoding="ascii") as stream:
        stream.write(canonical(value) + "\n")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _cube(origin=(0.0, 0.0, 0.0), side=1.0) -> TriangleMesh:
    x, y, z = origin
    s = float(side)
    vertices = tuple((x + dx * s, y + dy * s, z + dz * s)
                     for dz in (0.0, 1.0) for dy in (0.0, 1.0) for dx in (0.0, 1.0))
    # Index = dx + 2*dy + 4*dz.
    faces = (
        (0, 3, 1), (0, 2, 3), (4, 5, 7), (4, 7, 6),
        (0, 4, 6), (0, 6, 2), (1, 3, 7), (1, 7, 5),
        (0, 1, 5), (0, 5, 4), (2, 6, 7), (2, 7, 3),
    )
    return TriangleMesh(vertices, faces)


def _tetrahedron() -> TriangleMesh:
    vertices = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    faces = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))
    return TriangleMesh(vertices, faces)


def _matrix_add(*matrices):
    return tuple(tuple(sum(matrix[i][j] for matrix in matrices)
                       for j in range(3)) for i in range(3))


def _outer(vector):
    return tuple(tuple(vector[i] * vector[j] for j in range(3)) for i in range(3))


def _matrix_scale(matrix, factor):
    return tuple(tuple(factor * matrix[i][j] for j in range(3)) for i in range(3))


def _matrix_sub(left, right):
    return tuple(tuple(left[i][j] - right[i][j] for j in range(3)) for i in range(3))


def _inertia_from_centered(second):
    trace = sum(second[i][i] for i in range(3))
    return tuple(tuple((trace if i == j else 0.0) - second[i][j]
                       for j in range(3)) for i in range(3))


def _max_abs(left, right):
    return max(abs(float(a) - float(b)) for a, b in zip(left, right))


def _flat(matrix):
    return tuple(value for row in matrix for value in row)


def _cube_expected(origin, side, density):
    s = float(side)
    volume = s ** 3
    mass = float(density) * volume
    center = tuple(float(origin[i]) + s / 2.0 for i in range(3))
    # Keep row/column ordering explicit and auditable.
    centered_second = tuple(tuple(mass / 12.0 if i == j else 0.0
                                  for j in range(3)) for i in range(3))
    raw_second = _matrix_add(centered_second, _matrix_scale(_outer(center), mass))
    return {
        "volume": volume,
        "mass": mass,
        "first_moment": tuple(mass * value for value in center),
        "second_moment": raw_second,
        "center_of_mass": center,
        "central_inertia": _inertia_from_centered(centered_second),
    }


def _tetra_expected(density):
    d = float(density)
    volume = 1.0 / 6.0
    mass = d * volume
    first = (mass / 4.0, mass / 4.0, mass / 4.0)
    second = tuple(tuple(
        mass / 10.0 if i == j else mass / 20.0 for j in range(3)
    ) for i in range(3))
    center = (0.25, 0.25, 0.25)
    centered = _matrix_sub(second, _matrix_scale(_outer(center), mass))
    return {
        "volume": volume,
        "mass": mass,
        "first_moment": first,
        "second_moment": second,
        "center_of_mass": center,
        "central_inertia": _inertia_from_centered(centered),
    }


def _observe(mesh, density, expected):
    geometric = mass_moments.integrate_geometric_moments(mesh)
    properties = mass_moments.compute_mass_properties(mesh, density)
    observed = {
        "volume": geometric.volume,
        "mass": properties.mass,
        "first_moment": properties.first_moment,
        "second_moment": properties.second_moment,
        "center_of_mass": properties.center_of_mass,
        "central_inertia": properties.central_inertia,
    }
    errors = {
        key: _max_abs(_flat(value) if isinstance(value, tuple) and value and isinstance(value[0], tuple)
                      else (value,) if not isinstance(value, tuple) else value,
                      _flat(expected[key]) if isinstance(expected[key], tuple) and expected[key] and isinstance(expected[key][0], tuple)
                      else (expected[key],) if not isinstance(expected[key], tuple) else expected[key])
        for key, value in observed.items()
    }
    return observed, errors


def _aggregate_expected():
    a = _cube_expected((0.0, 0.0, 0.0), 1.0, 1.0)
    b = _cube_expected((2.0, 1.0, 3.0), 1.0, 2.0)
    mass = a["mass"] + b["mass"]
    first = tuple(a["first_moment"][i] + b["first_moment"][i] for i in range(3))
    center = tuple(value / mass for value in first)
    second = _matrix_add(a["second_moment"], b["second_moment"])
    centered = _matrix_sub(second, _matrix_scale(_outer(center), mass))
    return {
        "mass": mass,
        "first_moment": first,
        "second_moment": second,
        "center_of_mass": center,
        "central_inertia": _inertia_from_centered(centered),
    }


def run(output_directory: Path, *, dependency_roots: dict[str, Path]) -> dict[str, object]:
    output = output_directory.absolute()
    output.mkdir(parents=False, exist_ok=False)
    started, tick, begin, events = datetime.now(timezone.utc).isoformat(), time.perf_counter(), None, []
    try:
        inputs = parse(canonical({
            "schema": "shared-independent-moments-inputs-v1",
            "classification": "SYNTHETIC_REGRESSION",
            "fixtures": ["unit_cube_density_2", "unit_tetrahedron_density_3",
                         "two_translated_cubes_non_diagonal_composite"],
            "tolerance": _TOLERANCE,
            "integrator": "python/make_it_stand/verification/mass_moments.py",
        }))
        input_path = output / "inputs.json"
        _write(input_path, inputs)
        begin = begin_run(
            repository_root=ROOT,
            runner_files=[Path(__file__)],
            checker_files=[ROOT / "python/make_it_stand/verification/mass_moments.py"],
            asset_files=[input_path],
            build_configuration_files=[ROOT / "CMakeLists.txt"], parameters=inputs,
            seed_protocol={"seed": 0, "protocol": "deterministic analytic closed-form fixtures"},
            dependency_roots=dependency_roots, native_modules=[], pure_python=True,
            source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY,
        )
        _write(output / "begin.json", begin.as_dict())

        fixtures = []
        observed, errors = _observe(_cube(), 2.0, _cube_expected((0, 0, 0), 1, 2))
        fixtures.append({"name": "unit_cube_density_2", "observed": observed,
                        "maximum_absolute_error": max(errors.values()), "errors": errors})
        observed, errors = _observe(_tetrahedron(), 3.0, _tetra_expected(3))
        fixtures.append({"name": "unit_tetrahedron_density_3", "observed": observed,
                        "maximum_absolute_error": max(errors.values()), "errors": errors})

        regions = (DensityRegion(_cube((0, 0, 0)), 1.0),
                   DensityRegion(_cube((2, 1, 3)), 2.0))
        aggregate = mass_moments.aggregate_mass_properties(regions)
        expected = _aggregate_expected()
        aggregate_observed = {
            "mass": aggregate.mass, "first_moment": aggregate.first_moment,
            "second_moment": aggregate.second_moment,
            "center_of_mass": aggregate.center_of_mass,
            "central_inertia": aggregate.central_inertia,
        }
        aggregate_errors = {
            key: _max_abs(_flat(value) if isinstance(value, tuple) and value and isinstance(value[0], tuple)
                          else value if isinstance(value, tuple) else (value,),
                          _flat(expected[key]) if isinstance(expected[key], tuple) and expected[key] and isinstance(expected[key][0], tuple)
                          else expected[key] if isinstance(expected[key], tuple) else (expected[key],))
            for key, value in aggregate_observed.items()
        }
        fixtures.append({"name": "two_translated_cubes_non_diagonal_composite",
                         "observed": aggregate_observed,
                         "maximum_absolute_error": max(aggregate_errors.values()),
                         "errors": aggregate_errors})
        maximum_error = max(row["maximum_absolute_error"] for row in fixtures)
        if not math.isfinite(maximum_error) or maximum_error > _TOLERANCE:
            raise RuntimeError(f"independent mass-moment fixtures exceed tolerance: {maximum_error}")
        scientific = {
            "schema": "shared-independent-moments-scientific-v1",
            "classification": "SYNTHETIC_REGRESSION", "scope": _SCOPE,
            "integrator": "signed_tetrahedron_python_independent",
            "fixtures": fixtures, "maximum_absolute_error": maximum_error,
            "tolerance": _TOLERANCE,
            "density_applied_once": True, "non_diagonal_inertia_fixture": True,
        }
        scientific_path = output / "scientific.json"
        _write(scientific_path, scientific)
        events.append({"iteration": 0, "status": "ok", "data": {
            "kind": "SHARED_INDEPENDENT_MASS_MOMENTS", "maximum_absolute_error": maximum_error,
        }})
        evidence = finish_run(begin, events=events,
            output_files=sorted(path for path in output.rglob("*")
                                if path.is_file() and path.name != "begin.json"),
            measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick})
        _write(output / "evidence.json", evidence.as_dict())
        freshness = check_freshness(evidence)
        summary = {
            "status": "SYNTHETIC_INDEPENDENT_MASS_MOMENTS_VALIDATED",
            "classification": "SYNTHETIC_REGRESSION", "project_status": "PARTIAL_BLOCKED",
            "original_evidence_admission": "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS",
            "paper_case_admission": "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS",
            "freshness": {"status": freshness.status, "mismatches": list(freshness.mismatches)},
        }
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
    parser.add_argument("--dependency-root", action="append", required=True, metavar="NAME=PATH")
    args = parser.parse_args()
    roots: dict[str, Path] = {}
    for item in args.dependency_root:
        name, separator, value = item.partition("=")
        if not name or not separator or name in roots:
            parser.error("invalid or duplicate dependency root")
        roots[name] = Path(value).resolve()
    result = run(args.output_directory, dependency_roots=roots)
    print(canonical(result))
    return 0 if result["freshness"]["status"] == "CURRENT" else 1


if __name__ == "__main__":
    raise SystemExit(main())
