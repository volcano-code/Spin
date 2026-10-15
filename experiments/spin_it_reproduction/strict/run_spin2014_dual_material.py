#!/usr/bin/env python3
"""Capture the native Spin-It 2014 dual-material regression as strict evidence.

The fixture is deliberately synthetic.  Its native regression checks the exact
equal-density reduction, two distinct densities, material labels, and the sum
of the separately reported ten-moment contributions.  It is not a Sumone
asset, and this runner never upgrades it to an original-paper experiment.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
from dataclasses import asdict
import math
from pathlib import Path
import subprocess
import time
import traceback

from make_it_stand.verification.run_evidence import (
    EXPLICIT_EFFECTIVE_SOURCE_POLICY, begin_run, check_freshness, finish_run,
)
from make_it_stand.verification.spin_common import canonical, native, parse
from make_it_stand.verification.mass_moments import TriangleMesh, compute_mass_properties


ROOT = Path(__file__).resolve().parents[3]
_BINARY = "spin-it-2014-dual-material-regression"
_MARKERS = (
    "[PASS] scenario=spin-it-2014-dual-material mass=",
    "paper_givens=validated",
    "independent_label_moments=validated",
)


def _write(path: Path, value: object) -> None:
    with path.open("x", encoding="ascii") as stream:
        stream.write(canonical(value) + "\n")


def _label_moments(manifest: dict, label: int, density: float) -> list[float]:
    """Analytically integrate labeled cuboids without using solver moments."""
    dimensions = manifest["grid_dimensions"]
    minimum = manifest["grid_minimum"]
    maximum = manifest["grid_maximum"]
    labels = manifest["labels"]
    spacing = [(maximum[axis] - minimum[axis]) / dimensions[axis]
               for axis in range(3)]
    volume = spacing[0] * spacing[1] * spacing[2]
    result = [0.] * 10
    for index, value in enumerate(labels):
        if value != label:
            continue
        x = index % dimensions[0]
        y = (index // dimensions[0]) % dimensions[1]
        z = index // (dimensions[0] * dimensions[1])
        lower = [minimum[0] + x * spacing[0], minimum[1] + y * spacing[1],
                 minimum[2] + z * spacing[2]]
        upper = [lower[axis] + spacing[axis] for axis in range(3)]
        center = [(lower[axis] + upper[axis]) * .5 for axis in range(3)]
        second = [(lower[axis] ** 2 + lower[axis] * upper[axis] + upper[axis] ** 2) / 3.
                  for axis in range(3)]
        result[0] += density * volume
        for axis in range(3):
            result[1 + axis] += density * volume * center[axis]
            result[4 + axis] += density * volume * second[axis]
        result[7] += density * volume * center[1] * center[2]
        result[8] += density * volume * center[0] * center[2]
        result[9] += density * volume * center[0] * center[1]
    return result


def _constraint_values(manifest: dict, moments: list[float]) -> list[float]:
    """Independent Spin-It constraints for the exported target-frame identity cases."""
    frame = manifest["target_frame"]
    if (frame["origin"] != [0., 0., 0.] or
            frame["axes"] != [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]):
        raise ValueError("dual-material independent constraint checker only admits identity target frames")
    mass, lx, ly, lz, qxx, qyy, qzz, qyz, qxz, qxy = moments
    if not mass > 0.:
        raise ValueError("nonpositive exported material mass")
    if manifest["mode"] == "STATIC_BALANCE":
        return [lx, ly]
    if manifest["mode"] == "SPINNING_TOP":
        return [lx, ly, -qxz + lx * lz / mass, -qyz + ly * lz / mass]
    raise ValueError("unsupported dual-material mode for independent constraint checker")


def _validate_partition(manifest: dict) -> dict:
    labels = manifest["labels"]
    material = manifest["input_material"]
    carveable = manifest["input_carveable"]
    if len(material) != len(labels) or len(carveable) != len(labels):
        raise ValueError("input masks do not match actual material labels")
    if any(type(value) is not int or value not in (0, 1) for value in material + carveable):
        raise ValueError("input masks must be binary")
    counts = {label: labels.count(label) for label in (0, 1, 2)}
    fixed = removed = 0
    for index, label in enumerate(labels):
        if label and not material[index]:
            raise ValueError("retained material lies outside the input material domain")
        if label == 2 and not carveable[index]:
            raise ValueError("dense insert lies outside the carveable design domain")
        if material[index] and not carveable[index]:
            fixed += 1
            if label != 1:
                raise ValueError("fixed shell cell was removed or relabeled")
        if material[index] and carveable[index] and label == 0:
            removed += 1
    expected = manifest["partition"]
    observed = {
        "void_cells": counts[0], "boundary_cells": counts[1], "interior_cells": counts[2],
        "fixed_shell_cells": fixed, "removed_carveable_cells": removed,
    }
    if expected != observed:
        raise ValueError("manifest partition counts disagree with independently inspected labels")
    dimensions = manifest["grid_dimensions"]
    minimum, maximum = manifest["grid_minimum"], manifest["grid_maximum"]
    spacing = [(maximum[axis] - minimum[axis]) / dimensions[axis]
               for axis in range(3)]
    shared_interface_area = 0.0
    for z in range(dimensions[2]):
        for y in range(dimensions[1]):
            for x in range(dimensions[0]):
                index = (z * dimensions[1] + y) * dimensions[0] + x
                if labels[index] != 1:
                    continue
                if x + 1 < dimensions[0] and labels[index + 1] == 2:
                    shared_interface_area += spacing[1] * spacing[2]
                if y + 1 < dimensions[1] and labels[index + dimensions[0]] == 2:
                    shared_interface_area += spacing[0] * spacing[2]
                if z + 1 < dimensions[2] and labels[index + dimensions[0] * dimensions[1]] == 2:
                    shared_interface_area += spacing[0] * spacing[1]
    return {"label_partition_status": "DISJOINT_VOXEL_PARTITION_VALIDATED",
            "positive_volume_overlap": "EXCLUDED_BY_DISJOINT_LABEL_PARTITION",
            "shared_material_interface_area": shared_interface_area,
            **observed}


def _validate_region_mesh_matches_labels(mesh: TriangleMesh, manifest: dict, label: int) -> None:
    """Require each exported triangle to be exactly a labeled-voxel boundary face.

    The expected face multiset is independently enumerated from the material
    labels. Triangle winding and face order are ignored here; separate closed
    manifold and positive signed-volume checks remain in force.
    """
    dimensions = manifest["grid_dimensions"]
    minimum, maximum = manifest["grid_minimum"], manifest["grid_maximum"]
    spacing = [(maximum[axis] - minimum[axis]) / dimensions[axis]
               for axis in range(3)]
    if any(value <= 0.0 or not math.isfinite(value) for value in spacing):
        raise ValueError("invalid voxel spacing for label-derived boundary audit")
    labels = manifest["labels"]
    offsets = ((-1, 0, 0), (1, 0, 0), (0, -1, 0),
               (0, 1, 0), (0, 0, -1), (0, 0, 1))
    face_corners = (
        ((0,0,0),(0,0,1),(0,1,0),(0,1,1)),
        ((1,0,0),(1,1,0),(1,0,1),(1,1,1)),
        ((0,0,0),(1,0,0),(0,0,1),(1,0,1)),
        ((0,1,0),(0,1,1),(1,1,0),(1,1,1)),
        ((0,0,0),(0,1,0),(1,0,0),(1,1,0)),
        ((0,0,1),(1,0,1),(0,1,1),(1,1,1)),
    )

    def face_key(points):
        # Sorting vertices makes this a geometric (not OBJ-index or winding)
        # comparison while retaining each triangle's exact grid-corner set.
        return tuple(sorted(points))

    expected = Counter()
    for z in range(dimensions[2]):
        for y in range(dimensions[1]):
            for x in range(dimensions[0]):
                index = (z * dimensions[1] + y) * dimensions[0] + x
                if labels[index] != label:
                    continue
                for direction, (dx, dy, dz) in enumerate(offsets):
                    nx, ny, nz = x + dx, y + dy, z + dz
                    if (0 <= nx < dimensions[0] and
                            0 <= ny < dimensions[1] and
                            0 <= nz < dimensions[2]):
                        neighbor = (nz * dimensions[1] + ny) * dimensions[0] + nx
                        if labels[neighbor] == label:
                            continue
                    corners = face_corners[direction]
                    points = [
                        (x + corner[0], y + corner[1], z + corner[2])
                        for corner in corners
                    ]
                    expected[face_key((points[0], points[1], points[2]))] += 1
                    expected[face_key((points[2], points[1], points[3]))] += 1

    observed = Counter()
    grid_vertices = []
    for point in mesh.vertices:
        grid_point = []
        for axis in range(3):
            coordinate = (point[axis] - minimum[axis]) / spacing[axis]
            nearest = round(coordinate)
            if abs(coordinate - nearest) > 1e-9 or not 0 <= nearest <= dimensions[axis]:
                raise ValueError("exported material vertex is not on the voxel grid")
            grid_point.append(nearest)
        grid_vertices.append(tuple(grid_point))
    for face in mesh.faces:
        observed[face_key(tuple(grid_vertices[index] for index in face))] += 1
    if observed != expected:
        missing = sum((expected - observed).values())
        extra = sum((observed - expected).values())
        raise ValueError(
            f"exported label {label} mesh does not match label-derived voxel surface: "
            f"missing_triangles={missing} extra_triangles={extra}"
        )


def _validate_mesh(mesh: TriangleMesh, manifest: dict) -> None:
    if not mesh.vertices or not mesh.faces:
        raise ValueError("nonempty material label lacks a closed boundary mesh")
    minimum, maximum = manifest["grid_minimum"], manifest["grid_maximum"]
    edges = {}
    for face in mesh.faces:
        if len(face) != 3 or len(set(face)) != 3 or any(index < 0 or index >= len(mesh.vertices) for index in face):
            raise ValueError("invalid triangle in exported material boundary")
        for index in face:
            if not all(math.isfinite(value) for value in mesh.vertices[index]):
                raise ValueError("nonfinite exported material vertex")
            if any(mesh.vertices[index][axis] < minimum[axis] - 1e-12 or
                   mesh.vertices[index][axis] > maximum[axis] + 1e-12
                   for axis in range(3)):
                raise ValueError("exported material vertex lies outside input domain bounds")
        for a, b in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            edge = tuple(sorted((a, b)))
            edges[edge] = edges.get(edge, 0) + 1
    if any(count != 2 for count in edges.values()):
        raise ValueError("exported material boundary is not a closed two-manifold")


def check_exported_geometry(directory: Path) -> dict:
    """Reintegrate each actual exported region; no seal or optimization grant."""
    result = {}
    for case in ("dense_regression", "sumone_reconstructed"):
        root = Path(directory) / case
        manifest = parse((root / "manifest.json").read_text())
        labels = manifest["labels"]
        if len(labels) != math.prod(manifest["grid_dimensions"]) or any(
                type(label) is not int or label not in (0, 1, 2) for label in labels):
            raise ValueError("invalid actual material labels")
        partition = _validate_partition(manifest)
        combined = [0.] * 10
        regions = []
        materials = manifest["materials"]
        if [entry["label"] for entry in materials] != [1, 2]:
            raise ValueError("missing or duplicate material")
        for entry in materials:
            expected = entry["moments"]
            if len(expected) != 10 or any(type(v) not in (int, float) or not math.isfinite(v) for v in expected):
                raise ValueError("invalid material moment claim")
            if labels.count(entry["label"]) != entry["cell_count"]:
                raise ValueError("material count differs from actual label field")
            if entry["cell_count"] == 0:
                if entry["mesh"] is not None or max(map(abs, expected)) > 1e-10:
                    raise ValueError("empty material has a nonempty geometric contribution")
                regions.append({"label": entry["label"], "status": "EXPLICIT_EMPTY",
                                "native_roundoff_moments": expected})
                continue
            path = root / entry["mesh"]
            if path.resolve().parent != root.resolve():
                raise ValueError("material path escapes export directory")
            vertices, faces = [], []
            with path.open() as stream:
                for line in stream:
                    fields = line.split()
                    if fields and fields[0] == "v":
                        vertices.append(tuple(map(float, fields[1:])))
                    elif fields and fields[0] == "f":
                        faces.append(tuple(int(v) - 1 for v in fields[1:]))
            mesh = TriangleMesh(vertices, faces)
            _validate_mesh(mesh, manifest)
            _validate_region_mesh_matches_labels(mesh, manifest, entry["label"])
            properties = compute_mass_properties(mesh, entry["density"])
            q = properties.second_moment
            moments = [properties.mass, *properties.first_moment, q[0][0], q[1][1], q[2][2], q[1][2], q[0][2], q[0][1]]
            error = max(abs(a-b) for a,b in zip(moments, expected))
            if error > 1e-10:
                raise ValueError("exported material moments differ from solver")
            label_moments = _label_moments(manifest, entry["label"], entry["density"])
            label_error = max(abs(a-b) for a,b in zip(moments, label_moments))
            if label_error > 1e-10:
                raise ValueError("exported material geometry differs from independent label integration")
            combined = [a+b for a,b in zip(combined,moments)]
            regions.append({"label": entry["label"], "status": "INTEGRATED",
                            "label_surface_status": "EXACT_TRIANGLE_MULTISET_MATCH",
                            "label_surface_triangle_count": len(mesh.faces),
                            "moment_error": error, "label_moment_error": label_error,
                            "mass_properties": asdict(properties),
                            "mesh_sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        total = manifest["total_moments"]
        if len(total) != 10 or any(type(v) not in (int,float) or not math.isfinite(v) for v in total):
            raise ValueError("invalid combined moment claim")
        error = max(abs(a-b) for a,b in zip(combined,total))
        if error > 1e-10:
            raise ValueError("combined exported geometry differs from solver")
        constraints = _constraint_values(manifest, combined)
        claimed_constraints = manifest["binary_constraints"]
        constraint_error = max(abs(a-b) for a,b in zip(constraints, claimed_constraints))
        if len(constraints) != len(claimed_constraints) or constraint_error > 1e-10:
            raise ValueError("independent exported constraints differ from solver")
        constraint_inf = max(map(abs, constraints), default=0.)
        if abs(constraint_inf - manifest["binary_constraint_inf"]) > 1e-10:
            raise ValueError("constraint infinity norm differs from independently recomputed constraints")
        result[case] = {"status": "TWO_MATERIAL_GEOMETRY_INTEGRATED" if labels.count(2) else "NO_RETAINED_DENSE_INSERT",
                       "combined_moment_error": error, "partition": partition, "regions": regions,
                       "constraints": {"recomputation_error": constraint_error,
                                       "infinity_norm": constraint_inf,
                                       "current_1e-6_status": "PASS" if constraint_inf <= 1e-6 else "FAIL",
                                       "legacy_tolerance": manifest["binary_constraint_tolerance"]},
                       "scope": "independent boundary moments and manifoldness only; no original-case, full geometry, or seal grant"}
    return result


def run(
    output_directory: Path, *, build_directory: Path, dependency_roots: dict[str, Path],
    export_geometry: bool = False,
) -> dict[str, object]:
    output = output_directory.absolute()
    output.mkdir(parents=False, exist_ok=False)
    build = build_directory.resolve()
    executable = build / _BINARY
    if not executable.is_file() or not executable.stat().st_mode & 0o111:
        raise ValueError("required native dual-material regression executable is unavailable")

    started, tick, begin, events = (
        datetime.now(timezone.utc).isoformat(), time.perf_counter(), None, []
    )
    try:
        sumone_asset = (ROOT / "experiments/spin_it_reproduction/spin2014/paper-assets/"
                        "sumone-dual-density.vox.json").resolve(strict=True)
        inputs = parse(canonical({
            "schema": "spin2014-dual-material-inputs-v1",
            "classification": "SYNTHETIC_REGRESSION",
            "binary": {
                "path": str(executable),
                "sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
            },
            "sumone_asset": {
                "path": str(sumone_asset),
                "sha256": hashlib.sha256(sumone_asset.read_bytes()).hexdigest(),
            },
            "claims": [
                "equal_density_reduces_to_single_material",
                "separate_density_weighted_moments_sum_to_final_moments",
                "material_labels_are_in_0_1_2",
                "paper_givens_dual_material_path_is_exercised",
            ],
            "export_geometry": export_geometry,
        }))
        input_path = output / "inputs.json"
        _write(input_path, inputs)
        extension = native()
        if Path(extension.__file__).resolve().parents[1] != build:
            raise ValueError("loaded extension differs from declared build directory")
        begin = begin_run(
            repository_root=ROOT,
            runner_files=[Path(__file__)],
            checker_files=[
                ROOT / "tests/spin_it_2014_dual_material_regression.cpp",
                ROOT / "src/spin2014/DualMaterialSpinSolver.cpp",
                ROOT / "include/spin2014/DualMaterialSpinSolver.h",
                ROOT / "experiments/spin_it_reproduction/strict/dual_material_geometry_export.h",
                ROOT / "src/voxel/MaterialOccupancy.cpp",
                ROOT / "python/make_it_stand/verification/mass_moments.py",
                ROOT / "python/make_it_stand/verification/run_evidence.py",
                ROOT / "python/make_it_stand/verification/spin_common.py",
            ],
            asset_files=[input_path, executable, sumone_asset],
            build_configuration_files=[ROOT / "CMakeLists.txt", build / "CMakeCache.txt"],
            parameters=inputs,
            seed_protocol={"seed": 0, "protocol": "deterministic native dual-material regression"},
            dependency_roots=dependency_roots,
            native_modules=[extension],
            pure_python=False,
            source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY,
        )
        _write(output / "begin.json", begin.as_dict())
        completed = subprocess.run(
            [str(executable)] + (["--export-directory", str(output / "materials")] if export_geometry else []),
            cwd=ROOT, capture_output=True, text=True, check=False
        )
        transcript = completed.stdout + completed.stderr
        transcript_path = output / "dual-material.stdout"
        with transcript_path.open("x", encoding="utf-8") as stream:
            stream.write(transcript)
        if completed.returncode != 0 or any(marker not in transcript for marker in _MARKERS):
            raise RuntimeError("native dual-material regression failed")
        if export_geometry:
            _write(output / "independent_geometry.json", check_exported_geometry(output / "materials"))
        result = {
            "path": str(transcript_path),
            "sha256": hashlib.sha256(transcript.encode()).hexdigest(),
            "returncode": completed.returncode,
            "required_markers": list(_MARKERS),
        }
        scientific = {
            "schema": "spin2014-dual-material-scientific-v1",
            "classification": "SYNTHETIC_REGRESSION",
            "native_regression": result,
        }
        _write(output / "scientific.json", scientific)
        events.append({
            "iteration": 0,
            "status": "ok",
            "data": {"kind": "DUAL_MATERIAL_NATIVE_REGRESSION", **result},
        })
        evidence = finish_run(
            begin,
            events=events,
            output_files=sorted(
                path for path in output.rglob("*")
                if path.is_file() and path.name != "begin.json"
            ),
            measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick},
        )
        _write(output / "evidence.json", evidence.as_dict())
        freshness = check_freshness(evidence)
        summary = {
            "status": "SYNTHETIC_DUAL_MATERIAL_VALIDATED",
            "classification": "SYNTHETIC_REGRESSION",
            "project_status": "PARTIAL_BLOCKED",
            "original_evidence_admission": "BLOCKED_UNKNOWN",
            "paper_case_admission": "BLOCKED_MISSING_AUTHOR_SUMONE_ASSET_AND_PARAMETERS",
            "freshness": {"status": freshness.status, "mismatches": list(freshness.mismatches)},
        }
        _write(output / "summary.json", summary)
        return summary
    except Exception as error:
        if begin is not None and not (output / "evidence.json").exists():
            events.append({"iteration": len(events), "status": "failed", "data": {"error": str(error)}})
            _write(output / "evidence.json", finish_run(
                begin,
                events=events,
                output_files=sorted(
                    path for path in output.rglob("*")
                    if path.is_file() and path.name != "begin.json"
                ),
                measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick},
            ).as_dict())
        _write(output / "failure.json", {
            "status": "FAILED", "error": str(error), "traceback": traceback.format_exc(),
        })
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--build-directory", type=Path, required=True)
    parser.add_argument("--dependency-root", action="append", required=True, metavar="NAME=PATH")
    parser.add_argument("--export-geometry", action="store_true")
    args = parser.parse_args()
    roots: dict[str, Path] = {}
    for item in args.dependency_root:
        name, separator, value = item.partition("=")
        if not name or not separator or name in roots:
            parser.error("invalid or duplicate dependency root")
        roots[name] = Path(value).resolve()
    result = run(args.output_directory, build_directory=args.build_directory, dependency_roots=roots,
                 export_geometry=args.export_geometry)
    print(canonical(result))
    return 0 if result["freshness"]["status"] == "CURRENT" else 1


if __name__ == "__main__":
    raise SystemExit(main())
