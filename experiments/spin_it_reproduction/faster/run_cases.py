#!/usr/bin/env python3
"""Run the reconstructed complex-OBJ Spin-It Faster matrix twice."""

from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import time
from typing import Any

from make_it_stand.spin import (
    SpinItFasterRequest,
    SpinSolverConfig,
    evaluate_spin_it_faster_quadric,
    integrate_raw_mass_moments,
    solve_spin_it_faster,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_obj(path: Path):
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        fields = line.split()
        if not fields or fields[0] == "#":
            continue
        if fields[0] == "v" and len(fields) == 4:
            vertices.append(tuple(float(value) for value in fields[1:4]))
        elif fields[0] == "f" and len(fields) == 4:
            faces.append(tuple(int(value.split("/", 1)[0]) - 1 for value in fields[1:4]))
        else:
            raise ValueError(f"{path}:{line_number}: unsupported OBJ record")
    return tuple(vertices), tuple(faces)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def git(repository: Path, *arguments: str) -> str:
    return subprocess.run(
        ["git", *arguments], cwd=repository, check=True,
        capture_output=True, text=True,
    ).stdout.strip()


def run_case(case: dict[str, Any], root: Path, design, output: Path) -> dict[str, Any]:
    source = root / case["source"]
    if sha256(source) != case["source_sha256"]:
        raise RuntimeError(f"{case['id']}: source hash mismatch")
    outer_vertices, outer_faces = read_obj(source)
    design_vertices, design_faces = design
    outer_moments = integrate_raw_mass_moments(outer_vertices, outer_faces)
    design_moments = integrate_raw_mass_moments(design_vertices, design_faces)
    fixed = tuple(outer - variable for outer, variable in zip(outer_moments, design_moments))
    request = SpinItFasterRequest(
        outer_vertices=outer_vertices,
        outer_faces=outer_faces,
        design_vertices=design_vertices,
        design_faces=design_faces,
        fixed_shell_moments=fixed,
        target_origin=tuple(case["target_origin"]),
        mode=case["mode"],
        evidence="RECONSTRUCTED_ASSET",
        minimum_wall_thickness=float(case["minimum_wall_thickness"]),
        center_of_mass_weight=float(case.get("center_of_mass_weight", 0.0)),
        settings=SpinSolverConfig(maximum_restarts=1),
    )
    initial = tuple(float(value) for value in case["initial_quadric"])
    repeats = []
    for repeat in (1, 2):
        start = time.perf_counter()
        result = solve_spin_it_faster(request, initial_quadrics=(initial,))
        elapsed = time.perf_counter() - start
        repeat_root = output / case["id"] / f"repeat-{repeat}"
        repeat_root.mkdir(parents=True, exist_ok=True)
        result.save_outer_mesh(repeat_root / "outer.obj")
        result.save_material_mesh(repeat_root / "material.obj")
        result.save_cavity_mesh(repeat_root / "cavity.obj")
        result.save_diagnostics(repeat_root / "diagnostics.json")
        diagnostics = result.diagnostics
        interface = evaluate_spin_it_faster_quadric(
            request, diagnostics["quadric_coefficients"]
        ).interface_type
        repeats.append({
            "repeat": repeat,
            "elapsed_seconds": elapsed,
            "status": result.status,
            "evidence": result.evidence,
            "geometry_sealed": result.geometry_sealed,
            "input_hash": diagnostics["input_hash"],
            "geometry_hash": diagnostics["geometry_hash"],
            "trace_hash": diagnostics["trace_hash"],
            "output_hash": diagnostics["output_hash"],
            "diagnostics_sha256": sha256(repeat_root / "diagnostics.json"),
            "normalized_equation_residual": diagnostics["normalized_equation_residual"],
            "physical_relative_residual": diagnostics["physical_relative_residual"],
            "interface_type": interface,
            "restart_count": 1 + max(entry.restart for entry in result.trace),
            "function_evaluations": max(entry.function_evaluations for entry in result.trace),
            "geometry": diagnostics["geometry"],
        })
    deterministic = all(
        repeats[0][field] == repeats[1][field]
        for field in ("status", "input_hash", "geometry_hash", "trace_hash",
                      "output_hash", "diagnostics_sha256", "interface_type")
    )
    passed = all(
        repeat["status"] == "CONVERGED"
        and repeat["geometry_sealed"]
        and repeat["normalized_equation_residual"] <= 1e-12
        and repeat["physical_relative_residual"] <= 1e-6
        for repeat in repeats
    ) and deterministic
    return {
        "id": case["id"],
        "mode": case["mode"],
        "evidence": "RECONSTRUCTED_ASSET",
        "source": case["source"],
        "source_sha256": case["source_sha256"],
        "outer_vertex_count": len(outer_vertices),
        "outer_face_count": len(outer_faces),
        "fixed_shell_moments": fixed,
        "deterministic": deterministic,
        "passed": passed,
        "repeats": repeats,
    }


def report(summary: dict[str, Any]) -> str:
    lines = [
        "# Spin-It Faster reconstructed OBJ report",
        "",
        "All inputs are deterministic, license-free reconstructed assets. This is a",
        "`RECONSTRUCTED_ASSET` formula reproduction, not author-CAD or author-code equivalence.",
        "Fixed-shell moments are integrated as `outer - design domain`; the visible",
        "outer surface is unchanged and both material and cavity boundaries are audited.",
        "",
        "| Case | Mode | Faces | Interface | Eq. 14 residual | Physical residual | Repeat | Result |",
        "|---|---|---:|---|---:|---:|---|---|",
    ]
    for case in summary["cases"]:
        first = case["repeats"][0]
        lines.append(
            f"| {case['id']} | `{case['mode']}` | {case['outer_face_count']} | "
            f"`{first['interface_type']}` | {first['normalized_equation_residual']:.3e} | "
            f"{first['physical_relative_residual']:.3e} | "
            f"{'identical' if case['deterministic'] else 'DIFF'} | "
            f"{'PASS' if case['passed'] else 'FAIL'} |"
        )
    lines.extend((
        "",
        f"Release gate: **{'PASS' if summary['all_cases_passed'] else 'BLOCKED'}**.",
        "Wall-clock durations are recorded but excluded from deterministic hashes.",
        "",
    ))
    return "\n".join(lines)


def main() -> int:
    root = Path(__file__).resolve().parent
    repository = root.parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=root / "cases.json")
    parser.add_argument("--output", type=Path, default=root / "results")
    arguments = parser.parse_args()
    manifest = json.loads(arguments.manifest.read_text(encoding="utf-8"))
    design_path = root / manifest["design_domain"]
    if sha256(design_path) != manifest["design_domain_sha256"]:
        raise RuntimeError("design-domain hash mismatch")
    design = read_obj(design_path)
    worktree_status = git(repository, "status", "--porcelain=v1", "--untracked-files=all")
    module_spec = importlib.util.find_spec("_make_it_stand")
    if module_spec is None or module_spec.origin is None:
        raise RuntimeError("cannot resolve the native extension binary")
    cases = [run_case(case, root, design, arguments.output) for case in manifest["cases"]]
    summary = {
        "schema": "spin-it-faster-complex-obj-results-v1",
        "schema_version": 1,
        "claim_level": "RECONSTRUCTED_ASSET",
        "git_commit": git(repository, "rev-parse", "HEAD"),
        "worktree_clean_at_start": not bool(worktree_status),
        "worktree_status_sha256": hashlib.sha256(worktree_status.encode()).hexdigest(),
        "manifest_sha256": sha256(arguments.manifest),
        "native_extension_sha256": sha256(Path(module_spec.origin)),
        "all_cases_passed": all(case["passed"] for case in cases),
        "cases": cases,
    }
    write_json(arguments.output / "summary.json", summary)
    (root / "REPORT.md").write_text(report(summary), encoding="utf-8")
    gate = {
        "schema": "spin-it-faster-release-gate-v1",
        "schema_version": 1,
        "claim_level": "RECONSTRUCTED_ASSET",
        "passed": summary["all_cases_passed"] and summary["worktree_clean_at_start"],
        "requirements": {
            "four_modes_present": len(cases) == 4,
            "all_converged_and_sealed": all(case["passed"] for case in cases),
            "byte_identical_scientific_repeats": all(case["deterministic"] for case in cases),
            "clean_source_at_run_start": summary["worktree_clean_at_start"],
        },
        "summary_sha256": sha256(arguments.output / "summary.json"),
    }
    write_json(root / "release-gate.json", gate)
    print(report(summary), end="")
    return 0 if gate["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
