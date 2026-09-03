#!/usr/bin/env python3
"""Run the public Spin-It Faster Table 1 protocol on reconstructed assets."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import random
import statistics
import time
from typing import Any


ROOT = Path(__file__).resolve().parent
DEFAULT_CASES = ROOT / "paper_cases.json"
TIMING_KEYS = {"elapsed_seconds", "residual_seconds", "derivative_seconds"}


def _without_timings(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: _without_timings(item)
            for key, item in value.items()
            if key not in TIMING_KEYS
        }
    if isinstance(value, list):
        return [_without_timings(item) for item in value]
    return value


def scientific_hash(payload: Any) -> str:
    encoded = json.dumps(
        _without_timings(payload),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def generate_initial_quadrics(
    mode: str,
    seed: int,
    count: int,
    domain_points: tuple[tuple[float, float, float], ...] | None = None,
) -> list[list[float]]:
    if count <= 0:
        raise ValueError("initial-quadric count must be positive")
    if mode not in {"STATIC_BALANCE", "BUOYANT_STABILITY", "FIDGET_SPINNER", "YO_YO", "SPINNING_TOP"}:
        raise ValueError(f"unsupported mode: {mode}")
    generator = random.Random(seed)
    if domain_points:
        lower = [min(point[axis] for point in domain_points) for axis in range(3)]
        upper = [max(point[axis] for point in domain_points) for axis in range(3)]
    else:
        lower = [-1.0, -1.0, -1.0]
        upper = [1.0, 1.0, 1.0]
    center = [(a + b) * 0.5 for a, b in zip(lower, upper)]
    half = [max((b - a) * 0.5, 1e-6) for a, b in zip(lower, upper)]
    rows: list[list[float]] = []
    while len(rows) < count:
        values = [0.0] * 10
        if mode in {"STATIC_BALANCE", "BUOYANT_STABILITY"}:
            normal = [generator.gauss(0.0, 1.0) for _ in range(3)]
            if mode == "STATIC_BALANCE":
                normal[2] = 0.0
            normal_norm = math.sqrt(sum(value * value for value in normal))
            if normal_norm <= 1e-15:
                continue
            normal = [value / normal_norm for value in normal]
            projected_half_width = sum(abs(normal[i]) * half[i] for i in range(3))
            signed_offset = generator.uniform(-0.65, 0.65) * projected_half_width
            values[0] = signed_offset - sum(normal[i] * center[i] for i in range(3))
            values[1:4] = normal
        elif mode == "FIDGET_SPINNER":
            cx = center[0] + generator.uniform(-0.2, 0.2) * half[0]
            cy = center[1] + generator.uniform(-0.2, 0.2) * half[1]
            radius = generator.uniform(0.25, 0.75) * min(half[0], half[1])
            values[0] = radius * radius - cx * cx - cy * cy
            values[1] = 2.0 * cx
            values[2] = 2.0 * cy
            values[4] = -1.0
            values[5] = -1.0
        elif mode == "YO_YO":
            radial_scale = 1.0 / max(half[0] * half[0], half[1] * half[1])
            axial_scale = generator.uniform(0.2, 1.2) / (half[2] * half[2])
            values[3] = generator.uniform(-0.8, 0.8) / half[2]
            values[4] = radial_scale
            values[5] = radial_scale
            values[6] = -axial_scale
        elif mode == "SPINNING_TOP":
            cx = center[0] + generator.uniform(-0.15, 0.15) * half[0]
            cy = center[1] + generator.uniform(-0.15, 0.15) * half[1]
            cz = center[2] + generator.uniform(-0.15, 0.15) * half[2]
            radius = generator.uniform(0.25, 0.75) * min(half)
            values[0] = radius * radius - cx * cx - cy * cy - cz * cz
            values[1:4] = [2.0 * cx, 2.0 * cy, 2.0 * cz]
            values[4:7] = [-1.0, -1.0, -1.0]
        norm = math.sqrt(sum(value * value for value in values))
        if norm <= 1e-15:
            continue
        rows.append([value / norm for value in values])
    return rows


def summarize_attempts(attempts: list[dict[str, Any]], expected_runs: int) -> dict[str, Any]:
    successful: list[dict[str, Any]] = []
    restarts: list[int] = []
    for run in range(expected_runs):
        run_attempts = [item for item in attempts if item["run"] == run]
        converged = next((item for item in run_attempts if item["status"] == "CONVERGED"), None)
        if converged is not None:
            successful.append(converged)
            restarts.append(int(converged["restart"]))
    def median(key: str) -> float | None:
        return statistics.median(float(item[key]) for item in successful) if successful else None
    return {
        "expected_runs": expected_runs,
        "successful_runs": len(successful),
        "all_runs_converged": len(successful) == expected_runs,
        "average_restarts": sum(restarts) / len(restarts) if restarts else None,
        "median_newton_iterations": median("newton_iterations"),
        "median_function_evaluations": median("function_evaluations"),
        "median_residual_seconds": median("residual_seconds"),
        "median_derivative_seconds": median("derivative_seconds"),
    }


def read_obj(path: Path) -> tuple[tuple[tuple[float, float, float], ...], tuple[tuple[int, int, int], ...]]:
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if fields and fields[0] == "v":
            vertices.append(tuple(float(value) for value in fields[1:4]))
        elif fields and fields[0] == "f":
            faces.append(tuple(int(value) - 1 for value in fields[1:4]))
    return tuple(vertices), tuple(faces)


def build_request(case: dict[str, Any]):
    from make_it_stand.spin import SpinItFasterRequest, SpinSolverConfig, integrate_raw_mass_moments

    outer = read_obj(ROOT / case["outer"]["path"])
    design = read_obj(ROOT / case["design_domain"]["path"])
    density = float(case.get("material_density", 1.0))
    outer_moments = integrate_raw_mass_moments(*outer, density=density)
    design_moments = integrate_raw_mass_moments(*design, density=density)
    fixed = tuple(a - b for a, b in zip(outer_moments, design_moments))
    buoyancy = case.get("buoyancy", {})
    return SpinItFasterRequest(
        outer_vertices=outer[0],
        outer_faces=outer[1],
        design_vertices=design[0],
        design_faces=design[1],
        fixed_shell_moments=fixed,
        target_origin=tuple(case["target_origin"]),
        target_axes=tuple(tuple(row) for row in case["target_axes"]),
        mode=case["mode"],
        static_objective=case.get("static_objective", "MAXIMUM_MASS"),
        evidence="RECONSTRUCTED_ASSET",
        material_density=density,
        minimum_wall_thickness=0.0,
        center_of_mass_weight=float(case.get("center_of_mass_weight", 0.0)),
        liquid_density=float(buoyancy.get("liquid_density", 0.0)),
        displaced_liquid_volume=float(buoyancy.get("displaced_liquid_volume", 0.0)),
        center_of_buoyancy_target=tuple(buoyancy.get("center_of_buoyancy_target", (0.0, 0.0, 0.0))),
        waterline_area=float(buoyancy.get("waterline_area", 0.0)),
        waterline_central_second_moment=tuple(
            tuple(row) for row in buoyancy.get(
                "waterline_central_second_moment",
                ((0.0, 0.0), (0.0, 0.0)),
            )
        ),
        fixed_shell_volume=fixed[0] / density if case["mode"] == "BUOYANT_STABILITY" else 0.0,
        settings=SpinSolverConfig(maximum_restarts=1),
    )


def run_case(case: dict[str, Any], starts: int, max_restarts: int) -> dict[str, Any]:
    from make_it_stand.spin import solve_spin_it_faster

    request = build_request(case)
    design = read_obj(ROOT / case["design_domain"]["path"])
    origin = tuple(float(value) for value in case["target_origin"])
    axes = tuple(tuple(float(value) for value in row) for row in case["target_axes"])
    domain_points = tuple(
        tuple(
            sum(axes[row][column] * (point[row] - origin[row]) for row in range(3))
            for column in range(3)
        )
        for point in design[0]
    )
    attempts: list[dict[str, Any]] = []
    for run in range(starts):
        initials = generate_initial_quadrics(
            case["mode"],
            seed=2024 + 1009 * run + 7919 * sum(map(ord, case["id"])),
            count=max_restarts,
            domain_points=domain_points,
        )
        for restart, initial in enumerate(initials, 1):
            begin = time.perf_counter()
            result = solve_spin_it_faster(request, initial_quadrics=(initial,))
            elapsed = time.perf_counter() - begin
            trace = result.trace
            attempts.append({
                "run": run,
                "restart": restart,
                "status": result.status,
                "newton_iterations": max((entry.iteration for entry in trace), default=0),
                "function_evaluations": max((entry.function_evaluations for entry in trace), default=0),
                "normalized_equation_residual": result.diagnostics["normalized_equation_residual"],
                "physical_relative_residual": result.diagnostics["physical_relative_residual"],
                "output_hash": result.diagnostics["output_hash"],
                "initial_quadric": initial,
                "elapsed_seconds": elapsed,
                "residual_seconds": result.residual_evaluation_seconds,
                "derivative_seconds": result.derivative_evaluation_seconds,
                "trace": [
                    {
                        "iteration": entry.iteration,
                        "function_evaluations": entry.function_evaluations,
                        "normalized_equation_residual": entry.normalized_equation_residual,
                        "physical_relative_residual": entry.physical_relative_residual,
                        "accepted_geodesic_step": entry.accepted_geodesic_step,
                    }
                    for entry in trace
                ],
            })
            if result.status == "CONVERGED":
                break
    summary = summarize_attempts(attempts, starts)
    payload = {
        "id": case["id"],
        "mode": case["mode"],
        "asset_class": "RECONSTRUCTED_ASSET",
        "attempts": attempts,
        "summary": summary,
    }
    payload["scientific_hash"] = scientific_hash(payload)
    return payload


def plan_payload(cases: list[dict[str, Any]], starts: int) -> dict[str, Any]:
    return {
        "schema": "spin-it-faster-table1-plan-v1",
        "schema_version": 1,
        "starts_per_case": starts,
        "cases": [
            {
                "id": case["id"],
                "mode": case["mode"],
                "initial_quadrics": generate_initial_quadrics(
                    case["mode"],
                    seed=2024 + 7919 * sum(map(ord, case["id"])),
                    count=starts,
                ),
            }
            for case in cases
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--starts", type=int, default=50)
    parser.add_argument("--max-restarts", type=int, default=8)
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument("--plan-only", action="store_true")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    cases = [case for case in manifest["cases"] if not args.case or case["id"] in args.case]
    if args.case and len(cases) != len(set(args.case)):
        raise SystemExit("one or more requested case IDs are absent")
    if args.starts <= 0 or args.max_restarts <= 0:
        raise SystemExit("starts and max-restarts must be positive")
    if args.plan_only:
        payload = plan_payload(cases, args.starts)
    else:
        rows = [run_case(case, args.starts, args.max_restarts) for case in cases]
        payload = {
            "schema": "spin-it-faster-table1-results-v1",
            "schema_version": 1,
            "paper_doi": "10.1145/3658194",
            "starts_per_case": args.starts,
            "maximum_restarts_per_run": args.max_restarts,
            "asset_class": "RECONSTRUCTED_ASSET",
            "author_exact": False,
            "physical_validation": False,
            "cases": rows,
        }
        payload["scientific_hash"] = scientific_hash(payload)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
