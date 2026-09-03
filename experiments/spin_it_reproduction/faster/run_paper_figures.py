#!/usr/bin/env python3
"""Generate auditable data for Spin-It Faster Figures 10--13 and Sec. 6.6."""

from __future__ import annotations

import argparse
from dataclasses import replace
import importlib.util
import json
import math
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parent
TABLE_RUNNER = ROOT / "run_paper_table1.py"
APPLICATION_GROUP = {
    "STATIC_BALANCE": "static_stability",
    "BUOYANT_STABILITY": "buoyant_stability",
    "SPINNING_TOP": "spinning_tops",
    "YO_YO": "yo_yos",
}


def _table_module():
    spec = importlib.util.spec_from_file_location("spin_it_faster_table", TABLE_RUNNER)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load Faster Table 1 runner")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def quadratic_convergence_pairs(table: dict[str, Any]) -> dict[str, list[tuple[float, float]]]:
    grouped: dict[str, list[tuple[float, float]]] = {
        name: [] for name in APPLICATION_GROUP.values()
    }
    for case in table.get("cases", []):
        group = APPLICATION_GROUP.get(case.get("mode"))
        if group is None:
            continue
        for attempt in case.get("attempts", []):
            if attempt.get("status") != "CONVERGED":
                continue
            residuals = [
                float(row["normalized_equation_residual"])
                for row in attempt.get("trace", [])
            ]
            if not all(math.isfinite(value) and value >= 0.0 for value in residuals):
                raise ValueError("Figure 10 residuals must be finite and nonnegative")
            for current, following in zip(residuals, residuals[1:]):
                if following > current * (1.0 + 1e-12):
                    raise ValueError("accepted Newton residuals must be monotone")
                if current > 0.0 and following > 0.0:
                    grouped[group].append((current, following))
    return grouped


def pareto_weights(count: int) -> list[float]:
    if count < 2:
        raise ValueError("Pareto sweep needs at least two weights")
    return [index / (count - 1) for index in range(count)]


def validate_pareto_front(rows: list[dict[str, Any]]) -> dict[str, bool]:
    values = [
        (float(row["first_objective"]), float(row["center_height"]))
        for row in rows
    ]
    finite = bool(values) and all(math.isfinite(a) and math.isfinite(b) for a, b in values)
    nondominated = finite and all(
        not any(
            j != i and other_a <= a and other_b <= b and
            (other_a < a or other_b < b)
            for j, (other_a, other_b) in enumerate(values)
        )
        for i, (a, b) in enumerate(values)
    )
    return {"finite": finite, "nondominated": nondominated}


def binary_search_critical_distance(
    is_feasible: Callable[[float], bool],
    *,
    lower: float,
    upper: float,
    tolerance: float = 1e-9,
) -> float:
    if not (math.isfinite(lower) and math.isfinite(upper) and lower < upper):
        raise ValueError("critical-distance bracket must be finite and ordered")
    if tolerance <= 0.0 or not math.isfinite(tolerance):
        raise ValueError("critical-distance tolerance must be finite and positive")
    if not is_feasible(lower) or is_feasible(upper):
        raise ValueError("critical-distance bracket must be feasible/infeasible")
    while upper - lower > tolerance:
        midpoint = 0.5 * (lower + upper)
        if is_feasible(midpoint):
            lower = midpoint
        else:
            upper = midpoint
    return 0.5 * (lower + upper)


def feasibility_distances(critical: float) -> list[float]:
    if not math.isfinite(critical):
        raise ValueError("critical distance must be finite")
    return [critical - 10.0 ** (-exponent) for exponent in range(1, 10)]


def degenerate_one_dimensional_examples() -> list[dict[str, Any]]:
    return [
        {
            "domain": [0.0, 1.0],
            "objective": "(v-1/2)^2",
            "target_volume": 0.5,
            "optimizer_cardinality": "UNCOUNTABLY_INFINITE",
            "gradient_at_optimum": 0.0,
            "p_zero": True,
        },
        {
            "domain": [0.0, 0.5],
            "objective": "(v-1/2)^2",
            "target_volume": 0.5,
            "optimizer_cardinality": "ONE",
            "unique_fill": 1.0,
            "gradient_at_optimum": 0.0,
            "p_zero": True,
        },
    ]


def moment_metrics(raw: list[float] | tuple[float, ...]) -> dict[str, float]:
    if len(raw) != 10 or raw[0] <= 0.0 or not all(math.isfinite(v) for v in raw):
        raise ValueError("raw moment vector must have ten finite values and positive mass")
    v, l1, l2, l3, q11, q22, q33, _q23, _q13, q12 = raw
    i11 = q22 + q33 - (l2 * l2 + l3 * l3) / v
    i22 = q11 + q33 - (l1 * l1 + l3 * l3) / v
    i33 = q11 + q22 - (l1 * l1 + l2 * l2) / v
    i12 = -q12 + l1 * l2 / v
    if abs(i33) <= 1e-15:
        raise ValueError("axial inertia vanishes")
    return {
        "volume_or_mass": v,
        "center_height": l3 / v,
        "f_yo": (i11 * i11 + i22 * i22 + 2.0 * i12 * i12) / (i33 * i33),
    }


def run_pareto_sweep(case: dict[str, Any], weights: list[float]) -> dict[str, Any]:
    from make_it_stand.spin import evaluate_spin_it_faster_quadric, solve_spin_it_faster

    table = _table_module()
    base = table.build_request(case)
    design = table.read_obj(ROOT / case["design_domain"]["path"])
    target_points = table.target_domain_points(case, design[0])
    candidates = table.generate_initial_quadrics(
        case["mode"],
        2024 + 7919 * sum(map(ord, case["id"])),
        8,
        target_points,
    )
    reference_weight = 1.0 if case["mode"] == "STATIC_BALANCE" else float(
        case.get("center_of_mass_weight", 0.0)
    )

    def request_for(weight: float):
        if case["mode"] == "STATIC_BALANCE":
            return replace(
                base,
                static_objective="PARETO_TRADEOFF",
                static_pareto_weight=weight,
            )
        if case["mode"] == "SPINNING_TOP":
            return replace(base, center_of_mass_weight=weight)
        raise ValueError("Pareto sweep supports static balance or spinning tops")

    reference = None
    for candidate in candidates:
        attempted = solve_spin_it_faster(
            request_for(reference_weight),
            initial_quadrics=(candidate,),
        )
        if attempted.status == "CONVERGED":
            reference = attempted
            break
    if reference is None:
        return {
            "case_id": case["id"],
            "mode": case["mode"],
            "warm_started": True,
            "reference_weight": reference_weight,
            "reference_status": attempted.status,
            "rows": [],
            "front_audit": {"finite": False, "nondominated": False},
        }

    reference_coefficients = list(reference.diagnostics["quadric_coefficients"])
    if case["mode"] == "STATIC_BALANCE":
        paths = [sorted(set(weights), reverse=True)]
    else:
        paths = [
            sorted({weight for weight in weights if weight <= reference_weight}, reverse=True),
            sorted({weight for weight in weights if weight > reference_weight}),
        ]
    rows: list[dict[str, Any]] = []
    for path in paths:
        initial = reference_coefficients
        for weight in path:
            request = request_for(weight)
            result = solve_spin_it_faster(request, initial_quadrics=(initial,))
            row: dict[str, Any] = {"weight": weight, "status": result.status}
            if result.status == "CONVERGED":
                raw = list(result.diagnostics["final_raw_moments"])
                metrics = moment_metrics(raw)
                coefficients = list(result.diagnostics["quadric_coefficients"])
                evaluated = evaluate_spin_it_faster_quadric(request, coefficients)
                row.update(metrics)
                row["first_objective"] = (
                    -metrics["volume_or_mass"]
                    if case["mode"] == "STATIC_BALANCE"
                    else metrics["f_yo"]
                )
                row["interface_subtype"] = evaluated.interface_subtype
                row["quadric_coefficients"] = coefficients
                initial = coefficients
            rows.append(row)
    rows.sort(key=lambda row: row["weight"])
    successful = [row for row in rows if row["status"] == "CONVERGED"]
    return {
        "case_id": case["id"],
        "mode": case["mode"],
        "warm_started": True,
        "reference_weight": reference_weight,
        "reference_status": reference.status,
        "rows": rows,
        "front_audit": validate_pareto_front(successful),
    }


_BOX_FACES = (
    (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
    (0, 1, 5), (0, 5, 4), (3, 7, 6), (3, 6, 2),
    (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5),
)


def _tilted_box(
    distance: float,
    half_extents: tuple[float, float, float],
) -> tuple[tuple[tuple[float, float, float], ...], tuple[tuple[int, int, int], ...]]:
    angle = math.radians(18.0)
    cosine, sine = math.cos(angle), math.sin(angle)
    hx, hy, hz = half_extents
    local = (
        (-hx, -hy, -hz), (hx, -hy, -hz), (hx, hy, -hz), (-hx, hy, -hz),
        (-hx, -hy, hz), (hx, -hy, hz), (hx, hy, hz), (-hx, hy, hz),
    )
    vertices = tuple(
        (
            distance * cosine + cosine * x + sine * z,
            y,
            -distance * sine - sine * x + cosine * z + 1.0,
        )
        for x, y, z in local
    )
    return vertices, _BOX_FACES


def tilted_cube_request(distance: float):
    from make_it_stand.spin import (
        SpinItFasterRequest,
        SpinSolverConfig,
        integrate_raw_mass_moments,
    )

    if not math.isfinite(distance) or distance < 0.0:
        raise ValueError("tilted-cube distance must be finite and nonnegative")
    outer = _tilted_box(distance, (0.9, 0.75, 1.0))
    design = _tilted_box(distance, (0.8, 0.65, 0.9))
    outer_moments = integrate_raw_mass_moments(*outer)
    design_moments = integrate_raw_mass_moments(*design)
    fixed = tuple(a - b for a, b in zip(outer_moments, design_moments))
    angle = math.radians(18.0)
    cosine, sine = math.cos(angle), math.sin(angle)
    return SpinItFasterRequest(
        outer_vertices=outer[0],
        outer_faces=outer[1],
        design_vertices=design[0],
        design_faces=design[1],
        fixed_shell_moments=fixed,
        target_origin=(-sine, 0.0, 1.0 - cosine),
        target_axes=((cosine, 0.0, sine), (0.0, 1.0, 0.0), (-sine, 0.0, cosine)),
        mode="SPINNING_TOP",
        center_of_mass_weight=0.25,
        evidence="RECONSTRUCTED_ASSET",
        minimum_wall_thickness=0.09,
        settings=SpinSolverConfig(maximum_restarts=1),
    )


def run_tilted_cube_feasibility(
    *,
    tolerance: float = 1e-9,
    upper: float = 0.8,
) -> dict[str, Any]:
    from make_it_stand.spin import solve_spin_it_faster

    table = _table_module()
    base = tilted_cube_request(0.0)
    case = {
        "mode": "SPINNING_TOP",
        "target_origin": list(base.target_origin),
        "target_axes": [list(row) for row in base.target_axes],
    }
    target_points = table.target_domain_points(case, base.design_vertices)
    candidates = table.generate_initial_quadrics(
        "SPINNING_TOP", 999, 8, target_points
    )
    base_result = None
    for candidate in candidates:
        attempted = solve_spin_it_faster(base, initial_quadrics=(candidate,))
        if attempted.status == "CONVERGED":
            base_result = attempted
            break
    if base_result is None:
        raise RuntimeError("unperturbed tilted-cube reconstruction did not converge")
    initial = list(base_result.diagnostics["quadric_coefficients"])
    cache: dict[float, Any] = {0.0: base_result}

    def solve(distance: float):
        if distance not in cache:
            cache[distance] = solve_spin_it_faster(
                tilted_cube_request(distance),
                initial_quadrics=(initial,),
            )
        return cache[distance]

    critical = binary_search_critical_distance(
        lambda distance: solve(distance).status == "CONVERGED",
        lower=0.0,
        upper=upper,
        tolerance=tolerance,
    )
    rows = []
    for exponent, distance in enumerate(feasibility_distances(critical), 1):
        if distance < 0.0:
            rows.append({"m": exponent, "distance": distance, "status": "OUTSIDE_BRACKET"})
            continue
        result = solve(distance)
        rows.append({
            "m": exponent,
            "distance": distance,
            "status": result.status,
            "trace": [
                {
                    "iteration": entry.iteration,
                    "normalized_equation_residual": entry.normalized_equation_residual,
                }
                for entry in result.trace
            ],
        })
    return {
        "asset_class": "RECONSTRUCTED_ASSET",
        "critical_distance": critical,
        "binary_search_tolerance": tolerance,
        "unperturbed_status": base_result.status,
        "rows": rows,
    }


def build_payload(table: dict[str, Any] | None = None) -> dict[str, Any]:
    return {
        "schema": "spin-it-faster-paper-figures-v1",
        "schema_version": 1,
        "paper_doi": "10.1145/3658194",
        "author_exact": False,
        "physical_validation": False,
        "figure10": {
            "definition": "Eq.14 l2 error e_n versus e_(n+1), successful Table 1 runs",
            "groups": quadratic_convergence_pairs(table or {"cases": []}),
            "threshold": 1e-12,
        },
        "figure11": {
            "definition": "min -kappa*v+(1-kappa)*c3; l1=l2=0",
            "weights": pareto_weights(41),
            "status": "READY_RECONSTRUCTED_ASSET",
        },
        "figure12": {
            "definition": "min f_yo+w_c*c3; l1=l2=I13=I23=0",
            "weights": pareto_weights(41),
            "status": "READY_RECONSTRUCTED_ASSET",
        },
        "figure13": {
            "definition": "d_m=d_crit-10^-m, m=1..9",
            "binary_search_tolerance": 1e-9,
            "status": "REQUIRES_RECONSTRUCTED_TILTED_CUBE_FAMILY",
        },
        "section6_6": degenerate_one_dimensional_examples(),
    }


def finalize_payload(payload: dict[str, Any]) -> dict[str, Any]:
    result = dict(payload)
    result.pop("scientific_hash", None)
    result["scientific_hash"] = _table_module().scientific_hash(result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table1", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--run-pareto", action="store_true")
    parser.add_argument("--pareto-count", type=int, default=41)
    parser.add_argument("--run-figure13", action="store_true")
    parser.add_argument("--figure13-tolerance", type=float, default=1e-9)
    args = parser.parse_args()
    table = json.loads(args.table1.read_text(encoding="utf-8")) if args.table1 else None
    payload = build_payload(table)
    if args.run_pareto:
        manifest = json.loads((ROOT / "paper_cases.json").read_text(encoding="utf-8"))
        cases = {case["id"]: case for case in manifest["cases"]}
        weights = pareto_weights(args.pareto_count)
        payload["figure11"]["result"] = run_pareto_sweep(
            cases["static-tilted-vase"], weights
        )
        payload["figure12"]["result"] = run_pareto_sweep(
            cases["heart-spinning-top"], weights
        )
    if args.run_figure13:
        payload["figure13"]["result"] = run_tilted_cube_feasibility(
            tolerance=args.figure13_tolerance
        )
    payload = finalize_payload(payload)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
