#!/usr/bin/env python3
"""Run the public Spin-It Faster Table 1 protocol on reconstructed assets."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import random
import statistics
import sys
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
            if key not in TIMING_KEYS and not key.endswith("_seconds")
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
            # The material evaluator retains q>=0.  The paper's fidget void
            # is therefore initialized as the interior of a positive-curvature
            # circular cylinder, not its complement.
            values[0] = cx * cx + cy * cy - radius * radius
            values[1] = -2.0 * cx
            values[2] = -2.0 * cy
            values[4] = 1.0
            values[5] = 1.0
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
    if type(expected_runs) is not int or expected_runs <= 0:
        raise ValueError("statistics require a positive registered top-level population")
    successful: list[dict[str, Any]] = []
    restarts: list[int] = []
    attempt_counts: list[int] = []
    evaluation_counts: list[int] = []
    resource_limited_runs = 0
    resource_limited_attempts = 0
    runs_with_resource_limited_attempts = 0

    def budget_exhausted(attempt: dict[str, Any]) -> bool:
        # Dogleg uses its own terminal enum: do not require the Eq.14 wrapper's
        # RESOURCE_LIMIT spelling or misclassify geometry/line-search failure.
        if attempt.get("solver_family") == "PLANAR_DOGLEG_V1":
            return attempt.get("status") in {
                "MAXIMUM_FUNCTION_EVALUATIONS", "MAXIMUM_ITERATIONS"}
        if attempt.get("solver_family") == "STATIC_PLANE_FAMILY_V1":
            return attempt.get("status") == "ITERATION_LIMIT"
        return attempt.get("status") == "RESOURCE_LIMIT"

    for run in range(expected_runs):
        run_attempts = [item for item in attempts if item["run"] == run]
        attempt_counts.append(len(run_attempts))
        evaluation_counts.append(sum(item["function_evaluations"] for item in run_attempts))
        limited = sum(budget_exhausted(item) for item in run_attempts)
        resource_limited_attempts += limited
        runs_with_resource_limited_attempts += int(limited > 0)
        converged = next((item for item in run_attempts if item["status"] == "CONVERGED"), None)
        if converged is not None:
            successful.append(converged)
            restarts.append(int(converged["restart"]))
        elif run_attempts and budget_exhausted(run_attempts[-1]):
            resource_limited_runs += 1
    def median(key: str) -> float | None:
        values = [float(item[key]) for item in successful if item.get(key) is not None]
        return statistics.median(values) if values else None
    # Preserve historical native-status statistics, but never let those labels
    # serve as raw Eq.14 approval. Recompute from actual terminal values; a
    # stored numerical_acceptance/passed flag is not an input to this verdict.
    qualified: dict[int, dict[str, Any]] = {}
    rejected = []
    native_converged_but_unqualified = 0
    for index, attempt in enumerate(attempts):
        run = attempt.get('run')
        if type(run) is not int or not 0 <= run < expected_runs:
            raise ValueError("attempt does not belong to the registered top-level population")
        verdict = numerical_attempt_acceptance(attempt)
        if verdict['passed']:
            qualified.setdefault(run, attempt)
        else:
            rejected.append({'attempt_index':index,'run':run,
                             'restart':attempt.get('restart'),'reasons':verdict['reasons']})
            if attempt.get('status') == 'CONVERGED':
                native_converged_but_unqualified += 1
    strict_statistics = {
        'scope':'TERMINAL_NUMERICS_STATISTICS_NOT_MATERIAL_OR_PAPER_ACCEPTANCE',
        'expected_runs':expected_runs,
        'numerically_qualified_runs':len(qualified),
        'unqualified_runs':expected_runs-len(qualified),
        'qualified_run_ids':sorted(qualified),
        'all_runs_numerically_qualified':len(qualified)==expected_runs,
        'native_converged_but_unqualified_attempts':native_converged_but_unqualified,
        'rejected_attempts':rejected,
        'total_function_evaluations_all_attempts':sum(evaluation_counts),
        'average_function_evaluations_all_runs':sum(evaluation_counts)/expected_runs,
    }
    return {
        "expected_runs": expected_runs,
        "successful_runs": len(successful),
        "failed_runs": expected_runs - len(successful),
        "resource_limited_runs": resource_limited_runs,
        "resource_limited_runs_scope": "UNSUCCESSFUL_RUNS_WITH_BUDGET_EXHAUSTED_FINAL_ATTEMPT",
        "resource_limited_attempts": resource_limited_attempts,
        "runs_with_resource_limited_attempts": runs_with_resource_limited_attempts,
        "resource_limit_statistics_schema": "solver-family-budget-accounting-v1",
        "all_runs_converged": len(successful) == expected_runs,
        "average_attempts_all_runs": sum(attempt_counts) / expected_runs if expected_runs else None,
        "total_function_evaluations_all_attempts": sum(evaluation_counts),
        "function_evaluations_by_run": evaluation_counts,
        "average_function_evaluations_all_runs": sum(evaluation_counts) / expected_runs if expected_runs else None,
        "conditional_statistics_population": "SUCCESSFUL_TERMINAL_ATTEMPTS_ONLY",
        "average_restarts": sum(restarts) / len(restarts) if restarts else None,
        "median_newton_iterations": median("newton_iterations"),
        "median_function_evaluations": median("function_evaluations"),
        "median_residual_seconds": median("residual_seconds"),
        "median_derivative_seconds": median("derivative_seconds"),
        "native_statistics_scope": "NATIVE_STATUS_STATISTICS_NOT_STRICT_ACCEPTANCE",
        "strict_terminal_numerics": strict_statistics,
    }


def initial_run_plan(case: dict[str, Any], run: int, maximum_restarts: int,
                     domain_points, *, seed: int | None = None) -> dict[str, Any]:
    """One input-only initialization protocol shared by planning and execution."""
    if type(run) is not int or run < 0 or type(maximum_restarts) is not int or maximum_restarts <= 0:
        raise ValueError("run and restart budget must be nonnegative/positive integers")
    seed = run if seed is None else seed
    if type(seed) is not int or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    return {
        "run": run,
        "seed": seed,
        "initial_quadrics": generate_initial_quadrics(
            case["mode"], seed=seed, count=maximum_restarts,
            domain_points=domain_points),
    }


def registered_seeds(starts: int, initial_seeds=None) -> tuple[int, ...]:
    if type(starts) is not int or starts <= 0:
        raise ValueError("starts must be a positive integer")
    seeds = tuple(range(starts)) if initial_seeds is None else tuple(initial_seeds)
    if len(seeds) != starts or any(type(seed) is not int or seed < 0 for seed in seeds) or len(set(seeds)) != starts:
        raise ValueError("initial seed sequence must contain one distinct nonnegative integer per run")
    return seeds


def solver_event_payload(event: Any) -> dict[str, Any]:
    """Preserve every deterministic native solver event for one restart."""
    fields = (
        "type", "restart", "iteration", "evaluation_index",
        "restart_evaluation_index", "step", "outcome", "reason", "detail",
        "quadric_coefficients", "quadric_hash", "normalized_equation_residual",
        "physical_relative_residual", "objective", "multiplier_mu",
        "material_geometry_hash", "material_topology_hash", "equation_profile",
        "solver_coordinate_identity", "raw_equation_residual",
    )
    return {
        field: list(value) if field == "quadric_coefficients" and value is not None
        else value
        for field in fields
        if (value := getattr(event, field, None)) is not None
    }


def validate_attempt_ledger(
    attempts: list[dict[str, Any]], *, expected_runs: int, maximum_restarts: int
) -> None:
    """Fail closed unless each top-level run has a complete, budgeted restart ledger."""
    if type(expected_runs) is not int or type(maximum_restarts) is not int:
        raise ValueError("run/restart budgets must be integer counters")
    if expected_runs <= 0 or maximum_restarts <= 0:
        raise ValueError("expected runs and maximum restarts must be positive")
    grouped: dict[int, list[dict[str, Any]]] = {run: [] for run in range(expected_runs)}
    seen: set[tuple[int, int]] = set()
    for attempt in attempts:
        if any(type(attempt.get(key)) is not int
               for key in ("run", "restart", "function_evaluations")):
            raise ValueError("attempt identity and budget must be integer counters")
        run, restart = attempt["run"], attempt["restart"]
        key = (run, restart)
        if run not in grouped or restart < 1 or restart > maximum_restarts or key in seen:
            raise ValueError("attempt ledger has an out-of-range or duplicate restart")
        seen.add(key)
        grouped[run].append(attempt)
        evaluations = attempt["function_evaluations"]
        if attempt.get("solver_family") == "STATIC_PLANE_FAMILY_V1":
            validate_static_plane_ledger(attempt)
            continue
        if evaluations < 0 or evaluations > 60:
            raise ValueError("per-attempt function budget exceeded")
        if attempt.get("solver_family") == "PLANAR_DOGLEG_V1":
            validate_planar_evaluation_ledger(attempt)
            continue
        events = attempt.get("solver_events")
        if not isinstance(events, list) or not events:
            raise ValueError("solver event ledger is missing")
        starts = [event for event in events if event.get("type") == "attempt_start"]
        terminals = [event for event in events if event.get("type") == "attempt_terminal"]
        if len(starts) != 1 or len(terminals) != 1:
            raise ValueError("each solver attempt requires one start and one terminal event")
        if events[0] != starts[0] or events[-1] != terminals[0]:
            raise ValueError("solver attempt event order is inconsistent")
        if attempt.get("status") == "CONVERGED" and (
                terminals[0].get("outcome") != "terminated" or
                terminals[0].get("reason") != "converged"):
            raise ValueError("claimed convergence differs from native terminal event")
        native_restart = restart - 1 if attempt.get("solver_family") == "CHERCHI_RAW_EQ14_V1" else 0
        if attempt.get("solver_family") == "CHERCHI_RAW_EQ14_V1" and attempt.get("native_restart") != native_restart:
            raise ValueError("Cherchi native restart identity differs from table attempt")
        if any(type(event.get("restart")) is not int or event["restart"] != native_restart
               for event in events):
            raise ValueError("native solve has inconsistent restart identity")
        evaluation_events = [event for event in events if event.get("type") == "evaluation"]
        if len(evaluation_events) != evaluations:
            raise ValueError("solver event evaluation count differs from attempt count")
        # Each table attempt calls a fresh single-start native solve, so both
        # counters begin at one. Counting rows alone would permit a duplicated
        # accepted evaluation to replace a missing rejected trial.
        offset = attempt.get("native_evaluation_offset", 0) if native_restart else 0
        if type(offset) is not int or offset < 0:
            raise ValueError("invalid native evaluation offset")
        if attempt.get("solver_family") == "CHERCHI_RAW_EQ14_V1" and (
                starts[0].get("evaluation_index") != offset or
                starts[0].get("restart_evaluation_index") != 0 or
                terminals[0].get("evaluation_index") != offset + evaluations or
                terminals[0].get("restart_evaluation_index") != evaluations):
            raise ValueError("Cherchi terminal evaluation count differs from complete ledger")
        for ordinal, event in enumerate(evaluation_events, start=1):
            for field, expected in (("evaluation_index", offset + ordinal), ("restart_evaluation_index", ordinal)):
                if type(event.get(field)) is not int or event[field] != expected:
                    raise ValueError("solver evaluation ordinal sequence is incomplete or duplicated")
    for run, rows in grouped.items():
        rows.sort(key=lambda item: int(item["restart"]))
        if not rows:
            raise ValueError(f"run {run} has no recorded restart attempts")
        ordinals = [int(row["restart"]) for row in rows]
        if ordinals != list(range(1, len(rows) + 1)):
            raise ValueError(f"run {run} restart sequence is incomplete")
        successes = [i for i, row in enumerate(rows) if row.get("status") == "CONVERGED"]
        if len(successes) > 1 or (successes and successes[0] != len(rows) - 1):
            raise ValueError(f"run {run} contains attempts after its first convergence")
        expected_count = ordinals[-1] if successes else maximum_restarts
        if len(rows) != expected_count:
            raise ValueError(f"run {run} is missing failed or unattempted restart records")
        offset = 0
        for row in rows:
            if row.get("solver_family") == "CHERCHI_RAW_EQ14_V1":
                if row.get("native_evaluation_offset") != offset:
                    raise ValueError("Cherchi cumulative evaluation offsets are inconsistent")
                offset += row["function_evaluations"]


def read_obj(
    path: Path, data: bytes | None = None
) -> tuple[tuple[tuple[float, float, float], ...], tuple[tuple[int, int, int], ...]]:
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []
    for line in (path.read_bytes() if data is None else data).decode("utf-8").splitlines():
        fields = line.split()
        if fields and fields[0] == "v":
            vertices.append(tuple(float(value) for value in fields[1:4]))
        elif fields and fields[0] == "f":
            faces.append(tuple(int(value) - 1 for value in fields[1:4]))
    return tuple(vertices), tuple(faces)


def read_case_mesh(
    case: dict[str, Any], field: str, *, asset_root: Path = ROOT
) -> tuple[tuple[tuple[float, float, float], ...], tuple[tuple[int, int, int], ...]]:
    """Read exactly the bytes identified by the case manifest, inside its root."""
    entry = case[field]
    relative = Path(entry["path"])
    root = asset_root.resolve()
    path = (root / relative).resolve()
    if relative.is_absolute() or not path.is_relative_to(root):
        raise ValueError(f"{field} path must remain inside the manifest directory")
    data = path.read_bytes()
    expected = entry.get("sha256")
    if not isinstance(expected, str) or len(expected) != 64 or hashlib.sha256(data).hexdigest() != expected:
        raise ValueError(f"{field} SHA-256 does not match the manifest")
    return read_obj(path, data)


def target_domain_points(
    case: dict[str, Any],
    design_points: tuple[tuple[float, float, float], ...],
) -> tuple[tuple[float, float, float], ...]:
    origin = tuple(float(value) for value in case["target_origin"])
    axes = tuple(tuple(float(value) for value in row) for row in case["target_axes"])
    return tuple(
        tuple(
            sum(axes[row][column] * (point[row] - origin[row]) for row in range(3))
            for column in range(3)
        )
        for point in design_points
    )


def build_request(case: dict[str, Any], *, asset_root: Path = ROOT):
    from make_it_stand.spin import SpinItFasterRequest, SpinSolverConfig, integrate_raw_mass_moments

    outer = read_case_mesh(case, "outer", asset_root=asset_root)
    design = read_case_mesh(case, "design_domain", asset_root=asset_root)
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
        static_pareto_weight=float(case.get("static_pareto_weight", 0.0)),
        static_objective_measure=case.get("static_objective_measure","DENSITY_WEIGHTED_MASS_V1"),
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


def numerical_attempt_acceptance(attempt: dict[str, Any]) -> dict[str, Any]:
    if attempt.get("solver_family") == "STATIC_PLANE_FAMILY_V1":
        return static_plane_attempt_acceptance(attempt)
    """Check raw terminal numerics, never asset/material or paper completion.

    The caller must also validate the native event ledger. Keeping a legacy
    normalized solve available does not make its normalized stopping rule an
    Eq.14 acceptance certificate, even when a summary labels it CONVERGED.
    """
    if attempt.get("solver_family") == "PLANAR_DOGLEG_V1":
        return planar_attempt_acceptance(attempt)
    reasons = []
    if attempt.get('status') != 'CONVERGED':
        reasons.append('NATIVE_CONVERGENCE_REQUIRED')
    if attempt.get('equation_profile') != 'DUAL_MONOMIAL_RAW_V1':
        reasons.append('RAW_EQ14_PROFILE_REQUIRED')

    def finite(value):
        return type(value) in (int, float) and math.isfinite(value)

    raw = attempt.get('raw_equation_residual')
    physical = attempt.get('physical_relative_residual')
    multiplier = attempt.get('multiplier_mu')
    if not finite(raw) or not 0. <= raw <= 1e-12:
        reasons.append('RAW_EQ14_TOLERANCE')
    if not finite(physical) or not 0. <= physical <= 1e-6:
        reasons.append('PHYSICAL_CONSTRAINT_TOLERANCE')
    if not finite(multiplier) or multiplier <= 0.:
        reasons.append('POSITIVE_MULTIPLIER_REQUIRED')
    trace = attempt.get('trace')
    if not isinstance(trace, list) or not trace or not isinstance(trace[-1], dict):
        reasons.append('TERMINAL_TRACE_REQUIRED')
    else:
        terminal = trace[-1]
        if (terminal.get('equation_profile') != attempt.get('equation_profile') or
                not finite(terminal.get('raw_equation_residual')) or
                not finite(terminal.get('multiplier_mu')) or
                terminal.get('raw_equation_residual') != raw or
                terminal.get('multiplier_mu') != multiplier):
            reasons.append('TERMINAL_NUMERICS_MISMATCH')
    return {'passed':not reasons, 'reasons':reasons,
            'scope':'TERMINAL_NUMERICS_ONLY_NOT_MATERIAL_OR_PAPER_ACCEPTANCE'}


def planar_case_request(case: dict[str, Any]) -> dict[str, Any]:
    """The 2D method needs actual layer geometry, not a reconstructed 3D mesh."""
    problem = case.get("planar_problem")
    if not isinstance(problem, dict):
        raise ValueError("FIDGET_SPINNER requires an explicit planar_problem")
    if "initial_center" in problem:
        raise ValueError("initial_center must come from the registered seed protocol")
    allowed = {"design_region_loops", "fixed_material_loops", "fixed_areal_density",
               "design_areal_density", "maximum_iterations", "maximum_function_evaluations",
               "residual_tolerance", "initial_trust_radius", "maximum_trust_radius",
               "extrusion_interval", "extrusion_volume_density", "fixed_extrusion_layers"}
    if set(problem)-allowed:
        raise ValueError("unsupported planar_problem fields")
    if not problem.get("fixed_extrusion_layers") or "extrusion_interval" not in problem \
            or "extrusion_volume_density" not in problem:
        raise ValueError("planar table dispatch requires actual fixed layers and design extrusion")
    if case.get("target_origin", [0., 0., 0.]) != [0., 0., 0.] or \
            case.get("target_axes", [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]) != \
            [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]:
        raise ValueError("planar_problem must already use the target-axis frame")
    request = {"maximum_iterations": 100, "maximum_function_evaluations": 60,
               "residual_tolerance": 1e-10, "initial_trust_radius": .25,
               "maximum_trust_radius": 1., **problem}
    budget = request["maximum_function_evaluations"]
    tolerance = request["residual_tolerance"]
    if type(budget) is not int or not 1 <= budget <= 60:
        raise ValueError("planar function evaluation budget must be an integer in [1,60]")
    if type(tolerance) not in (int, float) or not math.isfinite(tolerance) or \
            not 0 < tolerance <= 1e-10:
        raise ValueError("planar residual tolerance must be positive and at most 1e-10")
    # Serialization freezes the caller's mutable geometry before solving.
    return json.loads(json.dumps(request, allow_nan=False))


def planar_initial_run_plan(request, run, seed, maximum_restarts):
    if type(maximum_restarts) is not int or maximum_restarts <= 0:
        raise ValueError("maximum_restarts must be a positive integer")
    points = [point for loop in request["design_region_loops"] for point in loop]
    if not points or any(len(point) != 2 or any(type(x) not in (int, float) or
            not math.isfinite(x) for x in point) for point in points):
        raise ValueError("planar initial bounds require finite two-dimensional geometry")
    lo = [min(p[i] for p in points) for i in (0, 1)]
    hi = [max(p[i] for p in points) for i in (0, 1)]
    if any(hi[i] <= lo[i] for i in (0, 1)):
        raise ValueError("planar initial bounds must have positive area")
    generator = random.Random(seed)
    centers = [[.5*(lo[i]+hi[i]) + generator.uniform(-.2, .2)*.5*(hi[i]-lo[i])
                for i in (0, 1)] for _ in range(maximum_restarts)]
    return {"run": run, "seed": seed, "initial_centers": centers,
            "protocol": "INPUT_BOUNDING_BOX_PLANAR_CENTERS_V1"}


def validate_planar_evaluation_ledger(attempt):
    ledger = attempt.get("evaluation_ledger")
    if not isinstance(ledger, dict) or ledger.get("schema") != "planar-dogleg-evaluation-ledger-v1":
        raise ValueError("actual planar evaluation ledger is missing")
    solver, auxiliary = ledger.get("solver_evaluations"), ledger.get("auxiliary_evaluations")
    if not isinstance(solver, list) or not isinstance(auxiliary, list):
        raise ValueError("planar call records are missing")
    expected = {"solver_function_evaluations": len(solver),
                "auxiliary_geometry_evaluations": len(auxiliary),
                "total_residual_evaluations": len(solver)+len(auxiliary)}
    if any(type(ledger.get(k)) is not int or ledger[k] != v for k, v in expected.items()) or \
            type(attempt.get("function_evaluations")) is not int or \
            attempt["function_evaluations"] != len(solver):
        raise ValueError("planar evaluation counts disagree with actual ledger")
    budget = attempt.get("function_evaluation_budget")
    if type(budget) is not int or not 1 <= budget <= 60 or len(solver) > budget:
        raise ValueError("planar evaluation budget is invalid")
    for records in (solver, auxiliary):
        for ordinal, event in enumerate(records, 1):
            if type(event.get("evaluation_index")) is not int or event["evaluation_index"] != ordinal:
                raise ValueError("planar evaluation ordinal sequence is incomplete or duplicated")
            if event.get("outcome") not in {"evaluated", "failed"} or \
                    type(event.get("accepted_trial")) is not bool:
                raise ValueError("planar call outcome is incomplete")
            if event["accepted_trial"] and event.get("purpose") != "trial":
                raise ValueError("non-trial planar call was marked accepted")
            if event["outcome"] == "failed" and (event.get("has_residual") is not False or
                    event.get("normalized_residual") is not None):
                raise ValueError("failed planar evaluation fabricated a residual")
    result = attempt.get("solver_result")
    if result is not None:
        if not solver or solver[0].get("purpose") != "initial":
            raise ValueError("returned planar solve has no initial call")
        for event in solver + auxiliary:
            if event.get("outcome") != "evaluated" or event.get("has_residual") is not True:
                raise ValueError("returned planar solve contains an exceptional or absent evaluation")
            if result.get("status") == "CONVERGED":
                vector = event.get("normalized_residual_vector")
                center = event.get("center")
                value = event.get("normalized_residual")
                if any(not isinstance(pair, list) or len(pair) != 2 or
                       any(type(x) not in (int, float) or not math.isfinite(x) for x in pair)
                       for pair in (vector, center)) or \
                        type(value) not in (int, float) or not math.isfinite(value) or \
                        not math.isclose(value, math.hypot(*vector), rel_tol=1e-14, abs_tol=0.):
                    raise ValueError("converged planar evaluation has invalid vector or norm")
        if result.get("status") != attempt.get("status") or \
                result.get("function_evaluations") != len(solver) or \
                result.get("accepted_steps") != sum(row["accepted_trial"] for row in solver):
            raise ValueError("planar result disagrees with actual calls or terminal status")
        if ledger.get("terminal_phase") != "complete" or len(auxiliary) != 2 or \
                [row.get("purpose") for row in auxiliary] != \
                ["binding_initial_assessment", "binding_final_assessment"]:
            raise ValueError("completed planar binding lost auxiliary assessments")
        terminal = auxiliary[-1]
        accepted = [event for event in solver if event["accepted_trial"]]
        last_state = accepted[-1] if accepted else solver[0]
        if result.get("center") != terminal.get("center") or \
                terminal.get("center") != last_state.get("center") or \
                terminal.get("normalized_residual_vector") != last_state.get("normalized_residual_vector"):
            raise ValueError("planar terminal assessment differs from the actual accepted state")
    elif attempt.get("status") != "EVALUATION_FAILED":
        raise ValueError("planar attempt lost its result or exception")


def planar_attempt_acceptance(attempt):
    reasons = []
    try:
        validate_planar_evaluation_ledger(attempt)
    except (KeyError, TypeError, ValueError) as error:
        return {"passed": False, "reasons": ["PLANAR_LEDGER_INVALID: " + str(error)],
                "scope": "PLANAR_TERMINAL_NUMERICS_NOT_MATERIAL_OR_PAPER_ACCEPTANCE"}
    result = attempt.get("solver_result")
    if attempt.get("status") != "CONVERGED" or result is None:
        reasons.append("PLANAR_DOGLEG_CONVERGENCE_REQUIRED")
    else:
        tolerance = attempt.get("residual_tolerance")
        residual = result.get("normalized_residual")
        valid_tolerance = type(tolerance) in (int, float) and \
            math.isfinite(tolerance) and 0 < tolerance <= 1e-10
        if not valid_tolerance or \
                type(residual) not in (int, float) or not math.isfinite(residual) or \
                not 0 <= residual <= tolerance:
            reasons.append("PLANAR_DOGLEG_RESIDUAL_TOLERANCE")
        terminal = attempt["evaluation_ledger"]["auxiliary_evaluations"][-1]
        value = terminal.get("normalized_residual")
        if not valid_tolerance or terminal.get("outcome") != "evaluated" or \
                type(value) not in (int, float) or \
                not math.isfinite(value) or not 0 <= value <= tolerance:
            reasons.append("PLANAR_FINAL_ASSESSMENT_REQUIRED")
        if residual != attempt.get("normalized_equation_residual") or residual != value:
            reasons.append("PLANAR_TERMINAL_NUMERICS_MISMATCH")
    return {"passed": not reasons, "reasons": reasons,
            "scope": "PLANAR_TERMINAL_NUMERICS_NOT_MATERIAL_OR_PAPER_ACCEPTANCE"}


def _planar_json_value(value):
    if hasattr(value, "tolist"):
        return _planar_json_value(value.tolist())
    if isinstance(value, dict):
        return {key: _planar_json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_planar_json_value(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return "nan" if math.isnan(value) else "+inf" if value > 0 else "-inf"
    return value


def run_planar_case(case, starts, max_restarts, *, initial_seeds=None):
    from make_it_stand.spin import solve_spin_it_faster_planar_fidget
    request = planar_case_request(case)
    seeds = registered_seeds(starts, initial_seeds)
    plans = [planar_initial_run_plan(request, run, seed, max_restarts)
             for run, seed in enumerate(seeds)]
    attempts = []
    for plan in plans:
        for restart, center in enumerate(plan["initial_centers"], 1):
            began = time.perf_counter()
            result, error_record = None, None
            try:
                result = _planar_json_value(solve_spin_it_faster_planar_fidget(
                    **request, initial_center=center))
                ledger = result["evaluation_ledger"]
            except Exception as error:
                ledger = getattr(error, "evaluation_ledger", None)
                if ledger is None:
                    # Invalid request/infrastructure is not a zero-cost solver
                    # attempt whose unknown evaluations can be invented.
                    raise
                error_record = {"type": type(error).__name__, "message": str(error)}
            row = {"run": plan["run"], "seed": plan["seed"], "restart": restart,
                   "solver_family": "PLANAR_DOGLEG_V1", "initial_center": center,
                   "status": result["status"] if result else "EVALUATION_FAILED",
                   "function_evaluations": ledger["solver_function_evaluations"],
                   "function_evaluation_budget": request["maximum_function_evaluations"],
                   "residual_tolerance": request["residual_tolerance"],
                   "normalized_equation_residual": result["normalized_residual"] if result else None,
                   "dogleg_iterations": result["iterations"] if result else None,
                   "evaluation_ledger": ledger, "solver_result": result,
                   "exception": error_record, "elapsed_seconds": time.perf_counter()-began}
            row["numerical_acceptance"] = planar_attempt_acceptance(row)
            attempts.append(row)
            if row["status"] == "CONVERGED":
                break
    validate_attempt_ledger(attempts, expected_runs=starts, maximum_restarts=max_restarts)
    summary = summarize_attempts(attempts, starts)
    payload = {"id": case["id"], "mode": case["mode"],
               "solver_family": "PLANAR_DOGLEG_V1", "asset_class": case.get("asset_class"),
               "initialization_protocol": "INPUT_BOUNDING_BOX_PLANAR_CENTERS_V1",
               "seed_sequence": list(seeds), "initialization_plan": plans,
               "planar_problem": request, "planar_problem_sha256": scientific_hash(request),
               "attempts": attempts, "summary": summary,
               "total_auxiliary_residual_evaluations": sum(
                   row["evaluation_ledger"]["auxiliary_geometry_evaluations"] for row in attempts),
               "total_residual_evaluations": sum(
                   row["evaluation_ledger"]["total_residual_evaluations"] for row in attempts),
               "passed": summary["strict_terminal_numerics"]["all_runs_numerically_qualified"],
               "acceptance_scope": "PLANAR_TERMINAL_NUMERICS_NOT_MATERIAL_OR_PAPER_ACCEPTANCE",
               "material_validation": "PENDING_INDEPENDENT_SAVED_GEOMETRY_CHECK",
               "formal_paper_admission": "NOT_GRANTED_BY_METHOD_RUNNER"}
    payload["scientific_hash"] = scientific_hash(payload)
    return payload


def static_plane_controls(case):
    controls = {"maximum_iterations": 60, "maximum_line_search_steps": 24,
                "constraint_tolerance": 1e-10}
    supplied = case.get("static_plane_settings", {})
    if not isinstance(supplied, dict) or set(supplied)-set(controls):
        raise ValueError("static_plane_settings supports only registered budgets and tolerance; no initial_coefficients")
    controls.update(supplied)
    for key, low, high in (("maximum_iterations",0,60), ("maximum_line_search_steps",1,24)):
        if type(controls[key]) is not int or not low <= controls[key] <= high:
            raise ValueError(f"invalid static plane {key}")
    tolerance = controls["constraint_tolerance"]
    if type(tolerance) not in (int,float) or not math.isfinite(tolerance) or not 0 < tolerance <= 1e-10:
        raise ValueError("static plane tolerance must be finite, positive and at most 1e-10")
    return controls


def static_plane_problem(case, *, asset_root=ROOT):
    outer = read_case_mesh(case, "outer", asset_root=asset_root)
    design = read_case_mesh(case, "design_domain", asset_root=asset_root)
    controls = static_plane_controls(case)
    objective = case.get("static_objective", "MAXIMUM_MASS")
    weight = case.get("static_pareto_weight", 0.)
    density = case.get("material_density", 1.)
    measure = case.get("static_objective_measure","DENSITY_WEIGHTED_MASS_V1")
    if measure not in ("DENSITY_WEIGHTED_MASS_V1","UNIFORM_DENSITY_VOLUME_V1"):
        raise ValueError("invalid static objective measure")
    if objective not in ("MAXIMUM_MASS", "LOWER_CENTER_OF_MASS", "PARETO_TRADEOFF"):
        raise ValueError("invalid static objective")
    if type(weight) not in (int,float) or not math.isfinite(weight) or not 0 <= weight <= 1 or (
            objective != "PARETO_TRADEOFF" and weight != 0):
        raise ValueError("invalid static Pareto weight")
    if type(density) not in (int,float) or not math.isfinite(density) or density <= 0:
        raise ValueError("invalid material density")
    points = target_domain_points(case, outer[0])
    length = max(max(p[i] for p in points)-min(p[i] for p in points) for i in range(3))
    if not math.isfinite(length) or length <= 0:
        raise ValueError("invalid static input normalization")
    return {"schema": "static-plane-input-v1", "outer_vertices": outer[0], "outer_faces": outer[1],
            "design_vertices": design[0], "design_faces": design[1],
            "fixed_material_definition": "OUTER_MINUS_DESIGN_SAME_DENSITY",
            "material_density": float(density), "target_origin": case["target_origin"],
            "target_axes": case["target_axes"], "target_length": length,
            "static_objective": objective, "static_pareto_weight": float(weight), "controls": controls,
            "static_objective_measure": measure}


def static_plane_initial_plan(problem, run, seed, max_restarts):
    if type(max_restarts) is not int or max_restarts <= 0:
        raise ValueError("static restart budget must be positive")
    points = target_domain_points(problem, problem["design_vertices"])
    length = problem["target_length"]
    center = [(min(p[i] for p in points)+max(p[i] for p in points))/(2*length) for i in range(3)]
    generator = random.Random(seed)
    maximum_mass = problem["static_objective"] == "MAXIMUM_MASS" or (
        problem["static_objective"] == "PARETO_TRADEOFF" and problem["static_pareto_weight"] == 1.)
    initials = []
    for _ in range(max_restarts):
        if maximum_mass:
            angle = generator.uniform(0., 2*math.pi)
            normal = [math.cos(angle), math.sin(angle), 0.]
            lower = min(sum(normal[i]*p[i] for i in range(3))/length for p in points)
            coefficients = [-lower/2 if lower < 0 else 1., *normal]
        else:
            normal = [generator.uniform(-.2,.2), generator.uniform(-.2,.2), -1.]
            coefficients = [-sum(normal[i]*center[i] for i in range(3)), *normal]
        norm = math.sqrt(sum(v*v for v in coefficients))
        initials.append([v/norm for v in coefficients])
    return {"run": run, "seed": seed, "initial_planes": initials,
            "initialization_protocol": "INPUT_BOUNDING_BOX_STATIC_PLANES_V1",
            "coefficient_coordinates": "DIMENSIONLESS_TARGET_X_OVER_L"}


def validate_static_plane_ledger(attempt):
    result = attempt.get("solver_result")
    if not isinstance(result, dict) or result.get("status") != attempt.get("status"):
        raise ValueError("static terminal status differs from native result")
    trace = result.get("trace")
    if not isinstance(trace, list) or not trace or len(trace) != attempt.get("function_evaluations"):
        raise ValueError("static native evaluation trace is missing or incomplete")
    controls = attempt["controls"]
    if len(trace) > 1+controls["maximum_iterations"]*controls["maximum_line_search_steps"]:
        raise ValueError("static iteration/line-search evaluation budget exceeded")
    for index, event in enumerate(trace):
        if type(event.get("evaluation")) is not int or event["evaluation"] != index:
            raise ValueError("static evaluation sequence is incomplete or duplicated")
        if type(event.get("accepted")) is not bool or type(event.get("evaluation_succeeded")) is not bool:
            raise ValueError("invalid static evaluation outcome")
        if event["accepted"] and not event["evaluation_succeeded"]:
            raise ValueError("failed static evaluation was accepted")
    if sum(row["accepted"] for row in trace) != result.get("accepted_steps"):
        raise ValueError("static accepted-step count differs from native trace")


def static_plane_attempt_acceptance(attempt):
    import numpy as np
    reasons = []
    try:
        validate_static_plane_ledger(attempt)
        result = attempt["solver_result"]
        static_plane_controls({"static_plane_settings": attempt["controls"]})
        tolerance = attempt["controls"]["constraint_tolerance"]
        if type(tolerance) not in (int,float) or not math.isfinite(tolerance) or not 0 < tolerance <= 1e-10:
            raise ValueError("invalid static acceptance tolerance")
        if attempt["status"] != "CONVERGED": reasons.append("NATIVE_CONVERGENCE_REQUIRED")
        total = np.asarray(result["total_raw_moments"], dtype=float)
        plane = np.asarray(result["plane_coefficients"], dtype=float)
        length, density = attempt["target_length"], attempt["material_density"]
        if total.shape != (10,) or plane.shape != (4,) or not np.isfinite(total).all() or not np.isfinite(plane).all():
            raise ValueError("nonfinite static terminal state")
        if (type(length) not in (int,float) or type(density) not in (int,float) or
                total[0] <= 0 or not math.isfinite(length*density) or min(length,density) <= 0):
            raise ValueError("invalid static normalization or mass")
        if abs(np.linalg.norm(plane)-1.) > 1e-12:
            raise ValueError("unnormalized static terminal plane")
        normalized = total[1:3]/(density*length**4)
        if np.linalg.norm(normalized) > tolerance: reasons.append("STATIC_BALANCE_TOLERANCE")
        terminal = next((row for row in reversed(result["trace"]) if row["accepted"]), result["trace"][0])
        traced_constraints = np.asarray(terminal["normalized_constraints"], dtype=float)
        if (traced_constraints.shape != (2,) or not np.isfinite(traced_constraints).all() or
                not terminal["evaluation_succeeded"] or not np.allclose(
                    normalized, traced_constraints, rtol=0, atol=1e-13)):
            reasons.append("STATIC_TERMINAL_TRACE_MISMATCH")
        objective, kappa = attempt["static_objective"], attempt["static_pareto_weight"]
        measure = attempt.get("static_objective_measure","DENSITY_WEIGHTED_MASS_V1")
        if measure not in ("DENSITY_WEIGHTED_MASS_V1","UNIFORM_DENSITY_VOLUME_V1"):
            raise ValueError("invalid static objective measure")
        mass_scale = 1/density if measure == "UNIFORM_DENSITY_VOLUME_V1" else 1.
        if (objective not in ("MAXIMUM_MASS", "LOWER_CENTER_OF_MASS", "PARETO_TRADEOFF") or
                type(kappa) not in (int,float) or not math.isfinite(kappa) or not 0 <= kappa <= 1 or
                (objective != "PARETO_TRADEOFF" and kappa != 0)):
            raise ValueError("invalid static objective or Pareto weight")
        maximum = objective == "MAXIMUM_MASS" or (objective == "PARETO_TRADEOFF" and kappa == 1.)
        traced_plane = np.asarray(terminal["plane_coefficients"], dtype=float)
        # The older maximum-mass trace stores physical lambda rather than the
        # four normalized coefficients; keep that historical representation.
        if maximum and np.linalg.norm(traced_plane) == 0.:
            traced_plane = np.array([1., *[v*length for v in terminal["lambda"]], 0.])
            traced_plane /= np.linalg.norm(traced_plane)
        if traced_plane.shape != (4,) or not np.allclose(plane, traced_plane, rtol=0, atol=1e-13):
            reasons.append("STATIC_TERMINAL_PLANE_MISMATCH")
        if maximum:
            if plane[0] <= 0 or abs(plane[3]) > 1e-12: reasons.append("STATIC_MAXIMUM_PLANE_KKT")
        else:
            if objective == "LOWER_CENTER_OF_MASS": kappa = 0.
            stationarity = abs(-plane[0]*(1-kappa)/total[0]-plane[3]/length*(
                kappa*mass_scale+(1-kappa)*total[3]/total[0]**2))
            witness = result.get("boundary_optimality")
            if result.get("material_region") == "STRICT_FULL" and witness is None:
                raise ValueError("typed full-domain state requires its boundary witness")
            if witness is not None:
                if (not isinstance(witness, dict) or result.get("material_region") != "STRICT_FULL" or
                    witness.get("kind") != "STRICT_FULL_DOMAIN_ZERO_EQUALITY_MULTIPLIERS_V1" or
                    witness.get("objective_measure") != measure):
                    raise ValueError("invalid full-domain witness identity")
                points = np.asarray(attempt["design_target_vertices"], dtype=float)
                if points.ndim != 2 or points.shape[1] != 3 or len(points) < 4 or not np.isfinite(points).all():
                    raise ValueError("missing full-domain geometry bounds")
                a, bz = kappa*mass_scale+(1-kappa)*total[3]/total[0]**2, -(1-kappa)/total[0]
                polynomial = a+bz*points[:,2]
                margin = 64*np.finfo(float).eps*np.max(abs(a)+abs(bz*points[:,2]))
                expected = np.array([a,0.,0.,length*bz]); expected /= np.linalg.norm(expected)
                recorded = np.asarray(witness["plane_coefficients"], dtype=float)
                minimum = witness["minimum_gradient_polynomial"]
                if (not np.isfinite(expected).all() or recorded.shape != (4,) or not np.isfinite(recorded).all() or
                    not np.allclose(recorded, expected, rtol=0, atol=1e-13) or
                    not np.min(polynomial) > margin or type(minimum) not in (int,float) or
                    not math.isfinite(minimum) or not math.isclose(minimum, float(np.min(polynomial)), rel_tol=1e-12, abs_tol=0.) or
                    not np.min(plane[0]+points@plane[1:]/length) > 1e-12):
                    reasons.append("STATIC_FULL_DOMAIN_OPTIMALITY")
            elif plane[3] >= -1e-12 or stationarity > tolerance: reasons.append("STATIC_OBJECTIVE_KKT")
    except (ValueError, KeyError, TypeError, IndexError, OverflowError) as error:
        reasons.append("INVALID_STATIC_EVIDENCE: "+str(error))
    return {"passed": not reasons, "reasons": reasons,
            "scope": "STATIC_PLANE_TERMINAL_NUMERICS_NOT_MATERIAL_OR_PAPER_ACCEPTANCE"}


def run_static_plane_case(case, starts, max_restarts, *, asset_root=ROOT, initial_seeds=None):
    from make_it_stand.spin import solve_spin_it_faster_static_plane_family
    problem = static_plane_problem(case, asset_root=asset_root)
    request = build_request(case, asset_root=asset_root)
    seeds = registered_seeds(starts, initial_seeds)
    plans = [static_plane_initial_plan(problem, run, seed, max_restarts) for run,seed in enumerate(seeds)]
    attempts = []
    for plan in plans:
        for restart, initial in enumerate(plan["initial_planes"], 1):
            began = time.perf_counter()
            result = _planar_json_value(solve_spin_it_faster_static_plane_family(
                request, initial_coefficients=initial, **problem["controls"]))
            row = {"run": plan["run"], "seed": plan["seed"], "restart": restart,
                   "solver_family": "STATIC_PLANE_FAMILY_V1", "initial_plane": initial,
                   "status": result["status"], "solver_result": result,
                   "controls": problem["controls"], "target_length": problem["target_length"],
                   "material_density": problem["material_density"],
                   "design_target_vertices": target_domain_points(problem, problem["design_vertices"]),
                   "static_objective": problem["static_objective"],
                   "static_pareto_weight": problem["static_pareto_weight"],
                   "static_objective_measure": problem["static_objective_measure"],
                   "function_evaluations": len(result["trace"]),
                   "evaluation_scope": "NATIVE_PLANE_EVALUATIONS_AND_ORIENTATION_PRECHECKS_NOT_EQ14_BUDGET",
                   "plane_iterations": result["iterations"], "elapsed_seconds": time.perf_counter()-began}
            row["numerical_acceptance"] = static_plane_attempt_acceptance(row)
            attempts.append(row)
            if row["status"] == "CONVERGED": break
    validate_attempt_ledger(attempts, expected_runs=starts, maximum_restarts=max_restarts)
    summary = summarize_attempts(attempts, starts)
    payload = {"id": case["id"], "mode": case["mode"], "solver_family": "STATIC_PLANE_FAMILY_V1",
               "asset_class": case.get("asset_class"), "seed_sequence": list(seeds),
               "initialization_protocol": "INPUT_BOUNDING_BOX_STATIC_PLANES_V1", "initialization_plan": plans,
               "static_problem": _planar_json_value(problem), "static_problem_sha256": scientific_hash(problem),
               "attempts": attempts, "summary": summary,
               "passed": summary["strict_terminal_numerics"]["all_runs_numerically_qualified"],
               "acceptance_scope": "STATIC_PLANE_TERMINAL_NUMERICS_NOT_MATERIAL_OR_PAPER_ACCEPTANCE",
               "material_validation": "PENDING_INDEPENDENT_SAVED_GEOMETRY_CHECK",
               "formal_paper_admission": "NOT_GRANTED_BY_METHOD_RUNNER"}
    payload["scientific_hash"] = scientific_hash(payload)
    return payload


def _cherchi_runner():
    """Reuse the audited geometry/solver/seal path, not another implementation."""
    path=ROOT.parent/"strict/run_faster_cherchi_newton_integration.py"
    name="run_faster_cherchi_newton_integration"
    if name in sys.modules and Path(sys.modules[name].__file__).resolve()==path.resolve():
        return sys.modules[name]
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    sys.path.insert(0,str(path.parent))
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.pop(0)
    sys.modules[name]=module
    return module


def cherchi_case_plan(case, starts, max_restarts, *, asset_root=ROOT, initial_seeds=None):
    """Bind real regions and native seeded starts before any numerical run."""
    seeds=registered_seeds(starts,initial_seeds)
    if type(max_restarts) is not int or not 1<=max_restarts<=2**31-1:
        raise ValueError("Cherchi restart budget must be a positive signed integer")
    if case.get("mode") not in {"SPINNING_TOP","YO_YO"}:
        raise ValueError("Cherchi Table route supports spinning tops and yo-yos")
    identity=case.get("id")
    if not isinstance(identity,str) or not identity or identity in {".",".."} or any(
            not (c.isascii() and (c.isalnum() or c in "_-")) for c in identity):
        raise ValueError("case id must be a safe local directory name")
    if any(key in case for key in ("initial_quadric","initial_quadrics","initial_coefficients",
                                   "initialization_mode","checkpoint","known_solution")):
        raise ValueError("Cherchi application initialization must be input-and-seed only")
    entry=case.get("cherchi_problem")
    if not isinstance(entry,dict) or set(entry)!={"path","sha256","assets"}:
        raise ValueError("Cherchi Table requires a hashed cherchi_problem and all material assets")
    root=Path(asset_root).resolve()
    def bound_file(item):
        if not isinstance(item,dict) or not {"path","sha256"}.issubset(item):
            raise ValueError("invalid Cherchi asset identity")
        relative=Path(item["path"])
        path=(root/relative).resolve()
        if relative.is_absolute() or not path.is_relative_to(root):
            raise ValueError("Cherchi asset path must stay inside manifest directory")
        data=path.read_bytes()
        if hashlib.sha256(data).hexdigest()!=item["sha256"]:
            raise ValueError("Cherchi asset SHA-256 mismatch")
        return path,data
    path,data=bound_file(entry)
    problem=json.loads(data)
    if problem.get("schema")!="faster-cherchi-problem-manifest-v1" or problem.get("target",{}).get("mode")!=case["mode"]:
        raise ValueError("Cherchi problem mode/schema differs from case")
    assets=entry["assets"]
    if not isinstance(assets,list):raise ValueError("Cherchi material assets must be a list")
    indexed={}
    for asset in assets:
        resolved,_=bound_file(asset)
        if resolved in indexed:raise ValueError("duplicate Cherchi material asset")
        indexed[resolved]=asset["sha256"]
    try:
        referenced={ (path.parent/mesh).resolve() for mesh in [problem["outer_mesh"],
            *(region["mesh"] for kind in ("fixed_regions","design_regions") for region in problem[kind])] }
    except (KeyError,TypeError) as error:
        raise ValueError("Cherchi problem omits actual material geometry") from error
    if referenced!=set(indexed):raise ValueError("Cherchi asset list must bind every and only referenced material mesh")
    for key in ("outer","design_domain","material_density","fixed_shell_moments","target_origin","target_axes"):
        if key in case:raise ValueError("Cherchi geometry/coordinates/densities must come only from locked problem")
    settings={"maximum_function_evaluations":60,"icosphere_subdivisions":5,
              "maximum_radius_doublings":6,"cherchi_timeout_seconds":180.}
    supplied=case.get("cherchi_settings",{})
    if not isinstance(supplied,dict) or not set(supplied)<=set(settings)|{"paper_projective_hyperboloid","algebraic_finite_map"}:
        raise ValueError("unknown Cherchi settings; normal discretization is fixed during each solve")
    if "paper_projective_hyperboloid" in supplied and type(supplied["paper_projective_hyperboloid"]) is not bool:
        raise ValueError("paper_projective_hyperboloid must be a boolean")
    if "algebraic_finite_map" in supplied and type(supplied["algebraic_finite_map"]) is not bool:
        raise ValueError("algebraic_finite_map must be a boolean")
    settings.update(supplied)
    for key,low,high in (("maximum_function_evaluations",1,60),("icosphere_subdivisions",0,7),
                          ("maximum_radius_doublings",0,16)):
        if type(settings[key]) is not int or not low<=settings[key]<=high:
            raise ValueError("invalid Cherchi "+key)
    timeout=settings["cherchi_timeout_seconds"]
    if type(timeout) not in (int,float) or not math.isfinite(timeout) or not 1<=timeout<=300:
        raise ValueError("invalid Cherchi timeout")
    weight=case.get("center_of_mass_weight",0.)
    if type(weight) not in (int,float) or not math.isfinite(weight) or weight<0 or (case["mode"]!="SPINNING_TOP" and weight!=0):
        raise ValueError("invalid Cherchi center_of_mass_weight")
    from make_it_stand.verification.faster_material import rotational_objective_definition
    selector=case.get("spinning_top_objective","INERTIA_AND_LOWERING")
    objective=rotational_objective_definition(case["mode"],weight,selector)
    return {"id":identity,"mode":case["mode"],"solver_family":"CHERCHI_RAW_EQ14_V1",
        "problem":problem,"problem_identity":entry,"settings":settings,
        "center_of_mass_weight":weight,
        "spinning_top_objective":selector,"objective_definition":objective,
        "initialization_protocol":"NATIVE_MODE_CONSTRAINED_SEEDS_V1",
        "runs":[{"run":i,"seed":seed,"initialization_mode":"deterministic_mode_seed",
                 "maximum_restarts":max_restarts} for i,seed in enumerate(seeds)]}


def cherchi_attempt_rows(record, *, run, seed, maximum_restarts, maximum_evaluations):
    """Split a complete native multi-restart ledger without changing its indices."""
    events=record.get("events"); trace=record.get("trace")
    if not isinstance(events,list) or not events or not isinstance(trace,list):
        raise ValueError("Cherchi native event/trace ledger is missing")
    count=record.get("restart_count")
    counts=record.get("restart_function_evaluations")
    if type(count) is not int or not 1<=count<=maximum_restarts or not isinstance(counts,list) or len(counts)!=count:
        raise ValueError("Cherchi restart count differs from complete budget ledger")
    if any(type(n) is not int or not 0<=n<=maximum_evaluations for n in counts) or sum(counts)!=record.get("function_evaluations"):
        raise ValueError("Cherchi function budget differs from complete ledger")
    if any(type(e.get("restart")) is not int or not 0<=e["restart"]<count for e in events):
        raise ValueError("Cherchi event has an invalid restart")
    if [e["restart"] for e in events]!=sorted(e["restart"] for e in events):
        raise ValueError("Cherchi event order regresses across restarts")
    rows=[]; offset=0
    for i,n in enumerate(counts):
        subset=[e for e in events if e["restart"]==i]
        history=[t for t in trace if t.get("restart")==i]
        terminal=history[-1] if history else {}
        reason=subset[-1].get("reason") if subset else None
        status=("CONVERGED" if reason=="converged" else
                "RESOURCE_LIMIT" if reason=="maximum_evaluations_reached" else
                "OPTIMIZATION_DID_NOT_CONVERGE")
        row={"run":run,"seed":seed,"restart":i+1,"native_restart":i,
            "native_evaluation_offset":offset,"solver_family":"CHERCHI_RAW_EQ14_V1",
            "status":status,"native_terminal_reason":subset[-1].get("reason") if subset else None,
            "function_evaluations":n,"maximum_function_evaluations":maximum_evaluations,
            "equation_profile":record.get("equation_profile"),
            "initial_quadric":subset[0].get("quadric_coefficients") if subset else None,
            "newton_iterations":max((t["iteration"] for t in history),default=0),
            "trace":history,"solver_events":subset,
            **{key:terminal.get(key) for key in ("raw_equation_residual","normalized_equation_residual",
                                               "physical_relative_residual","multiplier_mu")}}
        row["numerical_acceptance"]=numerical_attempt_acceptance(row)
        rows.append(row);offset+=n
    if any(type(t.get("restart")) is not int or not 0<=t["restart"]<count for t in trace):
        raise ValueError("Cherchi trace has an invalid restart")
    if (record.get("solver_status")=="CONVERGED")!=(rows[-1]["status"]=="CONVERGED"):
        raise ValueError("Cherchi terminal status contradicts event ledger")
    return rows


def run_cherchi_case(case,starts,max_restarts,bundle_root=None,asset_root=ROOT,*,
                     initial_seeds=None,cherchi_build=None,cherchi_probe=None):
    plan=cherchi_case_plan(case,starts,max_restarts,asset_root=asset_root,initial_seeds=initial_seeds)
    if bundle_root is None or cherchi_build is None or cherchi_probe is None:
        raise ValueError("Cherchi applications require explicit bundle_root, cherchi_build and cherchi_probe")
    runner=_cherchi_runner()
    root=Path(bundle_root).resolve()/case["id"]
    root.mkdir(parents=True,exist_ok=False)
    path=(Path(asset_root)/case["cherchi_problem"]["path"]).resolve()
    attempts=[];runs=[]
    for planned in plan["runs"]:
        directory=root/f"run-{planned['run']:06d}"
        record=runner.run(directory,build=Path(cherchi_build),probe=Path(cherchi_probe),
            problem_manifest=path,maximum_restarts=max_restarts,
            initialization_mode="deterministic_mode_seed",deterministic_seed=planned["seed"],
            center_of_mass_weight=plan["center_of_mass_weight"],
            spinning_top_objective=plan["spinning_top_objective"],**plan["settings"])
        # Preflight infeasibility is retained separately, never an invented
        # Newton attempt. Normal failures retain every attempted restart.
        if record["status"]!="INFEASIBLE_CERTIFIED":
            rows=cherchi_attempt_rows(record,run=planned["run"],seed=planned["seed"],
                maximum_restarts=max_restarts,maximum_evaluations=plan["settings"]["maximum_function_evaluations"])
            local=[{**row,"run":0} for row in rows]
            validate_attempt_ledger(local,expected_runs=1,maximum_restarts=max_restarts)
            attempts.extend(rows)
        runs.append({**planned,"status":record["status"],"solver_status":record["solver_status"],
            "relative_output":str(directory.relative_to(Path(bundle_root).resolve())),
            "function_evaluations":record["function_evaluations"],
            "geometry_callback_evaluations":len(record.get("evaluations",[])),
            "scientific_event_hash_v2":record.get("scientific_event_hash_v2"),
            "trace_hash":record.get("trace_hash"),
            "sealed_material_digest":record.get("sealed_artifact",{}).get("scientific_digest"),
            "locked_problem_digest":record.get("locked_problem_digest"),
            "preflight_input_hash":record.get("input_hash")})
    # Re-read input bindings after the complete run; never package stale bytes.
    if plan!=cherchi_case_plan(case,starts,max_restarts,asset_root=asset_root,initial_seeds=initial_seeds):
        raise ValueError("STALE_RESULT: Cherchi application input changed")
    summary=summarize_attempts(attempts,starts)
    payload={"id":case["id"],"mode":case["mode"],"solver_family":"CHERCHI_RAW_EQ14_V1",
        "asset_class":case.get("asset_class"),"initialization_protocol":plan["initialization_protocol"],
        "seed_sequence":[r["seed"] for r in plan["runs"]],"initialization_plan":plan,
        "attempts":attempts,"runs":runs,"summary":summary,
        "passed":summary["strict_terminal_numerics"]["all_runs_numerically_qualified"] and
                 all(r["status"]=="CONVERGED" and r["sealed_material_digest"] is not None for r in runs),
        "acceptance_scope":"CHERCHI_RUN_REPORTED_METHOD_STATUS_PENDING_INDEPENDENT_RELOAD",
        "material_validation":"SAVED_CHERCHI_RUN_RECORDS_REQUIRE_INDEPENDENT_RELOAD",
        "formal_paper_admission":"NOT_GRANTED_BY_METHOD_RUNNER"}
    payload["scientific_hash"]=scientific_hash(payload)
    return payload


def run_case(
    case: dict[str, Any],
    starts: int,
    max_restarts: int,
    bundle_root: Path | None = None,
    asset_root: Path = ROOT,
    *, initial_seeds: tuple[int, ...] | None = None,
    cherchi_build: Path | None = None, cherchi_probe: Path | None = None,
) -> dict[str, Any]:
    if case["mode"] == "BUOYANT_STABILITY":
        from make_it_stand.verification.buoyant_application import run_buoyant_case
        return run_buoyant_case(case, starts, max_restarts, bundle_root,
                               registered_seeds(starts, initial_seeds))
    if case["mode"] == "FIDGET_SPINNER":
        return run_planar_case(case, starts, max_restarts, initial_seeds=initial_seeds)
    if case["mode"] == "STATIC_BALANCE":
        return run_static_plane_case(case, starts, max_restarts, asset_root=asset_root, initial_seeds=initial_seeds)
    if case["mode"] in {"SPINNING_TOP","YO_YO"}:
        return run_cherchi_case(case,starts,max_restarts,bundle_root,asset_root,
            initial_seeds=initial_seeds,cherchi_build=cherchi_build,cherchi_probe=cherchi_probe)
    return run_legacy_quadric_case(case, starts, max_restarts, bundle_root, asset_root,
                                   initial_seeds=initial_seeds)


def run_legacy_quadric_case(case, starts, max_restarts, bundle_root=None, asset_root=ROOT,
                             *, initial_seeds=None):
    """Explicit historical regression entry; never dispatch static methods here."""
    from make_it_stand.spin import solve_spin_it_faster

    request = build_request(case, asset_root=asset_root)
    design = read_case_mesh(case, "design_domain", asset_root=asset_root)
    domain_points = target_domain_points(case, design[0])
    seeds = registered_seeds(starts,initial_seeds)
    attempts: list[dict[str, Any]] = []
    viewer_bundle: dict[str, str] | None = None
    for run in range(starts):
        initial_plan = initial_run_plan(case,run,max_restarts,domain_points,seed=seeds[run])
        initials = initial_plan["initial_quadrics"]
        for restart, initial in enumerate(initials, 1):
            begin = time.perf_counter()
            result = solve_spin_it_faster(request, initial_quadrics=(initial,))
            elapsed = time.perf_counter() - begin
            trace = result.trace
            status = result.status
            if (
                status != "CONVERGED"
                and trace
                and trace[-1].normalized_equation_residual <= request.settings.equation_tolerance
                and trace[-1].multiplier_mu <= 0.0
            ):
                status = "LOCAL_MAXIMUM"
            attempts.append({
                "run": run,
                "seed": initial_plan["seed"],
                "restart": restart,
                "status": status,
                "newton_iterations": max((entry.iteration for entry in trace), default=0),
                # Accepted-step traces omit rejected initial/trial evaluations.
                # Count the native solve, then cross-check the full event ledger.
                "function_evaluations": result.solve_function_evaluations,
                "normalized_equation_residual": result.diagnostics["normalized_equation_residual"],
                "raw_equation_residual": result.raw_equation_residual,
                "equation_profile": result.equation_profile,
                "multiplier_mu": trace[-1].multiplier_mu if trace else None,
                "physical_relative_residual": result.diagnostics["physical_relative_residual"],
                "output_hash": result.diagnostics["output_hash"],
                "initial_quadric": initial,
                "elapsed_seconds": elapsed,
                "residual_seconds": result.residual_evaluation_seconds,
                "derivative_seconds": result.derivative_evaluation_seconds,
                "trace": [
                    {
                        "iteration": entry.iteration,
                        "restart": entry.restart,
                        "function_evaluations": entry.function_evaluations,
                        "raw_equation_residual": entry.raw_equation_residual,
                        "multiplier_mu": entry.multiplier_mu,
                        "equation_profile": entry.equation_profile,
                        "solver_coordinate_identity": entry.solver_coordinate_identity,
                        "normalized_equation_residual": entry.normalized_equation_residual,
                        "physical_relative_residual": entry.physical_relative_residual,
                        "accepted_geodesic_step": entry.accepted_geodesic_step,
                    }
                    for entry in trace
                ],
                "solver_events": [solver_event_payload(event) for event in result.events],
            })
            if result.status == "CONVERGED":
                if bundle_root is not None and viewer_bundle is None:
                    directory = bundle_root / case["id"]
                    paths = result.save_viewer_bundle(directory)
                    viewer_bundle = {
                        Path(path).name: hashlib.sha256(Path(path).read_bytes()).hexdigest()
                        for path in paths
                    }
                break
    validate_attempt_ledger(attempts, expected_runs=starts, maximum_restarts=max_restarts)
    summary = summarize_attempts(attempts, starts)
    converged = [row for row in attempts if row["status"] == "CONVERGED"]
    for attempt in attempts:
        attempt['numerical_acceptance'] = numerical_attempt_acceptance(attempt)
    passed = (
        summary["all_runs_converged"]
        and len({int(row["run"]) for row in attempts}) == starts
        and all(row['numerical_acceptance']['passed'] for row in converged)
    )
    # Preserve the old numerical regression verdict separately. It cannot
    # become strict Eq.14/material/paper approval by renaming its fields.
    legacy_numerical_passed = (
        summary['all_runs_converged']
        and len({int(row['run']) for row in attempts}) == starts
        and all(row['normalized_equation_residual'] <= 1e-12 and
                row['physical_relative_residual'] <= 1e-6 for row in converged)
    )
    payload = {
        "id": case["id"],
        "mode": case["mode"],
        "asset_class": "RECONSTRUCTED_ASSET",
        "initialization_protocol": "MODE_CONSTRAINED_TOP_LEVEL_SEEDS_V1",
        "seed_sequence": list(seeds),
        "attempts": attempts,
        "summary": summary,
        "viewer_bundle": viewer_bundle,
        "passed": passed,
        "legacy_numerical_passed": legacy_numerical_passed,
        "legacy_acceptance_scope": "HISTORICAL_NORMALIZED_NUMERICS_NOT_PAPER_APPROVAL",
        "acceptance_scope": "TERMINAL_NUMERICS_ONLY_NOT_MATERIAL_OR_PAPER_ACCEPTANCE",
        "material_validation": "NOT_VERIFIED_BY_THIS_LEGACY_DISPATCH",
        "formal_paper_admission": "NOT_GRANTED_RECONSTRUCTED_ASSET",
    }
    payload["scientific_hash"] = scientific_hash(payload)
    return payload


def application_metadata(cases):
    protocols = {"INPUT_BOUNDING_BOX_PLANAR_CENTERS_V1" if case["mode"] == "FIDGET_SPINNER"
                 else "LOCKED_EXPLICIT_BUOYANT_PLANES_V1" if case["mode"] == "BUOYANT_STABILITY"
                 else "INPUT_BOUNDING_BOX_STATIC_PLANES_V1" if case["mode"] == "STATIC_BALANCE"
                 else "NATIVE_MODE_CONSTRAINED_SEEDS_V1" if case["mode"] in {"SPINNING_TOP","YO_YO"}
                 else "MODE_CONSTRAINED_TOP_LEVEL_SEEDS_V1" for case in cases}
    classes = {case.get("asset_class", "RECONSTRUCTED_ASSET") for case in cases}
    return {"initialization_protocol": next(iter(protocols)) if len(protocols) == 1
            else "PER_CASE_PROTOCOL_V1",
            "asset_class": next(iter(classes)) if len(classes) == 1 else "PER_CASE_ASSET_CLASS"}


def plan_payload(
    cases: list[dict[str, Any]], starts: int, *, asset_root: Path = ROOT,
    max_restarts: int = 8,
    initial_seeds: tuple[int, ...] | None = None,
) -> dict[str, Any]:
    seeds = registered_seeds(starts,initial_seeds)
    plans = []
    for case in cases:
        if case["mode"] == "BUOYANT_STABILITY":
            from make_it_stand.verification.buoyant_application import buoyant_case_plan
            if "buoyant_problem" not in case:
                plans.append({"id": case["id"], "mode": case["mode"],
                    "solver_family": "BUOYANT_PLANE_GEOMETRY_V1", "runs": [],
                    "status": "BLOCKED_MISSING_BUOYANT_GEOMETRY",
                    "blocker": "buoyant_problem with actual displacement and fixed geometry is required"})
            else:
                plans.append(buoyant_case_plan(case, starts, max_restarts, seeds))
            continue
        if case["mode"] in {"SPINNING_TOP","YO_YO"}:
            if "cherchi_problem" not in case:
                plans.append({"id":case["id"],"mode":case["mode"],"solver_family":"CHERCHI_RAW_EQ14_V1",
                    "status":"BLOCKED_MISSING_CHERCHI_PROBLEM","runs":[],
                    "blocker":"locked actual outer/fixed/design material geometry is required"})
            else:
                plans.append(cherchi_case_plan(case,starts,max_restarts,asset_root=asset_root,initial_seeds=initial_seeds))
            continue
        if case["mode"] == "STATIC_BALANCE":
            problem = static_plane_problem(case, asset_root=asset_root)
            plans.append({"id": case["id"], "mode": case["mode"], "solver_family": "STATIC_PLANE_FAMILY_V1",
                "static_problem": _planar_json_value(problem), "static_problem_sha256": scientific_hash(problem),
                "runs": [static_plane_initial_plan(problem, run, seed, max_restarts)
                         for run,seed in enumerate(seeds)]})
            continue
        if case["mode"] == "FIDGET_SPINNER":
            if "planar_problem" not in case:
                plans.append({"id": case["id"], "mode": case["mode"],
                              "solver_family": "PLANAR_DOGLEG_V1", "runs": [],
                              "status": "BLOCKED_MISSING_PLANAR_GEOMETRY",
                              "blocker": "planar_problem with actual layers and design contours is required"})
                continue
            request = planar_case_request(case)
            plans.append({"id": case["id"], "mode": case["mode"],
                          "solver_family": "PLANAR_DOGLEG_V1",
                          "planar_problem": request,
                          "planar_problem_sha256": scientific_hash(request),
                          "runs": [planar_initial_run_plan(request, run, seed, max_restarts)
                                   for run, seed in enumerate(seeds)]})
            continue
        read_case_mesh(case, "outer", asset_root=asset_root)
        domain_points = target_domain_points(
            case,read_case_mesh(case,"design_domain",asset_root=asset_root)[0])
        runs = [initial_run_plan(case,run,max_restarts,domain_points,seed=seeds[run])
                for run in range(starts)]
        plans.append({"id":case["id"],"mode":case["mode"],"runs":runs,
                      "initial_quadrics":[row["initial_quadrics"][0] for row in runs],
                      "initial_quadrics_scope":"FIRST_ATTEMPT_OF_EACH_TOP_LEVEL_RUN"})
    return {
        "schema": "spin-it-faster-table1-plan-v2",
        "schema_version": 2,
        "starts_per_case": starts,
        "maximum_restarts_per_run": max_restarts,
        "seed_sequence": list(seeds),
        **application_metadata(cases),
        "author_exact": False,
        "cases": plans,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--starts", type=int, default=50)
    parser.add_argument("--max-restarts", type=int, default=8)
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--bundle-root", type=Path)
    parser.add_argument("--cherchi-build", type=Path)
    parser.add_argument("--cherchi-probe", type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    cases = [case for case in manifest["cases"] if not args.case or case["id"] in args.case]
    if args.case and len(cases) != len(set(args.case)):
        raise SystemExit("one or more requested case IDs are absent")
    if args.starts <= 0 or args.max_restarts <= 0:
        raise SystemExit("starts and max-restarts must be positive")
    if args.plan_only:
        payload = plan_payload(cases, args.starts, asset_root=args.manifest.parent,
                               max_restarts=args.max_restarts)
    else:
        rows = [
            run_case(case, args.starts, args.max_restarts, args.bundle_root,
                     asset_root=args.manifest.parent,cherchi_build=args.cherchi_build,
                     cherchi_probe=args.cherchi_probe)
            for case in cases
        ]
        payload = {
            "schema": "spin-it-faster-table1-results-v1",
            "schema_version": 1,
            "paper_doi": "10.1145/3658194",
            "starts_per_case": args.starts,
            "maximum_restarts_per_run": args.max_restarts,
            **application_metadata(cases),
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
