#!/usr/bin/env python3
"""Generate auditable data for Spin-It Faster Figures 10--13 and Sec. 6.6."""

from __future__ import annotations

import argparse
from dataclasses import replace
import importlib.util
import hashlib
import json
import math
from pathlib import Path
import sys
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
            for residuals in _raw_residual_segments(attempt):
                for current, following in zip(residuals, residuals[1:]):
                    if following > current * (1.0 + 1e-12):
                        raise ValueError("accepted Newton residuals must be monotone")
                    if current > 0.0 and following > 0.0:
                        grouped[group].append((current, following))
    return grouped


def _raw_residuals(attempt: dict[str, Any]) -> list[float]:
    residuals = []
    for row in attempt.get("trace", []):
        value = row.get("raw_equation_residual")
        if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
            raise ValueError("Figure 10 requires finite nonnegative raw Eq.14 residuals; no normalized fallback")
        residuals.append(float(value))
    return residuals


def _raw_residual_segments(attempt: dict[str, Any]) -> list[list[float]]:
    values = _raw_residuals(attempt)
    segments = []
    previous_restart = None
    for row, value in zip(attempt.get("trace", []), values):
        restart = row.get("restart", 0)
        if type(restart) is not int or restart < 0:
            raise ValueError("Figure 10 restart identities must be nonnegative integers")
        if previous_restart is not None and restart < previous_restart:
            raise ValueError("Figure 10 restart identities cannot regress or be reused")
        if not segments or restart != previous_restart:
            segments.append([])
        segments[-1].append(value)
        previous_restart = restart
    return segments


def _fixed_branch_terminal_order(attempt: dict[str, Any], noise_floor: float) -> dict[str, Any]:
    """Report only a contiguous suffix on the endpoint's recorded stratum.

    Cherchi's material_topology_hash binds construction, CSG side, source type,
    discretization and cap branch. This is a journal diagnostic, not an audit of
    those underlying objects. Missing historical identities remain unknown.
    """
    rows = attempt.get("trace", [])
    for row in rows:
        identity = row.get("material_topology_hash")
        if identity is not None and (not isinstance(identity, str) or len(identity) != 64 or
                any(char not in "0123456789abcdef" for char in identity)):
            raise ValueError("Figure 10 material branch identity must be a SHA-256 or absent")
    endpoint = rows[-1] if rows else {}
    identity = endpoint.get("material_topology_hash")
    restart = endpoint.get("restart", 0) if rows else None
    suffix = []
    if identity is not None:
        for row in reversed(rows):
            if row.get("restart", 0) != restart or row.get("material_topology_hash") != identity:
                break
            suffix.append(row)
        suffix.reverse()
    # Select the branch BEFORE dropping noise samples. A switch even below
    # the noise floor cannot borrow convergence orders from the old branch.
    while suffix and suffix[-1]["raw_equation_residual"] <= noise_floor:
        suffix.pop()
    contiguous = suffix[-1:] if suffix else []
    for row in reversed(suffix[:-1]):
        if row["raw_equation_residual"] <= contiguous[0]["raw_equation_residual"]:
            break
        contiguous.insert(0, row)
        if len(contiguous) == 7:  # Same predeclared five-triple reporting window.
            break
    values = [float(row["raw_equation_residual"]) for row in contiguous]
    orders = [(math.log(c) - math.log(b)) / (math.log(b) - math.log(a))
              for a, b, c in zip(values, values[1:], values[2:])]
    # Keep the same finite-order convention without selecting isolated triples.
    if not all(math.isfinite(value) for value in orders):
        orders = []
    ordered = sorted(orders)
    count = len(orders)
    median = (ordered[count // 2] if count % 2 else
              (ordered[count // 2 - 1] + ordered[count // 2]) / 2) if count else None
    return {
        "schema": "faster-recorded-fixed-branch-order-v1",
        "status": "ESTIMATED" if count >= 2 else "INSUFFICIENT_EVIDENCE",
        "restart": restart, "material_topology_hash": identity,
        "identity_availability": "RECORDED" if identity is not None else "UNAVAILABLE",
        "selection": "LAST_FIVE_CONTIGUOUS_ABOVE_NOISE_TRIPLES_ON_ENDPOINT_RECORDED_STRATUM",
        "maximum_triples": 5, "usable_triples": count,
        "orders": orders, "median_order": median,
        "sample_evaluation_indices": [row.get("evaluation_index") for row in contiguous],
        "material_identity_attestation": "NOT_VERIFIED_BY_FIGURE_CHECKER",
        "scope": "RECORDED_BRANCH_NUMERICS_NOT_QUADRATIC_CONVERGENCE_OR_MATERIAL_CERTIFICATE",
    }


def convergence_order_evidence(
    table: dict[str, Any], *, noise_floor: float = 1e-12,
) -> dict[str, list[dict[str, Any]]]:
    """Local q-order estimates from contiguous triples above declared noise.

    An estimate is numerical evidence, not automatic paper reproduction or an
    asymptotic convergence proof. Missing/short traces remain insufficient.
    """
    if type(noise_floor) not in (int, float) or not math.isfinite(noise_floor) or noise_floor <= 0:
        raise ValueError("noise floor must be finite and positive")
    # Validate the same raw values and monotonicity as the plotted pairs.
    quadratic_convergence_pairs(table)
    grouped = {name: [] for name in APPLICATION_GROUP.values()}
    for case in table.get("cases", []):
        group = APPLICATION_GROUP.get(case.get("mode"))
        if group is None:
            continue
        for index, attempt in enumerate(case.get("attempts", [])):
            if attempt.get("status") != "CONVERGED":
                continue
            orders = []
            excluded = 0
            segments = _raw_residual_segments(attempt)
            terminal_orders = []
            # Fixed reporting window, not a fitted cutoff or convergence gate.
            # Earlier failed restarts must not establish the terminal q-order.
            maximum_terminal_triples = 5
            for segment_index, values in enumerate(segments):
                triples = list(zip(values, values[1:], values[2:]))
                # Drop the terminal numerical-noise suffix before selecting
                # the late window. Do not bridge gaps or cherry-pick orders.
                above_noise_end = len(triples)
                while above_noise_end and triples[above_noise_end - 1][2] <= noise_floor:
                    above_noise_end -= 1
                late_start = max(0, above_noise_end - maximum_terminal_triples)
                for triple_index, (a, b, c) in enumerate(triples):
                    in_terminal_window = (segment_index == len(segments) - 1 and
                        late_start <= triple_index < above_noise_end)
                    if not a > b > c > noise_floor:
                        excluded += 1
                        if in_terminal_window:
                            terminal_orders.clear()
                        continue
                    # Log differences avoid underflow in ratios near the noise floor.
                    order = (math.log(c) - math.log(b)) / (math.log(b) - math.log(a))
                    if math.isfinite(order):
                        orders.append(order)
                        if in_terminal_window:
                            terminal_orders.append(order)
                    else:
                        excluded += 1
                        if in_terminal_window:
                            terminal_orders.clear()
            ordered = sorted(orders)
            count = len(ordered)
            median = (ordered[count // 2] if count % 2 else
                      (ordered[count // 2 - 1] + ordered[count // 2]) / 2) if count else None
            terminal_sorted = sorted(terminal_orders)
            terminal_count = len(terminal_sorted)
            terminal_median = (terminal_sorted[terminal_count // 2] if terminal_count % 2 else
                (terminal_sorted[terminal_count // 2 - 1] + terminal_sorted[terminal_count // 2]) / 2
            ) if terminal_count else None
            grouped[group].append({
                "case_id": case.get("id"), "attempt_index": index,
                "status": "ESTIMATED" if count >= 2 else "INSUFFICIENT_EVIDENCE",
                "residual_kind": "RAW_EQ14", "noise_floor": noise_floor,
                "orders": orders, "median_order": median,
                "usable_triples": count, "excluded_triples": excluded,
                "fixed_branch_terminal_order_evidence": _fixed_branch_terminal_order(attempt, noise_floor),
                "terminal_restart_order_evidence": {
                    "restart": attempt["trace"][-1].get("restart", 0) if segments else None,
                    "status": "ESTIMATED" if terminal_count >= 2 else "INSUFFICIENT_EVIDENCE",
                    "selection": "LAST_FIVE_CONTIGUOUS_TRIPLES_ABOVE_NOISE_IN_TERMINAL_RESTART",
                    "maximum_triples": maximum_terminal_triples,
                    "orders": terminal_orders, "median_order": terminal_median,
                    "usable_triples": terminal_count,
                    "scope": "LOCAL_NUMERICAL_EVIDENCE_NOT_QUADRATIC_CONVERGENCE_CERTIFICATE",
                },
            })
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
    """Legacy helper for independently known monotone analytic predicates only.

    A Newton success/failure predicate is not an admissible feasibility oracle.
    This helper is not used by the Figure 13 experiment entry.
    """
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


def refine_feasibility_bracket(
    oracle: Callable[[float], Any], *, lower: float, upper: float,
    tolerance: float = 1e-11, maximum_steps: int = 128,
) -> dict[str, Any]:
    """Pointwise continuous proof protocol, not a global critical-boundary proof.

    The caller must separately establish family identity and monotonicity before
    using any returned interval as a Figure 13 critical boundary. UNKNOWN never
    updates an endpoint. Discrete LP proofs cannot certify a continuous family.
    """
    from make_it_stand.verification.feasibility import FeasibilityReport

    if not (math.isfinite(lower) and math.isfinite(upper) and lower < upper):
        raise ValueError("feasibility bracket must be finite and ordered")
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("feasibility tolerance must be finite and positive")
    if type(maximum_steps) is not int or maximum_steps <= 0:
        raise ValueError("maximum steps must be a positive integer")
    observations = []

    def observe(distance):
        report = oracle(distance)
        if not isinstance(report, FeasibilityReport):
            raise TypeError("oracle must return a checker-backed FeasibilityReport")
        observations.append({"distance": distance, "evidence": report.to_dict()})
        # Unknown reports may retain their default discrete scope; never promote.
        if report.scope not in {
            "CONTINUOUS_GEOMETRIC_MOMENT_CAPACITY",
            "CONTINUOUS_GEOMETRIC_SYMMETRY_WITNESS",
        }:
            return "UNKNOWN"
        return report.status

    def finish(status, reason):
        return {"status": status, "reason": reason, "lower": lower,
                "upper": upper, "interval_width": upper - lower,
                "critical_distance": None, "family_monotonicity": "NOT_VERIFIED",
                "scope": "POINTWISE_CONTINUOUS_PROOF_PROTOCOL",
                "observations": observations}

    low_status, high_status = observe(lower), observe(upper)
    if low_status != "FEASIBLE_WITNESS" or high_status != "INFEASIBLE_CERTIFIED":
        return finish("UNKNOWN", "CONTINUOUS_ENDPOINT_PROOFS_REQUIRED")
    for _ in range(maximum_steps):
        if upper - lower <= tolerance:
            return finish("POINTWISE_BRACKET_WIDTH_REACHED", "FAMILY_PROOF_STILL_REQUIRED")
        midpoint = lower + (upper - lower) / 2
        if midpoint == lower or midpoint == upper:
            return finish("UNKNOWN", "FLOAT_RESOLUTION_LIMIT")
        status = observe(midpoint)
        if status == "FEASIBLE_WITNESS":
            lower = midpoint
        elif status == "INFEASIBLE_CERTIFIED":
            upper = midpoint
        else:
            return finish("UNKNOWN", "MIDPOINT_PROOF_UNAVAILABLE")
    return finish("UNKNOWN", "REFINEMENT_BUDGET_EXHAUSTED")


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


def run_pareto_sweep(case: dict[str, Any], weights: list[float], *, asset_root=ROOT,
                     bundle_root=None, cherchi_build=None, cherchi_probe=None,
                     initial_seeds=(0, 1), max_restarts=8) -> dict[str, Any]:
    if case["mode"] != "STATIC_BALANCE":
        return run_cherchi_top_pareto_sweep(case, weights, asset_root=asset_root,
            bundle_root=bundle_root, cherchi_build=cherchi_build, cherchi_probe=cherchi_probe,
            initial_seeds=initial_seeds, max_restarts=max_restarts)
    if (not weights or len(set(weights)) != len(weights) or
            any(type(w) not in (int,float) or not math.isfinite(w) or not 0 <= w <= 1 for w in weights)):
        raise ValueError("static Pareto weights must be distinct finite values in [0,1]")
    table = _table_module()
    rows = []
    measure = case.get("static_objective_measure","UNIFORM_DENSITY_VOLUME_V1")
    # Independent pure endpoints; no result, successful start or previous
    # weight is supplied to the next optimization. All failures are retained.
    for weight in weights:
        objective = "LOWER_CENTER_OF_MASS" if weight == 0 else (
            "MAXIMUM_MASS" if weight == 1 else "PARETO_TRADEOFF")
        problem = {**case, "static_objective": objective,
                   "static_objective_measure": measure,
                   "static_pareto_weight": weight if objective == "PARETO_TRADEOFF" else 0.}
        application = table.run_static_plane_case(problem, 1, 8, asset_root=asset_root,
                                                   initial_seeds=(0,))
        terminal = application["attempts"][-1]
        row = {"weight": weight, "status": terminal["status"], "application": application}
        if application["passed"]:
            metrics = moment_metrics(terminal["solver_result"]["total_raw_moments"])
            row.update(metrics)
            row["material_mass"] = metrics["volume_or_mass"]
            row["material_volume"] = row["material_mass"]/application["static_problem"]["material_density"]
            row["first_objective"] = -(row["material_volume"] if measure == "UNIFORM_DENSITY_VOLUME_V1" else row["material_mass"])
            row["static_objective_measure"] = measure
            row["interface_subtype"] = "PLANE"
        rows.append(row)
    qualified = [row for row in rows if row["application"]["passed"]]
    return {"case_id": case["id"], "mode": case["mode"], "solver_family": "STATIC_PLANE_FAMILY_V1",
            "static_objective_measure": measure,
            "warm_started": False, "initialization_protocol": "INDEPENDENT_INPUT_SEED_0_PER_WEIGHT_V1",
            "rows": rows, "front_audit": validate_pareto_front(qualified) if qualified else
                {"finite": False, "nondominated": False},
            "cost_ledger": {"reference_attempts": 0, "weight_solves": len(rows),
                "function_evaluations_all_solves": sum(a["function_evaluations"] for row in rows
                                                       for a in row["application"]["attempts"])},
            "formal_paper_admission": "NOT_GRANTED_BY_METHOD_RUNNER"}


def run_cherchi_top_pareto_sweep(case, weights, *, asset_root=ROOT, bundle_root=None,
                                cherchi_build=None, cherchi_probe=None,
                                initial_seeds=(0, 1), max_restarts=8):
    """Finite lambda method runs; this does not implement the pure-height endpoint.

    Every weight starts from request + seed, never another terminal quadric.
    The preflight plan and each finished case are persisted before advancing.
    Native exceptions are not swallowed or relabeled as failed convergence.
    """
    if case.get("mode") != "SPINNING_TOP":
        raise ValueError("strict finite-weight Pareto requires SPINNING_TOP")
    if case.get("spinning_top_objective","INERTIA_AND_LOWERING") != "INERTIA_AND_LOWERING":
        raise ValueError("finite-weight objective sweep cannot substitute the pure-height endpoint")
    if (not weights or any(type(w) not in (int, float) or not math.isfinite(w) or w < 0 for w in weights)
            or len(set(weights)) != len(weights)):
        raise ValueError("top Pareto weights must be distinct finite nonnegative lambdas")
    table = _table_module()
    seeds = table.registered_seeds(len(initial_seeds), initial_seeds)
    problems = [{**case, "id":f"{case['id']}-pareto-{i:06d}", "center_of_mass_weight":weight}
                for i, weight in enumerate(weights)]
    plans = [table.cherchi_case_plan(p, len(seeds), max_restarts, asset_root=asset_root,
                                    initial_seeds=seeds) for p in problems]
    if bundle_root is None or cherchi_build is None or cherchi_probe is None:
        raise ValueError("strict top Pareto needs bundle_root, cherchi_build and cherchi_probe")
    root = Path(bundle_root)
    root.mkdir(parents=True, exist_ok=False)
    payload = {"schema":"cherchi-top-finite-pareto-v1", "case_id":case["id"],
        "mode":"SPINNING_TOP", "solver_family":"CHERCHI_RAW_EQ14_V1",
        "weights":list(weights), "initialization_plans":plans,
        "warm_started":False, "initialization_protocol":"INDEPENDENT_INPUT_SEEDS_PER_WEIGHT_V1",
        "endpoint_status":{"pure_inertia":"REQUESTED" if 0 in weights else "NOT_REQUESTED",
                           "pure_lowering":"BLOCKED_MISSING_GEOMETRY_ENDPOINT_VALIDATION"},
        "formal_paper_admission":"NOT_GRANTED_BY_METHOD_RUNNER"}
    (root/"pareto-plan.json").write_text(json.dumps(payload, sort_keys=True, indent=2, allow_nan=False)+"\n")
    cases = []
    for i, problem in enumerate(problems):
        result = table.run_cherchi_case(problem, len(seeds), max_restarts,
            asset_root=asset_root, bundle_root=root, initial_seeds=seeds,
            cherchi_build=cherchi_build, cherchi_probe=cherchi_probe)
        (root/f"pareto-case-{i:06d}.json").write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+"\n")
        cases.append(result)
    payload["cases"] = cases
    payload["cost_ledger"] = {"reference_attempts":0, "weight_solves":len(cases),
        "top_level_runs":sum(len(c["runs"]) for c in cases),
        "function_evaluations_all_solves":sum(a["function_evaluations"] for c in cases for a in c["attempts"])}
    payload["scientific_hash"] = table.scientific_hash(payload)
    (root/"pareto-result.json").write_text(json.dumps(payload, sort_keys=True, indent=2, allow_nan=False)+"\n")
    return payload


def run_legacy_quadric_pareto_sweep(case: dict[str, Any], weights: list[float]) -> dict[str, Any]:
    """Explicit historical path, retained without relabeling its old results."""
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

    def solve_record(result, initial):
        record = {
            'run':0, 'restart':1, 'status':result.status,
            'initial_quadric':list(initial),
            'function_evaluations':result.solve_function_evaluations,
            'raw_equation_residual':result.raw_equation_residual,
            'normalized_equation_residual':result.diagnostics['normalized_equation_residual'],
            'physical_relative_residual':result.diagnostics['physical_relative_residual'],
            'equation_profile':result.equation_profile,
            'multiplier_mu':result.trace[-1].multiplier_mu if result.trace else None,
            'trace':[{'iteration':entry.iteration, 'restart':entry.restart,
                      'raw_equation_residual':entry.raw_equation_residual,
                      'multiplier_mu':entry.multiplier_mu,
                      'equation_profile':entry.equation_profile}
                     for entry in result.trace],
            'solver_events':[table.solver_event_payload(event) for event in result.events],
        }
        # Do not emit invalid JSON or silently drop failed nonfinite attempts.
        # Explicit tagged values remain in the scientific ledger and cannot
        # satisfy numerical acceptance (which requires actual finite numbers).
        def tagged(value):
            if isinstance(value, dict):
                return {key:tagged(item) for key,item in value.items()}
            if isinstance(value, list):
                return [tagged(item) for item in value]
            if type(value) is float and not math.isfinite(value):
                return {'nonfinite':('NaN' if math.isnan(value) else
                        'POSITIVE_INFINITY' if value>0 else 'NEGATIVE_INFINITY')}
            return value
        record = tagged(record)
        table.validate_attempt_ledger([record], expected_runs=1, maximum_restarts=1)
        record['numerical_acceptance'] = table.numerical_attempt_acceptance(record)
        return record

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
    reference_attempts = []
    for candidate_index, candidate in enumerate(candidates):
        attempted = solve_spin_it_faster(
            request_for(reference_weight),
            initial_quadrics=(candidate,),
        )
        reference_attempts.append({'candidate_index':candidate_index,
                                   'solve_record':solve_record(attempted, candidate)})
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
            "reference_attempts": reference_attempts,
            "cost_ledger": {'reference_attempts':len(reference_attempts), 'weight_solves':0,
                'function_evaluations_all_solves':sum(row['solve_record']['function_evaluations']
                    for row in reference_attempts)},
            "formal_paper_admission": "NOT_GRANTED_LEGACY_RECONSTRUCTED_DISPATCH",
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
            row: dict[str, Any] = {"weight": weight, "status": result.status,
                                  "solve_record":solve_record(result, initial)}
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
        "reference_attempts": reference_attempts,
        "cost_ledger": {'reference_attempts':len(reference_attempts), 'weight_solves':len(rows),
            'function_evaluations_all_solves':sum(row['solve_record']['function_evaluations']
                for row in reference_attempts)+sum(row['solve_record']['function_evaluations'] for row in rows)},
        "formal_paper_admission": "NOT_GRANTED_LEGACY_RECONSTRUCTED_DISPATCH",
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
    # Existing box capacity/symmetry checkers do not support this rotated TOP
    # family. Running Newton cannot supply a continuous infeasibility certificate.
    if not math.isfinite(tolerance) or tolerance <= 0 or not math.isfinite(upper) or upper <= 0:
        raise ValueError("Figure 13 tolerance and upper distance must be finite and positive")
    return {
        "asset_class": "RECONSTRUCTED_ASSET",
        "status": "BLOCKED",
        "feasibility_status": "UNKNOWN",
        "blocker": "MISSING_CONTINUOUS_FEASIBILITY_ORACLE_FOR_ROTATED_TOP_FAMILY",
        "critical_distance": None,
        "critical_interval": None,
        "binary_search_tolerance": tolerance,
        "requested_upper_distance": upper,
        "unperturbed_status": "NOT_RUN",
        "rows": [],
        "historical_solver_status_search": "NOT_VALID_AS_FEASIBILITY_PROOF",
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
            "residual_kind": "RAW_EQ14",
            "convergence_order_evidence": convergence_order_evidence(table or {"cases": []}),
            "threshold": 1e-12,
        },
        "figure11": {
            "definition": "min -kappa*v+(1-kappa)*c3; l1=l2=0",
            "weights": pareto_weights(41),
            "status": "BLOCKED_MISSING_ORIGINAL_FIGURE_INPUT",
        },
        "figure12": {
            "definition": "min f_yo+w_c*c3; l1=l2=I13=I23=0",
            "weights": pareto_weights(41),
            "status": "BLOCKED_MISSING_ORIGINAL_FIGURE_INPUT",
        },
        "figure13": {
            "definition": "d_m=d_crit-10^-m, m=1..9",
            "binary_search_tolerance": 1e-9,
            "status": "BLOCKED_MISSING_CONTINUOUS_FEASIBILITY_ORACLE",
        },
        "section6_6": degenerate_one_dimensional_examples(),
    }


def formal_pareto_preflight(readiness_path: Path | None = None) -> dict[str, Any]:
    """Report original-figure input gaps, never substitute Table 1 method assets.

    Reuse the strict audit's exhaustive, source-bound readiness validation.
    A ready inventory row is not an authenticated locked solver request; that
    binding remains a separate implementation gate, not an inferred pass.
    """
    strict = ROOT.parent / "strict"
    path = Path(readiness_path) if readiness_path is not None else (
        ROOT.parent / "sources/case_asset_readiness.json")
    spec = importlib.util.spec_from_file_location("faster_figure_input_audit", strict / "audit.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load existing strict input audit")
    added = str(strict) not in sys.path
    if added:
        sys.path.insert(0, str(strict))
    validation = {"status": "VALIDATED_INVENTORY"}
    readiness = None
    source_binding = None
    try:
        audit = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(audit)
        try:
            record = audit.load_case_asset_readiness_record(path)
            readiness = record["cases"]
            source_binding = record["source_audit"]
        except (ValueError, OSError) as error:
            validation = {"status": "REJECTED", "exception_type": type(error).__name__,
                          "diagnostic": str(error)}
    finally:
        if added:
            sys.path.remove(str(strict))
    digests = {}
    read_errors = []
    for name, source in (("readiness_sha256", path),):
        try:
            digests[name] = hashlib.sha256(source.read_bytes()).hexdigest()
        except OSError as error:
            digests[name] = None
            read_errors.append({"path": str(source), "exception_type": type(error).__name__,
                                "diagnostic": str(error)})
            readiness = None
            source_binding = None
            if validation["status"] != "REJECTED":
                validation = {"status": "REJECTED", "exception_type": type(error).__name__,
                              "diagnostic": str(error)}
    figures = []
    for requirement in ("faster.figure.11", "faster.figure.12"):
        if readiness is None:
            figures.append({"requirement_id": requirement, "asset_status": "UNVERIFIED",
                            "blocker": "INVALID_CASE_ASSET_READINESS"})
            continue
        row = dict(readiness[requirement])
        missing = bool(row["missing_artifacts"] or row["missing_critical_parameters"])
        row["blocker"] = ("MISSING_ORIGINAL_FIGURE_INPUT" if missing else
                          "MISSING_AUTHENTICATED_LOCKED_FIGURE_REQUEST_BINDING")
        figures.append(row)
    return {"schema": "faster-formal-pareto-input-preflight-v1",
            "status": "BLOCKED", "solver_executed": False,
            "readiness_validation": validation,
            **digests, "input_read_errors": read_errors,
            "source_audit_binding": source_binding,
            "source_audit_sha256": source_binding["sha256"] if source_binding else None,
            "source_authenticity_verified": False, "figures": figures}


def _journal_application_binding(path: Path, run_config: Path) -> dict[str, Any]:
    """Bind the reporting group, not all run inputs, to the pinned request.

    Source-job receipts and the final material still require their own audit.
    In particular, this does not authenticate the journal's inputs_sha256.
    """
    config_bytes = Path(run_config).read_bytes()
    config = json.loads(config_bytes)
    if not isinstance(config, dict):
        raise ValueError("invalid run config")
    argv, identities = config.get('argv'), config.get('identity')
    if (not isinstance(argv, list) or not all(isinstance(arg, str) for arg in argv)
            or not isinstance(identities, dict)):
        raise ValueError("run config lacks argv or pinned identity")

    def argument_path(option):
        # The solver's argparse permits long-option abbreviations, including
        # later overrides. Require canonical arguments rather than guessing
        # which prefix was unique in the solver revision used by this run.
        matches = [i for i, arg in enumerate(argv)
                   if arg.startswith('--') and len(arg.split('=', 1)[0]) > 2
                   and option.startswith(arg.split('=', 1)[0])]
        if len(matches) != 1:
            raise ValueError(f"run config requires one unambiguous {option}")
        index = matches[0]
        if argv[index] != option or index+1 >= len(argv):
            raise ValueError(f"run config requires an explicit path after {option}")
        value = Path(argv[index+1])
        if not value.is_absolute():
            raise ValueError(f"run config requires an absolute {option} path")
        return value

    manifest = argument_path('--problem-manifest')
    output = argument_path('--output')
    if Path(path).resolve() != (output/'solver-events.jsonl').resolve():
        raise ValueError("solver journal does not belong to run config output")
    identity = identities.get(str(manifest))
    if (not isinstance(identity, dict) or identity.get('path') != str(manifest)):
        raise ValueError("problem manifest has no matching pinned identity")
    manifest_bytes = manifest.read_bytes()
    manifest_hash = hashlib.sha256(manifest_bytes).hexdigest()
    if (identity.get('sha256') != manifest_hash or type(identity.get('size')) is not int
            or identity['size'] != len(manifest_bytes)):
        raise ValueError("problem manifest identity has changed")
    problem = json.loads(manifest_bytes)
    if (not isinstance(problem, dict)
            or problem.get('schema') != 'faster-cherchi-problem-manifest-v1'
            or not isinstance(problem.get('target'), dict)
            or problem['target'].get('mode') not in APPLICATION_GROUP):
        raise ValueError("problem manifest lacks a supported target mode")
    return {
        'schema': 'faster-figure-journal-application-binding-v1',
        'scope': 'APPLICATION_GROUP_ONLY_NOT_FULL_RUN_ATTESTATION',
        'mode': problem['target']['mode'],
        'manifest_sha256': manifest_hash,
        'run_config_sha256': hashlib.sha256(config_bytes).hexdigest(),
    }


def build_payload_from_solver_journal(
    path: Path, *, case_id: str, mode: str | None = None, run_config: Path | None = None,
) -> dict[str, Any]:
    """Stream actual native events into the existing Figure 10 checks.

    A live prefix is explicitly incomplete. Journal identity and numerical
    qualification do not attest the input assets or final material geometry.
    """
    binding = None
    if run_config is not None:
        binding = _journal_application_binding(path, run_config)
        if mode is not None and mode != binding['mode']:
            raise ValueError("requested Figure 10 mode disagrees with the pinned request mode")
        mode = binding['mode']
    if not case_id or mode not in APPLICATION_GROUP:
        raise ValueError("a case identity and supported Figure 10 application are required")
    digest = hashlib.sha256()
    attempts: dict[int, dict[str, Any]] = {}
    counts: dict[str, int] = {}
    previous_restart = -1
    evaluation_count = 0
    header = None
    with Path(path).open('rb') as stream:
        for line_number, line in enumerate(stream, 1):
            digest.update(line)
            if not line.endswith(b'\n'):
                raise ValueError("truncated solver journal line")
            row = json.loads(line)
            if line_number == 1:
                if not isinstance(row, dict) or row.get('schema') not in {
                    'faster-streamed-solver-events-v1', 'faster-streamed-solver-events-v2',
                }:
                    raise ValueError("unsupported solver journal header")
                identity = row.get('inputs_sha256')
                if (not isinstance(identity, str) or len(identity) != 64 or
                        any(char not in '0123456789abcdef' for char in identity)):
                    raise ValueError("solver journal input identity is missing")
                header = row
                continue
            if not isinstance(row, dict) or not isinstance(row.get('record'), dict):
                raise ValueError("malformed solver journal record")
            if row.get('record_type') in ('geometry_evaluation_started', 'geometry_evaluation_finished'):
                continue
            if row.get('record_type') != 'solver_event':
                raise ValueError("unsupported solver journal record type")
            record = row['record']
            if record.get('type') not in ('evaluation', 'attempt_terminal'):
                continue
            restart = record.get('restart')
            if type(restart) is not int or restart < 0 or restart < previous_restart:
                raise ValueError("solver journal restart identity regressed")
            if restart > previous_restart and restart != previous_restart + 1:
                raise ValueError("solver journal has a missing restart attempt")
            previous_restart = restart
            attempt = attempts.setdefault(restart, {
                'restart': restart, 'status': 'INCOMPLETE', 'trace': [],
                'evaluations': [], 'function_evaluations': 0,
            })
            if attempt['status'] != 'INCOMPLETE':
                raise ValueError("solver journal continues an already terminated attempt")
            index = record.get('evaluation_index')
            if type(index) is not int:
                raise ValueError("solver journal evaluation identity is missing")
            if record['type'] == 'evaluation':
                if index != evaluation_count + 1:
                    raise ValueError("solver journal has a dropped or repeated evaluation")
                evaluation_count += 1
                outcome = record.get('outcome')
                if outcome not in ('accepted', 'rejected', 'failed'):
                    raise ValueError("unsupported solver evaluation outcome")
                counts[outcome] = counts.get(outcome, 0) + 1
                attempt['function_evaluations'] += 1
                attempt['evaluations'].append({key: record.get(key) for key in
                    ('evaluation_index', 'restart', 'iteration', 'outcome', 'reason')})
                if outcome == 'accepted':
                    # Do not expand large geometry/Jacobian payloads into memory.
                    accepted = {key: record.get(key) for key in (
                        'evaluation_index', 'restart', 'iteration', 'equation_profile',
                        'raw_equation_residual', 'normalized_equation_residual',
                        'multiplier_mu', 'physical_relative_residual', 'material_topology_hash')}
                    _raw_residuals({'trace': [accepted]})
                    attempt['trace'].append(accepted)
            else:
                if index != evaluation_count or record.get('outcome') != 'terminated':
                    raise ValueError("terminal event is not bound to the actual evaluation count")
                attempt['native_terminal_reason'] = record.get('reason')
                attempt['status'] = 'CONVERGED' if record.get('reason') == 'converged' else 'FAILED'
    if header is None:
        raise ValueError("empty solver journal")
    checker = _table_module().numerical_attempt_acceptance
    rows = list(attempts.values())
    for attempt in rows:
        if attempt['trace']:
            for key in ('equation_profile', 'raw_equation_residual',
                        'physical_relative_residual', 'multiplier_mu'):
                attempt[key] = attempt['trace'][-1][key]
        attempt['native_status'] = attempt['status']
        attempt['numerical_acceptance'] = checker(attempt)
        if attempt['status'] == 'CONVERGED' and not attempt['numerical_acceptance']['passed']:
            # Preserve the native terminal record but do not plot an invalid
            # stationary point as a qualified successful Eq.14 trajectory.
            attempt['status'] = 'NUMERICALLY_UNQUALIFIED'
    payload = build_payload({'cases': [{'id': case_id, 'mode': mode, 'attempts': rows}]})
    payload['figure10']['solver_journal'] = {
        'scope': 'JOURNAL_NUMERICS_NOT_MATERIAL_OR_ORIGINAL_PAPER_ACCEPTANCE',
        'sha256': digest.hexdigest(), 'inputs_sha256': header['inputs_sha256'],
        'figure_checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'numerical_checker_sha256': hashlib.sha256(TABLE_RUNNER.read_bytes()).hexdigest(),
        'input_identity_attestation': 'NOT_VERIFIED_BY_FIGURE_CHECKER',
        'evaluations': evaluation_count, 'evaluation_outcomes': counts,
        'failed_evaluations': counts.get('failed', 0),
        'terminal_attempts': sum(row['status'] != 'INCOMPLETE' for row in rows),
        'attempts': rows,
    }
    if binding is not None:
        payload['figure10']['solver_journal']['application_binding'] = binding
    return payload


def finalize_payload(payload: dict[str, Any]) -> dict[str, Any]:
    result = dict(payload)
    result.pop("scientific_hash", None)
    result["scientific_hash"] = _table_module().scientific_hash(result)
    return result


def attach_top_pareto_method_result(payload, result, bundle_root):
    """Do not inherit original-figure scope/weights for a finite method run.

    Feed <the explicit --pareto-bundle-root>/pareto-result.json, not this
    multi-figure report, to check_saved_cherchi_application.py --application.
    The display-only destination must not change the scientific hash.
    """
    payload["figure12"].update({"result":result,"weights":list(result["weights"]),
        "status":"PARTIAL_METHOD_ONLY", "scope":"FINITE_TOP_LAMBDAS_NOT_COMPLETE_FIGURE12",
        "endpoint_status":dict(result["endpoint_status"]),
        "checker_input":{"relative_path":"pareto-result.json","base":"--pareto-bundle-root"}})


def main() -> int:
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--table1", type=Path)
    source.add_argument("--solver-journal", type=Path)
    parser.add_argument("--journal-case-id")
    parser.add_argument("--journal-mode", choices=tuple(APPLICATION_GROUP))
    parser.add_argument("--run-config", type=Path,
                        help="source-job config: derive and verify the journal application group")
    parser.add_argument("--output", type=Path, required=True)
    pareto = parser.add_mutually_exclusive_group()
    pareto.add_argument("--run-pareto", action="store_true",
                        help="original Figures 11/12 only; blocked until their own inputs are admitted")
    parser.add_argument("--pareto-count", type=int, default=41)
    pareto.add_argument("--pareto-top-case", type=Path,
                        help="locked method case JSON; finite top lambdas, not a complete Figure 12")
    pareto.add_argument("--pareto-static-case", type=Path,
                        help="explicit static method case, never a substitute for Figure 11 eggs")
    parser.add_argument("--pareto-weights", nargs="+", type=float)
    parser.add_argument("--pareto-asset-root", type=Path)
    parser.add_argument("--pareto-bundle-root", type=Path)
    parser.add_argument("--cherchi-build", type=Path)
    parser.add_argument("--cherchi-probe", type=Path)
    parser.add_argument("--initial-seeds", nargs="+", type=int, default=[0,1])
    parser.add_argument("--max-restarts", type=int, default=8)
    parser.add_argument("--run-figure13", action="store_true")
    parser.add_argument("--figure13-tolerance", type=float, default=1e-9)
    args = parser.parse_args()
    if args.run_config and not args.solver_journal:
        parser.error("--run-config requires --solver-journal")
    if args.pareto_static_case and args.pareto_weights is None:
        parser.error("static method case requires explicit --pareto-weights")
    if args.pareto_top_case and any(x is None for x in
            (args.pareto_weights,args.pareto_bundle_root,args.cherchi_build,args.cherchi_probe)):
        parser.error("locked top sweep needs explicit weights, bundle root, build and probe")
    table = json.loads(args.table1.read_text(encoding="utf-8")) if args.table1 else None
    if args.solver_journal:
        if not args.journal_case_id or not (args.journal_mode or args.run_config):
            parser.error("solver journal requires --journal-case-id and --journal-mode or --run-config")
        payload = build_payload_from_solver_journal(
            args.solver_journal, case_id=args.journal_case_id, mode=args.journal_mode,
            run_config=args.run_config)
    else:
        payload = build_payload(table)
    if args.pareto_static_case:
        case = json.loads(args.pareto_static_case.read_text())
        if case.get("mode") != "STATIC_BALANCE":
            parser.error("static method case requires STATIC_BALANCE")
        payload["static_pareto_method"] = run_pareto_sweep(
            case, args.pareto_weights,
            asset_root=args.pareto_asset_root or args.pareto_static_case.resolve().parent)
    if args.pareto_top_case:
        case=json.loads(args.pareto_top_case.read_text())
        result=run_pareto_sweep(case,args.pareto_weights,
            asset_root=args.pareto_asset_root or args.pareto_top_case.resolve().parent,
            bundle_root=args.pareto_bundle_root,cherchi_build=args.cherchi_build,
            cherchi_probe=args.cherchi_probe,initial_seeds=tuple(args.initial_seeds),
            max_restarts=args.max_restarts)
        attach_top_pareto_method_result(payload,result,args.pareto_bundle_root)
    if args.run_pareto:
        payload["formal_pareto_preflight"] = formal_pareto_preflight()
        for row in payload["formal_pareto_preflight"]["figures"]:
            key = "figure" + row["requirement_id"].rsplit(".", 1)[-1]
            payload[key]["status"] = "BLOCKED_" + row["blocker"]
    if args.run_figure13:
        payload["figure13"]["result"] = run_tilted_cube_feasibility(
            tolerance=args.figure13_tolerance
        )
    payload = finalize_payload(payload)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 2 if args.run_pareto else 0


if __name__ == "__main__":
    raise SystemExit(main())
