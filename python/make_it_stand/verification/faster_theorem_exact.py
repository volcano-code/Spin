"""Exact rational instance checks for Spin-It Faster Theorems 1 and 2.

This checker validates exact applications of the published conditions.  It
does not prove the general infinite-dimensional theorems or authenticate an
author experiment.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Any


_SCHEMA = "faster-exact-theorem-instance-certificate-v1"
_MONOMIALS = (
    (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
    (2, 0, 0), (0, 2, 0), (0, 0, 2), (0, 1, 1), (1, 0, 1), (1, 1, 0),
)


def _q(value: Any, label: str) -> Fraction:
    if not isinstance(value, str):
        raise ValueError(f"{label} must be encoded as a rational string")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError(f"{label} is not a finite rational") from error


def _vector(value: Any, size: int, label: str) -> list[Fraction]:
    if not isinstance(value, list) or len(value) != size:
        raise ValueError(f"{label} must have exactly {size} entries")
    return [_q(item, f"{label}[{index}]") for index, item in enumerate(value)]


def _box(value: Any, label: str) -> tuple[list[Fraction], list[Fraction]]:
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a box object")
    lo = _vector(value.get("minimum"), 3, f"{label}.minimum")
    hi = _vector(value.get("maximum"), 3, f"{label}.maximum")
    if any(hi[axis] <= lo[axis] for axis in range(3)):
        raise ValueError(f"{label} must have positive extent on every axis")
    return lo, hi


def _integral(lo: list[Fraction], hi: list[Fraction], powers: tuple[int, int, int]) -> Fraction:
    result = Fraction(1)
    for axis, power in enumerate(powers):
        result *= (hi[axis] ** (power + 1) - lo[axis] ** (power + 1)) / (power + 1)
    return result


def _moments(lo: list[Fraction], hi: list[Fraction], density: Fraction = Fraction(1)) -> list[Fraction]:
    return [density * _integral(lo, hi, powers) for powers in _MONOMIALS]


def _fraction_vector(values: list[Fraction]) -> list[str]:
    return [str(value) for value in values]


def reference_certificate() -> dict[str, Any]:
    """Return two hand-specified exact instances covering Theorem 2 branches."""
    return {
        "schema": _SCHEMA,
        "instances": [
            {
                "case_id": "affine_cut_with_equality_constraint",
                "kind": "P_NONZERO_AFFINE_CUT",
                "domain_box": {"minimum": ["0", "0", "0"], "maximum": ["1", "1", "1"]},
                "beta_coefficients": ["3/5", "0", "0", "-4/5", "0", "0", "0", "0", "0", "0"],
                "mu": "1",
                "material_box": {"minimum": ["0", "0", "0"], "maximum": ["1", "1", "3/4"]},
                "raw_moments": ["3/4", "3/8", "3/8", "9/32", "1/4", "1/4", "9/64", "9/64", "9/64", "3/16"],
                "objective_kind": "LINEAR_RAW_MOMENTS",
                "objective_gradient": ["-3/5", "1", "0", "4/5", "0", "0", "0", "0", "0", "0"],
                "constraints": [
                    {
                        "name": "first_x_moment",
                        "gradient": ["0", "1", "0", "0", "0", "0", "0", "0", "0", "0"],
                        "target": "3/8",
                        "multiplier": "1",
                    }
                ],
            },
            {
                "case_id": "flat_objective_p_zero_branch",
                "kind": "P_ZERO_FLAT_OBJECTIVE",
                "domain_box": {"minimum": ["0", "0", "0"], "maximum": ["1", "1", "1"]},
                "characteristic_value": "1/2",
                "objective_kind": "LINEAR_RAW_MOMENTS",
                "objective_gradient": ["0"] * 10,
                "constraints": [],
                "raw_moments": ["1/2", "1/4", "1/4", "1/4", "1/6", "1/6", "1/6", "1/8", "1/8", "1/8"],
            },
        ],
    }


def _audit_nonzero_branch(instance: dict[str, Any]) -> dict[str, Any]:
    lo, hi = _box(instance.get("domain_box"), "domain_box")
    p = _vector(instance.get("beta_coefficients"), 10, "beta_coefficients")
    if sum(value * value for value in p) != 1:
        raise ValueError("Theorem 1 coefficient vector is not exactly unit norm")
    if any(p[index] != 0 for index in range(4, 10)):
        raise ValueError("exact affine-cut checker does not accept quadratic coefficients")
    active_axes = [axis for axis in range(3) if p[axis + 1] != 0]
    if len(active_axes) != 1:
        raise ValueError("exact affine-cut checker requires one nonzero linear coefficient")
    axis = active_axes[0]
    slope = p[axis + 1]
    root = -p[0] / slope
    if not lo[axis] < root < hi[axis]:
        raise ValueError("the affine zero plane must cross the interior of the domain")
    mu = _q(instance.get("mu"), "mu")
    if mu <= 0:
        raise ValueError("Theorem 1 multiplier mu must be exactly positive")

    material_lo, material_hi = lo.copy(), hi.copy()
    if slope > 0:
        material_lo[axis] = root
    else:
        material_hi[axis] = root
    submitted_material_lo, submitted_material_hi = _box(instance.get("material_box"), "material_box")
    if (submitted_material_lo, submitted_material_hi) != (material_lo, material_hi):
        raise ValueError("material box does not equal the exact beta-positive region")

    exact_moments = _moments(material_lo, material_hi)
    submitted_moments = _vector(instance.get("raw_moments"), 10, "raw_moments")
    if submitted_moments != exact_moments:
        raise ValueError("raw moments do not match exact box integration")

    if instance.get("objective_kind") != "LINEAR_RAW_MOMENTS":
        raise ValueError("exact instance must declare a linear raw-moment objective")
    objective_gradient = _vector(instance.get("objective_gradient"), 10, "objective_gradient")
    constraints = instance.get("constraints")
    if not isinstance(constraints, list):
        raise ValueError("constraints must be a list")
    total_constraint_gradient = [Fraction(0) for _ in range(10)]
    constraint_values = []
    for row in constraints:
        if not isinstance(row, dict):
            raise ValueError("constraint row must be an object")
        gradient = _vector(row.get("gradient"), 10, "constraint.gradient")
        target = _q(row.get("target"), "constraint.target")
        multiplier = _q(row.get("multiplier"), "constraint.multiplier")
        value = sum(a * b for a, b in zip(gradient, exact_moments)) - target
        if value != 0:
            raise ValueError("exact equality constraint is not satisfied")
        constraint_values.append((row.get("name"), value))
        for index in range(10):
            total_constraint_gradient[index] += multiplier * gradient[index]

    stationarity = [objective_gradient[index] + mu * p[index] - total_constraint_gradient[index]
                    for index in range(10)]
    if any(stationarity):
        raise ValueError("Eq. 7 stationarity is not exact")

    # A nonconstant affine polynomial has a plane as its zero set.  The plane
    # has 3-D Lebesgue measure zero; beta has strict sign on either side.
    # Combined with 0 <= chi <= 1, the exact pointwise sign argument yields a
    # strict objective increase for every distinct feasible characteristic
    # function, modulo equality almost everywhere.
    return {
        "case_id": instance.get("case_id"),
        "raw_moments": _fraction_vector(exact_moments),
        "p_squared_norm": "1",
        "mu": str(mu),
        "constraint_residuals": {name: str(value) for name, value in constraint_values},
        "eq7_stationarity_residual": _fraction_vector(stationarity),
        "beta_zero_set": "INTERIOR_PLANE",
        "zero_set_measure": "ZERO",
        "theorem_1": "EXACT_HYPOTHESES_AND_STRICT_SIGN_CONCLUSION",
        "theorem_2_branch": "P_NONZERO_STRICT_LOCAL_OPTIMUM",
        "exact_global_minimum_for_this_linear_instance": True,
        "scope": "EXACT_RATIONAL_AFFINE_BOX_INSTANCE_ONLY",
    }


def _audit_zero_branch(instance: dict[str, Any]) -> dict[str, Any]:
    lo, hi = _box(instance.get("domain_box"), "domain_box")
    if instance.get("objective_kind") != "LINEAR_RAW_MOMENTS":
        raise ValueError("p=0 branch instance must have a linear raw-moment objective")
    p = _vector(instance.get("beta_coefficients", ["0"] * 10), 10, "beta_coefficients")
    gradient = _vector(instance.get("objective_gradient"), 10, "objective_gradient")
    if any(p) or any(gradient) or instance.get("constraints") != []:
        raise ValueError("p=0 branch requires exactly zero objective gradient and no equalities")
    characteristic = _q(instance.get("characteristic_value"), "characteristic_value")
    if not 0 <= characteristic <= 1:
        raise ValueError("p=0 example characteristic function must lie in [0,1]")
    moments = [characteristic * value for value in _moments(lo, hi)]
    submitted = _vector(instance.get("raw_moments"), 10, "raw_moments")
    if submitted != moments:
        raise ValueError("p=0 branch moments do not match exact constant characteristic")
    return {
        "case_id": instance.get("case_id"),
        "raw_moments": _fraction_vector(moments),
        "p_squared_norm": "0",
        "relaxed_stationarity_residual": ["0"] * 10,
        "theorem_1": "NOT_APPLICABLE_P_ZERO",
        "theorem_2_branch": "P_ZERO_UNBOUNDED_FIRST_ORDER_STATIONARITY",
        "flat_objective_makes_this_instance_a_global_optimizer": True,
        "scope": "EXACT_RATIONAL_FLAT_OBJECTIVE_INSTANCE_ONLY",
    }


def verify_certificate(certificate: Any) -> dict[str, Any]:
    if not isinstance(certificate, dict) or certificate.get("schema") != _SCHEMA:
        raise ValueError("unsupported exact theorem-instance certificate schema")
    instances = certificate.get("instances")
    if not isinstance(instances, list) or len(instances) != 2 or any(not isinstance(row, dict) for row in instances):
        raise ValueError("exact certificate must contain both declared theorem instances")
    by_id = {row.get("case_id"): row for row in instances}
    if len(by_id) != len(instances) or set(by_id) != {
        "affine_cut_with_equality_constraint", "flat_objective_p_zero_branch",
    }:
        raise ValueError("exact theorem instance IDs are incomplete or duplicated")
    if by_id["affine_cut_with_equality_constraint"].get("kind") != "P_NONZERO_AFFINE_CUT":
        raise ValueError("nonzero-p theorem instance kind mismatch")
    if by_id["flat_objective_p_zero_branch"].get("kind") != "P_ZERO_FLAT_OBJECTIVE":
        raise ValueError("p=0 theorem instance kind mismatch")
    cases = [
        _audit_nonzero_branch(by_id["affine_cut_with_equality_constraint"]),
        _audit_zero_branch(by_id["flat_objective_p_zero_branch"]),
    ]
    return {
        "schema": "faster-exact-theorem-instance-result-v1",
        "status": "EXACT_INSTANCE_CONDITIONS_VERIFIED",
        "general_theorem_proof": "NOT_CLAIMED",
        "author_case_reproduction": "NOT_CLAIMED",
        "cases": cases,
    }
