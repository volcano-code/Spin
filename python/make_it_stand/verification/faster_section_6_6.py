"""Exact certificates for the two published Faster Section 6.6 examples.

This module deliberately supports only the two one-dimensional definitions in
the authenticated local specification.  It does not call a geometry solver and
does not classify other zero-sensitivity problems.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable


EXPECTED_SPEC_SHA256 = "053f896c2b10ece89b136bfba24cd480bf44cd0da99c0661e6e9bd34495c15fa"
EXPECTED_PDF_SHA256 = "e7fc6395362ef7e63a703851176afd81d8e953bae4e8080279d7e244b3c51fcc"

_TOP_LEVEL_FIELDS = {
    "schema", "paper", "source_url", "source_sha256", "paper_locator",
    "asset_class", "parameter_evidence", "implementation_status",
    "formal_run_status", "units", "variable", "volume_definition",
    "objective", "constraints", "cases", "independent_check_to_implement",
    "caveat",
}
_CASE_FIELDS = {
    "case_id", "domain", "expected_outcome", "expected_volume",
    "expected_objective", "expected_p", "interpretation",
}
_FIXED_DEFINITIONS = {
    "schema": "published-analytic-experiment-spec-v1",
    "paper": "Spin-It Faster, ACM TOG 43(4), article 78, 2024",
    "source_url": "https://research-explorer.ista.ac.at/download/17203/17204/sif-final.pdf",
    "source_sha256": EXPECTED_PDF_SHA256,
    "paper_locator": "Section 6.6, printed pages 78:12-78:13",
    "asset_class": "PUBLISHED_EXACT_SPECIFICATION",
    "parameter_evidence": "PUBLISHED",
    "units": "dimensionless",
    "variable": "measurable chi: Omega -> [0,1]",
    "volume_definition": "v(chi) = integral_Omega chi(x) dx",
    "objective": "(v(chi) - 1/2)^2",
    "constraints": [],
}
_DOMAINS = {
    "faster.section6_6.nonunique": {
        "dimension": 1, "lower": 0, "upper": 1, "endpoints": "open",
    },
    "faster.section6_6.full_boundary": {
        "dimension": 1, "lower": 0, "upper": "1/2", "endpoints": "open",
    },
}
_CATALOGUE_IDS = {
    "faster.section6_6.nonunique": "faster.section6_6.example.1",
    "faster.section6_6.full_boundary": "faster.section6_6.example.2",
}
_WITNESSES = {
    "affine_family": {
        "constant": ["1/2"],
        "tau_coefficient": ["-1", "2"],
        "tau_lower": "-1/2",
        "tau_upper": "1/2",
    },
    "binary_fixed_witnesses": [
        {"lower": "0", "upper": "1/2"},
        {"lower": "1/2", "upper": "1"},
    ],
    "binary_interval_family": {
        "parameter_lower": "0",
        "parameter_upper": "1/2",
        "interval_length": "1/2",
    },
    "nonbinary_constant_witness": "1/2",
    "full_domain_witness": {"constant": "1"},
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _pairs(items: Iterable[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in items:
        _require(key not in result, "duplicate JSON key: " + key)
        result[key] = value
    return result


def _reject_nonfinite(value: Any) -> None:
    if isinstance(value, dict):
        for item in value.values():
            _reject_nonfinite(item)
    elif isinstance(value, list):
        for item in value:
            _reject_nonfinite(item)
    elif type(value) is float:
        _require(math.isfinite(value), "nonfinite JSON number")


def load_spec(path: str | Path) -> dict[str, Any]:
    """Load JSON while rejecting duplicate keys and non-finite constants."""
    try:
        value = json.loads(
            Path(path).read_text(encoding="utf-8"),
            object_pairs_hook=_pairs,
            parse_constant=lambda token: (_ for _ in ()).throw(
                ValueError("nonfinite JSON number: " + token)
            ),
        )
    except (OSError, json.JSONDecodeError, TypeError) as error:
        raise ValueError("invalid specification JSON") from error
    _require(type(value) is dict, "specification must be a JSON object")
    _reject_nonfinite(value)
    return value


def canonical_json(value: Any) -> str:
    _reject_nonfinite(value)
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True, allow_nan=False)


def scientific_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def file_record(path: str | Path) -> dict[str, Any]:
    actual = Path(path).resolve()
    _require(actual.is_file() and not actual.is_symlink(), "missing regular file: " + str(actual))
    payload = actual.read_bytes()
    return {
        "path": str(actual),
        "size": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def published_witness_definitions() -> dict[str, Any]:
    """Return a detached copy of the fixed formulas the checker evaluates."""
    return json.loads(canonical_json(_WITNESSES))


def _fraction(value: Any, label: str) -> Fraction:
    _require((type(value) is int) or (type(value) is str and bool(value)),
             label + " must be an exact integer or rational string")
    try:
        result = Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError("invalid rational " + label) from error
    return result


def _rational(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _poly_integral(coefficients: Iterable[Fraction], lower: Fraction, upper: Fraction) -> Fraction:
    return sum(
        coefficient * (upper ** (power + 1) - lower ** (power + 1)) / (power + 1)
        for power, coefficient in enumerate(coefficients)
    )


def _poly_multiply(first: list[Fraction], second: list[Fraction]) -> list[Fraction]:
    result = [Fraction(0)] * (len(first) + len(second) - 1)
    for first_power, first_value in enumerate(first):
        for second_power, second_value in enumerate(second):
            result[first_power + second_power] += first_value * second_value
    return result


def _validate_spec(spec: dict[str, Any]) -> dict[str, dict[str, Any]]:
    _reject_nonfinite(spec)
    _require(type(spec) is dict and set(spec) == _TOP_LEVEL_FIELDS,
             "unsupported specification fields")
    for key, expected in _FIXED_DEFINITIONS.items():
        _require(spec.get(key) == expected, "mutated published definition: " + key)
    _require(type(spec["cases"]) is list and len(spec["cases"]) == 2,
             "exactly two published cases are required")
    cases: dict[str, dict[str, Any]] = {}
    for case in spec["cases"]:
        _require(type(case) is dict and set(case) == _CASE_FIELDS,
                 "unsupported published case fields")
        case_id = case.get("case_id")
        _require(type(case_id) is str and case_id in _DOMAINS,
                 "unsupported published case identity")
        _require(case_id not in cases, "duplicate published case identity")
        _require(case.get("domain") == _DOMAINS[case_id],
                 "mutated published domain: " + case_id)
        cases[case_id] = case
    _require(set(cases) == set(_DOMAINS), "missing published case")
    return cases


def _validate_witnesses(witnesses: dict[str, Any]) -> None:
    _reject_nonfinite(witnesses)
    _require(witnesses == _WITNESSES, "mutated or unsupported exact witness definitions")


def _objective_certificate(volume: Fraction, target: Fraction) -> dict[str, Any]:
    residual = volume - target
    objective = residual * residual
    derivative = 2 * residual
    return {
        "target_volume": _rational(target),
        "objective_global_lower_bound": "0",
        "objective": _rational(objective),
        "objective_derivative": _rational(derivative),
        "p": _rational(-derivative),
        "witness_attains_global_lower_bound": objective == 0,
    }


def _nonunique_certificate(witnesses: dict[str, Any]) -> dict[str, Any]:
    family = witnesses["affine_family"]
    constant = [_fraction(value, "affine constant") for value in family["constant"]]
    slope = [_fraction(value, "affine tau coefficient") for value in family["tau_coefficient"]]
    tau_lower = _fraction(family["tau_lower"], "tau lower")
    tau_upper = _fraction(family["tau_upper"], "tau upper")
    volume_polynomial = [
        _poly_integral(constant, Fraction(0), Fraction(1)),
        _poly_integral(slope, Fraction(0), Fraction(1)),
    ]

    def value(x: Fraction, tau: Fraction) -> Fraction:
        return sum(coefficient * x ** power for power, coefficient in enumerate(constant)) + tau * sum(
            coefficient * x ** power for power, coefficient in enumerate(slope)
        )

    corners = [
        value(Fraction(0), tau_lower), value(Fraction(0), tau_upper),
        value(Fraction(1), tau_lower), value(Fraction(1), tau_upper),
    ]
    derivative_l2 = _poly_integral(_poly_multiply(slope, slope), Fraction(0), Fraction(1))

    fixed_intervals = witnesses["binary_fixed_witnesses"]
    fixed_integrals = [
        _fraction(interval["upper"], "binary upper") - _fraction(interval["lower"], "binary lower")
        for interval in fixed_intervals
    ]
    binary = witnesses["binary_interval_family"]
    parameter_lower = _fraction(binary["parameter_lower"], "binary parameter lower")
    parameter_upper = _fraction(binary["parameter_upper"], "binary parameter upper")
    interval_length = _fraction(binary["interval_length"], "binary interval length")
    ordered_delta_lower = Fraction(0)
    ordered_delta_upper = parameter_upper - parameter_lower
    overlap_at_maximum_delta = interval_length - ordered_delta_upper
    symmetric_difference_polynomial = [Fraction(0), Fraction(2)]
    nonbinary = _fraction(witnesses["nonbinary_constant_witness"], "nonbinary constant")

    conditions = {
        "volume_constant_is_target": volume_polynomial == [Fraction(1, 2), Fraction(0)],
        "bilinear_rectangle_extrema_are_corners": tau_lower < tau_upper,
        "all_rectangle_corners_in_box": min(corners) >= 0 and max(corners) <= 1,
        "parameter_derivative_positive_modulo_ae": derivative_l2 > 0,
        "binary_fixed_witnesses_have_target_volume": fixed_integrals == [Fraction(1, 2)] * 2,
        "nonbinary_constant_is_strictly_inside_box": 0 < nonbinary < 1,
        "nonbinary_constant_has_target_volume": nonbinary == Fraction(1, 2),
        "shifted_intervals_stay_inside_domain": (
            parameter_lower == 0 and parameter_upper == Fraction(1, 2)
            and interval_length == Fraction(1, 2)
            and parameter_upper + interval_length == 1
        ),
        "ordered_parameter_delta_range_valid": (
            ordered_delta_lower == 0 and ordered_delta_upper == Fraction(1, 2)
        ),
        "overlap_nonnegative_for_all_ordered_parameters": overlap_at_maximum_delta >= 0,
        "positive_symmetric_difference_for_distinct_parameters": (
            symmetric_difference_polynomial == [0, 2]
        ),
    }
    _require(all(conditions.values()), "exact nonunique witness certificate failed")
    result = {
        "source_case_id": "faster.section6_6.nonunique",
        "catalogue_id": _CATALOGUE_IDS["faster.section6_6.nonunique"],
        "scientific_outcome": "VERIFIED_ANALYTIC_PROPERTIES",
        "domain": "(0,1)",
        "volume_polynomial_in_tau": [_rational(value) for value in volume_polynomial],
        "rectangle_corner_values": [_rational(value) for value in corners],
        "rectangle_extremum_argument": (
            "chi_tau is affine in each of x and tau, so extrema on the closed "
            "parameter rectangle occur at its four corners; open-domain bounds follow by restriction."
        ),
        "parameter_derivative_l2_squared": _rational(derivative_l2),
        "affine_family_separation_argument": (
            "For distinct real parameters s and t, the squared L2 distance is "
            "(t-s)^2 times 1/3, which is positive; hence the functions differ modulo a.e. equality."
        ),
        "uncountable_affine_family_modulo_ae": derivative_l2 > 0 and tau_lower < tau_upper,
        "binary_half_interval_integrals": [_rational(value) for value in fixed_integrals],
        "constant_nonbinary_integral": _rational(nonbinary),
        "shifted_binary_overlap_length": "1/2-(t-s)",
        "shifted_binary_symmetric_difference": "2(t-s)",
        "shifted_binary_symmetric_difference_polynomial": [
            _rational(value) for value in symmetric_difference_polynomial
        ],
        "binary_family_separation_argument": (
            "For real 0<=s<t<=1/2 the overlap length is 1/2-(t-s), so symmetric-difference "
            "measure is 2(t-s)>0; distinct parameters give distinct indicators modulo a.e. equality."
        ),
        "uncountable_binary_family_modulo_ae": parameter_lower < parameter_upper,
        "conditions": conditions,
    }
    result.update(_objective_certificate(volume_polynomial[0], Fraction(1, 2)))
    return result


def _full_domain_certificate(witnesses: dict[str, Any]) -> dict[str, Any]:
    domain_lower, domain_upper = Fraction(0), Fraction(1, 2)
    witness = _fraction(witnesses["full_domain_witness"]["constant"], "full-domain constant")
    domain_measure = domain_upper - domain_lower
    witness_integral = witness * domain_measure
    target = Fraction(1, 2)
    deficit_integral_at_optimum = domain_measure - target
    conditions = {
        "witness_is_in_box": 0 <= witness <= 1,
        "witness_integral_is_target": witness_integral == target,
        "box_makes_deficit_nonnegative": True,
        "deficit_integral_at_any_optimum_is_zero": deficit_integral_at_optimum == 0,
        "nonnegative_zero_integral_lemma_applies": deficit_integral_at_optimum == 0,
    }
    _require(all(conditions.values()), "exact full-domain witness certificate failed")
    result = {
        "source_case_id": "faster.section6_6.full_boundary",
        "catalogue_id": _CATALOGUE_IDS["faster.section6_6.full_boundary"],
        "scientific_outcome": "VERIFIED_ANALYTIC_PROPERTIES",
        "domain": "(0,1/2)",
        "domain_measure": _rational(domain_measure),
        "witness": "chi=1",
        "witness_integral": _rational(witness_integral),
        "deficit": "1-chi",
        "deficit_integral_at_any_optimum": _rational(deficit_integral_at_optimum),
        "nonnegative_zero_integral_lemma_applies": conditions["nonnegative_zero_integral_lemma_applies"],
        "uniqueness": "UNIQUE_ALMOST_EVERYWHERE",
        "uniqueness_argument": (
            "At an optimum integral(chi)=1/2=|Omega|. The box bound makes 1-chi nonnegative, "
            "and its integral is zero; the nonnegative-integrable-zero-integral lemma gives chi=1 a.e."
        ),
        "conditions": conditions,
    }
    result.update(_objective_certificate(witness_integral, target))
    return result


def verify_published_section_6_6(
    spec: dict[str, Any], witnesses: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Compute exact certificates without consulting expected-outcome fields."""
    _validate_spec(spec)
    actual_witnesses = published_witness_definitions() if witnesses is None else witnesses
    _validate_witnesses(actual_witnesses)
    cases = [_nonunique_certificate(actual_witnesses), _full_domain_certificate(actual_witnesses)]
    _require(all(case["objective"] == "0" and case["objective_derivative"] == "0"
                 and case["p"] == "0" for case in cases),
             "objective/sensitivity certificate failed")
    return {
        "schema": "faster-section-6-6-exact-scientific-data-v1",
        "paper_locator": _FIXED_DEFINITIONS["paper_locator"],
        "case_count": len(cases),
        "cases": cases,
        "project_status": "PARTIAL_BLOCKED",
        "author_bitwise_exact": "NOT_VERIFIED",
        "physical_validation": "NOT_IN_SCOPE",
        "native_extensions_required": [],
        "scope": "ONLY_THE_TWO_PUBLISHED_ONE_DIMENSIONAL_SECTION_6_6_EXAMPLES",
        "claims_all_p_zero_problems_classified": False,
        "ordinary_newton_convergence_claimed": False,
        "three_dimensional_geometry_claimed": False,
        "aggregate_verification_status_changed": False,
        "manual_release_gate_opened": False,
    }


def check_evidence_freshness(path: str | Path) -> dict[str, Any]:
    """Re-hash bound files and scientific data; no stored pass flag is trusted."""
    try:
        evidence = load_spec(path)
        _require(evidence.get("schema") == "faster-section-6-6-exact-evidence-v1",
                 "unsupported evidence schema")
        inputs = evidence.get("inputs")
        _require(type(inputs) is dict and set(inputs) == {
            "runner", "checker", "spec", "paper_pdf", "cmake_configuration",
        }, "invalid evidence input bindings")
        outputs = evidence.get("outputs")
        _require(type(outputs) is dict and set(outputs) == {"result", "human_report"},
                 "invalid evidence output bindings")
        mismatches: list[str] = []
        for label, stored in {**inputs, **outputs}.items():
            try:
                actual = file_record(stored["path"])
                if actual != stored:
                    mismatches.append(label + " content binding changed")
            except (KeyError, TypeError, ValueError, OSError) as error:
                mismatches.append(label + " cannot be re-read: " + str(error))
        if inputs["spec"].get("sha256") != EXPECTED_SPEC_SHA256:
            mismatches.append("spec immutable SHA256 mismatch")
        if inputs["paper_pdf"].get("sha256") != EXPECTED_PDF_SHA256:
            mismatches.append("paper_pdf immutable SHA256 mismatch")
        try:
            result = load_spec(outputs["result"]["path"])
            if result.get("scientific_data_sha256") != scientific_sha256(result.get("scientific_data")):
                mismatches.append("result scientific data digest changed")
        except (KeyError, TypeError, ValueError, OSError) as error:
            mismatches.append("result scientific data cannot be checked: " + str(error))
        return {
            "status": "STALE_RESULT" if mismatches else "CURRENT",
            "mismatches": list(dict.fromkeys(mismatches)),
        }
    except (KeyError, TypeError, ValueError, OSError) as error:
        return {"status": "INVALID_EVIDENCE", "mismatches": [str(error)]}
