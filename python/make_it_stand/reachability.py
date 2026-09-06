"""Legacy floating-relaxation schema; no physical infeasibility proof."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Callable, Iterable, Mapping

from .io import write_json_atomic, write_text_atomic


SCHEMA = "make-it-stand-carving-reachability-v1"
SCHEMA_VERSION = 1
STATUSES = frozenset({
    "ALREADY_SAFE",
    "PROVEN_UNREACHABLE",
    "RELAXATION_REACHABLE",
    "INDETERMINATE",
})
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def _finite(value: Any, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite number")
    return result


def _integer(value: Any, name: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")
    if value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return value


def _vector(value: Iterable[Any], name: str, length: int | None = None) -> tuple[float, ...]:
    try:
        result = tuple(
            _finite(item, f"{name}[{index}]")
            for index, item in enumerate(value)
        )
    except TypeError as error:
        raise ValueError(f"{name} must be an array") from error
    if length is not None and len(result) != length:
        raise ValueError(f"{name} must contain exactly {length} values")
    return result


@dataclass(frozen=True)
class FarkasCertificate:
    weights: tuple[float, ...]
    minimum_weight: float
    normalization_residual: float
    stationarity_infinity_norm: float
    contradiction: float
    valid: bool

    @classmethod
    def from_native(cls, native: Any) -> "FarkasCertificate":
        valid = native.valid
        if not isinstance(valid, bool):
            raise ValueError("certificate.valid must be a boolean")
        return cls(
            weights=_vector(native.weights, "certificate.weights"),
            minimum_weight=_finite(
                native.minimum_weight,
                "certificate.minimum_weight",
            ),
            normalization_residual=_finite(
                native.normalization_residual,
                "certificate.normalization_residual",
            ),
            stationarity_infinity_norm=_finite(
                native.stationarity_infinity_norm,
                "certificate.stationarity_infinity_norm",
            ),
            contradiction=_finite(
                native.contradiction,
                "certificate.contradiction",
            ),
            valid=valid,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "weights": list(self.weights),
            "minimum_weight": self.minimum_weight,
            "normalization_residual": self.normalization_residual,
            "stationarity_infinity_norm": self.stationarity_infinity_norm,
            "contradiction": self.contradiction,
            "valid": self.valid,
        }


@dataclass(frozen=True)
class ReachabilityResult:
    input_hash: str
    status: str
    input_cell_count: int
    variable_cell_count: int
    fixed_mass: float
    retained_fractions: tuple[float, ...]
    final_mass: float
    final_first_moment_2d: tuple[float, float]
    final_center_2d: tuple[float, float]
    maximum_primal_violation: float
    primal_solver_status: str
    certificate_solver_status: str
    certificate_system_nonzero_count: int
    certificate: FarkasCertificate

    @classmethod
    def from_native(cls, native: Any) -> "ReachabilityResult":
        input_hash = str(native.input_hash)
        if _SHA256.fullmatch(input_hash) is None:
            raise ValueError("input_hash must be a lower-case SHA-256 digest")
        status = str(native.status)
        if status not in STATUSES:
            raise ValueError(f"unsupported reachability status: {status!r}")
        input_cell_count = _integer(
            native.input_cell_count,
            "input_cell_count",
        )
        variable_cell_count = _integer(
            native.variable_cell_count,
            "variable_cell_count",
        )
        retained = _vector(native.retained_fractions, "retained_fractions")
        if retained and len(retained) != variable_cell_count:
            raise ValueError(
                "retained_fractions must match variable_cell_count"
            )
        if input_cell_count != variable_cell_count:
            raise ValueError(
                "input_cell_count must equal variable_cell_count in schema v1"
            )
        primal_status = str(native.primal_solver_status)
        certificate_status = str(native.certificate_solver_status)
        if not primal_status:
            raise ValueError("primal_solver_status must not be empty")
        return cls(
            input_hash=input_hash,
            status=status,
            input_cell_count=input_cell_count,
            variable_cell_count=variable_cell_count,
            fixed_mass=_finite(native.fixed_mass, "fixed_mass"),
            retained_fractions=retained,
            final_mass=_finite(native.final_mass, "final_mass"),
            final_first_moment_2d=_vector(
                native.final_first_moment_2d,
                "final_first_moment_2d",
                2,
            ),  # type: ignore[arg-type]
            final_center_2d=_vector(
                native.final_center_2d,
                "final_center_2d",
                2,
            ),  # type: ignore[arg-type]
            maximum_primal_violation=_finite(
                native.maximum_primal_violation,
                "maximum_primal_violation",
            ),
            primal_solver_status=primal_status,
            certificate_solver_status=certificate_status,
            certificate_system_nonzero_count=_integer(
                native.certificate_system_nonzero_count,
                "certificate_system_nonzero_count",
            ),
            certificate=FarkasCertificate.from_native(native.certificate),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "schema_version": SCHEMA_VERSION,
            "status_scope": "LEGACY_FLOATING_RELAXATION",
            "input_hash": self.input_hash,
            "status": self.status,
            "input_cell_count": self.input_cell_count,
            "variable_cell_count": self.variable_cell_count,
            "fixed_mass": self.fixed_mass,
            "retained_fractions": list(self.retained_fractions),
            "final_mass": self.final_mass,
            "final_first_moment_2d": list(self.final_first_moment_2d),
            "final_center_2d": list(self.final_center_2d),
            "maximum_primal_violation": self.maximum_primal_violation,
            "primal_solver_status": self.primal_solver_status,
            "certificate_solver_status": self.certificate_solver_status,
            "certificate_system_nonzero_count": (
                self.certificate_system_nonzero_count
            ),
            "certificate": self.certificate.to_dict(),
        }


def classify_result(status: str, observed_optimizer_outcome: str) -> str:
    """Label historical status names without upgrading floating evidence."""

    if status == "INDETERMINATE":
        return "LEGACY_FLOATING_RELAXATION_INDETERMINATE_REQUIRES_NUMERICAL_INVESTIGATION"
    if status == "PROVEN_UNREACHABLE":
        return "LEGACY_FLOATING_RELAXATION_REPORTED_UNREACHABLE"
    if observed_optimizer_outcome == "FAILED" and status in {
        "ALREADY_SAFE",
        "RELAXATION_REACHABLE",
    }:
        return "LEGACY_FLOATING_RELAXATION_CURRENT_PLANE_CARVER_SOLVER_GAP"
    return "LEGACY_FLOATING_RELAXATION_ALREADY_SAFE_OR_EXISTING_SUCCESS"


def _canonical_hash(value: Any) -> str:
    payload = json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _report_markdown(summary: Mapping[str, Any]) -> str:
    rows = [
        "# Carving reachability classification",
        "",
        "All statuses here are `LEGACY_FLOATING_RELAXATION` observations. ",
        "Retained fractions are continuous variables on discrete fixed cells; ",
        "they do not describe continuous geometric carving. `PROVEN_UNREACHABLE` ",
        "is a legacy solver label, not an independently exact or physical ",
        "infeasibility certificate. Geometry feasibility remains unknown.",
        "",
        "| Case | Prior optimizer outcome | Reachability | Cells | Certificate residual | Classification |",
        "|---|---|---|---:|---:|---|",
    ]
    for case in summary["cases"]:
        rows.append(
            "| {id} | {outcome} | {status} | {cells} | {residual:.3g} | {classification} |".format(
                id=case["id"],
                outcome=case["observed_optimizer_outcome"],
                status=case["result"]["status"],
                cells=case["result"]["variable_cell_count"],
                residual=case["result"]["certificate"][
                    "stationarity_infinity_norm"
                ],
                classification=case["classification"],
            )
        )
    rows.extend([
        "",
        f"Deterministic classification complete: `{str(summary['classification_complete']).lower()}`.",
        f"Legacy floating-relaxation release gate (not strict completion): `{str(summary['release_gate_passed']).lower()}`.",
        "",
    ])
    return "\n".join(rows)


def run_manifest(
    manifest: Mapping[str, Any],
    output_directory: str | Path,
    evaluator: Callable[
        [Mapping[str, Any]], ReachabilityResult | Mapping[str, Any]
    ],
) -> dict[str, Any]:
    """Evaluate a deterministic manifest and write auditable artifacts."""

    if manifest.get("schema") != "make-it-stand-reachability-cases-v1":
        raise ValueError("unsupported reachability case manifest schema")
    raw_cases = manifest.get("cases")
    if not isinstance(raw_cases, list) or not raw_cases:
        raise ValueError("manifest.cases must be a non-empty array")

    normalized: list[Mapping[str, Any]] = []
    identifiers: set[str] = set()
    for case in raw_cases:
        if not isinstance(case, Mapping):
            raise ValueError("every manifest case must be an object")
        identifier = case.get("id")
        outcome = case.get("observed_optimizer_outcome")
        if not isinstance(identifier, str) or not identifier:
            raise ValueError("every case requires a non-empty id")
        if identifier in identifiers:
            raise ValueError(f"duplicate case id: {identifier}")
        if outcome not in {"SUCCESS", "FAILED"}:
            raise ValueError(
                f"case {identifier}: observed_optimizer_outcome must be SUCCESS or FAILED"
            )
        identifiers.add(identifier)
        normalized.append(case)

    records: list[dict[str, Any]] = []
    for case in sorted(normalized, key=lambda value: str(value["id"])):
        evaluated = evaluator(case)
        result = (
            evaluated.to_dict()
            if isinstance(evaluated, ReachabilityResult)
            else dict(evaluated)
        )
        if result.get("schema") != SCHEMA or result.get("status") not in STATUSES:
            raise ValueError(f"case {case['id']}: evaluator returned an invalid result")
        certificate = result.get("certificate")
        if not isinstance(certificate, Mapping):
            raise ValueError(f"case {case['id']}: result certificate is invalid")
        weights = certificate.get("weights")
        if not isinstance(weights, list):
            raise ValueError(f"case {case['id']}: certificate weights are invalid")
        compact_certificate = dict(certificate)
        compact_certificate.pop("weights")
        compact_certificate["weights_count"] = len(weights)
        compact_certificate["weights_sha256"] = _canonical_hash(weights)
        compact_result = dict(result)
        compact_result["certificate"] = compact_certificate
        record = {
            "id": case["id"],
            "evidence_level": case.get("evidence_level", "RECONSTRUCTED_PARAMETER"),
            "observed_optimizer_outcome": case["observed_optimizer_outcome"],
            "case_definition": dict(case),
            "classification": classify_result(
                str(result["status"]),
                str(case["observed_optimizer_outcome"]),
            ),
            "case_definition_hash": _canonical_hash(case),
            "result_hash": _canonical_hash(result),
            "result": compact_result,
        }
        records.append(record)

    classifications = {record["classification"] for record in records}
    classification_complete = (
        "LEGACY_FLOATING_RELAXATION_INDETERMINATE_REQUIRES_NUMERICAL_INVESTIGATION" not in classifications
    )
    release_gate_passed = classification_complete and not any(
        record["classification"] == "LEGACY_FLOATING_RELAXATION_CURRENT_PLANE_CARVER_SOLVER_GAP"
        for record in records
    )
    summary = {
        "schema": "make-it-stand-reachability-report-v1",
        "schema_version": 1,
        "status_scope": "LEGACY_FLOATING_RELAXATION",
        "manifest_hash": _canonical_hash(manifest),
        "case_order": [record["id"] for record in records],
        "result_hashes": [record["result_hash"] for record in records],
        "classification_complete": classification_complete,
        "release_gate_passed": release_gate_passed,
        "cases": records,
    }
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    for record in records:
        write_json_atomic(output / f"{record['id']}.json", record)
    write_json_atomic(output / "summary.json", summary)
    write_text_atomic(output / "REPORT.md", _report_markdown(summary))
    return summary


__all__ = [
    "FarkasCertificate",
    "ReachabilityResult",
    "SCHEMA",
    "SCHEMA_VERSION",
    "STATUSES",
    "classify_result",
    "run_manifest",
]
