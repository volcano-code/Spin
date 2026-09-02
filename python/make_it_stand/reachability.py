"""Strict, deterministic carving-reachability result schema."""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Any, Iterable


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
            certificate=FarkasCertificate.from_native(native.certificate),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "schema_version": SCHEMA_VERSION,
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
            "certificate": self.certificate.to_dict(),
        }
