"""Independent, reproducible diagnostics for captured Faster Newton systems.

This module audits one *linear* Newton system.  It does not certify nonlinear
convergence, quadric feasibility, or validity of the material geometry.
"""

from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
from typing import Any

import numpy as np


_PAYLOAD_SCHEMA = "spin-newton-system-payload-v2"
_MAX_SYSTEM_SIZE = 32


def _auditor_identity() -> dict[str, str]:
    source_path = Path(__file__).resolve(strict=True)
    return {
        "source_path": str(source_path),
        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
    }


def _exact_binary64_rank(matrix: np.ndarray) -> int:
    """Compute the exact rational rank of the stored IEEE-754 binary64 values."""
    rows = [[Fraction.from_float(float(value)) for value in row] for row in matrix]
    row_count, column_count = matrix.shape
    pivot_row = 0
    for column in range(column_count):
        pivot = next((index for index in range(pivot_row, row_count)
                      if rows[index][column] != 0), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        pivot_value = rows[pivot_row][column]
        for index in range(pivot_row + 1, row_count):
            if rows[index][column] == 0:
                continue
            factor = rows[index][column] / pivot_value
            for offset in range(column, column_count):
                rows[index][offset] -= factor * rows[pivot_row][offset]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def _decimal_backward_error(matrix: np.ndarray, solution: np.ndarray,
                            rhs: np.ndarray) -> tuple[float, float]:
    """Evaluate the original-system residual from exact decimal forms of doubles."""
    with localcontext() as context:
        context.prec = 80
        a = [[Decimal.from_float(float(value)) for value in row] for row in matrix]
        x = [Decimal.from_float(float(value)) for value in solution]
        b = [Decimal.from_float(float(value)) for value in rhs]
        residuals = [sum((a[i][j] * x[j] for j in range(len(x))), Decimal(0)) - b[i]
                     for i in range(len(b))]
        residual_inf = max((abs(value) for value in residuals), default=Decimal(0))
        matrix_inf = max((sum((abs(value) for value in row), Decimal(0)) for row in a),
                         default=Decimal(0))
        solution_inf = max((abs(value) for value in x), default=Decimal(0))
        rhs_inf = max((abs(value) for value in b), default=Decimal(0))
        denominator = matrix_inf * solution_inf + rhs_inf
        backward_error = residual_inf / denominator if denominator else Decimal("Infinity")
        return float(residual_inf), float(backward_error)


def audit_newton_system_payload(payload: dict[str, Any]) -> dict[str, Any]:
    """Audit rank and a reversibly equilibrated solve for a captured Newton system.

    Exact rank is the rank over the rational numbers represented by the input
    binary64 values.  Numerical ranks use NumPy's standard SVD threshold and
    are diagnostic only.  The solve is accepted only when a high-precision
    residual check on the *unscaled original system* is small.
    """
    if not isinstance(payload, dict) or payload.get("schema") != _PAYLOAD_SCHEMA:
        raise ValueError(f"Newton payload schema must be {_PAYLOAD_SCHEMA}")
    try:
        matrix = np.asarray(payload["matrix"], dtype=np.float64)
        rhs = np.asarray(payload["right_hand_side"], dtype=np.float64)
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("Newton payload must contain a numeric matrix and right_hand_side") from error
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or matrix.shape[0] == 0:
        raise ValueError("Newton matrix must be nonempty and square")
    if matrix.shape[0] > _MAX_SYSTEM_SIZE:
        raise ValueError(f"Newton matrix size exceeds audit limit {_MAX_SYSTEM_SIZE}")
    if rhs.shape != (matrix.shape[0],):
        raise ValueError("Newton right_hand_side must match the square matrix dimension")
    if not np.isfinite(matrix).all() or not np.isfinite(rhs).all():
        raise ValueError("Newton matrix and right_hand_side must be finite")

    size = matrix.shape[0]
    exact_rank = _exact_binary64_rank(matrix)
    original_singular_values = np.linalg.svd(matrix, compute_uv=False)
    epsilon_threshold = np.finfo(np.float64).eps * size * original_singular_values[0]
    original_rank = int(np.count_nonzero(original_singular_values > epsilon_threshold))

    row_scales = np.max(np.abs(matrix), axis=1)
    if np.any(row_scales == 0.0):
        return {
            "schema": "spin-newton-system-audit-v1",
            "matrix_size": size,
            "exact_binary64_rank": exact_rank,
            "original_numerical_rank": original_rank,
            "scaled_numerical_rank": None,
            "interpretation": "EXACTLY_RANK_DEFICIENT" if exact_rank < size else "ZERO_ROW",
            "solution": None,
            "original_system_residual_inf": None,
            "original_system_backward_error": None,
        }
    row_scaled = matrix / row_scales[:, None]
    scaled_rhs = rhs / row_scales
    column_scales = np.max(np.abs(row_scaled), axis=0)
    if np.any(column_scales == 0.0):
        return {
            "schema": "spin-newton-system-audit-v1",
            "matrix_size": size,
            "exact_binary64_rank": exact_rank,
            "original_numerical_rank": original_rank,
            "scaled_numerical_rank": None,
            "interpretation": "EXACTLY_RANK_DEFICIENT" if exact_rank < size else "ZERO_COLUMN",
            "solution": None,
            "original_system_residual_inf": None,
            "original_system_backward_error": None,
        }
    equilibrated = row_scaled / column_scales[None, :]
    scaled_singular_values = np.linalg.svd(equilibrated, compute_uv=False)
    scaled_threshold = np.finfo(np.float64).eps * size * scaled_singular_values[0]
    scaled_rank = int(np.count_nonzero(scaled_singular_values > scaled_threshold))

    base = {
        "schema": "spin-newton-system-audit-v1",
        "matrix_size": size,
        "exact_binary64_rank": exact_rank,
        "original_numerical_rank": original_rank,
        "scaled_numerical_rank": scaled_rank,
        "original_min_singular_value": float(original_singular_values[-1]),
        "scaled_min_singular_value": float(scaled_singular_values[-1]),
        "original_rank_threshold": float(epsilon_threshold),
        "scaled_rank_threshold": float(scaled_threshold),
        "solution": None,
        "original_system_residual_inf": None,
        "original_system_backward_error": None,
    }
    if exact_rank < size:
        base["interpretation"] = "EXACTLY_RANK_DEFICIENT"
        return base
    if scaled_rank < size:
        base["interpretation"] = "NUMERICALLY_UNRESOLVED_AFTER_EQUILIBRATION"
        return base

    try:
        equilibrated_solution = np.linalg.solve(equilibrated, scaled_rhs)
        solution = equilibrated_solution / column_scales
    except np.linalg.LinAlgError:
        base["interpretation"] = "NUMERICAL_SOLVE_FAILED"
        return base
    if not np.isfinite(solution).all():
        base["interpretation"] = "NUMERICAL_SOLVE_NONFINITE"
        return base
    residual_inf, backward_error = _decimal_backward_error(matrix, solution, rhs)
    base["solution"] = [float(value) for value in solution]
    base["original_system_residual_inf"] = residual_inf
    base["original_system_backward_error"] = backward_error
    backward_error_limit = 128.0 * np.finfo(np.float64).eps * size
    base["backward_error_limit"] = float(backward_error_limit)
    if not math.isfinite(backward_error) or backward_error > backward_error_limit:
        base["interpretation"] = "ORIGINAL_SYSTEM_BACKWARD_ERROR_FAILED"
        return base
    base["interpretation"] = (
        "SCALE_SENSITIVE_FULL_RANK"
        if original_rank < size
        else "FULL_RANK_SOLVED_AND_VERIFIED"
    )
    return base


def audit_newton_event_journal(path: str | Path) -> dict[str, Any]:
    """Replay all captured Newton matrices in a streamed solver event journal."""
    journal_path = Path(path)
    digest = hashlib.sha256()
    events: list[dict[str, Any]] = []
    with journal_path.open("rb") as stream:
        line_number = 0
        for raw_line in stream:
            line_number += 1
            digest.update(raw_line)
            try:
                row = json.loads(raw_line)
            except (UnicodeDecodeError, json.JSONDecodeError) as error:
                raise ValueError(f"invalid JSON event journal line {line_number}") from error
            if not isinstance(row, dict) or row.get("record_type") != "solver_event":
                continue
            record = row.get("record")
            if not isinstance(record, dict):
                raise ValueError(f"solver event line {line_number} has no record object")
            payload = record.get("newton_system_payload")
            if payload is None:
                continue
            audit = audit_newton_system_payload(payload)
            events.append({
                "journal_line": line_number,
                "restart": record.get("restart"),
                "iteration": record.get("iteration"),
                "evaluation_index": record.get("evaluation_index"),
                "reason": record.get("reason"),
                "audit": audit,
            })
    return {
        "schema": "spin-newton-event-journal-audit-v1",
        "journal_path": str(journal_path.resolve()),
        "journal_sha256": digest.hexdigest(),
        "auditor": _auditor_identity(),
        "payload_count": len(events),
        "events": events,
        "scope_warning": (
            "linear-system rank and residual diagnostics only; this does not certify "
            "nonlinear Newton convergence or material geometry"
        ),
    }


def _main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("journal", type=Path, help="solver-events.jsonl to replay")
    parser.add_argument("--output", type=Path, help="optional JSON output path")
    args = parser.parse_args()
    report = audit_newton_event_journal(args.journal)
    encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is None:
        print(encoded, end="")
    else:
        args.output.write_text(encoded, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
