"""Shared CLI, diagnostics, reporting, and batch helpers for R9."""

from __future__ import annotations

from contextlib import contextmanager
from html import escape
import math
import os
from pathlib import Path
import sys
from typing import Any, Iterable, Mapping, Sequence

from .io import (
    DataFormatError,
    prepare_output_path,
    prepare_output_root,
    read_json,
    write_csv_atomic,
    write_json_atomic,
    write_text_atomic,
)
from .model import (
    BalanceConfig,
    ConfigurationError,
    DIAGNOSTICS_SCHEMA,
    run_config,
    validate_balance_config,
)


class DiagnosticsError(ValueError):
    """A diagnostics artifact is malformed or unsupported."""



def _validate_diagnostics(value: Any, source: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise DiagnosticsError(f"{source}: diagnostics must be a JSON object")
    if value.get("schema") != DIAGNOSTICS_SCHEMA:
        raise DiagnosticsError(
            f"{source}: diagnostics.schema must be {DIAGNOSTICS_SCHEMA!r}"
        )
    version = value.get("schema_version")
    if type(version) is not int or version != 1:
        raise DiagnosticsError(f"{source}: diagnostics.schema_version must be 1")
    status = value.get("status")
    if not isinstance(status, str) or not status:
        raise DiagnosticsError(f"{source}: diagnostics.status must be a non-empty string")
    if not isinstance(value.get("succeeded"), bool):
        raise DiagnosticsError(f"{source}: diagnostics.succeeded must be a boolean")

    top_ready = value.get("manufacturing_ready")
    if "manufacturing_ready" in value and not isinstance(top_ready, bool):
        raise DiagnosticsError(
            f"{source}: diagnostics.manufacturing_ready must be a boolean"
        )
    if "manufacturing" not in value:
        raise DiagnosticsError(
            f"{source}: diagnostics.manufacturing must be present"
        )
    manufacturing = value["manufacturing"]
    if not isinstance(manufacturing, Mapping):
        raise DiagnosticsError(
            f"{source}: diagnostics.manufacturing must be an object"
        )
    if "ready" not in manufacturing or not isinstance(manufacturing["ready"], bool):
        raise DiagnosticsError(
            f"{source}: diagnostics.manufacturing.ready must be a boolean"
        )
    nested_ready: Any = manufacturing["ready"]
    if (
        "manufacturing_ready" in value
        and "manufacturing" in value
        and top_ready is not nested_ready
    ):
        raise DiagnosticsError(
            f"{source}: top-level and nested manufacturing readiness disagree"
        )

    for field in ("carving_iterations", "deformation_iterations"):
        if field not in value or not isinstance(value[field], list):
            raise DiagnosticsError(f"{source}: diagnostics.{field} must be an array")

    def require_finite(number: Any, path: str) -> None:
        if (
            isinstance(number, bool)
            or not isinstance(number, (int, float))
            or not math.isfinite(float(number))
        ):
            raise DiagnosticsError(f"{source}: {path} must be a finite number")

    if "mass" in value:
        require_finite(value["mass"], "diagnostics.mass")
    if "mass_properties" not in value:
        raise DiagnosticsError(
            f"{source}: diagnostics.mass_properties must be present"
        )
    properties = value["mass_properties"]
    if not isinstance(properties, Mapping):
        raise DiagnosticsError(
            f"{source}: diagnostics.mass_properties must be an object"
        )
    if "mass" not in properties:
        raise DiagnosticsError(
            f"{source}: diagnostics.mass_properties.mass must be present"
        )
    for field in ("mass", "signed_mass"):
        if field in properties:
            require_finite(
                properties[field],
                f"diagnostics.mass_properties.{field}",
            )
    return value


def load_diagnostics(path: os.PathLike[str] | str) -> dict[str, Any]:
    return _validate_diagnostics(read_json(path), os.fspath(path))


def _finite_metric(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    result = float(value)
    return result if math.isfinite(result) else None


def _path_get(value: Mapping[str, Any], *paths: Sequence[str]) -> Any:
    for path in paths:
        current: Any = value
        for component in path:
            if not isinstance(current, Mapping) or component not in current:
                break
            current = current[component]
        else:
            return current
    return None


def _trace_count(value: Mapping[str, Any], name: str) -> int:
    trace = value.get(name)
    if isinstance(trace, list):
        return len(trace)
    count = _finite_metric(trace)
    return int(count) if count is not None and count >= 0 and count.is_integer() else 0


def summarize_diagnostics(records: Iterable[Mapping[str, Any]]) -> dict[str, Any]:
    """Summarize existing diagnostics without making paper-equivalence claims."""

    items = list(records)
    statuses: dict[str, int] = {}
    succeeded = 0
    manufacturing_ready = 0
    carving_iterations = 0
    deformation_iterations = 0
    masses: list[float] = []
    for index, item in enumerate(items):
        if not isinstance(item, dict):
            raise DiagnosticsError(f"record {index}: expected a JSON object")
        _validate_diagnostics(item, f"record {index}")
        status = item["status"]
        statuses[status] = statuses.get(status, 0) + 1
        succeeded += int(item.get("succeeded") is True)
        ready = _path_get(item, ("manufacturing_ready",), ("manufacturing", "ready"))
        manufacturing_ready += int(ready is True)
        carving_iterations += _trace_count(item, "carving_iterations")
        deformation_iterations += _trace_count(item, "deformation_iterations")
        mass = _finite_metric(
            _path_get(
                item,
                ("mass",),
                ("mass_properties", "mass"),
                ("mass_properties", "signed_mass"),
            )
        )
        if mass is not None:
            masses.append(mass)
    return {
        "schema": "make-it-stand-r9-diagnostics-summary-v1",
        "schema_version": 1,
        "record_count": len(items),
        "succeeded_count": succeeded,
        "manufacturing_ready_count": manufacturing_ready,
        "status_counts": dict(sorted(statuses.items())),
        "carving_iteration_count": carving_iterations,
        "deformation_iteration_count": deformation_iterations,
        "mean_mass": sum(masses) / len(masses) if masses else None,
    }


def compare_diagnostics(groups: Mapping[str, Sequence[Mapping[str, Any]]]) -> dict[str, Any]:
    """Compare named, pre-existing result groups; this does not run an ablation."""

    if not isinstance(groups, Mapping) or not groups:
        raise DiagnosticsError("at least one named diagnostics group is required")
    comparisons: dict[str, Any] = {}
    for name in sorted(groups):
        if not isinstance(name, str) or not name:
            raise DiagnosticsError("group names must be non-empty strings")
        records = groups[name]
        if isinstance(records, (str, bytes)) or not isinstance(records, Sequence):
            raise DiagnosticsError(f"group {name!r} must contain an array of records")
        comparisons[name] = summarize_diagnostics(records)
    return {
        "schema": "make-it-stand-r9-generic-diagnostics-comparison-v1",
        "schema_version": 1,
        "claim": (
            "Generic comparison of supplied diagnostics only; no R10 paper "
            "ablation or paper-equivalence experiment was run."
        ),
        "groups": comparisons,
    }


@contextmanager
def _exclusive_batch_claim(root: Path) -> Iterable[None]:
    """Hold a process-safe claim for one batch output root."""

    claim = root / ".r9-batch.claim"
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    flags |= getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(claim, flags, 0o600)
    except FileExistsError as error:
        raise FileExistsError(
            f"refusing active or stale batch output claim: {claim}"
        ) from error
    claimed_metadata = os.fstat(descriptor)
    try:
        yield
    finally:
        os.close(descriptor)
        try:
            current_metadata = os.lstat(claim)
        except FileNotFoundError:
            pass
        else:
            if (
                current_metadata.st_dev == claimed_metadata.st_dev
                and current_metadata.st_ino == claimed_metadata.st_ino
            ):
                os.unlink(claim)


def run_batch(
    configs: Sequence[BalanceConfig],
    *,
    output_directory: os.PathLike[str] | str,
    native_module: Any | None = None,
    keep_going: bool = False,
) -> dict[str, Any]:
    """Run validated cases in order and write a path-free batch summary."""

    validated_configs = tuple(validate_balance_config(config) for config in configs)
    root = prepare_output_root(output_directory)
    summary_path = root / "batch_summary.json"
    case_directories = [
        root / f"case-{index:04d}" for index in range(len(validated_configs))
    ]
    with _exclusive_batch_claim(root):
        if summary_path.exists() or summary_path.is_symlink():
            raise FileExistsError(
                f"refusing to overwrite existing batch artifact: {summary_path}"
            )
        for case_directory in case_directories:
            if case_directory.is_symlink() or case_directory.exists():
                raise FileExistsError(
                    f"refusing to reuse existing batch case path: {case_directory}"
                )
        cases: list[dict[str, Any]] = []
        for index, config in enumerate(validated_configs):
            case_directory = case_directories[index]
            try:
                os.mkdir(case_directory)
            except FileExistsError as error:
                raise FileExistsError(
                    f"refusing to reuse existing batch case path: {case_directory}"
                ) from error
            try:
                result = run_config(
                    config,
                    output_directory=case_directory,
                    native_module=native_module,
                )
                cases.append(
                    {
                        "index": index,
                        "completed": True,
                        "status": result.status,
                        "succeeded": result.succeeded,
                        "manufacturing_ready": result.manufacturing_ready,
                    }
                )
            except Exception as error:  # preserve completed cases in keep-going mode
                if not keep_going:
                    raise
                cases.append(
                    {
                        "index": index,
                        "completed": False,
                        "error_type": type(error).__name__,
                        "error": str(error),
                    }
                )
        summary = {
            "schema": "make-it-stand-r9-batch-summary-v1",
            "schema_version": 1,
            "case_count": len(validated_configs),
            "completed_count": sum(case.get("completed") is True for case in cases),
            "cases": cases,
        }
        write_json_atomic(
            prepare_output_path(root, "batch_summary.json"),
            summary,
        )
        return summary


def diagnostics_trace_rows(value: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Flatten numeric carving/deformation trace fields for CSV/SVG."""

    rows: list[dict[str, Any]] = []
    for kind, key in (("carving", "carving_iterations"), ("deformation", "deformation_iterations")):
        trace = value.get(key, [])
        if not isinstance(trace, list):
            continue
        for ordinal, record in enumerate(trace):
            if not isinstance(record, Mapping):
                continue
            iteration = record.get("iteration", ordinal)
            block = record.get("block", 0)
            for metric in sorted(record):
                if metric in {"block", "iteration"}:
                    continue
                number = _finite_metric(record[metric])
                if number is None:
                    continue
                rows.append(
                    {
                        "trace": kind,
                        "block": block,
                        "iteration": iteration,
                        "metric": metric,
                        "value": format(number, ".17g"),
                    }
                )
    return rows


def write_trace_csv(path: os.PathLike[str] | str, diagnostics: Mapping[str, Any]) -> None:
    write_csv_atomic(
        path,
        ("trace", "block", "iteration", "metric", "value"),
        diagnostics_trace_rows(diagnostics),
    )


def _svg_series(diagnostics: Mapping[str, Any]) -> tuple[str, list[float]]:
    preferences = (
        "refreshed_energy_after", "energy_after", "frozen_energy_after",
        "gradient_norm", "trial_energy", "energy_before",
    )
    rows = diagnostics_trace_rows(diagnostics)
    for preferred in preferences:
        values = [float(row["value"]) for row in rows if row["metric"] == preferred]
        if values:
            return preferred, values
    for row in rows:
        metric = str(row["metric"])
        values = [float(item["value"]) for item in rows if item["metric"] == metric]
        if values:
            return metric, values
    return "no numeric trace", []


def write_svg_report(
    path: os.PathLike[str] | str,
    diagnostics: Mapping[str, Any],
    *,
    title: str = "Make It Stand diagnostics",
) -> None:
    """Write a small deterministic SVG trace plot using only stdlib."""

    metric, values = _svg_series(diagnostics)
    width, height = 800, 360
    left, top, plot_width, plot_height = 70, 70, 680, 230
    points = ""
    if values:
        low, high = min(values), max(values)
        span = high - low
        coordinates: list[str] = []
        for index, value in enumerate(values):
            x = left + (plot_width * index / max(1, len(values) - 1))
            y = top + plot_height / 2 if span == 0.0 else top + plot_height * (high - value) / span
            coordinates.append(f"{x:.3f},{y:.3f}")
        points = (
            f'<polyline fill="none" stroke="#285f91" stroke-width="2" '
            f'points="{" ".join(coordinates)}" />'
        )
    status = escape(str(diagnostics.get("status", "unknown")))
    svg = f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="100%" height="100%" fill="#ffffff" />
  <text x="{left}" y="30" font-family="sans-serif" font-size="18">{escape(title)}</text>
  <text x="{left}" y="52" font-family="sans-serif" font-size="12">status: {status}; metric: {escape(metric)}</text>
  <line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#333333" />
  <line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}" stroke="#333333" />
  {points}
  <text x="{left}" y="340" font-family="sans-serif" font-size="11">R9 generic diagnostics visualization; not an R10 paper experiment.</text>
</svg>
'''
    write_text_atomic(path, svg)


def cli_error(program: str, error: BaseException) -> int:
    print(f"{program}: error: {error}", file=sys.stderr)
    return 2 if isinstance(error, (ConfigurationError, DataFormatError, DiagnosticsError, ValueError)) else 1


__all__ = [
    "DiagnosticsError", "cli_error", "compare_diagnostics",
    "diagnostics_trace_rows", "load_diagnostics", "run_batch",
    "summarize_diagnostics", "write_svg_report", "write_trace_csv",
]
