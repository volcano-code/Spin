#!/usr/bin/env python3
"""Validate an R10 evidence tree and render a portable local-evidence report.

The report is intentionally a presentation layer over the canonical R10
documents.  It revalidates every manifest, run record, aggregate, and bound
artifact before writing anything.  Paper-reference values remain in their own
reference-only section and are never plotted or pooled with local measurements.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from html import escape as xml_escape
import io
import math
import os
from pathlib import Path, PurePosixPath
import stat
from typing import Any, Mapping, Sequence

from make_it_stand.io import loads_json_strict
from make_it_stand.r10 import (
    aggregate_run_records,
    read_canonical_json,
    validate_aggregate,
    validate_manifest,
    validate_relative_path,
    validate_run_record,
)


REPORT_SCOPE = (
    "Local synthetic evidence from the modern implementation; not an exact "
    "paper reproduction."
)


class R10ReportError(ValueError):
    """Raised when an evidence tree cannot support a strict R10 report."""


def _fail(message: str) -> None:
    raise R10ReportError(message)


def _artifact_path(case_root: Path, relative_path: str) -> Path:
    relative = validate_relative_path(relative_path, "artifact.relative_path")
    anchor = case_root.resolve(strict=True)
    candidate = anchor
    for component in PurePosixPath(relative).parts:
        candidate /= component
        try:
            metadata = os.lstat(candidate)
        except FileNotFoundError:
            _fail(f"missing bound artifact: {relative}")
        if stat.S_ISLNK(metadata.st_mode):
            _fail(f"symbolic-link bound artifact is forbidden: {relative}")
    resolved = candidate.resolve(strict=True)
    try:
        resolved.relative_to(anchor)
    except ValueError:
        _fail(f"bound artifact escapes its case directory: {relative}")
    if not resolved.is_file():
        _fail(f"bound artifact is not an ordinary file: {relative}")
    return resolved


def _verify_binding(case_root: Path, binding: Mapping[str, Any]) -> Path:
    path = _artifact_path(case_root, str(binding["relative_path"]))
    digest = hashlib.sha256()
    byte_count = 0
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
            byte_count += len(chunk)
    if digest.hexdigest() != binding["sha256"]:
        _fail(f"SHA-256 mismatch for bound artifact: {binding['relative_path']}")
    if "byte_count" in binding and byte_count != binding["byte_count"]:
        _fail(f"byte-count mismatch for bound artifact: {binding['relative_path']}")
    return path


def _discover_case_roots(evidence_root: Path) -> tuple[Path, ...]:
    root = evidence_root.resolve(strict=True)
    if (root / "manifest.json").is_file():
        return (root,)
    cases = tuple(
        child
        for child in sorted(root.iterdir(), key=lambda item: item.name)
        if child.is_dir() and not child.is_symlink() and (child / "manifest.json").is_file()
    )
    if not cases:
        _fail("no R10 case directory containing manifest.json was found")
    return cases


def _metric_map(record: Mapping[str, Any]) -> dict[str, float]:
    return {
        str(metric["name"]): float(metric["value"])
        for metric in record["local_metrics"]
    }


def _read_convergence_curve(
    case_root: Path,
    records: Sequence[Mapping[str, Any]],
) -> tuple[tuple[int, float], ...]:
    """Read the first successful record's bound diagnostic event ledger."""

    for record in sorted(records, key=lambda item: int(item["repeat_index"])):
        if record["status"] != "succeeded":
            continue
        diagnostics_binding = next(
            (
                binding
                for binding in record["artifacts"]
                if binding["role"] == "diagnostics"
            ),
            None,
        )
        if diagnostics_binding is None:
            continue
        path = _verify_binding(case_root, diagnostics_binding)
        try:
            diagnostics = loads_json_strict(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, ValueError) as error:
            _fail(f"invalid bound diagnostics JSON: {diagnostics_binding['relative_path']} ({error})")
        if not isinstance(diagnostics, Mapping):
            _fail(f"bound diagnostics must be an object: {diagnostics_binding['relative_path']}")
        raw_events = diagnostics.get("r10_events", [])
        if not isinstance(raw_events, list):
            _fail(f"r10_events must be an array: {diagnostics_binding['relative_path']}")
        curve: list[tuple[int, float]] = []
        previous = -1
        for event in raw_events:
            if not isinstance(event, Mapping):
                _fail(f"r10_events contains a non-object: {diagnostics_binding['relative_path']}")
            ordinal = event.get("ordinal")
            energy = event.get("projected_target_energy")
            if isinstance(ordinal, bool) or not isinstance(ordinal, int) or ordinal <= previous:
                _fail(f"r10 event ordinals are not strictly increasing: {diagnostics_binding['relative_path']}")
            if isinstance(energy, bool) or not isinstance(energy, (int, float)):
                _fail(f"r10 event energy is not numeric: {diagnostics_binding['relative_path']}")
            measured = float(energy)
            if not math.isfinite(measured) or measured < 0.0:
                _fail(f"r10 event energy is invalid: {diagnostics_binding['relative_path']}")
            curve.append((ordinal, measured))
            previous = ordinal
        return tuple(curve)
    return ()


def load_evidence_tree(evidence_root: os.PathLike[str] | str) -> tuple[dict[str, Any], ...]:
    """Load and cross-check all complete cases below *evidence_root*."""

    cases: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    for case_root in _discover_case_roots(Path(evidence_root)):
        manifest = read_canonical_json(case_root / "manifest.json")
        validate_manifest(manifest)
        experiment_id = str(manifest["experiment_id"])
        if experiment_id in seen_ids:
            _fail(f"duplicate experiment_id in evidence tree: {experiment_id}")
        seen_ids.add(experiment_id)

        for binding in manifest["inputs"]:
            _verify_binding(case_root, binding)
        for binding in manifest["provenance"]["resources"]:
            _verify_binding(case_root, binding)

        records_directory = case_root / "records"
        if not records_directory.is_dir() or records_directory.is_symlink():
            _fail(f"missing ordinary records directory for case: {experiment_id}")
        record_paths = sorted(records_directory.glob("*.json"))
        if len(record_paths) != manifest["protocol"]["repeats"]:
            _fail(
                f"case {experiment_id} has {len(record_paths)} record files but "
                f"declares {manifest['protocol']['repeats']} repeats"
            )
        records: list[Mapping[str, Any]] = []
        for path in record_paths:
            if path.is_symlink():
                _fail(f"symbolic-link run record is forbidden in case: {experiment_id}")
            record = read_canonical_json(path)
            validate_run_record(record, manifest=manifest)
            for binding in record["artifacts"]:
                _verify_binding(case_root, binding)
            records.append(record)
        records.sort(key=lambda item: int(item["repeat_index"]))

        aggregate = read_canonical_json(case_root / "aggregate.json")
        validate_aggregate(aggregate, manifest=manifest)
        expected = aggregate_run_records(manifest, records)
        if aggregate != expected:
            _fail(f"aggregate is stale or inconsistent with records for case: {experiment_id}")

        local_values: dict[str, list[float]] = {}
        for record in records:
            if record["status"] != "succeeded":
                continue
            for name, value in _metric_map(record).items():
                local_values.setdefault(name, []).append(value)
        cases.append(
            {
                "manifest": manifest,
                "records": tuple(records),
                "aggregate": aggregate,
                "local_values": local_values,
                "convergence": _read_convergence_curve(case_root, records),
            }
        )
    cases.sort(key=lambda item: item["manifest"]["experiment_id"])
    return tuple(cases)


def _format_number(value: Any) -> str:
    if value is None or value == "":
        return ""
    return format(float(value), ".8g")


def _markdown(value: Any) -> str:
    return str(value).replace("\\", "\\\\").replace("|", "\\|").replace("\r", " ").replace("\n", " ")


def _summary_csv(cases: Sequence[Mapping[str, Any]]) -> str:
    stream = io.StringIO(newline="")
    columns = (
        "source", "experiment_id", "variant", "claim_level",
        "repeat_count", "successful_repeat_count", "metric_name", "unit",
        "sample_count", "median", "q1", "q3", "iqr", "minimum", "maximum",
        "reference_value", "paper_location", "comparability",
    )
    writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
    writer.writeheader()
    for case in cases:
        manifest = case["manifest"]
        aggregate = case["aggregate"]
        common = {
            "experiment_id": manifest["experiment_id"],
            "variant": manifest["protocol"]["variant"],
            "claim_level": manifest["claim"]["level"],
            "repeat_count": aggregate["repeat_count"],
            "successful_repeat_count": aggregate["successful_repeat_count"],
        }
        for metric in aggregate["local_metric_aggregates"]:
            writer.writerow(
                {
                    **common,
                    "source": "local_measurement",
                    "metric_name": metric["name"],
                    "unit": metric["unit"],
                    "sample_count": metric["sample_count"],
                    "median": _format_number(metric["median"]),
                    "q1": _format_number(metric["q1"]),
                    "q3": _format_number(metric["q3"]),
                    "iqr": _format_number(metric["iqr"]),
                    "minimum": _format_number(metric["minimum"]),
                    "maximum": _format_number(metric["maximum"]),
                }
            )
    for case in cases:
        manifest = case["manifest"]
        aggregate = case["aggregate"]
        common = {
            "experiment_id": manifest["experiment_id"],
            "variant": manifest["protocol"]["variant"],
            "claim_level": manifest["claim"]["level"],
            "repeat_count": aggregate["repeat_count"],
            "successful_repeat_count": aggregate["successful_repeat_count"],
        }
        for metric in aggregate["paper_reference_metrics"]:
            writer.writerow(
                {
                    **common,
                    "source": "paper_reference",
                    "metric_name": metric["name"],
                    "unit": metric["unit"],
                    "reference_value": _format_number(metric["value"]),
                    "paper_location": metric["paper_location"],
                    "comparability": "reference_only",
                }
            )
    return stream.getvalue()


def _summary_markdown(cases: Sequence[Mapping[str, Any]]) -> str:
    lines = [
        "# R10 evidence summary",
        "",
        f"> Scope: {REPORT_SCOPE}",
        "> Paper-reference values are reference-only and are never pooled with local measurements.",
        "",
        "## Integrity",
        "",
        f"Validated {len(cases)} case(s), their canonical manifests, complete run-record sets, canonical aggregates, and SHA-256-bound artifacts.",
        "",
        "## Local measurements",
        "",
        "All rows below are `local_measurement` values from repository-defined synthetic or local protocols.",
        "",
        "| Case | Variant | Runs | Algorithm success | Targets reached | Manufacturing ready | Target energy | Removed volume | RMS displacement | Total time (s) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for case in cases:
        manifest = case["manifest"]
        aggregate = case["aggregate"]
        values = case["local_values"]

        def rate(name: str) -> str:
            samples = values.get(name, [])
            return "n/a" if not samples else f"{sum(samples) / len(samples):.3f}"

        def median(name: str, *, timing: bool = False) -> str:
            source = aggregate["timing_aggregates"] if timing else aggregate["local_metric_aggregates"]
            item = next((candidate for candidate in source if candidate["name"] == name), None)
            return "n/a" if item is None else _format_number(item["median"])

        lines.append(
            "| "
            + " | ".join(
                (
                    _markdown(manifest["experiment_id"]),
                    _markdown(manifest["protocol"]["variant"]),
                    f"{aggregate['successful_repeat_count']}/{aggregate['repeat_count']}",
                    rate("algorithm_success"),
                    rate("targets_reached"),
                    rate("manufacturing_ready"),
                    median("projected_target_energy"),
                    median("removed_volume_fraction"),
                    median("normalized_rms_displacement"),
                    median("total", timing=True),
                )
            )
            + " |"
        )

    lines.extend(("", "## Paper-reference values", ""))
    references = [
        (case["manifest"], metric)
        for case in cases
        for metric in case["aggregate"]["paper_reference_metrics"]
    ]
    if references:
        lines.extend(
            (
                "These values are transcribed context only (`paper_reference`, `reference_only`). No local-to-paper parity claim is made.",
                "",
                "| Case | Metric | Reference value | Unit | Paper location | Comparability |",
                "|---|---:|---:|---:|---|---|",
            )
        )
        for manifest, metric in references:
            lines.append(
                "| "
                + " | ".join(
                    (
                        _markdown(manifest["experiment_id"]),
                        _markdown(metric["name"]),
                        _format_number(metric["value"]),
                        _markdown(metric["unit"]),
                        _markdown(metric["paper_location"]),
                        "reference_only",
                    )
                )
                + " |"
            )
    else:
        lines.append("No numerical paper-reference metric is declared by these manifests.")
    lines.extend(
        (
            "",
            "## Interpretation boundary",
            "",
            "This report supports implementation regression, ablation, and local robustness analysis only. It does not establish exact reproduction, object-class generality, or physical-print validation.",
            "",
        )
    )
    return "\n".join(lines)


def _summary_svg(cases: Sequence[Mapping[str, Any]]) -> str:
    width = 1120
    row_height = 28
    top = 116
    curve_count = sum(bool(case["convergence"]) for case in cases)
    legend_rows = (curve_count + 1) // 2
    left_x = 30
    label_width = 260
    cell_width = 88
    chart_x = 650
    chart_width = 430
    chart_top = top
    chart_height = max(150, row_height * len(cases))
    height = max(390, chart_top + chart_height + max(72, legend_rows * 17 + 82))
    colors = ("#2563eb", "#dc2626", "#059669", "#9333ea", "#ea580c", "#0891b2", "#4f46e5", "#be123c")

    pieces = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#172033}.title{font-size:21px;font-weight:700}.subtitle{font-size:12px;fill:#526071}.head{font-size:12px;font-weight:700}.label{font-size:11px}.axis{stroke:#9aa6b2;stroke-width:1}.grid{stroke:#d9e0e7;stroke-width:1}.curve{fill:none;stroke-width:2}</style>',
        '<text class="title" x="30" y="34">R10 local synthetic evidence</text>',
        f'<text class="subtitle" x="30" y="56">{xml_escape(REPORT_SCOPE)}</text>',
        '<text class="head" x="30" y="88">Outcome rates across successful executions</text>',
        '<text class="head" x="650" y="88">First successful run: committed-event convergence</text>',
        '<text class="subtitle" x="650" y="105">log10(max(projected target energy, 1e-16)); local measurements only</text>',
        f'<text class="head" x="{left_x}" y="{top - 10}">Case</text>',
    ]
    headings = (("Algorithm", 0), ("Targets", 1), ("Mfg.", 2))
    for label, offset in headings:
        pieces.append(
            f'<text class="head" text-anchor="middle" x="{left_x + label_width + (offset + 0.5) * cell_width}" y="{top - 10}">{label}</text>'
        )

    for index, case in enumerate(cases):
        y = top + index * row_height
        identifier = str(case["manifest"]["experiment_id"])
        display = identifier if len(identifier) <= 38 else identifier[:35] + "..."
        pieces.append(f'<text class="label" x="{left_x}" y="{y + 18}">{xml_escape(display)}</text>')
        for column, metric_name in enumerate(("algorithm_success", "targets_reached", "manufacturing_ready")):
            samples = case["local_values"].get(metric_name, [])
            rate = sum(samples) / len(samples) if samples else None
            fill = "#edf1f5" if rate is None else (
                "#d8f3e8" if rate >= 0.999 else "#fde7c7" if rate > 0.0 else "#f8d7da"
            )
            x = left_x + label_width + column * cell_width
            pieces.append(f'<rect x="{x}" y="{y + 2}" width="{cell_width - 6}" height="22" rx="3" fill="{fill}"/>')
            text = "n/a" if rate is None else f"{rate:.2f}"
            pieces.append(f'<text class="label" text-anchor="middle" x="{x + (cell_width - 6) / 2}" y="{y + 18}">{text}</text>')

    y_min = -16.0
    energies = [energy for case in cases for _, energy in case["convergence"]]
    positive_logs = [math.log10(max(value, 1.0e-16)) for value in energies]
    y_max = max(-1.0, max(positive_logs, default=-1.0))
    maximum_ordinal = max((ordinal for case in cases for ordinal, _ in case["convergence"]), default=1)
    maximum_ordinal = max(1, maximum_ordinal)

    pieces.extend(
        (
            f'<line class="axis" x1="{chart_x}" y1="{chart_top}" x2="{chart_x}" y2="{chart_top + chart_height}"/>',
            f'<line class="axis" x1="{chart_x}" y1="{chart_top + chart_height}" x2="{chart_x + chart_width}" y2="{chart_top + chart_height}"/>',
        )
    )
    for tick in range(-16, int(math.floor(y_max)) + 1, 4):
        y = chart_top + (y_max - tick) / (y_max - y_min) * chart_height
        pieces.append(f'<line class="grid" x1="{chart_x}" y1="{y:.2f}" x2="{chart_x + chart_width}" y2="{y:.2f}"/>')
        pieces.append(f'<text class="label" text-anchor="end" x="{chart_x - 7}" y="{y + 4:.2f}">{tick}</text>')
    pieces.append(f'<text class="label" text-anchor="middle" x="{chart_x + chart_width / 2}" y="{chart_top + chart_height + 28}">Committed event ordinal</text>')

    legend_y = chart_top + chart_height + 48
    legend_x = chart_x
    plotted_cases = tuple(
        (index, case) for index, case in enumerate(cases) if case["convergence"]
    )
    for plotted_index, (index, case) in enumerate(plotted_cases):
        curve = case["convergence"]
        color = colors[index % len(colors)]
        points = []
        for ordinal, energy in curve:
            x = chart_x + ordinal / maximum_ordinal * chart_width
            value = math.log10(max(energy, 1.0e-16))
            y = chart_top + (y_max - value) / (y_max - y_min) * chart_height
            points.append(f"{x:.2f},{y:.2f}")
        pieces.append(f'<polyline class="curve" stroke="{color}" points="{" ".join(points)}"/>')
        lx = legend_x + (plotted_index % 2) * 220
        ly = legend_y + (plotted_index // 2) * 17
        pieces.append(f'<line x1="{lx}" y1="{ly}" x2="{lx + 18}" y2="{ly}" stroke="{color}" stroke-width="2"/>')
        label = str(case["manifest"]["experiment_id"])
        if len(label) > 27:
            label = label[:24] + "..."
        pieces.append(f'<text class="label" x="{lx + 24}" y="{ly + 4}">{xml_escape(label)}</text>')
    pieces.append(f'<text class="subtitle" x="30" y="{height - 22}">Paper-reference values are excluded from this figure. Rates are derived only from succeeded local run records.</text>')
    pieces.append("</svg>")
    return "\n".join(pieces) + "\n"


def _write_text_atomic(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def generate_report(
    evidence_root: os.PathLike[str] | str,
    output_directory: os.PathLike[str] | str,
) -> tuple[str, ...]:
    """Validate *evidence_root* and write deterministic presentation files."""

    cases = load_evidence_tree(evidence_root)
    outputs = {
        "summary.csv": _summary_csv(cases),
        "summary.md": _summary_markdown(cases),
        "summary.svg": _summary_svg(cases),
    }
    root_text = os.fspath(Path(evidence_root).resolve(strict=True))
    for name, content in outputs.items():
        if root_text in content:
            _fail(f"generated {name} contains the absolute evidence-root path")
    destination = Path(output_directory)
    for name, content in outputs.items():
        _write_text_atomic(destination / name, content)
    return tuple(outputs)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    for name in generate_report(arguments.evidence_root, arguments.output_dir):
        print(name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
