"""Console-script entry points for the five R9 Python tools."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Sequence

from .cli_common import (
    cli_error,
    compare_diagnostics,
    load_diagnostics,
    run_batch,
    summarize_diagnostics,
    write_svg_report,
    write_trace_csv,
)
from .io import read_json, write_json_atomic
from .model import load_config, run_config


def _parser(program: str, description: str) -> argparse.ArgumentParser:
    return argparse.ArgumentParser(prog=program, description=description)


def run_balance_main(argv: Sequence[str] | None = None) -> int:
    parser = _parser("run_balance", "Run one validated R9 balance configuration.")
    parser.add_argument("config", type=Path, help="schema-v1 JSON configuration")
    parser.add_argument("--output-directory", type=Path, default=Path("."))
    args = parser.parse_args(argv)
    try:
        result = run_config(load_config(args.config), output_directory=args.output_directory)
    except Exception as error:
        return cli_error(parser.prog, error)
    print(
        json.dumps(
            {
                "status": result.status,
                "succeeded": result.succeeded,
                "manufacturing_ready": result.manufacturing_ready,
            },
            sort_keys=True,
        )
    )
    return 0 if result.succeeded and result.manufacturing_ready else 3


def _load_batch_manifest(path: Path) -> list[Any]:
    value = read_json(path)
    if not isinstance(value, dict):
        raise ValueError("batch manifest must be a JSON object")
    if set(value) != {"schema", "schema_version", "configs"}:
        raise ValueError("batch manifest fields must be schema, schema_version, configs")
    if value["schema"] != "make-it-stand-r9-batch-config-v1":
        raise ValueError("unsupported batch manifest schema")
    if type(value["schema_version"]) is not int or value["schema_version"] != 1:
        raise ValueError("batch manifest schema_version must be integer 1")
    configs = value["configs"]
    if not isinstance(configs, list) or not configs:
        raise ValueError("batch manifest configs must be a non-empty array")
    resolved = []
    for index, item in enumerate(configs):
        if not isinstance(item, str) or not item:
            raise ValueError(f"configs[{index}] must be a non-empty relative path")
        candidate = Path(item)
        if candidate.is_absolute() or ".." in candidate.parts:
            raise ValueError(f"configs[{index}] must not escape the manifest directory")
        resolved.append(load_config(path.parent / candidate))
    return resolved


def batch_experiment_main(argv: Sequence[str] | None = None) -> int:
    parser = _parser(
        "batch_experiment",
        "Run a user-supplied generic R9 batch (not an R10 paper experiment).",
    )
    parser.add_argument("manifest", type=Path, help="R9 batch manifest JSON")
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--keep-going", action="store_true")
    args = parser.parse_args(argv)
    try:
        summary = run_batch(
            _load_batch_manifest(args.manifest),
            output_directory=args.output_directory,
            keep_going=args.keep_going,
        )
    except Exception as error:
        return cli_error(parser.prog, error)
    print(json.dumps(summary, sort_keys=True))
    all_successful = (
        summary["completed_count"] == summary["case_count"]
        and all(
            case.get("succeeded") is True
            and case.get("manufacturing_ready") is True
            for case in summary["cases"]
        )
    )
    return 0 if all_successful else 3


def evaluate_main(argv: Sequence[str] | None = None) -> int:
    parser = _parser("evaluate", "Summarize existing R9 diagnostics artifacts.")
    parser.add_argument("diagnostics", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        summary = summarize_diagnostics([load_diagnostics(path) for path in args.diagnostics])
        write_json_atomic(args.output, summary)
    except Exception as error:
        return cli_error(parser.prog, error)
    return 0


def visualize_main(argv: Sequence[str] | None = None) -> int:
    parser = _parser(
        "visualize",
        "Export an existing diagnostics trace to deterministic CSV and SVG.",
    )
    parser.add_argument("diagnostics", type=Path)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--svg", type=Path, required=True)
    parser.add_argument("--title", default="Make It Stand diagnostics")
    args = parser.parse_args(argv)
    try:
        value = load_diagnostics(args.diagnostics)
        write_trace_csv(args.csv, value)
        write_svg_report(args.svg, value, title=args.title)
    except Exception as error:
        return cli_error(parser.prog, error)
    return 0


def _load_comparison_manifest(path: Path) -> dict[str, list[dict[str, Any]]]:
    value = read_json(path)
    if not isinstance(value, dict):
        raise ValueError("comparison manifest must be a JSON object")
    if set(value) != {"schema", "schema_version", "groups"}:
        raise ValueError("comparison manifest fields must be schema, schema_version, groups")
    if value["schema"] != "make-it-stand-r9-comparison-config-v1":
        raise ValueError("unsupported comparison manifest schema")
    if type(value["schema_version"]) is not int or value["schema_version"] != 1:
        raise ValueError("comparison manifest schema_version must be integer 1")
    groups = value["groups"]
    if not isinstance(groups, dict) or not groups:
        raise ValueError("comparison groups must be a non-empty object")
    loaded: dict[str, list[dict[str, Any]]] = {}
    for name in sorted(groups):
        paths = groups[name]
        if not isinstance(name, str) or not name or not isinstance(paths, list) or not paths:
            raise ValueError("each comparison group must be a named non-empty path array")
        loaded[name] = []
        for index, item in enumerate(paths):
            if not isinstance(item, str) or not item:
                raise ValueError(f"groups.{name}[{index}] must be a relative path")
            candidate = Path(item)
            if candidate.is_absolute() or ".." in candidate.parts:
                raise ValueError(f"groups.{name}[{index}] escapes the manifest directory")
            loaded[name].append(load_diagnostics(path.parent / candidate))
    return loaded


def compare_ablation_main(argv: Sequence[str] | None = None) -> int:
    parser = _parser(
        "compare_ablation",
        "Compare supplied diagnostics groups without running or claiming an R10 ablation.",
    )
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        result = compare_diagnostics(_load_comparison_manifest(args.manifest))
        write_json_atomic(args.output, result)
    except Exception as error:
        return cli_error(parser.prog, error)
    return 0


__all__ = [
    "batch_experiment_main", "compare_ablation_main", "evaluate_main",
    "run_balance_main", "visualize_main",
]
