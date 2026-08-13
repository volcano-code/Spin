#!/usr/bin/env python3
"""Validate and aggregate canonical R10 run records."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from make_it_stand.r10 import (
    R10RecordError,
    aggregate_run_records,
    read_canonical_json,
    validate_manifest,
    validate_run_record,
    write_canonical_json,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate same-manifest R10 local measurements. Paper-reported "
            "values remain labelled reference-only and are never pooled."
        )
    )
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--run-record", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="validate all bindings without writing an aggregate",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    try:
        manifest = validate_manifest(read_canonical_json(arguments.manifest))
        records = [
            validate_run_record(read_canonical_json(path), manifest=manifest)
            for path in arguments.run_record
        ]
        aggregate = aggregate_run_records(manifest, records)
        if not arguments.verify_only:
            write_canonical_json(arguments.output, aggregate)
    except (OSError, ValueError) as error:
        print(f"r10-aggregate: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
