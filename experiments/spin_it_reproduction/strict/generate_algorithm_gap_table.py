#!/usr/bin/env python3
"""Generate the current algorithm/case gap table from a strict audit.

The table is deliberately derived from the fail-closed audit rather than from
hand-maintained booleans.  It separates current verified method evidence from
unverified method coverage, original-case/asset blockers, and out-of-scope
physical rows.  The result is a planning artifact, not a completion claim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def _bucket(requirement_id: str, status: str) -> str:
    if status == "VERIFIED":
        return "verified_current_evidence"
    if status == "NOT_IN_SCOPE":
        return "postponed_out_of_scope"

    rid = requirement_id.lower()
    case_tokens = (
        ".case.",
        ".table1",
        ".figure",
        ".figures",
        ".original_",
        ".author_",
        ".cad.",
        ".determinism",
        ".simulation",
        ".physical",
    )
    if any(token in rid for token in case_tokens):
        return "asset_or_original_case_blocked"
    return "implemented_or_covered_but_unverified"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    audit_bytes = args.audit.read_bytes()
    audit = json.loads(audit_bytes)
    rows = audit.get("requirements", [])
    buckets: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        rid = str(row.get("requirement_id", ""))
        status = str(row.get("status", ""))
        bucket = _bucket(rid, status)
        buckets.setdefault(bucket, []).append(
            {
                "requirement_id": rid,
                "status": status,
                "evidence_description": row.get("evidence_description"),
            }
        )

    summary = {
        "schema": "spin-algorithm-gap-table-v1",
        "source_audit": str(args.audit.resolve()),
        "source_audit_sha256": hashlib.sha256(audit_bytes).hexdigest(),
        "project_status": audit.get("status"),
        "public_complete_reproduction": audit.get("public_complete_reproduction"),
        "counts": {
            key: len(value) for key, value in sorted(buckets.items())
        },
        "buckets": {key: buckets[key] for key in sorted(buckets)},
        "interpretation": {
            "verified_current_evidence": "Current, independently checked evidence; may be synthetic or method-level.",
            "implemented_or_covered_but_unverified": "Requires an independent check, integration run, or current evidence refresh.",
            "asset_or_original_case_blocked": "Original input, critical parameter, case matrix, or downstream evidence is missing or blocked.",
            "postponed_out_of_scope": "Explicitly outside this computational phase; never counted as passed.",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(args.output), "counts": summary["counts"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
