#!/usr/bin/env python3
"""Compatibility-only builder for the superseded legacy release gate.

The strict requirements audit in ``strict/audit.py`` is authoritative.  This
module remains importable for older callers, but caller-provided aggregate
booleans can no longer grant any strict reproduction claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


REQUIRED_CLAIMS = {
    "formula_coverage", "experiment_coverage", "asset_coverage",
    "author_exact", "physical_validation", "public_complete_reproduction",
}


def _all_cases_pass(payload: dict[str, Any]) -> bool:
    cases = payload.get("cases")
    return bool(cases) and all(row.get("passed") is True for row in cases)


def _figures_pass(payload: dict[str, Any]) -> bool:
    try:
        eleven = payload["figure11"]["result"]["front_audit"]
        twelve = payload["figure12"]["result"]["front_audit"]
        thirteen = payload["figure13"]["result"]["rows"]
    except (KeyError, TypeError):
        return False
    return (
        eleven.get("finite") is True
        and eleven.get("nondominated") is True
        and twelve.get("finite") is True
        and twelve.get("nondominated") is True
        and bool(thirteen)
        and all(row.get("status") == "CONVERGED" for row in thirteen)
    )


def aggregate(
    faster: dict[str, Any],
    spin2014: dict[str, Any],
    figures: dict[str, Any],
    *,
    tests_passed: bool,
) -> dict[str, Any]:
    faster_passed = _all_cases_pass(faster)
    spin2014_passed = _all_cases_pass(spin2014)
    paper_resolution = spin2014.get("paper_resolution_reproduced") is True
    paper_numerics = spin2014.get("paper_numerics_reproduced") is True
    cage_deformation = spin2014.get("cage_deformation_executed") is True
    figures_passed = _figures_pass(figures)
    # Retain the legacy observations for diagnostics only.  They are not
    # independent requirement evidence and must never be promoted to claims.
    formula_coverage = False
    experiment_coverage = False
    gate = {
        "schema": "spin-it-legacy-release-gate-v2",
        "schema_version": 2,
        "status": "SUPERSEDED_BLOCKED",
        "formula_coverage": formula_coverage,
        "experiment_coverage": experiment_coverage,
        "asset_coverage": "RECONSTRUCTED_COMPLETE",
        "author_exact": False,
        "physical_validation": False,
        "public_complete_reproduction": False,
        "requirements": {
            "caller_reported_tests_passed": bool(tests_passed),
            "spin_it_faster_cases_passed": faster_passed,
            "spin_it_2014_cases_passed": spin2014_passed,
            "spin_it_2014_paper_resolution": paper_resolution,
            "spin_it_2014_paper_numerics": paper_numerics,
            "spin_it_2014_cage_deformation": cage_deformation,
            "spin_it_faster_figures_passed": figures_passed,
            "author_code_and_cad_public": False,
            "fabricated_object_validation_available": False,
        },
        "limitations": [
            "superseded by the fixed strict requirements audit",
            "legacy aggregate inputs are untrusted diagnostics and cannot grant formula or experiment coverage",
            "author code and author CAD/voxel assets were not publicly released",
            "paper examples use deterministic reconstructed structural counterparts",
            "Spin-It 2014 level-9 octree resolution and cage deformation must be executed, not inferred",
            "no fabricated-object spin or balance measurements are claimed",
        ],
    }
    validate_gate(gate)
    encoded = json.dumps(gate, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    gate["scientific_hash"] = hashlib.sha256(encoded).hexdigest()
    return gate


def validate_gate(gate: dict[str, Any]) -> None:
    missing = REQUIRED_CLAIMS - gate.keys()
    if missing:
        raise ValueError(f"release gate lacks claim dimensions: {sorted(missing)}")
    for key in ("formula_coverage", "experiment_coverage", "author_exact", "physical_validation", "public_complete_reproduction"):
        if not isinstance(gate[key], bool):
            raise ValueError(f"{key} must be boolean")
    if gate["asset_coverage"] not in {"NONE", "PARTIAL", "RECONSTRUCTED_COMPLETE", "AUTHOR_EXACT"}:
        raise ValueError("unsupported asset coverage")
    if gate["public_complete_reproduction"] and not (
        gate["formula_coverage"] and gate["experiment_coverage"]
    ):
        raise ValueError("public-complete claim exceeds formula/experiment evidence")
    if gate["author_exact"] and gate["asset_coverage"] != "AUTHOR_EXACT":
        raise ValueError("author-exact requires author-exact assets")
    if gate.get("schema") == "spin-it-legacy-release-gate-v2" and any(
        gate[key]
        for key in (
            "formula_coverage",
            "experiment_coverage",
            "public_complete_reproduction",
        )
    ):
        raise ValueError("superseded legacy gates cannot grant strict claims")


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def report(gate: dict[str, Any]) -> str:
    verdict = "PASS" if gate["public_complete_reproduction"] else "BLOCKED"
    return "\n".join((
        "# Spin-It superseded legacy compatibility report",
        "",
        f"Computational public-complete gate: **{verdict}**.",
        "",
        f"- Formula coverage: `{str(gate['formula_coverage']).lower()}`",
        f"- Experiment coverage: `{str(gate['experiment_coverage']).lower()}`",
        f"- Asset coverage: `{gate['asset_coverage']}`",
        "- Author-exact: `false`",
        "- Physical validation: `false`",
        "",
        "This compatibility-only aggregate records legacy diagnostics but cannot grant strict reproduction claims.",
        "Use strict/audit.py for the fixed requirements inventory and fail-closed input preflight.",
        "Neither path currently proves source authenticity, independent verification,",
        "or the papers' fabricated-object measurements.",
        "",
    ))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--faster", type=Path, required=True)
    parser.add_argument("--spin2014", type=Path, required=True)
    parser.add_argument("--figures", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    gate = aggregate(
        _load(args.faster), _load(args.spin2014), _load(args.figures),
        tests_passed=False,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(gate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(report(gate), encoding="utf-8")
    return 0 if gate["public_complete_reproduction"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
