#!/usr/bin/env python3
"""Fail-closed audit for the fixed Spin-It computational requirements.

This M0/M1 foundation validates catalogue completeness, evidence shape, and
asset/parameter readiness.  It intentionally contains no independent
requirement checkers, so neither ready assets nor caller claims can produce a
verification or publication grant.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, NoReturn

from requirements_catalogue import (
    CATALOGUE_SCHEMA,
    CATALOGUE_VERSION,
    REQUIREMENTS,
    REQUIRED_REQUIREMENT_IDS,
    catalogue,
)


EVIDENCE_SCHEMA = "spin-it-strict-evidence-v1"
EVIDENCE_VERSION = 1
AUDIT_SCHEMA = "spin-it-strict-audit-v1"
AUDIT_VERSION = 1

_CATALOGUE_FIELDS = {
    "requirement_id",
    "paper_locator",
    "requirement_kind",
    "expected_outcome",
    "dependencies",
    "required_evidence",
}
_EVIDENCE_FIELDS = {
    "requirement_id",
    "status",
    "evidence_description",
    "details",
}
_ASSET_FIELDS = {
    "asset_id",
    "asset_class",
    "path",
    "sha256",
    "source_reference",
    "required_critical_parameters",
    "parameter_evidence",
}
_PARAMETER_FIELDS = {
    "parameter",
    "status",
    "evidence_description",
    "value",
}
_BLOCKED_REQUIREMENT_STATUSES = {"NOT_VERIFIED", "NOT_IMPLEMENTED", "BLOCKED"}
_ALLOWED_ASSET_CLASSES = {"AUTHOR_ORIGINAL", "PUBLISHED_EXACT_SPECIFICATION"}
_ALLOWED_PARAMETER_STATUSES = {"PUBLISHED", "AUTHOR_PROVIDED", "EXACTLY_DERIVED"}


def _raise_nonfinite(constant: str) -> NoReturn:
    raise ValueError(f"nonfinite JSON value is not allowed: {constant}")


def _reject_nonfinite(value: Any, location: str) -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"nonfinite value at {location}")
    if isinstance(value, dict):
        for key, child in value.items():
            _reject_nonfinite(child, f"{location}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _reject_nonfinite(child, f"{location}[{index}]")


def _mapping(value: Any, location: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"malformed {location}: expected an object")
    return value


def _list(value: Any, location: str) -> list[Any]:
    if not isinstance(value, list):
        raise ValueError(f"malformed {location}: expected an array")
    return value


def _nonempty_text(value: Any, location: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"empty evidence or malformed text at {location}")
    return value


def _exact_fields(row: dict[str, Any], expected: set[str], location: str) -> None:
    missing = sorted(expected - row.keys())
    unknown = sorted(row.keys() - expected)
    if missing or unknown:
        parts = []
        if missing:
            parts.append(f"missing fields {missing}")
        if unknown:
            parts.append(f"unknown fields {unknown}")
        raise ValueError(f"malformed {location}: {'; '.join(parts)}")


def load_catalogue() -> dict[str, Any]:
    """Return and validate a fresh copy of the versioned fixed catalogue."""

    document = catalogue()
    if document.get("schema") != CATALOGUE_SCHEMA:
        raise ValueError("malformed fixed catalogue schema")
    if document.get("schema_version") != CATALOGUE_VERSION:
        raise ValueError("malformed fixed catalogue version")
    validate_catalogue(document.get("requirements"))
    return document


def validate_catalogue(rows: Any) -> None:
    """Reject any incomplete, duplicated, unknown, or altered catalogue."""

    entries = _list(rows, "catalogue requirements")
    ids: list[str] = []
    for index, raw_row in enumerate(entries):
        row = _mapping(raw_row, f"catalogue requirement {index}")
        _exact_fields(row, _CATALOGUE_FIELDS, f"catalogue requirement {index}")
        requirement_id = _nonempty_text(
            row["requirement_id"], f"catalogue requirement {index}.requirement_id"
        )
        ids.append(requirement_id)

    duplicates = sorted({item for item in ids if ids.count(item) > 1})
    if duplicates:
        raise ValueError(f"duplicate catalogue requirement IDs: {duplicates}")

    required_ids = set(REQUIRED_REQUIREMENT_IDS)
    actual_ids = set(ids)
    unknown = sorted(actual_ids - required_ids)
    missing = sorted(required_ids - actual_ids)
    if unknown:
        raise ValueError(f"unknown catalogue requirement IDs: {unknown}; missing: {missing}")
    if missing:
        raise ValueError(f"missing catalogue requirement IDs: {missing}")

    canonical = {row["requirement_id"]: row for row in REQUIREMENTS}
    for row in entries:
        requirement_id = row["requirement_id"]
        if row != canonical[requirement_id]:
            raise ValueError(
                f"catalogue requirement {requirement_id} differs from the fixed catalogue"
            )


def _validated_parameter_evidence(
    asset_id: str,
    declared_parameters: list[str],
    raw_evidence: Any,
) -> list[dict[str, Any]]:
    evidence_rows = _list(raw_evidence, f"asset {asset_id} parameter_evidence")
    evidence: list[dict[str, Any]] = []
    names: list[str] = []
    for index, raw_row in enumerate(evidence_rows):
        row = _mapping(raw_row, f"asset {asset_id} parameter evidence {index}")
        _exact_fields(
            row,
            _PARAMETER_FIELDS,
            f"asset {asset_id} parameter evidence {index}",
        )
        name = _nonempty_text(
            row["parameter"], f"asset {asset_id} parameter evidence {index}.parameter"
        )
        names.append(name)
        status = _nonempty_text(
            row["status"], f"asset {asset_id} parameter {name}.status"
        )
        if status == "UNKNOWN_CRITICAL":
            raise ValueError(f"asset {asset_id} parameter {name} is UNKNOWN_CRITICAL")
        if status not in _ALLOWED_PARAMETER_STATUSES:
            raise ValueError(
                f"asset {asset_id} parameter {name} has unsupported evidence status {status}"
            )
        _nonempty_text(
            row["evidence_description"],
            f"asset {asset_id} parameter {name}.evidence_description",
        )
        _reject_nonfinite(row["value"], f"asset {asset_id} parameter {name}.value")
        evidence.append(dict(row))

    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"duplicate parameter evidence for asset {asset_id}: {duplicates}")
    declared = set(declared_parameters)
    actual = set(names)
    unknown = sorted(actual - declared)
    if unknown:
        raise ValueError(f"unknown parameter evidence for asset {asset_id}: {unknown}")
    missing = sorted(declared - actual)
    if missing:
        raise ValueError(f"missing critical parameter evidence for asset {asset_id}: {missing}")
    return evidence


def preflight_assets(assets: Any, asset_root: Path) -> dict[str, Any]:
    """Validate exact asset provenance, bytes, and critical parameters.

    ``READY`` means only that inputs are ready to be consumed by future
    independent checkers.  It is not a verification result.
    """

    rows = _list(assets, "assets")
    root = Path(asset_root).resolve()
    if not rows:
        return {
            "schema": "spin-it-strict-asset-preflight-v1",
            "schema_version": 1,
            "status": "BLOCKED",
            "ready": False,
            "source_authenticity_verified": False,
            "strict_requirement_verification": False,
            "assets": [],
            "blockers": ["no exact assets were declared"],
        }

    validated: list[dict[str, Any]] = []
    asset_ids: list[str] = []
    for index, raw_asset in enumerate(rows):
        asset = _mapping(raw_asset, f"asset {index}")
        _exact_fields(asset, _ASSET_FIELDS, f"asset {index}")
        asset_id = _nonempty_text(asset["asset_id"], f"asset {index}.asset_id")
        asset_ids.append(asset_id)

        asset_class = _nonempty_text(asset["asset_class"], f"asset {asset_id}.asset_class")
        if asset_class == "RECONSTRUCTED_ASSET":
            raise ValueError(f"reconstructed asset {asset_id} is not admissible")
        if asset_class not in _ALLOWED_ASSET_CLASSES:
            raise ValueError(f"unsupported asset class {asset_class} for asset {asset_id}")

        source_reference = _nonempty_text(
            asset["source_reference"], f"asset {asset_id}.source_reference"
        )
        declared_raw = _list(
            asset["required_critical_parameters"],
            f"asset {asset_id}.required_critical_parameters",
        )
        if not declared_raw:
            raise ValueError(f"asset {asset_id} declares no required critical parameters")
        declared = [
            _nonempty_text(value, f"asset {asset_id}.required_critical_parameters[{index}]")
            for index, value in enumerate(declared_raw)
        ]
        duplicates = sorted({name for name in declared if declared.count(name) > 1})
        if duplicates:
            raise ValueError(f"duplicate critical parameters for asset {asset_id}: {duplicates}")
        _validated_parameter_evidence(asset_id, declared, asset["parameter_evidence"])

        relative_path = Path(_nonempty_text(asset["path"], f"asset {asset_id}.path"))
        if relative_path.is_absolute():
            raise ValueError(f"asset {asset_id} path must be relative to the asset root")
        resolved = (root / relative_path).resolve()
        try:
            resolved.relative_to(root)
        except ValueError as error:
            raise ValueError(f"asset {asset_id} path escapes the asset root") from error
        if not resolved.is_file():
            raise ValueError(f"asset {asset_id} file does not exist: {relative_path}")

        declared_digest = _nonempty_text(asset["sha256"], f"asset {asset_id}.sha256")
        if len(declared_digest) != 64:
            raise ValueError(f"malformed SHA-256 digest for asset {asset_id}")
        try:
            int(declared_digest, 16)
        except ValueError as error:
            raise ValueError(f"malformed SHA-256 digest for asset {asset_id}") from error
        actual_digest = hashlib.sha256(resolved.read_bytes()).hexdigest()
        if actual_digest != declared_digest.lower():
            raise ValueError(
                f"digest mismatch for asset {asset_id}: expected {declared_digest.lower()}, "
                f"got {actual_digest}"
            )

        validated.append(
            {
                "asset_id": asset_id,
                "asset_class": asset_class,
                "path": relative_path.as_posix(),
                "sha256": actual_digest,
                "source_reference": source_reference,
                "critical_parameters": declared,
                "status": "READY",
            }
        )

    duplicate_assets = sorted(
        {asset_id for asset_id in asset_ids if asset_ids.count(asset_id) > 1}
    )
    if duplicate_assets:
        raise ValueError(f"duplicate asset IDs: {duplicate_assets}")

    return {
        "schema": "spin-it-strict-asset-preflight-v1",
        "schema_version": 1,
        "status": "READY",
        "ready": True,
        "source_authenticity_verified": False,
        "strict_requirement_verification": False,
        "assets": validated,
        "blockers": [],
        "readiness_scope": (
            "syntax, declared provenance, and checked file-byte integrity only; "
            "not source authenticity or strict requirement verification"
        ),
    }


def _validate_requirement_evidence(raw_requirements: Any) -> list[dict[str, Any]]:
    rows = _list(raw_requirements, "requirement evidence")
    ids: list[str] = []
    validated: list[dict[str, Any]] = []
    requirement_kinds = {
        row["requirement_id"]: row["requirement_kind"] for row in REQUIREMENTS
    }
    for index, raw_row in enumerate(rows):
        row = _mapping(raw_row, f"requirement evidence {index}")
        _exact_fields(row, _EVIDENCE_FIELDS, f"requirement evidence {index}")
        requirement_id = _nonempty_text(
            row["requirement_id"], f"requirement evidence {index}.requirement_id"
        )
        ids.append(requirement_id)
        if requirement_id not in requirement_kinds:
            continue
        status = _nonempty_text(
            row["status"], f"requirement {requirement_id}.status"
        )
        if status == "VERIFIED":
            raise ValueError(
                f"VERIFIED is not admissible for {requirement_id}: independent checkers are not implemented"
            )
        if requirement_kinds[requirement_id] == "PHYSICAL_MEASUREMENT":
            if status != "NOT_IN_SCOPE":
                raise ValueError(
                    f"physical requirement {requirement_id} must remain NOT_IN_SCOPE"
                )
        elif status not in _BLOCKED_REQUIREMENT_STATUSES:
            raise ValueError(f"unsupported blocked status {status} for {requirement_id}")
        _nonempty_text(
            row["evidence_description"],
            f"requirement {requirement_id}.evidence_description",
        )
        details = _mapping(row["details"], f"requirement {requirement_id}.details")
        _reject_nonfinite(details, f"requirement {requirement_id}.details")
        validated.append(
            {
                "requirement_id": requirement_id,
                "status": status,
                "evidence_description": row["evidence_description"],
                "details": dict(details),
            }
        )

    duplicates = sorted({item for item in ids if ids.count(item) > 1})
    if duplicates:
        raise ValueError(f"duplicate requirement evidence IDs: {duplicates}")
    required = set(REQUIRED_REQUIREMENT_IDS)
    actual = set(ids)
    unknown = sorted(actual - required)
    missing = sorted(required - actual)
    if unknown:
        raise ValueError(f"unknown requirement evidence IDs: {unknown}; missing: {missing}")
    if missing:
        raise ValueError(f"missing requirement evidence IDs: {missing}")
    return validated


def audit(payload: Any, *, asset_root: Path) -> dict[str, Any]:
    """Validate submitted evidence and return a deliberately blocked audit."""

    document = _mapping(payload, "evidence document")
    _reject_nonfinite(document, "evidence")
    _exact_fields(
        document,
        {"schema", "schema_version", "requirements", "assets"},
        "evidence document",
    )
    if document["schema"] != EVIDENCE_SCHEMA:
        raise ValueError(f"unsupported evidence schema: {document['schema']!r}")
    if document["schema_version"] != EVIDENCE_VERSION:
        raise ValueError(f"unsupported evidence schema version: {document['schema_version']!r}")

    catalogue_document = load_catalogue()
    requirements = _validate_requirement_evidence(document["requirements"])
    asset_preflight = preflight_assets(document["assets"], Path(asset_root))

    physical_count = sum(row["status"] == "NOT_IN_SCOPE" for row in requirements)
    blocked_requirements = [
        row for row in requirements if row["status"] != "NOT_IN_SCOPE"
    ]
    blockers = [
        "independent requirement checkers are not implemented in this M0/M1 foundation; input readiness cannot grant verification",
        "formula coverage remains blocked until each formula requirement has an independent checker result",
        "experiment coverage remains blocked until each original experiment has independently checked evidence",
    ]
    blockers.extend(
        f"{row['requirement_id']}: {row['status']} — {row['evidence_description']}"
        for row in blocked_requirements
    )
    blockers.extend(asset_preflight["blockers"])

    return {
        "schema": AUDIT_SCHEMA,
        "schema_version": AUDIT_VERSION,
        "catalogue_schema": catalogue_document["schema"],
        "catalogue_version": catalogue_document["schema_version"],
        "status": "PARTIAL_BLOCKED",
        "public_complete_reproduction": False,
        "formula_coverage": False,
        "experiment_coverage": False,
        "asset_preflight": asset_preflight,
        "requirements": requirements,
        "requirement_summary": {
            "required": len(REQUIRED_REQUIREMENT_IDS),
            "blocked_or_not_implemented": len(blocked_requirements),
            "physical_measurements_out_of_scope": physical_count,
            "verified": 0,
        },
        "blockers": blockers,
        "limitations": [
            "This is only the M0/M1 strict-audit foundation, not full M1 completion.",
            "No independent formula, experiment, Boolean, simulation, or visualization checker is implemented here.",
            "Physical fabricated-object measurements are retained explicitly and remain out of scope.",
            "READY asset preflight covers syntax, declared provenance, and checked file-byte integrity only; it does not prove source authenticity, strict verification, or publication readiness.",
        ],
    }


def _load_json(path: Path) -> Any:
    return json.loads(
        path.read_text(encoding="utf-8"),
        parse_constant=_raise_nonfinite,
    )


def _write_json(path: Path, document: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate strict Spin-It evidence without granting publication claims."
    )
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--asset-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    try:
        result = audit(_load_json(args.evidence), asset_root=args.asset_root)
    except (OSError, json.JSONDecodeError, ValueError) as error:
        result = {
            "schema": AUDIT_SCHEMA,
            "schema_version": AUDIT_VERSION,
            "status": "INVALID_BLOCKED",
            "public_complete_reproduction": False,
            "formula_coverage": False,
            "experiment_coverage": False,
            "blockers": [str(error)],
        }
    _write_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
