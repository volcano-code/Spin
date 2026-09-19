"""Read-only lineage audit for frozen adaptive-source stage snapshots.

The input carrier records authenticated partition snapshots, not a transaction
log for individual split/merge calls.  Consequently this module never
reconstructs or asserts an operation order from a set difference: such a
relationship is explicitly marked ``INFERRED_NET_PARTITION_CHANGE``.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
import math
from typing import Any


class LineageError(ValueError):
    """Malformed, incomplete, or unauthenticated lineage evidence."""


def _fail(message: str) -> None:
    raise LineageError(message)


def cell_key(cell: dict[str, Any]) -> str:
    if not isinstance(cell, dict) or type(cell.get("level")) is not int or cell["level"] < 0:
        _fail("MALFORMED_CELL_KEY")
    coordinates = cell.get("coordinates")
    if (not isinstance(coordinates, list) or len(coordinates) != 3 or
            any(type(value) is not int or value < 0 or value >= 2 ** cell["level"]
                for value in coordinates)):
        _fail("MALFORMED_CELL_KEY")
    return f"{cell['level']}:{coordinates[0]},{coordinates[1]},{coordinates[2]}"


def _cell_from_key(key: str) -> tuple[int, tuple[int, int, int]]:
    if not isinstance(key, str) or key.count(":") != 1:
        _fail("MALFORMED_CELL_KEY")
    level_text, coordinates_text = key.split(":")
    try:
        level = int(level_text)
        coordinates = tuple(int(value) for value in coordinates_text.split(","))
    except ValueError:
        _fail("MALFORMED_CELL_KEY")
    if (level < 0 or len(coordinates) != 3 or any(value < 0 or value >= 2 ** level
                                                   for value in coordinates) or
            key != f"{level}:{coordinates[0]},{coordinates[1]},{coordinates[2]}"):
        _fail("MALFORMED_CELL_KEY")
    return level, coordinates


def _validated_cells(cells: Any) -> dict[str, dict[str, Any]]:
    if not isinstance(cells, list) or not cells:
        _fail("INVALID_PARTITION")
    result: dict[str, dict[str, Any]] = {}
    for cell in cells:
        key = cell_key(cell)
        moments = cell.get("moments") if isinstance(cell, dict) else None
        if (not isinstance(moments, list) or len(moments) != 10 or
                any(isinstance(value, bool) or not isinstance(value, (int, float)) or
                    not math.isfinite(float(value)) for value in moments)):
            _fail("NONFINITE_OR_INVALID_MOMENTS")
        if key in result:
            _fail("DUPLICATE_CELL")
        result[key] = {"level": cell["level"], "coordinates": list(cell["coordinates"]),
                       "moments": [float(value) for value in moments]}
    return result


def partition_hash(cells: list[dict[str, Any]]) -> str:
    """Canonical test/audit digest; source-native digests are retained separately."""
    normalized = _validated_cells(cells)
    payload = [dict(key=key, **normalized[key]) for key in sorted(normalized)]
    data = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(data).hexdigest()


def _hash(value: Any) -> str:
    if not isinstance(value, str) or len(value) != 64:
        _fail("MISMATCHED_PARTITION_HASH")
    try:
        int(value, 16)
    except ValueError:
        _fail("MISMATCHED_PARTITION_HASH")
    return value


def _validate_identity(record: dict[str, Any], field: str, cells: list[dict[str, Any]]) -> dict[str, str]:
    digest = _hash(record.get(f"{field}_hash"))
    canonical = partition_hash(cells)
    scheme = record.get("hash_scheme", "CANONICAL_AUDIT_SHA256")
    if scheme == "CANONICAL_AUDIT_SHA256":
        if digest != canonical:
            _fail("MISMATCHED_PARTITION_HASH")
        verification = "RECOMPUTED_CANONICAL_SHA256"
    elif scheme == "RECORDED_NATIVE_HASH":
        verification = "AUTHENTICATED_RECORD_ONLY_NATIVE_HASH"
    else:
        _fail("UNKNOWN_PARTITION_HASH_SCHEME")
    return {"recorded_hash": digest, "canonical_audit_hash": canonical,
            "verification": verification}


def _beta(values: Any, expected: int) -> list[float]:
    if not isinstance(values, list) or len(values) != expected:
        _fail("INVALID_FILL_VECTOR")
    answer = [float(value) for value in values]
    if any(not math.isfinite(value) or value < 0.0 or value > 1.0 for value in answer):
        _fail("INVALID_FILL_VECTOR")
    return answer


def fractional_counts_by_level(cells: list[dict[str, Any]], beta: list[float], epsilon: float) -> dict[str, dict[str, int]]:
    normalized = _validated_cells(cells)
    values = _beta(beta, len(normalized))
    if not isinstance(epsilon, (int, float)) or not math.isfinite(epsilon) or not 0.0 < epsilon < 0.5:
        _fail("INVALID_FRACTIONAL_EPSILON")
    counts: dict[str, dict[str, int]] = {}
    for cell, value in zip(cells, values, strict=True):
        level = str(cell["level"])
        row = counts.setdefault(level, {"cells": 0, "fractional_cells": 0})
        row["cells"] += 1
        # This is deliberately inclusive, matching the frozen source predicate.
        if epsilon <= value <= 1.0 - epsilon:
            row["fractional_cells"] += 1
    return dict(sorted(counts.items(), key=lambda item: int(item[0])))


def _children_of(parent: str) -> list[str]:
    level, (x, y, z) = _cell_from_key(parent)
    return [f"{level + 1}:{2 * x + dx},{2 * y + dy},{2 * z + dz}"
            for dx in range(2) for dy in range(2) for dz in range(2)]


def _delta(before: dict[str, dict[str, Any]], after: dict[str, dict[str, Any]]) -> list[float]:
    result = []
    for index in range(10):
        result.append(sum(cell["moments"][index] for cell in after.values()) -
                      sum(cell["moments"][index] for cell in before.values()))
    scale = max(1.0, *(abs(value) for cell in before.values() for value in cell["moments"]),
                *(abs(value) for cell in after.values() for value in cell["moments"]))
    if any(abs(value) > 1e-12 * scale for value in result):
        _fail("INCORRECT_MOMENT_CONSERVATION")
    return [0.0 if value == 0.0 else value for value in result]


def _exact_operations(recorded: Any, before: dict[str, dict[str, Any]], after: dict[str, dict[str, Any]],
                      before_beta: list[float], after_beta: list[float]) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    if not isinstance(recorded, dict) or set(recorded) != {"split", "merge"}:
        _fail("INVALID_RECORDED_TRANSITION")
    split, merge = recorded["split"], recorded["merge"]
    if not isinstance(split, list) or not isinstance(merge, list):
        _fail("INVALID_RECORDED_TRANSITION")
    splits: dict[str, list[str]] = {}
    merges: dict[str, list[str]] = {}
    for operation in split:
        if not isinstance(operation, dict) or set(operation) != {"parent", "children"}:
            _fail("INVALID_RECORDED_TRANSITION")
        parent, children = operation["parent"], operation["children"]
        _cell_from_key(parent)
        if (parent not in before or not isinstance(children, list) or children != _children_of(parent) or
                len(set(children)) != 8 or any(child not in after for child in children)):
            _fail("INVALID_RECORDED_TRANSITION")
        splits[parent] = children
    for operation in merge:
        if not isinstance(operation, dict) or set(operation) != {"children", "parent"}:
            _fail("INVALID_RECORDED_TRANSITION")
        children, parent = operation["children"], operation["parent"]
        _cell_from_key(parent)
        if (parent not in after or not isinstance(children, list) or children != _children_of(parent) or
                len(set(children)) != 8 or any(child not in before for child in children)):
            _fail("INVALID_RECORDED_TRANSITION")
        indices = [list(before).index(child) for child in children]
        values = [before_beta[index] for index in indices]
        if any(value not in (0.0, 1.0) for value in values) or len(set(values)) != 1:
            _fail("UNEQUAL_OR_NONBINARY_SIBLING_MERGE")
        merges[parent] = children
    return splits, merges


def _inferred_operations(before: dict[str, dict[str, Any]], after: dict[str, dict[str, Any]]) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    splits = {key: _children_of(key) for key in before
              if all(child in after for child in _children_of(key))}
    merges = {key: _children_of(key) for key in after
              if all(child in before for child in _children_of(key))}
    return splits, merges


def analyze_transition(record: dict[str, Any], *, epsilon: float, before_beta: list[float], after_beta: list[float]) -> dict[str, Any]:
    """Audit one complete boundary without inventing unrecorded operation order."""
    if not isinstance(record, dict):
        _fail("INVALID_STAGE_RECORD")
    if record.get("stage") == 20 or record.get("complete") is False:
        _fail("STAGE_OUT_OF_SCOPE_OR_PARTIAL")
    before_cells = record.get("before")
    after_cells = record.get("after")
    before, after = _validated_cells(before_cells), _validated_cells(after_cells)
    before_values, after_values = _beta(before_beta, len(before)), _beta(after_beta, len(after))
    identities = {"input": _validate_identity(record, "before", before_cells),
                  "output": _validate_identity(record, "after", after_cells)}
    delta = _delta(before, after)
    exact = "recorded_transition" in record
    if exact:
        splits, merges = _exact_operations(record["recorded_transition"], before, after,
                                            before_values, after_values)
        evidence, order = "EXACT_RECORDED", "RECORDED"
    else:
        splits, merges = _inferred_operations(before, after)
        evidence = "INFERRED_NET_PARTITION_CHANGE" if before != after else "NO_PARTITION_CHANGE"
        order = "UNAVAILABLE"
    return {
        "input_partition": identities["input"], "output_partition": identities["output"],
        "transition_evidence": evidence,
        "split": {"parent_to_children": splits, "operation_order": order},
        "merge": {"parent_to_children": merges, "operation_order": order},
        "fractional_counts_by_level": {"input": fractional_counts_by_level(before_cells, before_values, epsilon),
                                        "output": fractional_counts_by_level(after_cells, after_values, epsilon)},
        "moment_conservation_delta": delta,
        "historical_subtimers": "UNAVAILABLE",
    }


def analyze_sequence(rows: list[tuple[dict[str, Any], list[float], list[float]]], *, epsilon: float) -> list[dict[str, Any]]:
    """Audit ordered transitions and surface cells refined then coarsened again."""
    result = []
    events: dict[str, list[str]] = defaultdict(list)
    for record, before_beta, after_beta in rows:
        stage = analyze_transition(record, epsilon=epsilon, before_beta=before_beta, after_beta=after_beta)
        for key in stage["split"]["parent_to_children"]:
            events[key].append("SPLIT")
        for key in stage["merge"]["parent_to_children"]:
            events[key].append("MERGE")
        stage["repeated_lineage"] = []
        result.append(stage)
    repeated = [{"cell": key, "events": value} for key, value in sorted(events.items())
                if "SPLIT" in value and "MERGE" in value]
    if result:
        result[-1]["repeated_lineage"] = repeated
    return result
