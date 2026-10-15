"""Read-only lineage and moment-conservation audit for sealed adaptive snapshots."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


MOMENT_COUNT = 10


def is_fractional(beta: float, epsilon: float) -> bool:
    """Match the source solver's inclusive fractional interval exactly."""
    return epsilon <= beta <= 1.0 - epsilon


def _binary(beta: float, epsilon: float) -> bool:
    return beta < epsilon or beta > 1.0 - epsilon


def _key(leaf: dict[str, Any]) -> tuple[int, int, int, int]:
    level = leaf.get("level")
    coordinates = leaf.get("coordinates")
    if type(level) is not int or level < 0 or not isinstance(coordinates, list) or len(coordinates) != 3:
        raise ValueError("invalid leaf key")
    if any(type(value) is not int or value < 0 or value >= (1 << level) for value in coordinates):
        raise ValueError("invalid leaf coordinates")
    return level, coordinates[0], coordinates[1], coordinates[2]


def _label(key: tuple[int, int, int, int]) -> str:
    level, x, y, z = key
    return f"{level}:{x},{y},{z}"


def _parent(key: tuple[int, int, int, int]) -> tuple[int, int, int, int] | None:
    level, x, y, z = key
    if level == 0:
        return None
    return level - 1, x // 2, y // 2, z // 2


def _ancestors(key: tuple[int, int, int, int]):
    current = _parent(key)
    while current is not None:
        yield current
        current = _parent(current)


def _vector(leaf: dict[str, Any]) -> tuple[float, ...]:
    values = leaf.get("moments")
    if not isinstance(values, list) or len(values) != MOMENT_COUNT:
        raise ValueError("leaf must contain ten moments")
    out = tuple(float(value) for value in values)
    if not all(math.isfinite(value) for value in out):
        raise ValueError("nonfinite leaf moment")
    return out


def _snapshot(snapshot: dict[str, Any], epsilon: float) -> dict[str, Any]:
    stage = snapshot.get("stage")
    leaves = snapshot.get("leaves")
    fills = snapshot.get("fill")
    if type(stage) is not int or stage < 0 or not isinstance(leaves, list) or not isinstance(fills, list):
        raise ValueError("invalid snapshot")
    if len(leaves) != len(fills):
        raise ValueError("fill/partition length mismatch")
    keys = [_key(leaf) for leaf in leaves]
    if len(set(keys)) != len(keys):
        raise ValueError("duplicate partition leaf")
    key_set = set(keys)
    for key in keys:
        if any(parent in key_set for parent in _ancestors(key)):
            raise ValueError("overlapping partition")
    vectors = [_vector(leaf) for leaf in leaves]
    beta = [float(value) for value in fills]
    if any(not math.isfinite(value) or value < 0.0 or value > 1.0 for value in beta):
        raise ValueError("fill outside [0,1] or nonfinite")
    moment_by_key = dict(zip(keys, vectors))
    beta_by_key = dict(zip(keys, beta))
    level_counts = Counter(key[0] for key in keys)
    fractional_by_level = Counter(key[0] for key, value in beta_by_key.items() if is_fractional(value, epsilon))
    total_moments = [math.fsum(vector[i] for vector in vectors) for i in range(MOMENT_COUNT)]
    weighted_moments = [math.fsum(vector[i] * beta[j] for j, vector in enumerate(vectors))
                        for i in range(MOMENT_COUNT)]
    return {
        "stage": stage,
        "keys": keys,
        "key_set": key_set,
        "moment_by_key": moment_by_key,
        "beta_by_key": beta_by_key,
        "total_moments": total_moments,
        "weighted_moments": weighted_moments,
        "partition_hash": snapshot.get("partition_hash"),
        "source_sha256": snapshot.get("source_sha256"),
        "leaf_count": len(keys),
        "fractional_count": sum(is_fractional(value, epsilon) for value in beta),
        "level_counts": {str(level): level_counts[level] for level in sorted(level_counts)},
        "fractional_by_level": {str(level): fractional_by_level[level] for level in sorted(level_counts)},
    }


def _linf_delta(left: list[float], right: list[float]) -> tuple[float, float]:
    absolute = max((abs(a - b) for a, b in zip(left, right)), default=0.0)
    scale = max(max((abs(value) for value in left), default=0.0),
                max((abs(value) for value in right), default=0.0), 1e-30)
    return absolute, absolute / scale


def analyze_snapshots(snapshots: list[dict[str, Any]], epsilon: float,
                      moment_tolerance: float = 1e-8) -> dict[str, Any]:
    if not math.isfinite(epsilon) or not 0.0 < epsilon < 0.5:
        raise ValueError("invalid fractional epsilon")
    if not math.isfinite(moment_tolerance) or moment_tolerance < 0.0:
        raise ValueError("invalid moment tolerance")
    states = [_snapshot(snapshot, epsilon) for snapshot in snapshots]
    if not states:
        raise ValueError("no snapshots")
    for left, right in zip(states, states[1:]):
        if right["stage"] != left["stage"] + 1:
            raise ValueError("snapshot stages must be consecutive")

    split_events: Counter[str] = Counter()
    merge_events: Counter[str] = Counter()
    split_seen: set[tuple[int, int, int, int]] = set()
    merge_seen: set[tuple[int, int, int, int]] = set()
    merge_then_resplit: set[tuple[int, int, int, int]] = set()
    refractionalized: set[tuple[int, int, int, int]] = set()
    transitions = []

    for before, after in zip(states, states[1:]):
        after_descendants: dict[tuple[int, int, int, int], list[tuple[int, int, int, int]]] = defaultdict(list)
        for key in after["keys"]:
            after_descendants[key].append(key)
            for ancestor in _ancestors(key):
                after_descendants[ancestor].append(key)

        before_descendants: dict[tuple[int, int, int, int], list[tuple[int, int, int, int]]] = defaultdict(list)
        for key in before["keys"]:
            before_descendants[key].append(key)
            for ancestor in _ancestors(key):
                before_descendants[ancestor].append(key)

        split_parents = []
        mapped_before_beta: dict[tuple[int, int, int, int], float] = {}
        for key in before["keys"]:
            if key in after["key_set"]:
                mapped_before_beta[key] = after["beta_by_key"][key]
                continue
            ancestor = next((item for item in _ancestors(key) if item in after["key_set"]), None)
            if ancestor is not None:
                mapped_before_beta[key] = after["beta_by_key"][ancestor]
                continue
            descendants = after_descendants.get(key, [])
            if not descendants:
                raise ValueError("partition lineage has an uncovered cell: " + _label(key))
            values = [after["beta_by_key"][item] for item in descendants]
            if max(values) - min(values) > 1e-12:
                raise ValueError("split descendants do not inherit one fill value: " + _label(key))
            mapped_before_beta[key] = values[0]
            split_parents.append(key)

        merge_parents = []
        for key in after["keys"]:
            if key in before["key_set"]:
                continue
            descendants = before_descendants.get(key, [])
            if len(descendants) > 1:
                merge_parents.append(key)

        for key in split_parents:
            label = _label(key)
            split_events[label] += 1
            if key in merge_seen:
                merge_then_resplit.add(key)
            split_seen.add(key)
        for key in merge_parents:
            merge_events[_label(key)] += 1
            merge_seen.add(key)

        refractionalized_this_transition = set()
        for key in before["key_set"] & after["key_set"]:
            old = before["beta_by_key"][key]
            new = after["beta_by_key"][key]
            if _binary(old, epsilon) and is_fractional(new, epsilon):
                refractionalized.add(key)
                refractionalized_this_transition.add(key)

        event_moment_errors = []
        for key in split_parents:
            descendants = after_descendants[key]
            pre = [before["moment_by_key"][key][i] * mapped_before_beta[key]
                   for i in range(MOMENT_COUNT)]
            post = [math.fsum(after["moment_by_key"][item][i] * after["beta_by_key"][item]
                              for item in descendants) for i in range(MOMENT_COUNT)]
            absolute, relative = _linf_delta(pre, post)
            event_moment_errors.append({"kind": "split", "cell": _label(key),
                                        "linf_absolute_delta": absolute,
                                        "linf_relative_delta": relative})
        for key in merge_parents:
            descendants = before_descendants[key]
            beta = after["beta_by_key"][key]
            pre = [math.fsum(before["moment_by_key"][item][i] * beta for item in descendants)
                   for i in range(MOMENT_COUNT)]
            post = [after["moment_by_key"][key][i] * beta for i in range(MOMENT_COUNT)]
            absolute, relative = _linf_delta(pre, post)
            event_moment_errors.append({"kind": "merge", "cell": _label(key),
                                        "linf_absolute_delta": absolute,
                                        "linf_relative_delta": relative})

        preweighted = [math.fsum(before["moment_by_key"][key][i] * mapped_before_beta[key]
                                  for key in before["keys"]) for i in range(MOMENT_COUNT)]
        before_geometric_abs, before_geometric_rel = _linf_delta(before["total_moments"], after["total_moments"])
        weighted_abs, weighted_rel = _linf_delta(preweighted, after["weighted_moments"])
        event_relative = max((item["linf_relative_delta"] for item in event_moment_errors), default=0.0)
        transitions.append({
            "from_stage": before["stage"], "to_stage": after["stage"],
            "before_leaf_count": before["leaf_count"], "after_leaf_count": after["leaf_count"],
            "split_count": len(split_parents), "merge_count": len(merge_parents),
            "split_parents": [_label(key) for key in split_parents],
            "merged_parents": [_label(key) for key in merge_parents],
            "refractionalized_cells": sorted(_label(key) for key in refractionalized_this_transition),
            "before_partition_moments": before["total_moments"],
            "after_partition_moments": after["total_moments"],
            "weighted_moments_before_transition": preweighted,
            "weighted_moments_after_transition": after["weighted_moments"],
            "partition_moment_linf_absolute_delta": before_geometric_abs,
            "partition_moment_linf_relative_delta": before_geometric_rel,
            "weighted_moment_linf_absolute_delta": weighted_abs,
            "weighted_moment_linf_relative_delta": weighted_rel,
            "transition_moment_event_count": len(event_moment_errors),
            "transition_moment_event_max_relative_delta": event_relative,
            "transition_moment_event_worst": max(
                event_moment_errors, key=lambda item: item["linf_relative_delta"], default=None),
            "moment_check_passed": max(before_geometric_rel, weighted_rel, event_relative) <= moment_tolerance,
        })

    return {
        "schema": "spin2014-adaptive-partition-lineage-audit-v1",
        "scope": "READ_ONLY_RECONSTRUCTION_FROM_PUBLISHED_CHECKPOINT_HEADERS",
        "fractional_rule": "epsilon <= beta <= 1-epsilon",
        "epsilon": epsilon,
        "moment_tolerance_relative_linf": moment_tolerance,
        "snapshots": [{key: state[key] for key in (
            "stage", "leaf_count", "fractional_count", "level_counts", "fractional_by_level",
            "partition_hash", "source_sha256")} for state in states],
        "transitions": transitions,
        "genealogy": {
            "split_events_by_cell": dict(sorted(split_events.items())),
            "merge_events_by_cell": dict(sorted(merge_events.items())),
            "repeatedly_split_cells": sorted(key for key, count in split_events.items() if count > 1),
            "merge_then_resplit_cells": sorted(_label(key) for key in merge_then_resplit),
            "refractionalized_cells": sorted(_label(key) for key in refractionalized),
        },
        "all_transition_moment_checks_passed": all(item["moment_check_passed"] for item in transitions),
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_attempt_snapshots(batch: Path, first_stage: int, last_stage: int, repeat: int = 0):
    if first_stage < 1 or last_stage < first_stage:
        raise ValueError("first stage needs its predecessor snapshot")
    snapshots = []
    for stage in range(first_stage - 1, last_stage + 1):
        path = batch / f"repeat-{repeat}/payload/attempts/attempt-{stage:08d}/payload/header.json"
        if not path.is_file():
            raise FileNotFoundError(path)
        header = json.loads(path.read_text())
        if len(header.get("trace_counts", [])) != stage + 1:
            raise ValueError(f"attempt {stage} does not contain the expected cumulative pass count")
        solve = header.get("payload", {}).get("solve", {})
        leaves = solve.get("final_partition")
        fill = solve.get("relaxed_fill")
        if not isinstance(leaves, list) or not isinstance(fill, list):
            raise ValueError(f"attempt {stage} has no published final partition/fill")
        snapshots.append({"stage": stage, "leaves": leaves, "fill": fill,
                          "partition_hash": solve.get("final_partition_hash"),
                          "source_sha256": _sha256(path)})
    return snapshots


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", type=Path, required=True)
    parser.add_argument("--first-stage", type=int, required=True)
    parser.add_argument("--last-stage", type=int, required=True)
    parser.add_argument("--repeat", type=int, default=0)
    parser.add_argument("--epsilon", type=float, required=True)
    parser.add_argument("--moment-tolerance", type=float, default=1e-8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    snapshots = load_attempt_snapshots(args.batch, args.first_stage, args.last_stage, args.repeat)
    report = analyze_snapshots(snapshots, args.epsilon, args.moment_tolerance)
    report["source_attempt_range"] = [args.first_stage - 1, args.last_stage]
    rendered = json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
