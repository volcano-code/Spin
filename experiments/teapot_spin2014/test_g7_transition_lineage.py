import copy
import hashlib
import json

import pytest

from g7_transition_lineage import (
    LineageError,
    analyze_sequence,
    analyze_transition,
    fractional_counts_by_level,
    partition_hash,
)


def cell(level, coordinates, *, mass=1.0):
    return {"level": level, "coordinates": list(coordinates),
            "moments": [mass] + [0.0] * 9}


def children(parent, *, mass=0.125):
    level = parent["level"] + 1
    x, y, z = parent["coordinates"]
    return [cell(level, (2 * x + dx, 2 * y + dy, 2 * z + dz), mass=mass)
            for dx in range(2) for dy in range(2) for dz in range(2)]


def exact(before, after, split=None, merge=None, *, stage=1, complete=True):
    return {
        "stage": stage,
        "complete": complete,
        "before": before,
        "after": after,
        "before_hash": partition_hash(before),
        "after_hash": partition_hash(after),
        "recorded_transition": {"split": split or [], "merge": merge or []},
    }


def test_fractional_classification_is_inclusive_at_source_boundaries():
    cells = [cell(2, (i, 0, 0)) for i in range(4)]
    result = fractional_counts_by_level(cells, [0.05, 0.95, 0.049999, 0.950001], 0.05)
    assert result == {"2": {"cells": 4, "fractional_cells": 2}}


def test_exact_parent_to_eight_children_split_conserves_all_ten_moments():
    parent = cell(2, (1, 2, 3))
    after = children(parent)
    record = exact([parent], after, split=[{"parent": "2:1,2,3",
                                             "children": [f"3:{c['coordinates'][0]},{c['coordinates'][1]},{c['coordinates'][2]}"
                                                          for c in after]}])
    result = analyze_transition(record, epsilon=0.05, before_beta=[0.5], after_beta=[0.5] * 8)
    assert result["transition_evidence"] == "EXACT_RECORDED"
    assert result["split"]["parent_to_children"]["2:1,2,3"] == [
        f"3:{c['coordinates'][0]},{c['coordinates'][1]},{c['coordinates'][2]}" for c in after]
    assert result["moment_conservation_delta"] == [0.0] * 10


def test_equal_binary_siblings_merge_and_repeated_ancestry_are_recorded():
    parent = cell(1, (0, 0, 0))
    leaves = children(parent)
    split = exact([parent], leaves, split=[{"parent": "1:0,0,0",
                                             "children": [f"2:{c['coordinates'][0]},{c['coordinates'][1]},{c['coordinates'][2]}"
                                                          for c in leaves]}])
    merge = exact(leaves, [parent], merge=[{"children": [f"2:{c['coordinates'][0]},{c['coordinates'][1]},{c['coordinates'][2]}"
                                                           for c in leaves], "parent": "1:0,0,0"}], stage=2)
    result = analyze_sequence([(split, [1.0], [1.0] * 8), (merge, [1.0] * 8, [1.0])], epsilon=0.05)
    assert result[1]["repeated_lineage"] == [{"cell": "1:0,0,0", "events": [
        {"event": "SPLIT", "evidence": "EXACT_TRANSACTION_EVIDENCE"},
        {"event": "MERGE", "evidence": "EXACT_TRANSACTION_EVIDENCE"},
    ]}]


@pytest.mark.parametrize("mutator", [
    lambda r: r["before"].__setitem__(0, cell(1, (0, 0, 0))),
    lambda r: r["before"].append(copy.deepcopy(r["before"][0])),
    lambda r: r.__setitem__("before_hash", "0" * 64),
    lambda r: r["before"][0]["moments"].__setitem__(0, float("nan")),
    lambda r: r["after"][0]["moments"].__setitem__(0, 2.0),
])
def test_rejects_malformed_or_unauthenticated_partition_evidence(mutator):
    record = exact([cell(0, (0, 0, 0))], [cell(0, (0, 0, 0))])
    mutator(record)
    with pytest.raises(LineageError):
        analyze_transition(record, epsilon=0.05, before_beta=[0.5], after_beta=[0.5])


def test_unrecorded_partition_difference_is_inferred_not_exact_and_history_is_unavailable():
    parent = cell(0, (0, 0, 0))
    result = analyze_transition({key: value for key, value in exact([parent], children(parent)).items()
                                 if key != "recorded_transition"}, epsilon=0.05,
                                before_beta=[0.5], after_beta=[0.5] * 8)
    assert result["transition_evidence"] == "INFERRED_NET_PARTITION_CHANGE"
    assert result["split"]["operation_order"] == "UNAVAILABLE"
    assert result["historical_subtimers"] == "UNAVAILABLE"


@pytest.mark.parametrize("stage, complete, remove", [
    (None, True, "stage"),
    (True, True, None),
    (0, True, None),
    (20, True, None),
    (21, True, None),
    (1, False, None),
    (1, None, "complete"),
])
def test_rejects_missing_or_out_of_scope_stage_metadata(stage, complete, remove):
    record = exact([cell(0, (0, 0, 0))], [cell(0, (0, 0, 0))])
    if remove:
        record.pop(remove)
    else:
        record["stage"] = stage
        record["complete"] = complete
    with pytest.raises(LineageError, match="STAGE_OUT_OF_SCOPE_OR_PARTIAL"):
        analyze_transition(record, epsilon=0.05, before_beta=[0.5], after_beta=[0.5])


def test_rejects_syntactically_valid_unverified_native_hash_scheme():
    record = exact([cell(0, (0, 0, 0))], [cell(0, (0, 0, 0))])
    record["hash_scheme"] = "RECORDED_NATIVE_HASH"
    record["before_hash"] = "a" * 64
    record["after_hash"] = "b" * 64
    with pytest.raises(LineageError, match="UNVERIFIED_NATIVE_HASH_SCHEME"):
        analyze_transition(record, epsilon=0.05, before_beta=[0.5], after_beta=[0.5])
