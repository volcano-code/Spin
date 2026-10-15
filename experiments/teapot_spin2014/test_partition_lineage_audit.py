import pytest

from partition_lineage_audit import analyze_snapshots, is_fractional


def leaf(level, coordinates, moments):
    return {"level": level, "coordinates": list(coordinates), "moments": list(moments)}


def test_fractional_interval_includes_both_epsilon_endpoints():
    assert is_fractional(0.05, 0.05)
    assert is_fractional(0.95, 0.05)
    assert not is_fractional(0.049999, 0.05)
    assert not is_fractional(0.950001, 0.05)


def test_split_then_merge_preserves_moments_and_tracks_repeated_cell():
    parent_moments = [float(i + 1) for i in range(10)]
    child_moments = [[value / 8.0 for value in parent_moments] for _ in range(8)]
    parent = leaf(0, (0, 0, 0), parent_moments)
    children = [leaf(1, (x, y, z), child_moments[x * 4 + y * 2 + z])
                for x in range(2) for y in range(2) for z in range(2)]

    report = analyze_snapshots([
        {"stage": 0, "leaves": [parent], "fill": [0.5]},
        {"stage": 1, "leaves": children, "fill": [0.5] * 8},
        {"stage": 2, "leaves": [parent], "fill": [0.0]},
        {"stage": 3, "leaves": [parent], "fill": [0.05]},
    ], epsilon=0.05)

    assert [item["split_count"] for item in report["transitions"]] == [1, 0, 0]
    assert [item["merge_count"] for item in report["transitions"]] == [0, 1, 0]
    assert report["genealogy"]["split_events_by_cell"]["0:0,0,0"] == 1
    assert report["genealogy"]["merge_events_by_cell"]["0:0,0,0"] == 1
    assert report["genealogy"]["refractionalized_cells"] == ["0:0,0,0"]
    assert [len(item["refractionalized_cells"]) for item in report["transitions"]] == [0, 0, 1]
    assert all(item["weighted_moment_linf_relative_delta"] <= 1e-12 for item in report["transitions"])
    assert all(item["transition_moment_event_max_relative_delta"] <= 1e-12 for item in report["transitions"])


def test_rejects_non_antichain_or_wrong_fill_cardinality():
    parent = leaf(0, (0, 0, 0), [1.0] * 10)
    child = leaf(1, (0, 0, 0), [0.125] * 10)
    with pytest.raises(ValueError, match="overlapping partition"):
        analyze_snapshots([{"stage": 0, "leaves": [parent, child], "fill": [0.5, 0.5]}], epsilon=0.05)
    with pytest.raises(ValueError, match="fill/partition length"):
        analyze_snapshots([{"stage": 0, "leaves": [parent], "fill": []}], epsilon=0.05)


def test_per_cell_moment_check_catches_cancelling_partition_errors():
    parents = [leaf(1, (x, 0, 0), [1.0] * 10) for x in (0, 1)]
    children = []
    for parent_index, parent in enumerate(parents):
        for xbit in range(2):
            for ybit in range(2):
                for zbit in range(2):
                    moments = [0.125] * 10
                    moments[0] += (0.0125 if parent_index == 0 else -0.0125)
                    children.append(leaf(2, (parent["coordinates"][0] * 2 + xbit, ybit, zbit), moments))
    report = analyze_snapshots([
        {"stage": 0, "leaves": parents, "fill": [0.5, 0.5]},
        {"stage": 1, "leaves": children, "fill": [0.5] * 16},
    ], epsilon=0.05)
    transition = report["transitions"][0]
    assert transition["partition_moment_linf_relative_delta"] < 1e-12
    assert transition["transition_moment_event_max_relative_delta"] > 1e-3
    assert not transition["moment_check_passed"]
