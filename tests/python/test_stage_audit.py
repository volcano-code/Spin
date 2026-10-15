"""Stage diagnostics are algebra checks, not solver or geometry admission."""
import copy
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "teapot_spin2014"))
sys.path.insert(0, str(ROOT / "python"))


def _cell(level, coordinates, moments):
    return {"level": level, "coordinates": coordinates, "moments": moments}


def fixture():
    # Deliberately asymmetric, hand-checkable moment stencils.
    outer = [20.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]
    a = [2.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
    b = [1.0, 4.0, 3.0, 2.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    a0 = [v / 2.0 for v in a]
    a1 = [v / 2.0 for v in a]
    before = [_cell(1, [0, 0, 0], a), _cell(2, [1, 0, 0], b)]
    after = [
        _cell(2, [0, 0, 0], a0),
        _cell(2, [1, 0, 0], a1),
        _cell(2, [2, 0, 0], b),
    ]
    removed = [1.25, 3.25, 2.75, 2.25, 1.75, 2.75, 3.75, 4.75, 5.75, 6.75]
    material = [37.5, 13.5, 34.5, 55.5, 76.5, 94.5, 112.5, 130.5, 148.5, 166.5]
    phi = math.pi / 6.0
    before_hash, after_hash = "b" * 64, "c" * 64
    kkt = {
        "valid": True,
        "status": "WITNESS_LP_OPTIMAL_FIRST_ORDER",
        "full": 5e-8,
        "equality_primal": 4e-9,
        "stationarity": 3e-8,
        "complementarity": 2e-8,
        "dual_sign": 1e-9,
        "bound_primal": 0.0,
        "lambda": [1.0, 2.0, 3.0, 4.0, 5.0],
        "alpha": [0.0, 0.0],
        "eta": [0.0, 0.0],
    }
    pass_payload = {
        "index": 0,
        "status": "CONVERGED",
        "evaluation_state": "AVAILABLE",
        "before": before,
        "after": after,
        "before_hash": before_hash,
        "after_hash": after_hash,
        "kkt_partition_hash": before_hash,
        "beta": [0.25, 0.75],
        "output_beta": [0.25, 0.25, 0.75],
        "final_state": [0.25, 0.75, phi],
        "final_phi": phi,
        "fractional_count": 2,
        "output_fractional_count": 3,
        "level_counts": [0, 1, 1],
        "fractional_level_counts": [0, 1, 1],
        "output_level_counts": [0, 0, 3],
        "output_fractional_level_counts": [0, 0, 3],
        "paper_resolution_complete": True,
        "final_kkt": kkt,
        "evaluation": {
            "final_moments": material,
            "constraints": [13.5, 34.5, 148.5, 130.5, 75.4557713659401],
        },
        "split": {
            "before_hash": before_hash,
            "after_hash": after_hash,
            "removed_before": removed,
            "removed_after": removed,
            "changes": [{"parent": {"level": 1, "coordinates": [0, 0, 0]}}],
        },
        "merge": {
            "before_hash": after_hash,
            "after_hash": after_hash,
            "removed_before": removed,
            "removed_after": removed,
            "changes": [],
        },
    }
    terminal = {
        "final_state": True,
        "state": [0.25, 0.75, phi],
        "phi": phi,
        "constraints": pass_payload["evaluation"]["constraints"],
        "kkt": kkt,
    }
    binary_moments = [38.0, 12.0, 34.0, 56.0, 78.0, 96.0, 114.0, 132.0, 150.0, 168.0]
    header = {
        "schema": "source-spin2014-native-material-v1",
        "output_digest": "d" * 64,
        "settings": {
            "density": 2.0,
            "source_settings": {"root_side": 2.0},
            "solver_settings": {"objective": "TOP_INERTIA_ONLY", "fractional_epsilon": 0.05},
        },
        "source": {
            "outer_surface_moments": outer,
            "source_hash": "e" * 64,
            "config_hash": "f" * 64,
            "partition_hash": "1" * 64,
        },
        "solve": {
            "status": "BINARY_VALIDATION_FAILED",
            "diagnostic": "fixture rounding failure",
            "final_phi": phi,
            "final_partition": after,
            "final_partition_hash": after_hash,
            "relaxed_fill": [0.25, 0.25, 0.75],
            "binary_fill": [0.0, 0.0, 1.0],
            "binary_moments": binary_moments,
            "rounding_delta": [-0.25, -0.25, 0.25],
            "terminal_relaxation_kkt_applies": False,
            "input_hash": "2" * 64,
            "output_hash": "3" * 64,
            "trace_hash": "4" * 64,
            "binary_check": {
                "passed": False,
                "constraints": [12.0, 34.0, 150.0, 132.0, 76.2057713659401],
                "normalized_constraint_residual": 2.34375,
            },
        },
        "geometry": {"state": "UNAVAILABLE", "status": "NOT_ATTEMPTED_UNKNOWN"},
    }
    return header, pass_payload, terminal


def failed_pass_fixture():
    header, stage, terminal = fixture()
    stage["status"] = "NUMERICAL_UNKNOWN"
    stage["paper_resolution_complete"] = False
    stage["after"] = copy.deepcopy(stage["before"])
    stage["after_hash"] = stage["before_hash"]
    stage["output_beta"] = list(stage["beta"])
    stage["output_level_counts"] = list(stage["level_counts"])
    stage["output_fractional_level_counts"] = list(stage["fractional_level_counts"])
    stage["output_fractional_count"] = stage["fractional_count"]
    for name in ("split", "merge"):
        stage[name] = dict(before_hash="", after_hash="", changes=[],
                           removed_before=[0.0] * 10, removed_after=[0.0] * 10)
    header["solve"].update(status="RELAXATION_UNKNOWN", final_partition=stage["after"],
                           final_partition_hash=stage["after_hash"], relaxed_fill=stage["output_beta"],
                           binary_fill=[], binary_moments=None, binary_check=None, rounding_delta=[])
    return header, stage, terminal


def test_failed_solver_without_adaptation_is_auditable_but_not_admitted():
    from stage_audit import audit_components

    header, stage, terminal = failed_pass_fixture()
    result = audit_components(header, [(stage, terminal)])
    item = result["passes"][0]
    assert item["solver_status"] == "NUMERICAL_UNKNOWN"
    assert item["transition"]["state"] == "NOT_EXECUTED_SOLVER_NOT_CONVERGED"
    assert item["split_merge_moment_delta"] == [0.0] * 10
    assert item["kkt_evidence"]["terminal_partition"] == "NOT_APPLICABLE_SOLVER_NOT_CONVERGED"
    assert result["terminal_rounding"]["available"] is False
    assert result["geometry_export_stage"] == "NOT_AVAILABLE"
    assert "passed" not in result


@pytest.mark.parametrize("mutation", ["partition", "beta", "hash", "changes", "moments",
                                     "partial_hash", "converged", "resolution_complete"])
def test_failed_no_transition_rejects_hidden_changes_or_false_success(mutation):
    from stage_audit import pass_diagnostics

    header, stage, terminal = failed_pass_fixture()
    if mutation == "partition":
        stage["after"][0]["coordinates"][0] = 9
    elif mutation == "beta":
        stage["output_beta"][0] += 1e-13
    elif mutation == "hash":
        stage["after_hash"] = "a" * 64
    elif mutation == "changes":
        stage["split"]["changes"] = [{"parent": {"level": 1}}]
    elif mutation == "moments":
        stage["merge"]["removed_after"][0] = 1e-13
    elif mutation == "partial_hash":
        stage["split"]["before_hash"] = stage["before_hash"]
    elif mutation == "converged":
        stage["status"] = "CONVERGED"
    else:
        stage["paper_resolution_complete"] = True
    with pytest.raises(ValueError):
        pass_diagnostics(header, stage, terminal)


def test_asymmetric_nonzero_phi_density_and_scale_are_hand_derived():
    from stage_audit import pass_diagnostics

    header, stage, terminal = fixture()
    result = pass_diagnostics(header, stage, terminal)
    assert result["pre_transition"]["material_moments"] == pytest.approx(
        [37.5, 13.5, 34.5, 55.5, 76.5, 94.5, 112.5, 130.5, 148.5, 166.5]
    )
    assert result["pre_transition"]["constraints_input_normalized"] == pytest.approx(
        [0.421875, 1.078125, 2.3203125, 2.0390625, 1.178996427592814]
    )


def test_split_preserves_moments_and_reports_actual_fractional_volumes_by_level():
    from stage_audit import pass_diagnostics

    header, stage, terminal = fixture()
    result = pass_diagnostics(header, stage, terminal)
    assert result["split_merge_moment_delta"] == pytest.approx([0.0] * 10, abs=1e-14)
    pre = result["pre_transition"]["fractional_distribution"]
    assert pre["fractional_cell_count"] == 2
    assert pre["fractional_cell_volume"] == pytest.approx(3.0)
    assert pre["fractional_removed_volume"] == pytest.approx(1.25)
    assert pre["by_level"] == {
        "1": {"cells": 1, "fractional_cells": 1, "cell_volume": 2.0,
              "fractional_cell_volume": 2.0, "fractional_removed_volume": 0.5},
        "2": {"cells": 1, "fractional_cells": 1, "cell_volume": 1.0,
              "fractional_cell_volume": 1.0, "fractional_removed_volume": 0.75},
    }


def test_fractional_epsilon_boundaries_are_inclusive_and_empty_case_is_explicit():
    from stage_audit import _fractional_distribution

    moments = __import__("numpy").asarray([
        [2.0] + [0.0] * 9,
        [3.0] + [0.0] * 9,
    ])
    levels = __import__("numpy").asarray([1, 2])
    boundary, mask = _fractional_distribution(
        __import__("numpy").asarray([0.05, 0.95]), moments, levels, 0.05
    )
    assert mask.tolist() == [True, True]
    assert boundary["fractional_cell_count"] == 2
    empty, mask = _fractional_distribution(
        __import__("numpy").asarray([0.0, 1.0]), moments, levels, 0.05
    )
    assert mask.tolist() == [False, False]
    assert empty["fractional_cell_count"] == 0
    assert empty["fractional_cell_volume"] == 0.0


def test_zero_fractional_epsilon_is_rejected_like_native_solver():
    from stage_audit import pass_diagnostics

    header, stage, terminal = fixture()
    header["settings"]["solver_settings"]["fractional_epsilon"] = 0.0
    with pytest.raises(ValueError, match="INVALID_FRACTIONAL_EPSILON"):
        pass_diagnostics(header, stage, terminal)


def test_partition_change_limits_kkt_and_final_rounding_is_existing_comparison():
    from stage_audit import audit_components

    header, stage, terminal = fixture()
    result = audit_components(header, [(stage, terminal)])
    item = result["passes"][0]
    assert item["kkt_evidence"]["pre_transition"] == "AVAILABLE_FIRST_ORDER_ONLY"
    assert item["kkt_evidence"]["terminal_partition"] == "NOT_APPLICABLE_PARTITION_CHANGED"
    assert set(item["kkt_evidence"]["pre_transition_record"]) == {
        "valid", "status", "full", "equality_primal", "stationarity",
        "complementarity", "dual_sign", "bound_primal", "lambda",
    }
    assert result["terminal_rounding"]["raw_moment_delta"] == pytest.approx(
        [0.5, -1.5, -0.5, 0.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5]
    )
    assert result["terminal_rounding"]["existing_binary_check"] == header["solve"]["binary_check"]
    assert result["geometry_export_stage"] == "NOT_AVAILABLE"
    assert "passed" not in result and "success" not in result


@pytest.mark.parametrize("mutation, reason", [
    ("invalid_kkt", "NOT_APPLICABLE_INVALID_PRE_TRANSITION_KKT"),
    ("changed_beta", "NOT_APPLICABLE_BETA_CHANGED"),
])
def test_same_partition_kkt_requires_valid_converged_witness_and_unchanged_beta(mutation, reason):
    from stage_audit import audit_components, pass_diagnostics

    header, stage, terminal = fixture()
    stage["after"] = copy.deepcopy(stage["before"])
    stage["after_hash"] = stage["before_hash"]
    stage["output_beta"] = list(stage["beta"])
    stage["output_level_counts"] = list(stage["level_counts"])
    stage["output_fractional_level_counts"] = list(stage["fractional_level_counts"])
    stage["output_fractional_count"] = stage["fractional_count"]
    stage["split"]["after_hash"] = stage["before_hash"]
    stage["merge"]["before_hash"] = stage["before_hash"]
    stage["merge"]["after_hash"] = stage["before_hash"]
    if mutation == "invalid_kkt":
        stage["final_kkt"]["valid"] = False
        terminal["kkt"]["valid"] = False
    else:
        stage["output_beta"] = [0.2500000000001, 0.75]
        removed = [0.2500000000001 * a + 0.75 * b for a, b in zip(
            stage["before"][0]["moments"], stage["before"][1]["moments"]
        )]
        stage["split"]["removed_after"] = removed
        stage["merge"]["removed_before"] = removed
        stage["merge"]["removed_after"] = removed
    diagnostic = pass_diagnostics(header, stage, terminal)
    assert diagnostic["kkt_evidence"]["terminal_partition"] == reason

    header["solve"]["final_partition"] = copy.deepcopy(stage["after"])
    header["solve"]["final_partition_hash"] = stage["after_hash"]
    header["solve"]["relaxed_fill"] = list(stage["output_beta"])
    header["solve"]["terminal_relaxation_kkt_applies"] = False
    header["solve"]["binary_fill"] = [0.0, 1.0]
    header["solve"]["rounding_delta"] = [
        -stage["output_beta"][0], 1.0 - stage["output_beta"][1]
    ]
    if mutation == "changed_beta":
        assert stage["output_beta"][0] - stage["beta"][0] == pytest.approx(1e-13)
    result = audit_components(header, [(stage, terminal)])
    assert result["passes"][0]["kkt_evidence"]["terminal_partition"] == reason


@pytest.mark.parametrize("mutation", ["bad_moments", "bad_state", "nan", "wrong_terminal_beta"])
def test_bad_dimensions_nonfinite_and_inconsistent_terminal_state_are_rejected(mutation):
    from stage_audit import pass_diagnostics

    header, stage, terminal = fixture()
    if mutation == "bad_moments":
        stage["before"][0]["moments"].pop()
    elif mutation == "bad_state":
        terminal["state"].pop()
    elif mutation == "nan":
        stage["beta"][0] = float("nan")
    else:
        terminal["state"][0] = 0.5
    with pytest.raises(ValueError):
        pass_diagnostics(header, stage, terminal)


def _events():
    header, stage, terminal = fixture()
    yield {"kind": "header", "trace_counts": [1], "payload": header}
    yield {"kind": "pass", "pass_index": 0, "payload": stage}
    yield {"kind": "iteration", "pass_index": 0, "trace_index": 0, "payload": terminal}


def test_archive_must_be_complete_and_cli_refuses_overwrite(tmp_path):
    from make_it_stand.verification.source_trace_archive import write_archive
    from stage_audit import audit_archive

    archive = tmp_path / "trace.jsonl"
    output = tmp_path / "stage.json"
    write_archive(archive, _events())
    data = archive.read_bytes()
    archive.write_bytes(data.rsplit(b"\n", 2)[0] + b"\n")
    with pytest.raises(ValueError):
        audit_archive(archive)

    archive.unlink()
    write_archive(archive, _events())
    command = [sys.executable, str(ROOT / "experiments" / "teapot_spin2014" / "stage_audit.py"),
               str(archive), str(output)]
    first = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    assert first.returncode == 0, first.stderr
    report = json.loads(output.read_text())
    assert report["input_archive_sha256"] == hashlib.sha256(archive.read_bytes()).hexdigest()
    second = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    assert second.returncode != 0
    assert "FileExistsError" in second.stderr


def test_header_terminal_partition_must_match_last_transition():
    from stage_audit import audit_components

    header, stage, terminal = fixture()
    header["solve"]["relaxed_fill"][0] = 0.5
    with pytest.raises(ValueError, match="TERMINAL_PARTITION_MISMATCH"):
        audit_components(header, [(stage, terminal)])
