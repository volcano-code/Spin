"""Bounded staged diagnostics for source-trace-archive-v1 evidence.

This module checks recorded algebra and provenance.  It does not admit a
solver result, certify geometry, or establish a local/global minimum.
"""
import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
PYTHON = ROOT / "python"
if str(PYTHON) not in sys.path:
    sys.path.insert(0, str(PYTHON))

from make_it_stand.verification.source_trace_archive import read_archive
from resolution_rounding_audit import rounding_diagnostics


def _array(value, shape, reason):
    result = np.asarray(value, dtype=float)
    if result.shape != shape:
        raise ValueError(reason)
    if not np.isfinite(result).all():
        raise ValueError("NONFINITE_STAGE_STATE")
    return result


def _scalar(value, reason, *, positive=False):
    if isinstance(value, bool):
        raise ValueError(reason)
    result = float(value)
    if not math.isfinite(result) or (positive and result <= 0.0):
        raise ValueError(reason)
    return result


def _same(left, right, reason, *, atol=1e-12):
    if not np.allclose(left, right, rtol=1e-11, atol=atol):
        raise ValueError(reason)


def _cells(value, reason):
    if not isinstance(value, list) or not value:
        raise ValueError(reason)
    moments = np.empty((len(value), 10), dtype=float)
    levels = np.empty(len(value), dtype=int)
    for index, cell in enumerate(value):
        if not isinstance(cell, dict) or type(cell.get("level")) is not int or cell["level"] < 0:
            raise ValueError(reason)
        moments[index] = _array(cell.get("moments"), (10,), reason)
        levels[index] = cell["level"]
    if np.any(moments[:, 0] < 0.0):
        raise ValueError("NEGATIVE_CELL_VOLUME")
    return moments, levels


def _transform(phi):
    c, s = math.cos(phi), math.sin(phi)
    transform = np.zeros((5, 10), dtype=float)
    transform[0, 1] = transform[1, 2] = transform[2, 8] = transform[3, 7] = 1.0
    transform[4, 4] = c * s
    transform[4, 5] = -c * s
    transform[4, 9] = c * c - s * s
    return transform


def _constraints(moments, phi, density, side):
    raw = _transform(phi) @ moments
    scale = density * np.asarray([side ** 4, side ** 4, side ** 5, side ** 5, side ** 5])
    if not np.isfinite(raw).all() or not np.isfinite(scale).all() or np.any(scale <= 0.0):
        raise ValueError("INVALID_NORMALIZATION")
    return raw, raw / scale


def _counts(levels):
    size = int(levels.max()) + 1
    return np.bincount(levels, minlength=size).astype(int)


def _recorded_counts(recorded, computed, reason):
    if not isinstance(recorded, list) or any(type(value) is not int or value < 0 for value in recorded):
        raise ValueError(reason)
    padded = np.pad(computed, (0, max(0, len(recorded) - len(computed))))
    if len(recorded) < len(computed) or not np.array_equal(padded, np.asarray(recorded, dtype=int)):
        raise ValueError(reason)


def _fractional_distribution(beta, moments, levels, epsilon):
    fractional = (beta >= epsilon) & (beta <= 1.0 - epsilon)
    by_level = {}
    for level in sorted(set(levels.tolist())):
        mask = levels == level
        fm = mask & fractional
        by_level[str(level)] = {
            "cells": int(np.count_nonzero(mask)),
            "fractional_cells": int(np.count_nonzero(fm)),
            "cell_volume": float(moments[mask, 0].sum()),
            "fractional_cell_volume": float(moments[fm, 0].sum()),
            "fractional_removed_volume": float((beta[fm] * moments[fm, 0]).sum()),
        }
    return {
        "fractional_epsilon": epsilon,
        "fractional_cell_count": int(np.count_nonzero(fractional)),
        "fractional_cell_volume": float(moments[fractional, 0].sum()),
        "fractional_removed_volume": float((beta[fractional] * moments[fractional, 0]).sum()),
        "by_level": by_level,
    }, fractional


def _partition(header, cells, beta, levels, phi, density, side, epsilon):
    outer = _array(header["source"].get("outer_surface_moments"), (10,), "BAD_OUTER_MOMENTS")
    removed = beta @ cells
    material = density * (outer - removed)
    if not np.isfinite(material).all():
        raise ValueError("NONFINITE_STAGE_STATE")
    raw, normalized = _constraints(material, phi, density, side)
    distribution, fractional = _fractional_distribution(beta, cells, levels, epsilon)
    return {
        "cell_count": len(beta),
        "removed_moments": removed.tolist(),
        "material_moments": material.tolist(),
        "constraints_at_actual_phi": raw.tolist(),
        "constraints_input_normalized": normalized.tolist(),
        "normalized_constraint_residual": float(np.max(np.abs(normalized))),
        "fractional_distribution": distribution,
    }, removed, material, fractional


def pass_diagnostics(header, stage, terminal_iteration):
    """Diagnose one pass from its authoritative pass state and final marker."""
    if not isinstance(header, dict) or not isinstance(stage, dict) or not isinstance(terminal_iteration, dict):
        raise ValueError("INVALID_STAGE_RECORD")
    before, before_levels = _cells(stage.get("before"), "BAD_BEFORE_PARTITION")
    after, after_levels = _cells(stage.get("after"), "BAD_AFTER_PARTITION")
    beta = _array(stage.get("beta"), (len(before),), "BAD_BETA_DIMENSIONS")
    output_beta = _array(stage.get("output_beta"), (len(after),), "BAD_OUTPUT_BETA_DIMENSIONS")
    if np.any((beta < 0.0) | (beta > 1.0)) or np.any((output_beta < 0.0) | (output_beta > 1.0)):
        raise ValueError("INVALID_FILL")

    density = _scalar(header["settings"].get("density"), "INVALID_NORMALIZATION", positive=True)
    side = _scalar(header["settings"]["source_settings"].get("root_side"),
                   "INVALID_NORMALIZATION", positive=True)
    epsilon = _scalar(header["settings"]["solver_settings"].get("fractional_epsilon"),
                      "INVALID_FRACTIONAL_EPSILON")
    if epsilon <= 0.0 or epsilon >= 0.5:
        raise ValueError("INVALID_FRACTIONAL_EPSILON")
    phi = _scalar(stage.get("final_phi"), "INVALID_FINAL_PHI")

    if terminal_iteration.get("final_state") is not True:
        raise ValueError("MISSING_ACCEPTED_TERMINAL_STATE")
    terminal_state = _array(terminal_iteration.get("state"), (len(before) + 1,), "BAD_TERMINAL_DIMENSIONS")
    pass_state = _array(stage.get("final_state"), (len(before) + 1,), "BAD_FINAL_STATE_DIMENSIONS")
    _same(terminal_state, pass_state, "TERMINAL_STATE_MISMATCH")
    _same(terminal_state[:-1], beta, "TERMINAL_BETA_MISMATCH")
    terminal_phi = _scalar(terminal_iteration.get("phi"), "INVALID_TERMINAL_PHI")
    if not math.isclose(terminal_state[-1], phi, rel_tol=1e-12, abs_tol=1e-12) or not math.isclose(
            terminal_phi, phi, rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("TERMINAL_PHI_MISMATCH")

    pre, removed_before, material_before, fractional_before = _partition(
        header, before, beta, before_levels, phi, density, side, epsilon)
    post, removed_after, material_after, fractional_after = _partition(
        header, after, output_beta, after_levels, phi, density, side, epsilon)

    evaluation = stage.get("evaluation")
    if stage.get("evaluation_state") == "AVAILABLE":
        if not isinstance(evaluation, dict):
            raise ValueError("MISSING_PASS_EVALUATION")
        _same(_array(evaluation.get("final_moments"), (10,), "BAD_EVALUATION_MOMENTS"),
              material_before, "EVALUATION_MOMENT_MISMATCH")
        _same(_array(evaluation.get("constraints"), (5,), "BAD_EVALUATION_CONSTRAINTS"),
              np.asarray(pre["constraints_at_actual_phi"]), "EVALUATION_CONSTRAINT_MISMATCH")

    terminal_constraints = _array(terminal_iteration.get("constraints"), (5,), "BAD_TERMINAL_CONSTRAINTS")
    _same(terminal_constraints, np.asarray(pre["constraints_at_actual_phi"]),
          "TERMINAL_CONSTRAINT_MISMATCH")

    split, merge = stage.get("split"), stage.get("merge")
    if not isinstance(split, dict) or not isinstance(merge, dict):
        raise ValueError("MISSING_TRANSITION_METADATA")
    before_hash, after_hash = stage.get("before_hash"), stage.get("after_hash")
    _digest(before_hash, "INVALID_PARTITION_DIGEST")
    _digest(after_hash, "INVALID_PARTITION_DIGEST")
    if stage.get("kkt_partition_hash") != before_hash:
        raise ValueError("TRANSITION_HASH_MISMATCH")
    transition_state = "EXECUTED"
    if stage.get("status") != "CONVERGED":
        # The native solver exits before adaptation on an incomplete relaxation.
        # Empty transition hashes and zero moments are unexecuted placeholders,
        # not measured removals. Verify exact immutability; never invent a chain.
        if (before_hash != after_hash or stage["before"] != stage["after"] or
                not np.array_equal(beta, output_beta) or
                stage.get("paper_resolution_complete") is not False):
            raise ValueError("FAILED_STAGE_CHANGED_PARTITION_OR_FILL")
        for operation in (split, merge):
            if (operation.get("before_hash") != "" or operation.get("after_hash") != "" or
                    operation.get("changes") != [] or
                    np.any(_array(operation.get("removed_before"), (10,), "BAD_UNUSED_TRANSITION") != 0.0) or
                    np.any(_array(operation.get("removed_after"), (10,), "BAD_UNUSED_TRANSITION") != 0.0)):
                raise ValueError("FAILED_STAGE_HAS_TRANSITION")
        transition_state = "NOT_EXECUTED_SOLVER_NOT_CONVERGED"
    else:
        if (split.get("before_hash") != before_hash or split.get("after_hash") != merge.get("before_hash") or
                merge.get("after_hash") != after_hash):
            raise ValueError("TRANSITION_HASH_MISMATCH")
        _same(_array(split.get("removed_before"), (10,), "BAD_SPLIT_MOMENTS"), removed_before,
              "SPLIT_INPUT_MISMATCH")
        split_after = _array(split.get("removed_after"), (10,), "BAD_SPLIT_MOMENTS")
        merge_before = _array(merge.get("removed_before"), (10,), "BAD_MERGE_MOMENTS")
        _same(split_after, merge_before, "SPLIT_MERGE_CHAIN_MISMATCH")
        _same(_array(merge.get("removed_after"), (10,), "BAD_MERGE_MOMENTS"), removed_after,
              "MERGE_OUTPUT_MISMATCH")

    _recorded_counts(stage.get("level_counts"), _counts(before_levels), "LEVEL_COUNT_MISMATCH")
    _recorded_counts(stage.get("output_level_counts"), _counts(after_levels), "LEVEL_COUNT_MISMATCH")
    _recorded_counts(stage.get("fractional_level_counts"), _counts(before_levels[fractional_before])
                     if np.any(fractional_before) else np.zeros(0, dtype=int), "FRACTIONAL_COUNT_MISMATCH")
    _recorded_counts(stage.get("output_fractional_level_counts"), _counts(after_levels[fractional_after])
                     if np.any(fractional_after) else np.zeros(0, dtype=int), "FRACTIONAL_COUNT_MISMATCH")
    if (stage.get("fractional_count") != int(np.count_nonzero(fractional_before)) or
            stage.get("output_fractional_count") != int(np.count_nonzero(fractional_after))):
        raise ValueError("FRACTIONAL_COUNT_MISMATCH")

    final_kkt = stage.get("final_kkt")
    terminal_kkt = terminal_iteration.get("kkt")
    if final_kkt != terminal_kkt:
        raise ValueError("TERMINAL_KKT_MISMATCH")
    if isinstance(final_kkt, dict) and final_kkt.get("valid") is True:
        pre_kkt = "AVAILABLE_FIRST_ORDER_ONLY"
    else:
        pre_kkt = "UNAVAILABLE"
    terminal_applicability = _terminal_kkt_applicability(stage, beta, output_beta)
    kkt_summary = _kkt_summary(final_kkt)
    raw_delta, normalized_delta = _constraints(material_after - material_before, phi, density, side)
    return {
        "pass_index": stage.get("index"),
        "solver_status": stage.get("status", "UNKNOWN"),
        "actual_final_phi": phi,
        "terminal_beta_source": "RECORDED_FINAL_ITERATION_ACCEPTED_TERMINAL_STATE",
        "pre_transition": pre,
        "post_split_merge": post,
        "split_merge_moment_delta": (material_after - material_before).tolist(),
        "split_merge_constraint_delta": raw_delta.tolist(),
        "split_merge_constraint_delta_input_normalized": normalized_delta.tolist(),
        "transition": {
            "state": transition_state,
            "before_partition_hash": before_hash,
            "after_partition_hash": after_hash,
            "split_changes": len(split.get("changes", [])),
            "merge_changes": len(merge.get("changes", [])),
            "paper_resolution_complete_recorded": stage.get("paper_resolution_complete"),
        },
        "kkt_evidence": {
            "pre_transition": pre_kkt,
            "pre_transition_record": kkt_summary,
            "terminal_partition": terminal_applicability,
        },
    }


def _kkt_summary(kkt):
    if not isinstance(kkt, dict) or type(kkt.get("valid")) is not bool or not isinstance(kkt.get("status"), str):
        raise ValueError("INVALID_KKT_RECORD")
    summary = {"valid": kkt["valid"], "status": kkt["status"]}
    for key in ("full", "equality_primal", "stationarity", "complementarity", "dual_sign", "bound_primal"):
        summary[key] = _scalar(kkt.get(key), "INVALID_KKT_RECORD")
    summary["lambda"] = _array(kkt.get("lambda"), (5,), "INVALID_KKT_RECORD").tolist()
    return summary


def _terminal_kkt_applicability(stage, beta=None, output_beta=None):
    if stage.get("status") != "CONVERGED":
        return "NOT_APPLICABLE_SOLVER_NOT_CONVERGED"
    kkt = stage.get("final_kkt")
    if not isinstance(kkt, dict) or kkt.get("valid") is not True:
        return "NOT_APPLICABLE_INVALID_PRE_TRANSITION_KKT"
    if stage.get("before_hash") != stage.get("after_hash"):
        return "NOT_APPLICABLE_PARTITION_CHANGED"
    if beta is None:
        beta = np.asarray(stage.get("beta"), dtype=float)
    if output_beta is None:
        output_beta = np.asarray(stage.get("output_beta"), dtype=float)
    if beta.shape != output_beta.shape or not np.array_equal(beta, output_beta):
        return "NOT_APPLICABLE_BETA_CHANGED"
    return "APPLICABLE_UNCHANGED_TERMINAL_STATE"


def _digest(value, reason):
    if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise ValueError(reason)
    return value


def _source_digests(header):
    source, solve = header.get("source"), header.get("solve")
    if not isinstance(source, dict) or not isinstance(solve, dict):
        raise ValueError("INVALID_HEADER_PAYLOAD")
    fields = {
        "native_output_digest": header.get("output_digest"),
        "source_hash": source.get("source_hash"),
        "source_config_hash": source.get("config_hash"),
        "source_partition_hash": source.get("partition_hash"),
        "solve_input_hash": solve.get("input_hash"),
        "solve_output_hash": solve.get("output_hash"),
        "solve_trace_hash": solve.get("trace_hash"),
    }
    return {key: _digest(value, "INVALID_SOURCE_DIGEST") for key, value in fields.items()}


def _rounding(header):
    result = rounding_diagnostics(header)
    result["existing_binary_check"] = header["solve"].get("binary_check")
    if not result.get("available"):
        return result
    solve = header["solve"]
    relaxed = _array(solve.get("relaxed_fill"), (result["total_cells"],), "BAD_ROUNDING_DIMENSIONS")
    binary = _array(solve.get("binary_fill"), (result["total_cells"],), "BAD_ROUNDING_DIMENSIONS")
    _same(_array(solve.get("rounding_delta"), (result["total_cells"],), "BAD_ROUNDING_DIMENSIONS"),
          binary - relaxed, "ROUNDING_DELTA_MISMATCH")
    check = solve.get("binary_check")
    if not isinstance(check, dict):
        raise ValueError("MISSING_BINARY_CHECK")
    density = _scalar(header["settings"].get("density"), "INVALID_NORMALIZATION", positive=True)
    side = _scalar(header["settings"]["source_settings"].get("root_side"),
                   "INVALID_NORMALIZATION", positive=True)
    phi = _scalar(solve.get("final_phi"), "INVALID_FINAL_PHI")
    raw, normalized = _constraints(np.asarray(result["binary_moments"]), phi, density, side)
    _same(_array(check.get("constraints"), (5,), "BAD_BINARY_CHECK"), raw,
          "BINARY_CHECK_CONSTRAINT_MISMATCH")
    residual = _scalar(check.get("normalized_constraint_residual"), "BAD_BINARY_CHECK")
    if not math.isclose(residual, float(np.max(np.abs(normalized))), rel_tol=1e-10, abs_tol=1e-12):
        raise ValueError("BINARY_CHECK_RESIDUAL_MISMATCH")
    return result


def audit_components(header, stages):
    """Audit an iterable of (pass payload, terminal iteration payload) pairs."""
    if not isinstance(header, dict) or header.get("schema") != "source-spin2014-native-material-v1":
        raise ValueError("INVALID_HEADER_PAYLOAD")
    source_digests = _source_digests(header)
    diagnostics = []
    previous = None
    last_stage = None
    for expected_index, pair in enumerate(stages):
        try:
            stage, terminal = pair
        except (TypeError, ValueError):
            raise ValueError("INVALID_STAGE_RECORD") from None
        if stage.get("index") != expected_index:
            raise ValueError("PASS_INDEX_MISMATCH")
        if previous is not None:
            if previous.get("after_hash") != stage.get("before_hash"):
                raise ValueError("PARTITION_CHAIN_MISMATCH")
            if "input_beta" in stage:
                _same(_array(stage["input_beta"], (len(previous["output_beta"]),), "BAD_INPUT_BETA"),
                      np.asarray(previous["output_beta"], dtype=float), "INPUT_BETA_CHAIN_MISMATCH")
        diagnostics.append(pass_diagnostics(header, stage, terminal))
        previous, last_stage = stage, stage
    if last_stage is None:
        raise ValueError("MISSING_PASSES")

    solve = header.get("solve", {})
    final_partition = solve.get("final_partition")
    if (solve.get("final_partition_hash") != last_stage.get("after_hash") or
            final_partition != last_stage.get("after")):
        raise ValueError("TERMINAL_PARTITION_MISMATCH")
    final_beta = _array(solve.get("relaxed_fill"), (len(last_stage["output_beta"]),),
                        "TERMINAL_PARTITION_MISMATCH")
    _same(final_beta, np.asarray(last_stage["output_beta"], dtype=float), "TERMINAL_PARTITION_MISMATCH")
    if not math.isclose(_scalar(solve.get("final_phi"), "INVALID_FINAL_PHI"),
                        _scalar(last_stage.get("final_phi"), "INVALID_FINAL_PHI"),
                        rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("TERMINAL_PHI_MISMATCH")
    recorded_applies = solve.get("terminal_relaxation_kkt_applies")
    expected_applies = _terminal_kkt_applicability(last_stage) == "APPLICABLE_UNCHANGED_TERMINAL_STATE"
    if type(recorded_applies) is not bool or recorded_applies != expected_applies:
        raise ValueError("TERMINAL_KKT_APPLICABILITY_MISMATCH")

    geometry = header.get("geometry", {})
    geometry_stage = ("NOT_AVAILABLE" if solve.get("status") != "BINARY_MOMENT_VALIDATED" else
                      geometry.get("state", "NOT_AVAILABLE"))
    return {
        "schema": "spin2014-stage-audit-v1",
        "scope": "DIAGNOSTIC_ALGEBRA_NOT_SOLVER_ADMISSION_OR_GEOMETRY_CERTIFICATION",
        "source_digests": source_digests,
        "passes": diagnostics,
        "terminal_rounding": _rounding(header),
        "geometry_export_stage": geometry_stage,
    }


def _stage_pairs(events):
    current = None
    terminal = None
    last_trace_index = None
    for event in events:
        kind = event["kind"]
        if kind == "header":
            raise ValueError("DUPLICATE_HEADER")
        if kind == "pass":
            if current is not None:
                if terminal is None or terminal[0] != last_trace_index:
                    raise ValueError("MISSING_ACCEPTED_TERMINAL_STATE")
                yield current, terminal[1]
            current = event["payload"]
            terminal = None
            last_trace_index = None
        elif kind == "iteration":
            if current is None:
                raise ValueError("ITERATION_WITHOUT_PASS")
            last_trace_index = event["trace_index"]
            if event["payload"].get("final_state") is True:
                if terminal is not None:
                    raise ValueError("MULTIPLE_TERMINAL_STATES")
                terminal = (event["trace_index"], event["payload"])
    if current is not None:
        if terminal is None or terminal[0] != last_trace_index:
            raise ValueError("MISSING_ACCEPTED_TERMINAL_STATE")
        yield current, terminal[1]


def _file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _footer(path):
    with Path(path).open("rb") as stream:
        stream.seek(0, 2)
        size = stream.tell()
        stream.seek(max(0, size - 4096))
        lines = stream.read().splitlines()
    if not lines:
        raise ValueError("MISSING_COMPLETION_FOOTER")
    value = json.loads(lines[-1])
    if not isinstance(value, dict) or value.get("kind") != "complete":
        raise ValueError("MISSING_COMPLETION_FOOTER")
    return value


def audit_archive(path):
    """Fully exhaust and verify an archive before returning diagnostics."""
    reader = read_archive(path)
    try:
        first = next(reader)
    except StopIteration:
        raise ValueError("MISSING_HEADER") from None
    if first.get("kind") != "header":
        raise ValueError("MISSING_HEADER")
    result = audit_components(first["payload"], _stage_pairs(reader))
    # audit_components exhausts the pairs, which exhausts read_archive and verifies
    # its sequence, completion footer, payload checksum, and trailing EOF.
    footer = _footer(path)
    result["archive_integrity"] = "VERIFIED_COMPLETE_SOURCE_TRACE_ARCHIVE_V1"
    result["input_archive_sha256"] = _file_sha256(path)
    result["input_payload_sha256"] = _digest(footer.get("payload_sha256"), "INVALID_ARCHIVE_DIGEST")
    result["checker_digests"] = {
        "stage_audit_sha256": _file_sha256(Path(__file__)),
        "rounding_diagnostics_sha256": _file_sha256(Path(__file__).with_name("resolution_rounding_audit.py")),
        "archive_reader_sha256": _file_sha256(PYTHON / "make_it_stand" / "verification" / "source_trace_archive.py"),
    }
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args(argv)
    result = audit_archive(args.archive)
    encoded = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    with args.output.open("x", encoding="utf-8") as stream:
        stream.write(encoded)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
