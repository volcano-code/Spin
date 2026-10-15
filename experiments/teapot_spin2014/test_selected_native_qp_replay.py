"""Command-line integration checks for a provenance-bound fixed-QP replay.

This diagnostic is invoked explicitly with a replay executable and a sealed
native trace; it is not part of the default test suite because paper traces
are external, multi-gigabyte run artifacts.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path


def replay_selected(binary: Path, trace: Path, *, bounded: int,
                    selected_iteration: int, history_capacity: int,
                    maximum_faces: int, maximum_iterations: int) -> dict:
    completed = subprocess.run(
        [str(binary), str(trace), str(bounded), str(history_capacity),
         str(maximum_faces), str(maximum_iterations), str(selected_iteration)],
        check=True, capture_output=True, text=True, timeout=180,
    )
    return json.loads(completed.stdout)


def assert_selected_replay(binary: Path, trace: Path, *, bounded: int,
                           selected_iteration: int, history_capacity: int,
                           maximum_faces: int, maximum_iterations: int) -> dict:
    result = replay_selected(binary, trace, bounded=bounded,
                             selected_iteration=selected_iteration,
                             history_capacity=history_capacity,
                             maximum_faces=maximum_faces,
                             maximum_iterations=maximum_iterations)
    assert result["selected_iteration"] == selected_iteration
    assert result["records"] > selected_iteration
    assert result["scope"] == "FIXED_QP_REPLAY_NOT_NONLINEAR_OR_BINARY_VALIDATION"
    assert result["terminal"]["iteration"] == selected_iteration
    return result
