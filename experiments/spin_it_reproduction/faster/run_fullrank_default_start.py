#!/usr/bin/env python3
"""Run one fixed synthetic full-rank Faster case with native default starts.

This is a method diagnostic only. It is not an author case or a paper
experiment. No solution quadric is supplied to the optimizer.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import time
from typing import Any

from make_it_stand.spin import (
    SpinEvaluatorConfig,
    SpinItFasterRequest,
    SpinSolverConfig,
    integrate_raw_mass_moments,
    solve_spin_it_faster_fullrank_quadric,
)


def cube(low: float, high: float, z_offset: float):
    vertices = (
        (low, low, low + z_offset), (high, low, low + z_offset),
        (high, high, low + z_offset), (low, high, low + z_offset),
        (low, low, high + z_offset), (high, low, high + z_offset),
        (high, high, high + z_offset), (low, high, high + z_offset),
    )
    faces = (
        (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4), (3, 7, 6), (3, 6, 2),
        (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5),
    )
    return vertices, faces


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite_json(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return "NaN" if math.isnan(value) else ("+Infinity" if value > 0 else "-Infinity")
    if isinstance(value, dict):
        return {str(key): finite_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [finite_json(item) for item in value]
    return value


def git(root: Path, *arguments: str) -> str:
    return subprocess.run(
        ["git", *arguments], cwd=root, check=True,
        capture_output=True, text=True,
    ).stdout.strip()


def run(seed: int, output: Path) -> int:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing run evidence: {output}")

    outer_vertices, outer_faces = cube(-1.0, 1.0, 2.0)
    design_vertices, design_faces = cube(-0.8, 0.8, 2.0)
    outer_moments = integrate_raw_mass_moments(outer_vertices, outer_faces)
    design_moments = integrate_raw_mass_moments(design_vertices, design_faces)
    fixed_shell = tuple(a - b for a, b in zip(outer_moments, design_moments))
    request = SpinItFasterRequest(
        outer_vertices=outer_vertices,
        outer_faces=outer_faces,
        design_vertices=design_vertices,
        design_faces=design_faces,
        fixed_shell_moments=fixed_shell,
        mode="SPINNING_TOP",
        center_of_mass_weight=1.0,
        minimum_wall_thickness=0.15,
        evidence="RECONSTRUCTED_ASSET",
        settings=SpinSolverConfig(
            equation_profile="DUAL_MONOMIAL_RAW_V1",
            maximum_iterations=60,
            maximum_restarts=3,
            maximum_function_evaluations_per_restart=60,
            equation_tolerance=1e-12,
            physical_relative_tolerance=1e-6,
            deterministic_seed=seed,
        ),
    )
    evaluator = SpinEvaluatorConfig(icosphere_subdivisions=1, truncation_radius=8.0)
    module = importlib.util.find_spec("_make_it_stand")
    if module is None or module.origin is None:
        raise RuntimeError("native extension could not be resolved")

    started = time.perf_counter()
    # Deliberately pass no initial quadrics: the solver's deterministic
    # mode-aware restart sequence is the independent initialization protocol.
    result = solve_spin_it_faster_fullrank_quadric(
        request, evaluator_settings=evaluator
    )
    elapsed = time.perf_counter() - started
    repository = Path(__file__).resolve().parents[3]
    status = git(repository, "status", "--porcelain=v1", "--untracked-files=all")
    payload = {
        "schema": "spin-faster-fullrank-default-start-diagnostic-v1",
        "claim_level": "SYNTHETIC_REGRESSION",
        "not_a_paper_case": True,
        "initialization_protocol": "native_deterministic_mode_aware_restarts_no_user_quadric",
        "initial_quadrics_supplied": [],
        "seed": seed,
        "request": {
            "outer_cube": [-1.0, 1.0, 2.0],
            "design_cube": [-0.8, 0.8, 2.0],
            "mode": request.mode,
            "center_of_mass_weight": request.center_of_mass_weight,
            "minimum_wall_thickness": request.minimum_wall_thickness,
            "equation_profile": request.settings.equation_profile,
            "maximum_iterations": request.settings.maximum_iterations,
            "maximum_restarts": request.settings.maximum_restarts,
            "maximum_function_evaluations_per_restart": request.settings.maximum_function_evaluations_per_restart,
            "equation_tolerance": request.settings.equation_tolerance,
            "physical_relative_tolerance": request.settings.physical_relative_tolerance,
            "icosphere_subdivisions": evaluator.icosphere_subdivisions,
            "truncation_radius": evaluator.truncation_radius,
        },
        "provenance": {
            "git_commit": git(repository, "rev-parse", "HEAD"),
            "worktree_clean_at_start": not bool(status),
            "worktree_status_sha256": hashlib.sha256(status.encode()).hexdigest(),
            "runner_sha256": sha256(Path(__file__)),
            "native_extension_sha256": sha256(Path(module.origin)),
        },
        "elapsed_seconds_non_scientific": elapsed,
        "result": finite_json(asdict(result)),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "output": str(output),
        "status": result.status,
        "diagnostic": result.diagnostic,
        "final_geometry_verified": result.final_geometry_verified,
        "restart_count": result.restart_count,
        "restart_function_evaluations": result.restart_function_evaluations,
        "function_evaluations": result.function_evaluations,
        "event_count": len(result.events),
        "event_hash": result.event_hash,
        "output_hash": result.output_hash,
        "elapsed_seconds": elapsed,
    }, sort_keys=True))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    return run(arguments.seed, arguments.output)


if __name__ == "__main__":
    raise SystemExit(main())
