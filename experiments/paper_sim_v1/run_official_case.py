#!/usr/bin/env python3
"""Run one exact-cache official case under an explicit execution profile."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import resource
import shutil
import signal
import time
from typing import Any

from make_it_stand.official_beta import ReleasedBetaArchive
from make_it_stand.paper_compat import (
    SupportFlattening,
    paper_run_config,
    released_beta_run_config,
)
from make_it_stand.paper_evidence import (
    FreshnessStatus,
    capture_execution_fingerprint,
    classify_evidence,
)
from make_it_stand.paper_trace import TraceProfile, build_iteration_trace_v2
from make_it_stand.sealed_artifact import SealedBalanceArtifact


SUPPORTED_CASES = ("spheres", "armadillo1", "armadillo2", "dinosaur", "horse")
PROFILE_BUILDERS = {
    "paper": paper_run_config,
    "released-beta": released_beta_run_config,
}
TRACE_PROFILES = {
    "paper": TraceProfile.PAPER_PROFILE,
    "released-beta": TraceProfile.RELEASED_BETA_PROFILE,
}


class WallClockLimit(RuntimeError):
    pass


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _alarm_handler(_signum: int, _frame: Any) -> None:
    raise WallClockLimit("official case exceeded its wall-clock limit")


def _input_files(root: Path, case: str, archive: ReleasedBetaArchive) -> dict[str, Path]:
    configs = root / "data" / "configs"
    meshes = root / "data" / "meshes"
    stem = f"{case}_{archive.voxel_resolution}"
    return {
        "mis": configs / f"{case}.mis",
        "obj": meshes / archive.specification.mesh,
        "vox": configs / f"{stem}.vox",
        "bbw": configs / f"{stem}.bbw",
    }


def _source_paths(repository_root: Path) -> list[Path]:
    return [
        repository_root / "CMakeLists.txt",
        repository_root / "CMakePresets.json",
        repository_root / "bindings",
        repository_root / "include",
        repository_root / "src",
        repository_root / "python",
        Path(__file__).resolve(),
    ]


def _executable_paths() -> list[Path]:
    paths = [Path(__import__("sys").executable).resolve()]
    for module_name in ("_make_it_stand", "make_it_stand._make_it_stand"):
        try:
            specification = importlib.util.find_spec(module_name)
        except (ImportError, ModuleNotFoundError, ValueError):
            continue
        if specification is not None and specification.origin:
            paths.append(Path(specification.origin).resolve())
            break
    return paths


def _flattening_reference_vertices(archive: ReleasedBetaArchive):
    """Use the exact outer surface optimized by the current archive object."""

    return archive.normalized_vertices


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--beta-root", required=True, type=Path)
    parser.add_argument("--case", required=True, choices=SUPPORTED_CASES)
    parser.add_argument("--profile", required=True, choices=tuple(PROFILE_BUILDERS))
    parser.add_argument("--maximum-carving-blocks", type=int, default=400)
    parser.add_argument("--maximum-deformation-iterations", type=int, default=400)
    parser.add_argument("--shape-scale-lambda", type=float)
    parser.add_argument(
        "--enable-carving",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument(
        "--enable-deformation",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument(
        "--enable-laplacian",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--replacement-normalized-surface", type=Path)
    parser.add_argument("--time-limit-seconds", type=int, default=3600)
    parser.add_argument("--memory-limit-gib", type=float, default=24.0)
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    if arguments.maximum_carving_blocks <= 0:
        parser.error("--maximum-carving-blocks must be positive")
    if arguments.maximum_deformation_iterations <= 0:
        parser.error("--maximum-deformation-iterations must be positive")
    if (
        arguments.shape_scale_lambda is not None and
        (
            not math.isfinite(arguments.shape_scale_lambda) or
            arguments.shape_scale_lambda < 0.0
        )
    ):
        parser.error("--shape-scale-lambda must be finite and nonnegative")
    if arguments.time_limit_seconds <= 0:
        parser.error("--time-limit-seconds must be positive")
    if not math.isfinite(arguments.memory_limit_gib) or arguments.memory_limit_gib <= 0:
        parser.error("--memory-limit-gib must be finite and positive")
    if arguments.repeat <= 0:
        parser.error("--repeat must be positive")

    beta_root = arguments.beta_root.expanduser().resolve()
    output = arguments.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    archive = ReleasedBetaArchive.load(
        beta_root / "data" / "configs" / f"{arguments.case}.mis",
        beta_root=beta_root,
    )
    replacement_surface = None
    if arguments.replacement_normalized_surface is not None:
        replacement_surface = (
            arguments.replacement_normalized_surface.expanduser().resolve()
        )
        archive = archive.with_reconstructed_normalized_surface(
            replacement_surface
        )
    if arguments.case == "horse" and replacement_surface is None:
        parser.error(
            "horse's untouched OBJ fails the compatibility boundary gate; "
            "provide --replacement-normalized-surface for a labelled reconstruction"
        )
    if arguments.case != "horse" and replacement_surface is not None:
        parser.error(
            "--replacement-normalized-surface is currently restricted to horse"
        )
    inputs = _input_files(beta_root, arguments.case, archive)
    if replacement_surface is not None:
        inputs["replacement_normalized_surface"] = replacement_surface
    input_hashes = {name: _sha256(path) for name, path in inputs.items()}
    evidence_class = classify_evidence(
        complete_cache=True,
        parameter_override=arguments.shape_scale_lambda is not None,
        reconstructed_surface=replacement_surface is not None,
    )

    memory_bytes = int(arguments.memory_limit_gib * 1024**3)
    resource.setrlimit(resource.RLIMIT_AS, (memory_bytes, memory_bytes))
    signal.signal(signal.SIGALRM, _alarm_handler)

    overrides = {
        "maximum_carving_blocks": arguments.maximum_carving_blocks,
        "maximum_deformation_iterations_per_block": (
            arguments.maximum_deformation_iterations
        ),
        "enable_carving": arguments.enable_carving,
        "enable_deformation": arguments.enable_deformation,
        "enable_laplacian": arguments.enable_laplacian,
    }
    if arguments.shape_scale_lambda is not None:
        overrides["shape_scale_lambda"] = arguments.shape_scale_lambda
    config = PROFILE_BUILDERS[arguments.profile](**overrides)
    repository_root = Path(__file__).resolve().parents[2]
    execution_fingerprint = capture_execution_fingerprint(
        repository_root,
        source_paths=[*_source_paths(repository_root), *inputs.values()],
        executable_paths=_executable_paths(),
        configuration={
            "case": arguments.case,
            "profile": arguments.profile,
            "run_config": asdict(config),
            "limits": {
                "memory_gib": arguments.memory_limit_gib,
                "wall_clock_seconds": arguments.time_limit_seconds,
            },
        },
    )

    summaries: list[dict[str, Any]] = []
    exit_code = 0
    for repetition in range(arguments.repeat):
        run_output = output / f"repeat-{repetition + 1}"
        run_output.mkdir(parents=True, exist_ok=True)
        model = archive.build_model()
        started = time.monotonic()
        signal.alarm(arguments.time_limit_seconds)
        try:
            result = model.optimize(config)
        except WallClockLimit as error:
            elapsed = time.monotonic() - started
            summary = {
                "completion": "timeout",
                "error": str(error),
                "elapsed_seconds": elapsed,
            }
            exit_code = 124
        except Exception as error:
            elapsed = time.monotonic() - started
            summary = {
                "completion": "error",
                "error_type": type(error).__name__,
                "error": str(error),
                "elapsed_seconds": elapsed,
            }
            exit_code = 1
        else:
            elapsed = time.monotonic() - started
            diagnostics_path = run_output / "diagnostics.json"
            outer_path = run_output / "outer.obj"
            result.save_diagnostics(diagnostics_path)
            raw_outer_path = run_output / "outer-optimizer-raw.obj"
            result.save_outer_mesh(raw_outer_path)
            raw_inner_path = run_output / "inner-optimizer-raw.obj"
            if result.native.inner_faces.size:
                # Preserve the expensive optimizer output before sealing so
                # an audit failure can be diagnosed without re-optimizing.
                result.save_inner_mesh(raw_inner_path)
            diagnostics = result.diagnostics
            trace = build_iteration_trace_v2(
                diagnostics,
                TRACE_PROFILES[arguments.profile],
            )
            trace_path = run_output / "iteration-trace-v2.json"
            _write_json(trace_path, trace)
            flattening = SupportFlattening.from_vertices(
                _flattening_reference_vertices(archive),
                archive.gravity,
                archive.specification.flatten_threshold,
                1.0e-4,
            )
            final_event = diagnostics["r10_events"][-1]
            support_polygon = diagnostics["problem"]["standing_constraint"][
                "support_polygon"
            ]
            sealed = SealedBalanceArtifact.create(
                outer_vertices=result.outer_vertices,
                outer_faces=result.outer_faces,
                inner_vertices=result.inner_vertices,
                inner_faces=result.inner_faces,
                flattenings=[flattening],
                occupancy_hash=str(final_event["alpha_hash"]),
                source_geometry_hash=str(diagnostics["geometry_hash"]),
                gravity=archive.gravity,
                support_polygon=support_polygon,
                evidence_class=evidence_class,
                native_module=model._native_module,
                structured_grid_dimensions=tuple(
                    diagnostics["r10_metrics"]["grid_dimensions"]
                ),
                repair_nonmanifold_cavity=True,
                safe_region_gamma_radians=float(
                    diagnostics["problem"]["standing_constraint"][
                        "gamma_radians"
                    ]
                ),
            )
            sealed_root = run_output / "sealed"
            sealed_manifest = sealed.write_directory(sealed_root)
            shutil.copyfile(sealed_root / "outer-flattened.obj", outer_path)
            output_hashes = {
                "diagnostics": _sha256(diagnostics_path),
                "outer": _sha256(outer_path),
                "outer_optimizer_raw": _sha256(raw_outer_path),
                "iteration_trace_v2": _sha256(trace_path),
                "sealed_artifact": _sha256(sealed_manifest),
            }
            if raw_inner_path.is_file():
                output_hashes["inner_optimizer_raw"] = _sha256(raw_inner_path)
            inner_path = run_output / "inner.obj"
            if result.native.inner_faces.size:
                result.save_inner_mesh(inner_path)
                output_hashes["inner"] = _sha256(inner_path)
            paper_iterations = {
                item["paper_outer_iteration"]
                for item in trace["events"]
                if item["paper_outer_iteration"] is not None
            }
            beta_calls = {
                item["beta_optimize_call"]
                for item in trace["deformation_trials"]
                if item["beta_optimize_call"] is not None
            }
            summary = {
                "completion": "stable" if result.succeeded else "not_converged",
                "elapsed_seconds": elapsed,
                "status": result.status,
                "succeeded": result.succeeded,
                "output_sha256": output_hashes,
                "paper_outer_iterations": len(paper_iterations),
                "beta_optimize_calls": len(beta_calls),
                "committed_events": len(trace["events"]),
                "carving_commits": len(diagnostics["carving_iterations"]),
                "deformation_trials": len(diagnostics["deformation_iterations"]),
                "first_standing_feasible_event": diagnostics["r10_metrics"][
                    "first_standing_feasible_event"
                ],
                "first_safe_event": diagnostics["r10_metrics"]["first_safe_event"],
                "final_center_of_mass": diagnostics["center_of_mass"],
                "final_projected_target_energy": diagnostics["r10_metrics"][
                    "final_projected_target_energy"
                ],
                "sealed_artifact_status": (
                    "SEALED" if sealed.ready else "REJECTED"
                ),
                "sealed_artifact_failures": list(sealed.failures),
                "sealed_geometry_sha256": sealed.geometry_sha256,
                "sealed_mass_properties_sha256": (
                    sealed.mass_properties_sha256
                ),
                "sealed_safe_region_state": sealed.safe_region_state,
                "sealed_safe_region_residual": sealed.safe_region_residual,
                "sealed_bounding_box_diagonal": sealed.bounding_box_diagonal,
                "sealed_source_occupancy_hash": sealed.source_occupancy_hash,
                "sealed_occupancy_hash": sealed.occupancy_hash,
                "sealed_refilled_cavity_cells": len(
                    sealed.occupancy_repair_audit.get(
                        "refilled_cells_xyz",
                        [],
                    )
                ),
                "sealed_outer_flattening_repair_status": (
                    sealed.outer_flattening_repair_audit.get("status")
                ),
                "sealed_outer_flattening_repair_patches": len(
                    sealed.outer_flattening_repair_audit.get("patches", [])
                ),
                "sealed_wall_thickness_status": (
                    sealed.wall_thickness_audit.get("status")
                ),
                "sealed_measured_minimum_wall_clearance": (
                    sealed.wall_thickness_audit.get(
                        "measured_minimum_clearance"
                    )
                ),
                "sealed_containment_status": sealed.containment_audit.get(
                    "status"
                ),
            }
            if not result.succeeded:
                exit_code = 2
        finally:
            signal.alarm(0)
        _write_json(run_output / "summary.json", summary)
        summaries.append(summary)
        if exit_code != 0:
            break

    deterministic = None
    if len(summaries) == arguments.repeat and arguments.repeat > 1:
        deterministic = all(
            item.get("output_sha256") == summaries[0].get("output_sha256")
            for item in summaries[1:]
        )
        if not deterministic:
            exit_code = 3

    manifest = {
        "schema": "make-it-stand-paper-sim-official-case-v1",
        "schema_version": 1,
        "case": arguments.case,
        "profile": arguments.profile,
        "classification": evidence_class.value,
        "freshness": {
            "status": FreshnessStatus.CURRENT.value,
            "mismatches": [],
        },
        "execution_fingerprint": execution_fingerprint,
        "surface_provenance": (
            "untouched_official_obj"
            if replacement_surface is None
            else "reconstructed_from_official_occupancy_not_original_surface"
        ),
        "parameter_provenance": {
            "shape_scale_lambda": (
                "profile_default"
                if arguments.shape_scale_lambda is None
                else "explicit_reconstructed_experiment_not_archived_by_authors"
            )
        },
        "voxel_resolution": archive.voxel_resolution,
        "occupied_cells": archive.occupied_cell_count,
        "active_nodes": archive.active_node_count,
        "input_sha256": input_hashes,
        "configuration": asdict(config),
        "limits": {
            "memory_gib": arguments.memory_limit_gib,
            "wall_clock_seconds": arguments.time_limit_seconds,
        },
        "repeat_count_requested": arguments.repeat,
        "deterministic_output_hashes": deterministic,
        "runs": summaries,
        "claim_boundary": (
            "A timeout or not_converged run is evidence of executable wiring only; "
            "it is not a successful paper-case reproduction."
        ),
    }
    _write_json(output / "manifest.json", manifest)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
