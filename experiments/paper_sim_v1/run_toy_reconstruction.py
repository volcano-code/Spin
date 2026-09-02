#!/usr/bin/env python3
"""Run an explicitly reconstructed rocking-toy case from the only public OBJ.

The beta archive contains no MIS, occupancy, BBW, pose, flatten threshold, or
handle history for this model.  This driver therefore regenerates every
missing item deterministically and permanently labels the result
``RECONSTRUCTED_ASSET``.  It cannot produce ``EXACT_CACHE`` evidence.
"""

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
from typing import Any, Sequence

import numpy as np

from make_it_stand.model import ModelConfig
from make_it_stand.paper_compat import (
    PaperCompatibilityModel,
    SupportFlattening,
    paper_run_config,
    released_beta_run_config,
)
from make_it_stand.paper_evidence import (
    EvidenceClass,
    FreshnessStatus,
    capture_execution_fingerprint,
)
from make_it_stand.paper_trace import TraceProfile, build_iteration_trace_v2
from make_it_stand.sealed_artifact import SealedBalanceArtifact


class WallClockLimit(RuntimeError):
    pass


def _alarm_handler(_signum: int, _frame: Any) -> None:
    raise WallClockLimit("toy reconstruction exceeded its wall-clock limit")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _write_obj(path: Path, vertices: np.ndarray, faces: np.ndarray) -> None:
    lines = [
        f"v {point[0]:.17g} {point[1]:.17g} {point[2]:.17g}"
        for point in vertices
    ]
    lines.extend(
        f"f {face[0] + 1} {face[1] + 1} {face[2] + 1}" for face in faces
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def reconstruct_toy_surface(
    path: Path,
    *,
    resolution: int,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Deterministically close the public 29-component open triangle model."""

    if isinstance(resolution, bool) or not isinstance(resolution, int) or resolution < 16:
        raise ValueError("surface reconstruction resolution must be at least 16")
    try:
        import trimesh
    except ImportError as error:
        raise RuntimeError(
            "rocking-toy reconstruction requires the optional trimesh dependency"
        ) from error
    source = trimesh.load_mesh(path, process=True, validate=True)
    if not isinstance(source, trimesh.Trimesh):
        raise ValueError("rocking-toy OBJ did not load as one triangle mesh")
    longest = float(np.max(source.extents))
    if not math.isfinite(longest) or longest <= 0.0:
        raise ValueError("rocking-toy OBJ has no positive extent")
    pitch = longest / float(resolution)
    voxelized = source.voxelized(pitch).fill()
    reconstructed_index = voxelized.marching_cubes
    reconstructed_vertices = trimesh.transform_points(
        reconstructed_index.vertices,
        voxelized.transform,
    )
    reconstructed = trimesh.Trimesh(
        vertices=reconstructed_vertices,
        faces=reconstructed_index.faces,
        process=True,
        validate=True,
    )
    if reconstructed.volume < 0.0:
        reconstructed.invert()
    components = reconstructed.split(only_watertight=False)
    if (
        not reconstructed.is_watertight
        or not reconstructed.is_winding_consistent
        or len(components) != 1
        or reconstructed.volume <= 0.0
    ):
        raise ValueError("voxel-union reconstruction did not produce one closed shell")
    metadata = {
        "algorithm": "trimesh_surface_voxelize_fill_marching_cubes_v1",
        "trimesh_version": str(trimesh.__version__),
        "surface_resolution_longest_axis": resolution,
        "surface_pitch_source_units": pitch,
        "source_vertices_after_processing": int(len(source.vertices)),
        "source_faces": int(len(source.faces)),
        "source_watertight": bool(source.is_watertight),
        "source_components": int(len(source.split(only_watertight=False))),
        "reconstructed_vertices": int(len(reconstructed.vertices)),
        "reconstructed_faces": int(len(reconstructed.faces)),
        "reconstructed_watertight": bool(reconstructed.is_watertight),
        "reconstructed_components": int(len(components)),
        "reconstructed_volume_source_units": float(reconstructed.volume),
    }
    return (
        np.ascontiguousarray(reconstructed.vertices, dtype=np.float64),
        np.ascontiguousarray(reconstructed.faces, dtype=np.int32),
        metadata,
    )


def select_reconstructed_handles(
    vertices: Sequence[Sequence[float]],
    *,
    gravity: Sequence[float],
    support_vertex_indices: Sequence[int],
    count: int,
) -> tuple[int, ...]:
    """Deterministic farthest-surface sampling outside the support band."""

    points = np.asarray(vertices, dtype=np.float64)
    direction = np.asarray(gravity, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError("vertices must be a finite N by 3 array")
    if direction.shape != (3,) or not np.isfinite(direction).all():
        raise ValueError("gravity must contain three finite values")
    norm = float(np.linalg.norm(direction))
    if norm <= 1.0e-15:
        raise ValueError("gravity must be nonzero")
    if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
        raise ValueError("count must be a positive integer")
    direction /= norm
    blocked = {int(index) for index in support_vertex_indices}
    candidates = np.asarray(
        [index for index in range(points.shape[0]) if index not in blocked],
        dtype=np.int64,
    )
    if candidates.size < count:
        raise ValueError("not enough non-support vertices for reconstructed handles")

    # Start at the greatest height opposite gravity; np.argmax supplies the
    # lowest vertex-index tie break.  Subsequent samples maximize distance to
    # the selected set, again with a stable vertex-index tie break.
    heights = -(points[candidates] @ direction)
    chosen = [int(candidates[int(np.argmax(heights))])]
    while len(chosen) < count:
        remaining = np.asarray(
            [index for index in candidates if int(index) not in chosen],
            dtype=np.int64,
        )
        differences = (
            points[remaining, np.newaxis, :]
            - points[np.asarray(chosen, dtype=np.int64), :][np.newaxis, :, :]
        )
        minimum_squared_distance = np.min(
            np.einsum("mki,mki->mk", differences, differences),
            axis=1,
        )
        chosen.append(int(remaining[int(np.argmax(minimum_squared_distance))]))
    return tuple(chosen)


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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--beta-root", required=True, type=Path)
    parser.add_argument("--profile", required=True, choices=("paper", "released-beta"))
    parser.add_argument("--voxel-resolution", type=int, default=50)
    parser.add_argument("--surface-reconstruction-resolution", type=int, default=80)
    parser.add_argument("--flatten-threshold", type=float, default=0.005)
    parser.add_argument("--handle-count", type=int, default=4)
    parser.add_argument("--shape-scale-lambda", type=float)
    parser.add_argument("--maximum-carving-blocks", type=int, default=400)
    parser.add_argument("--maximum-deformation-iterations", type=int, default=400)
    parser.add_argument("--time-limit-seconds", type=int, default=3600)
    parser.add_argument("--memory-limit-gib", type=float, default=24.0)
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    if arguments.voxel_resolution < 2:
        parser.error("--voxel-resolution must be at least two")
    if arguments.surface_reconstruction_resolution < 16:
        parser.error("--surface-reconstruction-resolution must be at least 16")
    if arguments.flatten_threshold < 0.0 or not math.isfinite(arguments.flatten_threshold):
        parser.error("--flatten-threshold must be finite and nonnegative")
    if arguments.handle_count <= 0:
        parser.error("--handle-count must be positive")
    if arguments.repeat <= 0:
        parser.error("--repeat must be positive")

    beta_root = arguments.beta_root.expanduser().resolve()
    mesh = beta_root / "data" / "meshes" / "toyrocking.obj"
    if not mesh.is_file():
        parser.error(f"missing rocking-toy OBJ: {mesh}")
    output = arguments.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    memory_bytes = int(arguments.memory_limit_gib * 1024**3)
    resource.setrlimit(resource.RLIMIT_AS, (memory_bytes, memory_bytes))
    signal.signal(signal.SIGALRM, _alarm_handler)
    reconstructed_vertices, reconstructed_faces, reconstruction_audit = (
        reconstruct_toy_surface(
            mesh,
            resolution=arguments.surface_reconstruction_resolution,
        )
    )
    reconstructed_input = output / "reconstructed-input.obj"
    _write_obj(reconstructed_input, reconstructed_vertices, reconstructed_faces)
    gravity = (0.0, -1.0, 0.0)
    gamma = math.radians(5.0)
    overrides: dict[str, Any] = {
        "maximum_carving_blocks": arguments.maximum_carving_blocks,
        "maximum_deformation_iterations_per_block": (
            arguments.maximum_deformation_iterations
        ),
    }
    if arguments.shape_scale_lambda is not None:
        overrides["shape_scale_lambda"] = arguments.shape_scale_lambda
    config = (
        paper_run_config(**overrides)
        if arguments.profile == "paper"
        else released_beta_run_config(**overrides)
    )
    trace_profile = (
        TraceProfile.PAPER_PROFILE
        if arguments.profile == "paper"
        else TraceProfile.RELEASED_BETA_PROFILE
    )
    repository_root = Path(__file__).resolve().parents[2]
    reconstruction = {
        "voxel_resolution": arguments.voxel_resolution,
        "surface_reconstruction": reconstruction_audit,
        "gravity": list(gravity),
        "mesh_rotation_wxyz": [1.0, 0.0, 0.0, 0.0],
        "flatten_threshold": arguments.flatten_threshold,
        "gamma_radians": gamma,
        "beta_minimum_stability_inset": True,
        "handle_policy": "deterministic_farthest_surface_v1",
        "handle_count": arguments.handle_count,
        "missing_author_items": [
            "mis", "vox", "bbw", "handles", "rotation", "flatten_threshold"
        ],
    }
    fingerprint = capture_execution_fingerprint(
        repository_root,
        source_paths=[
            repository_root / "CMakeLists.txt",
            repository_root / "CMakePresets.json",
            repository_root / "bindings",
            repository_root / "include",
            repository_root / "src",
            repository_root / "python",
            Path(__file__).resolve(),
            mesh,
            reconstructed_input,
        ],
        executable_paths=_executable_paths(),
        configuration={
            "case": "toyrocking",
            "profile": arguments.profile,
            "run_config": asdict(config),
            "reconstruction": reconstruction,
        },
    )

    summaries: list[dict[str, Any]] = []
    exit_code = 0
    selected_handles: tuple[int, ...] | None = None
    for repetition in range(1, arguments.repeat + 1):
        run_root = output / f"repeat-{repetition}"
        run_root.mkdir(parents=True, exist_ok=True)
        started = time.monotonic()
        signal.alarm(arguments.time_limit_seconds)
        try:
            model = PaperCompatibilityModel(
                reconstructed_vertices,
                reconstructed_faces,
                config=ModelConfig(
                    path="<reconstructed-toy-memory>",
                    voxel_resolution=arguments.voxel_resolution,
                    minimum_wall_thickness=0.0,
                ),
            )
            model.set_gravity(gravity)
            model.set_support_from_mesh(
                arguments.flatten_threshold,
                gamma,
                beta_minimum_inset=True,
            )
            flattening = SupportFlattening.from_vertices(
                model.normalized_vertices,
                gravity,
                arguments.flatten_threshold,
                1.0e-4,
            )
            handles = select_reconstructed_handles(
                model.normalized_vertices,
                gravity=gravity,
                support_vertex_indices=flattening.vertex_indices,
                count=arguments.handle_count,
            )
            if selected_handles is None:
                selected_handles = handles
            elif handles != selected_handles:
                raise RuntimeError("reconstructed handle selection is nondeterministic")
            for index in handles:
                model.add_normalized_handle(model.normalized_vertices[index])
            result = model.optimize(config)
            elapsed = time.monotonic() - started
            native_result = result.native_result
            diagnostics = native_result.diagnostics
            diagnostics_path = run_root / "diagnostics.json"
            native_result.save_diagnostics(diagnostics_path)
            raw_outer = run_root / "outer-optimizer-raw.obj"
            native_result.save_outer_mesh(raw_outer)
            raw_inner = run_root / "inner-optimizer-raw.obj"
            if native_result.inner_faces.size:
                native_result.save_inner_mesh(raw_inner)
            trace = build_iteration_trace_v2(diagnostics, trace_profile)
            trace_path = run_root / "iteration-trace-v2.json"
            _write_json(trace_path, trace)
            final_event = diagnostics["r10_events"][-1]
            sealed = SealedBalanceArtifact.create(
                outer_vertices=native_result.outer_vertices,
                outer_faces=native_result.outer_faces,
                inner_vertices=native_result.inner_vertices,
                inner_faces=native_result.inner_faces,
                flattenings=result.flattenings,
                occupancy_hash=str(final_event["alpha_hash"]),
                source_geometry_hash=str(diagnostics["geometry_hash"]),
                gravity=gravity,
                support_polygon=diagnostics["problem"]["standing_constraint"][
                    "support_polygon"
                ],
                evidence_class=EvidenceClass.RECONSTRUCTED_ASSET,
                native_module=model._native_module,
                structured_grid_dimensions=tuple(
                    diagnostics["r10_metrics"]["grid_dimensions"]
                ),
                repair_nonmanifold_cavity=True,
                safe_region_gamma_radians=gamma,
            )
            sealed_manifest = sealed.write_directory(run_root / "sealed")
            outer_path = run_root / "outer.obj"
            shutil.copyfile(run_root / "sealed" / "outer-flattened.obj", outer_path)
            hashes = {
                "diagnostics": _sha256(diagnostics_path),
                "outer": _sha256(outer_path),
                "outer_optimizer_raw": _sha256(raw_outer),
                "iteration_trace_v2": _sha256(trace_path),
                "sealed_artifact": _sha256(sealed_manifest),
            }
            if raw_inner.is_file():
                hashes["inner_optimizer_raw"] = _sha256(raw_inner)
            inner_path = run_root / "inner.obj"
            if native_result.inner_faces.size:
                native_result.save_inner_mesh(inner_path)
                hashes["inner"] = _sha256(inner_path)
            summary = {
                "completion": "stable" if result.succeeded else "not_converged",
                "elapsed_seconds": elapsed,
                "status": result.status,
                "succeeded": result.succeeded,
                "output_sha256": hashes,
                "paper_outer_iterations": len({
                    event["paper_outer_iteration"] for event in trace["events"]
                    if event["paper_outer_iteration"] is not None
                }),
                "beta_optimize_calls": len({
                    event["beta_optimize_call"]
                    for event in trace["deformation_trials"]
                    if event["beta_optimize_call"] is not None
                }),
                "committed_events": len(trace["events"]),
                "first_standing_feasible_event": diagnostics["r10_metrics"][
                    "first_standing_feasible_event"
                ],
                "first_safe_event": diagnostics["r10_metrics"]["first_safe_event"],
                "final_projected_target_energy": diagnostics["r10_metrics"][
                    "final_projected_target_energy"
                ],
                "sealed_artifact_status": "SEALED" if sealed.ready else "REJECTED",
                "sealed_artifact_failures": list(sealed.failures),
                "sealed_geometry_sha256": sealed.geometry_sha256,
                "sealed_mass_properties_sha256": sealed.mass_properties_sha256,
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
        except WallClockLimit as error:
            summary = {
                "completion": "timeout",
                "error": str(error),
                "elapsed_seconds": time.monotonic() - started,
            }
            exit_code = 124
        except Exception as error:
            summary = {
                "completion": "error",
                "error_type": type(error).__name__,
                "error": str(error),
                "elapsed_seconds": time.monotonic() - started,
            }
            exit_code = 1
        finally:
            signal.alarm(0)
        _write_json(run_root / "summary.json", summary)
        summaries.append(summary)
        if exit_code:
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
        "case": "toyrocking",
        "profile": arguments.profile,
        "classification": EvidenceClass.RECONSTRUCTED_ASSET.value,
        "freshness": {"status": FreshnessStatus.CURRENT.value, "mismatches": []},
        "execution_fingerprint": fingerprint,
        "surface_provenance": "official_obj_with_reconstructed_voxels_bbw_pose_and_handles",
        "parameter_provenance": {
            "shape_scale_lambda": (
                "profile_default" if arguments.shape_scale_lambda is None
                else "explicit_reconstructed_experiment_not_archived_by_authors"
            )
        },
        "voxel_resolution": arguments.voxel_resolution,
        "input_sha256": {
            "official_open_obj": _sha256(mesh),
            "reconstructed_closed_obj": _sha256(reconstructed_input),
        },
        "reconstruction": {
            **reconstruction,
            "selected_handle_vertex_indices": list(selected_handles or ()),
        },
        "configuration": asdict(config),
        "limits": {
            "memory_gib": arguments.memory_limit_gib,
            "wall_clock_seconds": arguments.time_limit_seconds,
        },
        "repeat_count_requested": arguments.repeat,
        "deterministic_output_hashes": deterministic,
        "runs": summaries,
        "claim_boundary": (
            "The author cache and interaction history are unavailable. This is "
            "a deterministic reconstructed asset and can never be EXACT_CACHE."
        ),
    }
    _write_json(output / "manifest.json", manifest)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
