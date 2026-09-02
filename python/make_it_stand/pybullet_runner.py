#!/usr/bin/env python3
"""Deterministic PyBullet standing tests for the ShadowArt-10 meshes.

The default discovery mode consumes the closed 64^3 physics proxies produced
by the R11 audit.  A custom JSON manifest can point at later Make-It-Stand
outputs without changing this runner.  Dynamic triangle meshes are represented
by cached V-HACD convex decompositions because Bullet does not support a
general concave triangle mesh as a moving rigid body.

Examples
--------
Headless ten-case regression with TinyRenderer snapshots::

    make-it-stand-pybullet batch --manifest optimized_cases.json

Replay one result in the interactive PyBullet GUI::

    make-it-stand-pybullet replay \
      --manifest optimized_cases.json \
      --case fig1_A_gear

Test optimized meshes by supplying a manifest (schema documented in README)::

    make-it-stand-pybullet batch \
      --manifest /path/to/optimized_cases.json --output results/optimized
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


SCRIPT = Path(__file__).resolve()
REPOSITORY_ROOT = SCRIPT.parents[2]
DEFAULT_R11_ROOT = REPOSITORY_ROOT / "experiments" / "r11" / "shadowart10"
DEFAULT_OUTPUT = Path.cwd() / "make-it-stand-pybullet-results"


@dataclass(frozen=True)
class Case:
    name: str
    collision_mesh: Path
    visual_mesh: Path
    expected_stable: bool | None = None
    expected_margin_cells: float | None = None
    provenance: str = "custom_manifest"
    inner_void_mesh: Path | None = None
    center_of_mass_native: tuple[float, float, float] | None = None
    fixture_mesh: Path | None = None
    occupancy_hash: str | None = None
    collision_mesh_sha256: str | None = None
    inner_void_mesh_sha256: str | None = None
    geometry_hash: str | None = None
    material_mass_native: float | None = None
    inertia_tensor_about_center_native: (
        tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]
        | None
    ) = None
    mass_properties_sha256: str | None = None
    gravity_direction_native: tuple[float, float, float] = (0.0, 0.0, -1.0)


@dataclass(frozen=True)
class SimulationConfig:
    scene: str = "ground"
    preset: str = "default"
    target_extent_m: float = 0.10
    table_width_m: float = 0.50
    table_depth_m: float = 0.38
    table_height_m: float = 0.72
    table_thickness_m: float = 0.035
    camera_distance_m: float = 0.50
    camera_yaw_deg: float = 38.0
    camera_pitch_deg: float = -12.0
    mass_kg: float = 0.20
    initial_tilt_deg: float = 0.0
    initial_tilt_azimuth_deg: float = 0.0
    perturbation_protocol: str = "reseeded_initial_pose_drop"
    duration_s: float = 8.0
    time_step_s: float = 1.0 / 240.0
    sample_hz: float = 30.0
    video_fps: float = 12.0
    gravity_m_s2: float = 9.81
    object_friction: float = 0.80
    ground_friction: float = 0.80
    rolling_friction: float = 0.001
    spinning_friction: float = 0.001
    restitution: float = 0.0
    linear_damping: float = 0.04
    angular_damping: float = 0.04
    solver_iterations: int = 150
    initial_clearance_m: float = 0.0002
    stable_tilt_deg: float = 5.0
    maximum_excursion_deg: float = 15.0
    fallen_tilt_deg: float = 45.0
    maximum_drift_m: float = 0.01
    settled_linear_speed_m_s: float = 0.005
    settled_angular_speed_rad_s: float = 0.10
    settle_window_s: float = 0.50
    vhacd_resolution: int = 100_000
    vhacd_concavity: float = 0.0025
    vhacd_max_vertices: int = 64
    vhacd_timeout_s: float = 180.0


@dataclass(frozen=True)
class SupportScene:
    body: int
    height_m: float
    camera_target_z_m: float


VHACD_PRESETS: dict[str, dict[str, float | int]] = {
    "default": {"resolution": 100_000, "concavity": 0.0025},
    "robust": {"resolution": 1_000_000, "concavity": 0.001},
}


def _require_modules() -> tuple[Any, Any]:
    try:
        import pybullet as pybullet
    except ImportError as error:
        raise SystemExit(
            "PyBullet is required. Install it in the active environment with "
            "`python -m pip install pybullet`."
        ) from error
    try:
        import trimesh
    except ImportError as error:
        raise SystemExit(
            "trimesh is required. Install it with `python -m pip install trimesh`."
        ) from error
    return pybullet, trimesh


def _json_dump(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _resolve_manifest_path(value: str, base: Path) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def discover_r11_cases(root: Path = DEFAULT_R11_ROOT) -> list[Case]:
    cases: list[Case] = []
    for audit_path in sorted(root.glob("*/audit.json")):
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        mesh = audit_path.parent / "physics_proxy64.obj"
        if not mesh.is_file():
            continue
        canonical = audit.get("canonical_z_min_support", {})
        cases.append(
            Case(
                name=audit_path.parent.name,
                collision_mesh=mesh.resolve(),
                visual_mesh=mesh.resolve(),
                expected_stable=canonical.get("stable"),
                expected_margin_cells=canonical.get("signed_margin"),
                provenance="r11_closed_proxy64_fixed_z_min",
            )
        )
    if len(cases) != 10:
        raise RuntimeError(
            f"Expected 10 R11 cases below {root}, discovered {len(cases)}. "
            "Run experiments/r11/run_shadowart10_batch.py first or pass --manifest."
        )
    return cases


def load_manifest(path: Path) -> list[Case]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    entries = payload["cases"] if isinstance(payload, dict) else payload
    base = path.resolve().parent
    cases: list[Case] = []
    for entry in entries:
        collision = _resolve_manifest_path(entry["collision_mesh"], base)
        visual = _resolve_manifest_path(entry.get("visual_mesh", str(collision)), base)
        inner_void = (
            _resolve_manifest_path(entry["inner_void_mesh"], base)
            if entry.get("inner_void_mesh")
            else None
        )
        center_override = entry.get("center_of_mass_native")
        if center_override is not None and len(center_override) != 3:
            raise ValueError(
                f"center_of_mass_native for {entry['name']} must contain three values"
            )
        inertia_override = entry.get("inertia_tensor_about_center_native")
        if inertia_override is not None:
            inertia_array = np.asarray(inertia_override, dtype=float)
            if inertia_array.shape != (3, 3) or not np.all(np.isfinite(inertia_array)):
                raise ValueError(
                    f"inertia_tensor_about_center_native for {entry['name']} "
                    "must be a finite 3x3 matrix"
                )
        gravity_direction = entry.get("gravity_direction_native", [0.0, 0.0, -1.0])
        gravity_array = np.asarray(gravity_direction, dtype=float)
        if (
            gravity_array.shape != (3,)
            or not np.all(np.isfinite(gravity_array))
            or float(np.linalg.norm(gravity_array)) <= 1.0e-12
        ):
            raise ValueError(
                f"gravity_direction_native for {entry['name']} must be a "
                "finite nonzero 3-vector"
            )
        gravity_array /= np.linalg.norm(gravity_array)
        analytic_fields = (
            entry.get("material_mass_native"),
            center_override,
            inertia_override,
            entry.get("mass_properties_sha256"),
        )
        if any(value is not None for value in analytic_fields):
            if any(value is None for value in analytic_fields):
                raise ValueError(
                    f"Analytic mass properties for {entry['name']} are incomplete"
                )
            if entry.get("mass_properties_sha256") != _mass_properties_fingerprint(entry):
                raise ValueError(f"Analytic mass-properties hash mismatch for {entry['name']}")
        cases.append(
            Case(
                name=str(entry["name"]),
                collision_mesh=collision,
                visual_mesh=visual,
                expected_stable=entry.get("expected_stable"),
                expected_margin_cells=entry.get("expected_margin_cells"),
                provenance=str(entry.get("provenance", path.name)),
                inner_void_mesh=inner_void,
                center_of_mass_native=(
                    tuple(float(value) for value in center_override)
                    if center_override is not None
                    else None
                ),
                fixture_mesh=(
                    _resolve_manifest_path(entry["fixture_mesh"], base)
                    if entry.get("fixture_mesh")
                    else None
                ),
                occupancy_hash=(
                    str(entry["occupancy_hash"])
                    if entry.get("occupancy_hash") is not None
                    else None
                ),
                collision_mesh_sha256=(
                    str(entry["collision_mesh_sha256"])
                    if entry.get("collision_mesh_sha256") is not None
                    else None
                ),
                inner_void_mesh_sha256=(
                    str(entry["inner_void_mesh_sha256"])
                    if entry.get("inner_void_mesh_sha256") is not None
                    else None
                ),
                geometry_hash=(
                    str(entry["geometry_hash"])
                    if entry.get("geometry_hash") is not None
                    else None
                ),
                material_mass_native=(
                    float(entry["material_mass_native"])
                    if entry.get("material_mass_native") is not None
                    else None
                ),
                inertia_tensor_about_center_native=(
                    tuple(tuple(float(value) for value in row) for row in inertia_override)
                    if inertia_override is not None
                    else None
                ),
                mass_properties_sha256=(
                    str(entry["mass_properties_sha256"])
                    if entry.get("mass_properties_sha256") is not None
                    else None
                ),
                gravity_direction_native=tuple(float(value) for value in gravity_array),
            )
        )
    if not cases:
        raise ValueError(f"Manifest {path} contains no cases")
    names = [case.name for case in cases]
    if len(names) != len(set(names)):
        raise ValueError(f"Manifest {path} contains duplicate case names")
    for case in cases:
        for mesh in (
            case.collision_mesh,
            case.visual_mesh,
            case.inner_void_mesh,
            case.fixture_mesh,
        ):
            if mesh is None:
                continue
            if not mesh.is_file():
                raise FileNotFoundError(f"Missing mesh for {case.name}: {mesh}")
        if case.occupancy_hash is not None:
            analytic = case.material_mass_native is not None
            if (
                not case.occupancy_hash or
                case.collision_mesh_sha256 is None or
                (
                    analytic and
                    (
                        case.geometry_hash is None or
                        case.center_of_mass_native is None or
                        case.inertia_tensor_about_center_native is None or
                        case.mass_properties_sha256 is None
                    )
                ) or
                (
                    not analytic and
                    (
                        case.inner_void_mesh is None or
                        case.inner_void_mesh_sha256 is None
                    )
                )
            ):
                raise ValueError(
                    f"Occupancy-bound case {case.name} has incomplete provenance"
                )
            if _file_sha256(case.collision_mesh) != case.collision_mesh_sha256:
                raise ValueError(f"Outer mesh hash mismatch for {case.name}")
            if (
                case.inner_void_mesh is not None and
                case.inner_void_mesh_sha256 is not None and
                _file_sha256(case.inner_void_mesh) != case.inner_void_mesh_sha256
            ):
                raise ValueError(f"Inner mesh hash mismatch for {case.name}")
            if not analytic and case.center_of_mass_native is not None:
                raise ValueError(
                    f"Occupancy-bound case {case.name} cannot override its CoM"
                )
    return cases


def _mass_properties_fingerprint(entry: Mapping[str, Any]) -> str:
    """Bind analytic material moments to occupancy, geometry, and collision input."""
    payload = {
        "schema": "make-it-stand-analytic-mass-properties-v1",
        "occupancy_hash": str(entry["occupancy_hash"]),
        "geometry_hash": str(entry["geometry_hash"]),
        "collision_mesh_sha256": str(entry["collision_mesh_sha256"]),
        "material_mass_native": float(entry["material_mass_native"]),
        "center_of_mass_native": [
            float(value) for value in entry["center_of_mass_native"]
        ],
        "inertia_tensor_about_center_native": [
            [float(value) for value in row]
            for row in entry["inertia_tensor_about_center_native"]
        ],
        "gravity_direction_native": [
            float(value)
            for value in entry.get("gravity_direction_native", [0.0, 0.0, -1.0])
        ],
    }
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _mesh_fingerprint(path: Path, config: SimulationConfig) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    digest.update(
        json.dumps(
            {
                "resolution": config.vhacd_resolution,
                "concavity": config.vhacd_concavity,
                "max_vertices": config.vhacd_max_vertices,
            },
            sort_keys=True,
        ).encode("utf-8")
    )
    return digest.hexdigest()


def ensure_vhacd(
    case: Case, cache_root: Path, config: SimulationConfig, force: bool = False
) -> tuple[Path, dict[str, Any]]:
    """Create/reuse a V-HACD mesh in a crash-isolated subprocess."""
    case_root = cache_root / case.name
    case_root.mkdir(parents=True, exist_ok=True)
    output = case_root / "collision_vhacd.obj"
    metadata_path = case_root / "collision_vhacd.json"
    fingerprint = _mesh_fingerprint(case.collision_mesh, config)
    if not force and output.is_file() and metadata_path.is_file():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("fingerprint") == fingerprint and output.stat().st_size > 0:
            return output, metadata

    vhacd_input = case.collision_mesh
    if vhacd_input.suffix.lower() != ".obj":
        _, trimesh = _require_modules()
        geometry = trimesh.load(vhacd_input, force="mesh", process=False)
        if isinstance(geometry, trimesh.Scene):
            geometry = geometry.to_geometry()
        if geometry.is_empty:
            raise ValueError(f"Empty V-HACD source mesh: {vhacd_input}")
        vhacd_input = case_root / "source_mesh.obj"
        geometry.export(vhacd_input, file_type="obj")

    process_log = case_root / "vhacd_process.log"
    command = [
        sys.executable,
        str(SCRIPT),
        "vhacd",
        "--input",
        str(vhacd_input),
        "--output",
        str(output),
        "--log",
        str(case_root / "vhacd.log"),
        "--resolution",
        str(config.vhacd_resolution),
        "--concavity",
        str(config.vhacd_concavity),
        "--max-vertices",
        str(config.vhacd_max_vertices),
    ]
    started = time.monotonic()
    with process_log.open("w", encoding="utf-8") as stream:
        completed = subprocess.run(
            command,
            stdout=stream,
            stderr=subprocess.STDOUT,
            timeout=config.vhacd_timeout_s,
            check=False,
        )
    elapsed = time.monotonic() - started
    if completed.returncode != 0 or not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(
            f"V-HACD failed for {case.name} with exit {completed.returncode}; "
            f"see {process_log}"
        )
    text = output.read_text(encoding="utf-8", errors="ignore")
    metadata = {
        "schema": "make-it-stand-pybullet-vhacd-v1",
        "case": case.name,
        "source_mesh": str(case.collision_mesh),
        "vhacd_input_mesh": str(vhacd_input),
        "fingerprint": fingerprint,
        "elapsed_seconds": elapsed,
        "output_bytes": output.stat().st_size,
        "object_groups": sum(
            1 for line in text.splitlines() if line.startswith(("o ", "g "))
        ),
        "parameters": {
            "resolution": config.vhacd_resolution,
            "concavity": config.vhacd_concavity,
            "max_vertices_per_convex_hull": config.vhacd_max_vertices,
        },
    }
    _json_dump(metadata_path, metadata)
    return output, metadata


def _translated_inertia(
    inertia_at_center: np.ndarray,
    mass: float,
    source_center: np.ndarray,
    target_center: np.ndarray,
) -> np.ndarray:
    displacement = source_center - target_center
    return inertia_at_center + mass * (
        float(displacement @ displacement) * np.eye(3)
        - np.outer(displacement, displacement)
    )


def _rotation_matrix_to_quaternion_xyzw(rotation: np.ndarray) -> np.ndarray:
    """Convert a proper 3x3 rotation matrix to an xyzw quaternion."""
    matrix = np.eye(4)
    matrix[:3, :3] = rotation
    # trimesh uses the Hamilton wxyz convention.
    _, trimesh = _require_modules()
    quaternion_wxyz = np.asarray(
        trimesh.transformations.quaternion_from_matrix(matrix), dtype=float
    )
    quaternion_xyzw = quaternion_wxyz[[1, 2, 3, 0]]
    if quaternion_xyzw[3] < 0.0:
        quaternion_xyzw *= -1.0
    return quaternion_xyzw


def _principal_inertia(
    inertia_native: np.ndarray,
    material_volume: float,
    mass_kg: float,
    scale: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return principal moments, principal-to-object rotation, and quaternion."""
    symmetric = 0.5 * (inertia_native + inertia_native.T)
    moments_native, principal_to_object = np.linalg.eigh(symmetric)
    if not np.all(np.isfinite(moments_native)) or np.any(moments_native <= 0.0):
        raise ValueError(
            "Mass properties produced a non-positive or non-finite principal "
            "moment; refusing to replace it with an invented inertia"
        )
    if moments_native[0] + moments_native[1] < moments_native[2] - (
        1.0e-10 * moments_native[2]
    ):
        raise ValueError(
            "Mass properties violate the principal-moment triangle inequality"
        )
    if np.linalg.det(principal_to_object) < 0.0:
        principal_to_object[:, -1] *= -1.0
    moments = moments_native / material_volume * mass_kg * scale * scale
    quaternion = _rotation_matrix_to_quaternion_xyzw(principal_to_object)
    return moments, principal_to_object, quaternion


def _load_mesh_properties(
    case: Case,
    target_extent_m: float,
    mass_kg: float,
    *,
    scale_override: float | None = None,
) -> dict[str, Any]:
    _, trimesh = _require_modules()
    path = case.collision_mesh
    loaded = trimesh.load(path, force="mesh", process=False)
    if isinstance(loaded, trimesh.Scene):
        loaded = loaded.to_geometry()
    if loaded.is_empty:
        raise ValueError(f"Empty mesh: {path}")
    bounds = np.asarray(loaded.bounds, dtype=float)
    extents = bounds[1] - bounds[0]
    maximum_extent = float(np.max(extents))
    if not np.isfinite(maximum_extent) or maximum_extent <= 0.0:
        raise ValueError(f"Degenerate mesh bounds: {path}")
    if scale_override is not None and (
        not math.isfinite(scale_override) or scale_override <= 0.0
    ):
        raise ValueError("scale_override must be finite and positive")
    scale = (
        float(scale_override)
        if scale_override is not None
        else target_extent_m / maximum_extent
    )
    if case.material_mass_native is not None:
        if (
            not math.isfinite(case.material_mass_native) or
            case.material_mass_native <= 0.0 or
            case.center_of_mass_native is None or
            case.inertia_tensor_about_center_native is None
        ):
            raise ValueError(
                f"Analytic mass properties are invalid for {case.name}"
            )
        center_mass = np.asarray(case.center_of_mass_native, dtype=float)
        inertia_native = np.asarray(
            case.inertia_tensor_about_center_native,
            dtype=float,
        )
        inertia_diagonal, principal_to_object, principal_quaternion = (
            _principal_inertia(
                inertia_native,
                case.material_mass_native,
                mass_kg,
                scale,
            )
        )
        inertia_tensor = (
            0.5 * (inertia_native + inertia_native.T)
            / case.material_mass_native
            * mass_kg
            * scale
            * scale
        )
        return {
            "bounds_native": bounds,
            "extent_native": extents,
            "scale": scale,
            "center_mass_native": center_mass,
            "mass_properties_source": (
                "occupancy_hash_bound_analytic_volume_moments"
            ),
            "material_volume_native": case.material_mass_native,
            "local_inertia_diagonal_kg_m2": inertia_diagonal,
            "inertia_tensor_object_frame_kg_m2": inertia_tensor,
            "principal_to_object_rotation": principal_to_object,
            "local_inertial_frame_orientation_xyzw": principal_quaternion,
            "inner_void": None,
            "occupancy_hash": case.occupancy_hash,
            "watertight": bool(loaded.is_watertight),
            "winding_consistent": bool(loaded.is_winding_consistent),
            "vertices": int(len(loaded.vertices)),
            "faces": int(len(loaded.faces)),
            "volume_native": abs(float(loaded.volume)),
        }
    outer_center = np.asarray(loaded.center_mass, dtype=float)
    if outer_center.shape != (3,) or not np.all(np.isfinite(outer_center)):
        outer_center = np.asarray(loaded.centroid, dtype=float)
    outer_signed_volume = float(loaded.volume)
    outer_volume = abs(outer_signed_volume)
    outer_inertia = np.asarray(loaded.moment_inertia, dtype=float) * (
        1.0 if outer_signed_volume >= 0.0 else -1.0
    )
    center_mass = outer_center
    material_volume = outer_volume
    inertia_native = outer_inertia
    mass_properties_source = "closed_collision_mesh"
    inner_report: dict[str, Any] | None = None
    if case.occupancy_hash is not None and (
        not loaded.is_watertight or
        not loaded.is_winding_consistent or
        not outer_signed_volume > 0.0
    ):
        raise ValueError(
            f"Occupancy-bound case {case.name} has no closed material boundary"
        )
    if case.inner_void_mesh is not None:
        inner = trimesh.load(case.inner_void_mesh, force="mesh", process=False)
        if isinstance(inner, trimesh.Scene):
            inner = inner.to_geometry()
        inner_signed_volume = float(inner.volume)
        if case.occupancy_hash is not None and (
            not inner.is_watertight or
            not inner.is_winding_consistent or
            not inner_signed_volume < 0.0
        ):
            raise ValueError(
                f"Occupancy-bound case {case.name} has no closed material boundary"
            )
        inner_volume = abs(inner_signed_volume)
        inner_center = np.asarray(inner.center_mass, dtype=float)
        material_volume = outer_volume - inner_volume
        if material_volume <= max(1.0e-12, outer_volume * 1.0e-9):
            raise ValueError(
                f"Inner void volume is not smaller than outer volume for {case.name}"
            )
        center_mass = (
            outer_volume * outer_center - inner_volume * inner_center
        ) / material_volume
        inertia_native = _translated_inertia(
            outer_inertia, outer_volume, outer_center, center_mass
        ) - _translated_inertia(
            np.asarray(inner.moment_inertia, dtype=float)
            * (1.0 if inner_signed_volume >= 0.0 else -1.0),
            inner_volume,
            inner_center,
            center_mass,
        )
        mass_properties_source = "outer_mesh_minus_inner_void_mesh"
        inner_report = {
            "path": str(case.inner_void_mesh),
            "volume_native": inner_volume,
            "center_of_mass_native": inner_center.tolist(),
            "watertight": bool(inner.is_watertight),
        }
        if case.occupancy_hash is not None:
            mass_properties_source = (
                "occupancy_hash_bound_closed_outer_minus_inner"
            )
    if case.center_of_mass_native is not None:
        center_override = np.asarray(case.center_of_mass_native, dtype=float)
        inertia_native = _translated_inertia(
            inertia_native, material_volume, center_mass, center_override
        )
        center_mass = center_override
        mass_properties_source += "+manifest_com_override"
    inertia_diagonal, principal_to_object, principal_quaternion = _principal_inertia(
        inertia_native, material_volume, mass_kg, scale
    )
    inertia_tensor = (
        0.5
        * (inertia_native + inertia_native.T)
        / material_volume
        * mass_kg
        * scale
        * scale
    )
    return {
        "bounds_native": bounds,
        "extent_native": extents,
        "scale": scale,
        "center_mass_native": center_mass,
        "mass_properties_source": mass_properties_source,
        "material_volume_native": material_volume,
        "local_inertia_diagonal_kg_m2": inertia_diagonal,
        "inertia_tensor_object_frame_kg_m2": inertia_tensor,
        "principal_to_object_rotation": principal_to_object,
        "local_inertial_frame_orientation_xyzw": principal_quaternion,
        "inner_void": inner_report,
        "occupancy_hash": case.occupancy_hash,
        "watertight": bool(loaded.is_watertight),
        "winding_consistent": bool(loaded.is_winding_consistent),
        "vertices": int(len(loaded.vertices)),
        "faces": int(len(loaded.faces)),
        "volume_native": outer_volume,
    }


def _quat_angle_deg(quaternion: Sequence[float]) -> float:
    scalar = min(1.0, max(-1.0, abs(float(quaternion[3]))))
    return math.degrees(2.0 * math.acos(scalar))


def _tilt_deg(pybullet: Any, quaternion: Sequence[float]) -> float:
    rotation = np.asarray(pybullet.getMatrixFromQuaternion(quaternion)).reshape(3, 3)
    cosine = min(1.0, max(-1.0, float(rotation[2, 2])))
    return math.degrees(math.acos(cosine))


def _alignment_rotation(gravity_direction_native: Sequence[float]) -> np.ndarray:
    """Return a proper rotation mapping native gravity to world -Z."""
    gravity = np.asarray(gravity_direction_native, dtype=float)
    if gravity.shape != (3,) or not np.all(np.isfinite(gravity)):
        raise ValueError("gravity_direction_native must be a finite 3-vector")
    norm = float(np.linalg.norm(gravity))
    if norm <= 1.0e-12:
        raise ValueError("gravity_direction_native must be nonzero")
    source = -gravity / norm
    target = np.asarray([0.0, 0.0, 1.0])
    cosine = float(np.clip(np.dot(source, target), -1.0, 1.0))
    if cosine >= 1.0 - 1.0e-14:
        return np.eye(3)
    if cosine <= -1.0 + 1.0e-14:
        basis = np.zeros(3)
        basis[int(np.argmin(np.abs(source)))] = 1.0
        axis = np.cross(source, basis)
        axis /= np.linalg.norm(axis)
        return 2.0 * np.outer(axis, axis) - np.eye(3)
    cross = np.cross(source, target)
    skew = np.asarray([
        [0.0, -cross[2], cross[1]],
        [cross[2], 0.0, -cross[0]],
        [-cross[1], cross[0], 0.0],
    ])
    return np.eye(3) + skew + (skew @ skew) / (1.0 + cosine)


def _matrix_to_quaternion_xyzw(rotation: np.ndarray) -> list[float]:
    """Convert a proper 3x3 rotation matrix to a normalized xyzw quaternion."""
    matrix = np.asarray(rotation, dtype=float)
    trace = float(np.trace(matrix))
    if trace > 0.0:
        scale = math.sqrt(trace + 1.0) * 2.0
        quaternion = np.asarray([
            (matrix[2, 1] - matrix[1, 2]) / scale,
            (matrix[0, 2] - matrix[2, 0]) / scale,
            (matrix[1, 0] - matrix[0, 1]) / scale,
            0.25 * scale,
        ])
    else:
        index = int(np.argmax(np.diag(matrix)))
        j = (index + 1) % 3
        k = (index + 2) % 3
        scale = math.sqrt(
            max(0.0, 1.0 + matrix[index, index] - matrix[j, j] - matrix[k, k])
        ) * 2.0
        quaternion = np.zeros(4)
        quaternion[index] = 0.25 * scale
        quaternion[j] = (matrix[j, index] + matrix[index, j]) / scale
        quaternion[k] = (matrix[k, index] + matrix[index, k]) / scale
        quaternion[3] = (matrix[k, j] - matrix[j, k]) / scale
    quaternion /= np.linalg.norm(quaternion)
    return [float(value) for value in quaternion]


def _tilt_from_native_gravity_deg(
    object_to_world_rotation: np.ndarray,
    gravity_direction_native: Sequence[float],
) -> float:
    gravity = np.asarray(gravity_direction_native, dtype=float)
    native_up = -gravity / np.linalg.norm(gravity)
    world_up = np.asarray(object_to_world_rotation, dtype=float) @ native_up
    cosine = float(np.clip(world_up[2] / np.linalg.norm(world_up), -1.0, 1.0))
    return math.degrees(math.acos(cosine))


def _camera_rgb(
    pybullet: Any,
    config: SimulationConfig,
    width: int = 720,
    height: int = 544,
    target_z: float = 0.045,
) -> np.ndarray:
    view = pybullet.computeViewMatrixFromYawPitchRoll(
        cameraTargetPosition=[0.0, 0.0, target_z],
        distance=config.camera_distance_m,
        yaw=config.camera_yaw_deg,
        pitch=config.camera_pitch_deg,
        roll=0.0,
        upAxisIndex=2,
    )
    projection = pybullet.computeProjectionMatrixFOV(
        fov=39.0, aspect=width / height, nearVal=0.01, farVal=3.0
    )
    image = pybullet.getCameraImage(
        width,
        height,
        viewMatrix=view,
        projectionMatrix=projection,
        renderer=pybullet.ER_TINY_RENDERER,
        shadow=1,
        lightDirection=[-1.8, -2.5, -4.0],
        lightColor=[1.0, 0.96, 0.90],
        lightDistance=2.0,
        lightAmbientCoeff=0.55,
        lightDiffuseCoeff=0.75,
        lightSpecularCoeff=0.12,
    )
    rgba = np.asarray(image[2], dtype=np.uint8).reshape(height, width, 4)
    return rgba[:, :, :3]


def _render(
    pybullet: Any,
    config: SimulationConfig,
    output: Path,
    target_z: float = 0.045,
) -> None:
    try:
        from PIL import Image
    except ImportError:
        return
    rgb = _camera_rgb(pybullet, config, target_z=target_z)
    output.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(rgb, mode="RGB").save(output)


def _display_font(size: int) -> Any:
    from PIL import ImageFont

    for candidate in (
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf",
    ):
        if Path(candidate).is_file():
            return ImageFont.truetype(candidate, size=size)
    return ImageFont.load_default()


def _video_frame(
    pybullet: Any,
    label: str,
    simulation_time: float,
    target_z: float,
    config: SimulationConfig,
) -> np.ndarray:
    from PIL import Image, ImageDraw

    frame = Image.fromarray(
        _camera_rgb(pybullet, config, target_z=target_z), mode="RGB"
    )
    draw = ImageDraw.Draw(frame)
    font = _display_font(20)
    text = f"{label}   t={simulation_time:05.2f} s"
    box = draw.textbbox((0, 0), text, font=font)
    draw.rectangle((8, 8, box[2] + 20, box[3] + 18), fill=(255, 255, 255))
    draw.text((14, 12), text, fill=(20, 20, 20), font=font)
    return np.asarray(frame)


def _static_box(
    pybullet: Any,
    client: int,
    half_extents: Sequence[float],
    position: Sequence[float],
    color: Sequence[float],
) -> int:
    collision = pybullet.createCollisionShape(
        pybullet.GEOM_BOX,
        halfExtents=half_extents,
        physicsClientId=client,
    )
    visual = pybullet.createVisualShape(
        pybullet.GEOM_BOX,
        halfExtents=half_extents,
        rgbaColor=color,
        specularColor=[0.10, 0.08, 0.06],
        physicsClientId=client,
    )
    return pybullet.createMultiBody(
        baseMass=0.0,
        baseCollisionShapeIndex=collision,
        baseVisualShapeIndex=visual,
        basePosition=position,
        physicsClientId=client,
    )


def _configure_world(
    pybullet: Any, client: int, config: SimulationConfig
) -> SupportScene:
    pybullet.resetSimulation(physicsClientId=client)
    pybullet.setPhysicsEngineParameter(
        fixedTimeStep=config.time_step_s,
        numSolverIterations=config.solver_iterations,
        deterministicOverlappingPairs=1,
        enableConeFriction=1,
        physicsClientId=client,
    )
    pybullet.setTimeStep(config.time_step_s, physicsClientId=client)
    pybullet.setGravity(0.0, 0.0, -config.gravity_m_s2, physicsClientId=client)
    ground = _static_box(
        pybullet,
        client,
        [0.65, 0.65, 0.012],
        [0.0, 0.0, -0.012],
        [0.86, 0.86, 0.84, 1.0],
    )
    pybullet.changeDynamics(
        ground,
        -1,
        lateralFriction=config.ground_friction,
        restitution=config.restitution,
        physicsClientId=client,
    )
    if config.scene == "ground":
        return SupportScene(body=ground, height_m=0.0, camera_target_z_m=0.05)

    tabletop_height = config.table_height_m
    tabletop_thickness = config.table_thickness_m
    table_width = config.table_width_m
    table_depth = config.table_depth_m
    tabletop = _static_box(
        pybullet,
        client,
        [table_width / 2.0, table_depth / 2.0, tabletop_thickness / 2.0],
        [0.0, 0.0, tabletop_height - tabletop_thickness / 2.0],
        [0.47, 0.25, 0.105, 1.0],
    )
    leg_width = min(0.04, 0.09 * min(table_width, table_depth))
    leg_half = [
        leg_width / 2.0,
        leg_width / 2.0,
        (tabletop_height - tabletop_thickness) / 2.0,
    ]
    leg_z = leg_half[2]
    leg_x = table_width / 2.0 - 1.35 * leg_width
    leg_y = table_depth / 2.0 - 1.35 * leg_width
    for x in (-leg_x, leg_x):
        for y in (-leg_y, leg_y):
            leg = _static_box(
                pybullet,
                client,
                leg_half,
                [x, y, leg_z],
                [0.34, 0.17, 0.07, 1.0],
            )
            pybullet.changeDynamics(
                leg,
                -1,
                lateralFriction=config.ground_friction,
                restitution=config.restitution,
                physicsClientId=client,
            )
    pybullet.changeDynamics(
        tabletop,
        -1,
        lateralFriction=config.ground_friction,
        restitution=config.restitution,
        physicsClientId=client,
    )
    return SupportScene(
        body=tabletop,
        height_m=tabletop_height,
        camera_target_z_m=tabletop_height + 0.18 * config.target_extent_m,
    )


def _create_body(
    pybullet: Any,
    client: int,
    case: Case,
    vhacd_mesh: Path,
    properties: dict[str, Any],
    config: SimulationConfig,
    support_height_m: float,
    horizontal_offset_m: Sequence[float] = (0.0, 0.0),
    visual_color: Sequence[float] = (0.76, 0.53, 0.22, 1.0),
    world_translation_m: Sequence[float] | None = None,
    body_mass_kg: float | None = None,
) -> tuple[int, list[float]]:
    _, trimesh = _require_modules()
    scale = float(properties["scale"])
    center_mass = np.asarray(properties["center_mass_native"], dtype=float)
    collision_geometry = trimesh.load(vhacd_mesh, force="mesh", process=False)
    if isinstance(collision_geometry, trimesh.Scene):
        collision_geometry = collision_geometry.to_geometry()
    if collision_geometry.is_empty:
        raise ValueError(f"Empty V-HACD collision mesh: {vhacd_mesh}")
    collision_minimum = np.asarray(collision_geometry.bounds[0], dtype=float)
    frame_offset = (-center_mass * scale).tolist()
    inertial_orientation = properties["local_inertial_frame_orientation_xyzw"].tolist()
    mesh_scale = [scale, scale, scale]
    collision = pybullet.createCollisionShape(
        pybullet.GEOM_MESH,
        fileName=str(vhacd_mesh),
        meshScale=mesh_scale,
        collisionFramePosition=frame_offset,
        physicsClientId=client,
    )
    visual = pybullet.createVisualShape(
        pybullet.GEOM_MESH,
        fileName=str(case.visual_mesh),
        meshScale=mesh_scale,
        visualFramePosition=frame_offset,
        rgbaColor=visual_color,
        specularColor=[0.16, 0.16, 0.16],
        physicsClientId=client,
    )
    nominal_orientation = _matrix_to_quaternion_xyzw(
        _alignment_rotation(case.gravity_direction_native)
    )
    initial_orientation = nominal_orientation
    if config.initial_tilt_deg != 0.0:
        azimuth = math.radians(config.initial_tilt_azimuth_deg)
        perturbation = list(pybullet.getQuaternionFromAxisAngle(
            [math.cos(azimuth), math.sin(azimuth), 0.0],
            math.radians(config.initial_tilt_deg),
        ))
        _, composed = pybullet.multiplyTransforms(
            [0.0, 0.0, 0.0],
            perturbation,
            [0.0, 0.0, 0.0],
            nominal_orientation,
        )
        initial_orientation = list(composed)
    if world_translation_m is None:
        rotation = np.asarray(
            pybullet.getMatrixFromQuaternion(initial_orientation),
            dtype=float,
        ).reshape(3, 3)
        centered_vertices = (
            np.asarray(collision_geometry.vertices, dtype=float) - center_mass
        ) * scale
        rotated_vertices = centered_vertices @ rotation.T
        lowest_relative_z = float(rotated_vertices[:, 2].min())
        initial_position = [
            float(horizontal_offset_m[0]),
            float(horizontal_offset_m[1]),
            float(support_height_m + config.initial_clearance_m - lowest_relative_z),
        ]
    else:
        if len(world_translation_m) != 3:
            raise ValueError("world_translation_m must contain three values")
        initial_position = [
            float(center_mass[index] * scale + world_translation_m[index])
            for index in range(3)
        ]
    mass = config.mass_kg if body_mass_kg is None else float(body_mass_kg)
    if not math.isfinite(mass) or mass <= 0.0:
        raise ValueError("body_mass_kg must be finite and positive")
    body = pybullet.createMultiBody(
        baseMass=mass,
        baseCollisionShapeIndex=collision,
        baseVisualShapeIndex=visual,
        basePosition=initial_position,
        baseOrientation=initial_orientation,
        baseInertialFramePosition=[0.0, 0.0, 0.0],
        baseInertialFrameOrientation=inertial_orientation,
        physicsClientId=client,
    )
    pybullet.changeDynamics(
        body,
        -1,
        lateralFriction=config.object_friction,
        rollingFriction=config.rolling_friction,
        spinningFriction=config.spinning_friction,
        restitution=config.restitution,
        linearDamping=config.linear_damping,
        angularDamping=config.angular_damping,
        localInertiaDiagonal=properties["local_inertia_diagonal_kg_m2"].tolist(),
        physicsClientId=client,
    )
    properties["collision_bounds_native"] = np.asarray(
        collision_geometry.bounds, dtype=float
    )
    properties["gravity_direction_native"] = np.asarray(
        case.gravity_direction_native, dtype=float
    )
    properties["nominal_orientation_xyzw"] = np.asarray(
        nominal_orientation, dtype=float
    )
    return body, initial_position


def _combined_mass_properties_audit(
    object_properties: dict[str, Any],
    object_mass_kg: float,
    fixture_properties: dict[str, Any],
    fixture_mass_kg: float,
) -> dict[str, Any]:
    object_center = np.asarray(object_properties["center_mass_native"], dtype=float)
    fixture_center = np.asarray(fixture_properties["center_mass_native"], dtype=float)
    total_mass = object_mass_kg + fixture_mass_kg
    combined = (
        object_mass_kg * object_center + fixture_mass_kg * fixture_center
    ) / total_mass
    valid = bool(
        total_mass > 0.0
        and np.all(np.isfinite(combined))
        and object_properties["material_volume_native"] > 0.0
        and fixture_properties["material_volume_native"] > 0.0
    )
    return {
        "valid": valid,
        "object_mass_kg": float(object_mass_kg),
        "fixture_mass_kg": float(fixture_mass_kg),
        "combined_mass_kg": float(total_mass),
        "combined_center_of_mass_native": combined.tolist(),
        "same_uniform_material_density": True,
    }


def _contact_records(contacts: Iterable[Sequence[Any]]) -> list[dict[str, Any]]:
    """Preserve contact geometry needed to audit the effective support hull."""

    records: list[dict[str, Any]] = []
    for contact in contacts:
        if len(contact) < 10:
            raise ValueError("PyBullet contact record is incomplete")
        records.append({
            "position_on_object_world_m": [float(value) for value in contact[5]],
            "position_on_support_world_m": [float(value) for value in contact[6]],
            "support_normal_world": [float(value) for value in contact[7]],
            "contact_distance_m": float(contact[8]),
            "normal_force_n": float(contact[9]),
        })
    return records


def _sample_state(
    pybullet: Any,
    client: int,
    body: int,
    support_body: int,
    simulation_time: float,
    local_inertial_orientation: Sequence[float],
    gravity_direction_native: Sequence[float],
    nominal_orientation_xyzw: Sequence[float],
) -> dict[str, Any]:
    position, inertial_orientation = pybullet.getBasePositionAndOrientation(
        body, physicsClientId=client
    )
    _, inverse_local_orientation = pybullet.invertTransform(
        [0.0, 0.0, 0.0], local_inertial_orientation
    )
    _, object_orientation = pybullet.multiplyTransforms(
        [0.0, 0.0, 0.0],
        inertial_orientation,
        [0.0, 0.0, 0.0],
        inverse_local_orientation,
    )
    _, inverse_nominal = pybullet.invertTransform(
        [0.0, 0.0, 0.0], nominal_orientation_xyzw
    )
    _, relative_orientation = pybullet.multiplyTransforms(
        [0.0, 0.0, 0.0],
        object_orientation,
        [0.0, 0.0, 0.0],
        inverse_nominal,
    )
    object_rotation = np.asarray(
        pybullet.getMatrixFromQuaternion(object_orientation), dtype=float
    ).reshape(3, 3)
    linear_velocity, angular_velocity = pybullet.getBaseVelocity(
        body, physicsClientId=client
    )
    contacts = pybullet.getContactPoints(
        bodyA=body, bodyB=support_body, physicsClientId=client
    )
    return {
        "time_s": simulation_time,
        "position_m": [float(value) for value in position],
        "orientation_xyzw": [float(value) for value in object_orientation],
        "inertial_orientation_xyzw": [float(value) for value in inertial_orientation],
        "tilt_deg": _tilt_from_native_gravity_deg(
            object_rotation, gravity_direction_native
        ),
        "orientation_delta_deg": _quat_angle_deg(relative_orientation),
        "linear_speed_m_s": float(np.linalg.norm(linear_velocity)),
        "angular_speed_rad_s": float(np.linalg.norm(angular_velocity)),
        "ground_contacts": len(contacts),
        "normal_force_n": float(sum(contact[9] for contact in contacts)),
        "contacts": _contact_records(contacts),
    }


def _classify(
    trajectory: list[dict[str, Any]],
    initial_position: Sequence[float],
    config: SimulationConfig,
) -> dict[str, Any]:
    final = trajectory[-1]
    window_start = max(0.0, config.duration_s - config.settle_window_s)
    window = [sample for sample in trajectory if sample["time_s"] >= window_start]
    max_tilt = max(sample["tilt_deg"] for sample in trajectory)
    max_orientation_delta = max(
        sample["orientation_delta_deg"] for sample in trajectory
    )
    maximum_drift = max(
        math.hypot(
            sample["position_m"][0] - initial_position[0],
            sample["position_m"][1] - initial_position[1],
        )
        for sample in trajectory
    )
    mean_linear_speed = float(
        np.mean([sample["linear_speed_m_s"] for sample in window])
    )
    mean_angular_speed = float(
        np.mean([sample["angular_speed_rad_s"] for sample in window])
    )
    contact_fraction = float(
        np.mean([sample["ground_contacts"] > 0 for sample in window])
    )
    settled = (
        mean_linear_speed <= config.settled_linear_speed_m_s
        and mean_angular_speed <= config.settled_angular_speed_rad_s
    )
    fallen = (
        max_tilt >= config.fallen_tilt_deg
        or final["tilt_deg"] >= config.fallen_tilt_deg
    )
    target_pose_retained = (
        final["tilt_deg"] <= config.stable_tilt_deg
        and max_tilt <= config.maximum_excursion_deg
        and maximum_drift <= config.maximum_drift_m
    )
    stable = (
        target_pose_retained and settled and contact_fraction >= 0.80 and not fallen
    )
    reasons: list[str] = []
    if fallen:
        reasons.append("fallen")
    if final["tilt_deg"] > config.stable_tilt_deg:
        reasons.append("final_tilt_exceeds_limit")
    if max_tilt > config.maximum_excursion_deg:
        reasons.append("tilt_excursion_exceeds_limit")
    if maximum_drift > config.maximum_drift_m:
        reasons.append("horizontal_drift_exceeds_limit")
    if not settled:
        reasons.append("not_settled")
    if contact_fraction < 0.80:
        reasons.append("insufficient_final_ground_contact")
    return {
        "stable": stable,
        "fallen": fallen,
        "target_pose_retained": target_pose_retained,
        "failure_reasons": reasons,
        "final_tilt_deg": final["tilt_deg"],
        "maximum_tilt_deg": max_tilt,
        "maximum_orientation_delta_deg": max_orientation_delta,
        "maximum_horizontal_drift_m": maximum_drift,
        "final_com_height_m": final["position_m"][2],
        "settled": settled,
        "settle_window_mean_linear_speed_m_s": mean_linear_speed,
        "settle_window_mean_angular_speed_rad_s": mean_angular_speed,
        "settle_window_contact_fraction": contact_fraction,
    }


def simulate_case(
    pybullet: Any,
    client: int,
    case: Case,
    vhacd_mesh: Path,
    vhacd_metadata: dict[str, Any],
    output_root: Path,
    config: SimulationConfig,
    render: bool,
    record_video: bool = False,
    realtime: bool = False,
) -> dict[str, Any]:
    scene = _configure_world(pybullet, client, config)
    properties = _load_mesh_properties(case, config.target_extent_m, config.mass_kg)
    body, initial_position = _create_body(
        pybullet,
        client,
        case,
        vhacd_mesh,
        properties,
        config,
        scene.height_m,
    )
    local_inertial_orientation = properties[
        "local_inertial_frame_orientation_xyzw"
    ].tolist()
    case_root = output_root / case.name
    case_root.mkdir(parents=True, exist_ok=True)
    if render:
        _render(
            pybullet,
            config,
            case_root / "initial.png",
            target_z=scene.camera_target_z_m,
        )
    video_writer = None
    if record_video:
        try:
            import imageio.v2 as imageio
        except ImportError as error:
            raise RuntimeError(
                "--record-video requires imageio and an ffmpeg backend"
            ) from error
        video_writer = imageio.get_writer(
            case_root / "simulation.mp4",
            fps=config.video_fps,
            codec="libx264",
            quality=7,
            macro_block_size=16,
        )
        video_writer.append_data(
            _video_frame(
                pybullet,
                case.name,
                0.0,
                scene.camera_target_z_m,
                config,
            )
        )
    steps = max(1, int(round(config.duration_s / config.time_step_s)))
    sample_interval = max(1, int(round(1.0 / (config.sample_hz * config.time_step_s))))
    video_interval = max(1, int(round(1.0 / (config.video_fps * config.time_step_s))))
    trajectory: list[dict[str, Any]] = [
        _sample_state(
            pybullet,
            client,
            body,
            scene.body,
            0.0,
            local_inertial_orientation,
            properties["gravity_direction_native"],
            properties["nominal_orientation_xyzw"],
        )
    ]
    try:
        for step in range(1, steps + 1):
            pybullet.stepSimulation(physicsClientId=client)
            if realtime:
                time.sleep(config.time_step_s)
            if step % sample_interval == 0 or step == steps:
                trajectory.append(
                    _sample_state(
                        pybullet,
                        client,
                        body,
                        scene.body,
                        step * config.time_step_s,
                        local_inertial_orientation,
                        properties["gravity_direction_native"],
                        properties["nominal_orientation_xyzw"],
                    )
                )
            if video_writer is not None and (
                step % video_interval == 0 or step == steps
            ):
                video_writer.append_data(
                    _video_frame(
                        pybullet,
                        case.name,
                        step * config.time_step_s,
                        scene.camera_target_z_m,
                        config,
                    )
                )
    finally:
        if video_writer is not None:
            video_writer.close()
    if render:
        _render(
            pybullet,
            config,
            case_root / "final.png",
            target_z=scene.camera_target_z_m,
        )
    classification = _classify(trajectory, initial_position, config)
    result = {
        "schema": "make-it-stand-pybullet-case-v1",
        "perturbation_protocol": config.perturbation_protocol,
        "case": case.name,
        "provenance": case.provenance,
        "collision_mesh": str(case.collision_mesh),
        "visual_mesh": str(case.visual_mesh),
        "expected_static_audit": {
            "stable": case.expected_stable,
            "signed_margin_cells": case.expected_margin_cells,
        },
        "mesh": {
            "vertices": properties["vertices"],
            "faces": properties["faces"],
            "watertight": properties["watertight"],
            "winding_consistent": properties["winding_consistent"],
            "volume_native": properties["volume_native"],
            "bounds_native": properties["bounds_native"].tolist(),
            "collision_bounds_native": properties["collision_bounds_native"].tolist(),
            "center_of_mass_native": properties["center_mass_native"].tolist(),
            "mass_properties_source": properties["mass_properties_source"],
            "occupancy_hash": properties["occupancy_hash"],
            "material_volume_native": properties["material_volume_native"],
            "local_inertia_diagonal_kg_m2": properties[
                "local_inertia_diagonal_kg_m2"
            ].tolist(),
            "inertia_tensor_object_frame_kg_m2": properties[
                "inertia_tensor_object_frame_kg_m2"
            ].tolist(),
            "principal_to_object_rotation": properties[
                "principal_to_object_rotation"
            ].tolist(),
            "local_inertial_frame_orientation_xyzw": properties[
                "local_inertial_frame_orientation_xyzw"
            ].tolist(),
            "gravity_direction_native": properties[
                "gravity_direction_native"
            ].tolist(),
            "nominal_orientation_xyzw": properties[
                "nominal_orientation_xyzw"
            ].tolist(),
            "inner_void": properties["inner_void"],
            "scale_m_per_native_unit": properties["scale"],
            "target_maximum_extent_m": config.target_extent_m,
        },
        "vhacd": vhacd_metadata,
        "initial_position_m": initial_position,
        "support_scene": {
            "type": config.scene,
            "support_height_m": scene.height_m,
            "support_body": "tabletop" if config.scene == "table" else "ground",
        },
        "classification": classification,
        "trajectory": trajectory,
    }
    _json_dump(case_root / "result.json", result)
    return result


def simulate_fixture_case(
    pybullet: Any,
    client: int,
    case: Case,
    object_vhacd_mesh: Path,
    object_vhacd_metadata: dict[str, Any],
    fixture_vhacd_mesh: Path,
    fixture_vhacd_metadata: dict[str, Any],
    output_root: Path,
    config: SimulationConfig,
    render: bool,
    record_video: bool = False,
    realtime: bool = False,
) -> dict[str, Any]:
    """Simulate a fixture delivery as two independent dynamic rigid bodies."""

    if case.fixture_mesh is None:
        raise ValueError("fixture simulation requires case.fixture_mesh")
    scene = _configure_world(pybullet, client, config)
    object_properties = _load_mesh_properties(
        case, config.target_extent_m, config.mass_kg
    )
    fixture_case = Case(
        name=f"{case.name}__fixture",
        collision_mesh=case.fixture_mesh,
        visual_mesh=case.fixture_mesh,
        provenance=case.provenance,
        gravity_direction_native=case.gravity_direction_native,
    )
    fixture_unit_properties = _load_mesh_properties(
        fixture_case,
        config.target_extent_m,
        1.0,
        scale_override=float(object_properties["scale"]),
    )
    fixture_mass_kg = float(
        config.mass_kg
        * fixture_unit_properties["material_volume_native"]
        / object_properties["material_volume_native"]
    )
    fixture_properties = _load_mesh_properties(
        fixture_case,
        config.target_extent_m,
        fixture_mass_kg,
        scale_override=float(object_properties["scale"]),
    )
    mass_audit = _combined_mass_properties_audit(
        object_properties,
        config.mass_kg,
        fixture_properties,
        fixture_mass_kg,
    )
    if not mass_audit["valid"]:
        raise ValueError("combined object-fixture mass properties are invalid")

    scale = float(object_properties["scale"])
    object_bounds = np.asarray(object_properties["bounds_native"], dtype=float)
    fixture_bounds = np.asarray(fixture_properties["bounds_native"], dtype=float)
    assembly_lower = np.minimum(object_bounds[0], fixture_bounds[0])
    assembly_upper = np.maximum(object_bounds[1], fixture_bounds[1])
    horizontal_center = 0.5 * (assembly_lower[:2] + assembly_upper[:2])
    translation = np.asarray(
        [
            -horizontal_center[0] * scale,
            -horizontal_center[1] * scale,
            scene.height_m + config.initial_clearance_m - fixture_bounds[0, 2] * scale,
        ],
        dtype=float,
    )
    fixture_body, fixture_initial = _create_body(
        pybullet,
        client,
        fixture_case,
        fixture_vhacd_mesh,
        fixture_properties,
        config,
        scene.height_m,
        visual_color=(0.24, 0.42, 0.62, 1.0),
        world_translation_m=translation,
        body_mass_kg=fixture_mass_kg,
    )
    object_body, object_initial = _create_body(
        pybullet,
        client,
        case,
        object_vhacd_mesh,
        object_properties,
        config,
        scene.height_m,
        visual_color=(0.76, 0.53, 0.22, 1.0),
        world_translation_m=translation,
        body_mass_kg=config.mass_kg,
    )
    object_local_inertia = object_properties[
        "local_inertial_frame_orientation_xyzw"
    ].tolist()
    fixture_local_inertia = fixture_properties[
        "local_inertial_frame_orientation_xyzw"
    ].tolist()
    case_root = output_root / case.name
    case_root.mkdir(parents=True, exist_ok=True)
    if render:
        _render(
            pybullet,
            config,
            case_root / "initial.png",
            target_z=scene.camera_target_z_m,
        )
    video_writer = None
    if record_video:
        try:
            import imageio.v2 as imageio
        except ImportError as error:
            raise RuntimeError(
                "--record-video requires imageio and an ffmpeg backend"
            ) from error
        video_writer = imageio.get_writer(
            case_root / "simulation.mp4",
            fps=config.video_fps,
            codec="libx264",
            quality=7,
            macro_block_size=16,
        )

    steps = max(1, int(round(config.duration_s / config.time_step_s)))
    sample_interval = max(1, int(round(1.0 / (config.sample_hz * config.time_step_s))))
    video_interval = max(1, int(round(1.0 / (config.video_fps * config.time_step_s))))

    def sample(simulation_time: float) -> tuple[dict[str, Any], dict[str, Any]]:
        return (
            _sample_state(
                pybullet,
                client,
                object_body,
                fixture_body,
                simulation_time,
                object_local_inertia,
                object_properties["gravity_direction_native"],
                object_properties["nominal_orientation_xyzw"],
            ),
            _sample_state(
                pybullet,
                client,
                fixture_body,
                scene.body,
                simulation_time,
                fixture_local_inertia,
                fixture_properties["gravity_direction_native"],
                fixture_properties["nominal_orientation_xyzw"],
            ),
        )

    first_object, first_fixture = sample(0.0)
    object_trajectory = [first_object]
    fixture_trajectory = [first_fixture]
    try:
        for step in range(1, steps + 1):
            pybullet.stepSimulation(physicsClientId=client)
            if realtime:
                time.sleep(config.time_step_s)
            if step % sample_interval == 0 or step == steps:
                object_state, fixture_state = sample(step * config.time_step_s)
                object_trajectory.append(object_state)
                fixture_trajectory.append(fixture_state)
            if video_writer is not None and (
                step % video_interval == 0 or step == steps
            ):
                video_writer.append_data(
                    _video_frame(
                        pybullet,
                        case.name,
                        step * config.time_step_s,
                        scene.camera_target_z_m,
                        config,
                    )
                )
    finally:
        if video_writer is not None:
            video_writer.close()
    if render:
        _render(
            pybullet,
            config,
            case_root / "final.png",
            target_z=scene.camera_target_z_m,
        )

    object_classification = _classify(object_trajectory, object_initial, config)
    fixture_classification = _classify(fixture_trajectory, fixture_initial, config)
    contact_observed = any(
        sample_state["ground_contacts"] > 0 for sample_state in object_trajectory
    )
    stable = bool(
        object_classification["stable"]
        and fixture_classification["stable"]
        and contact_observed
    )
    failure_reasons = [
        f"object:{reason}" for reason in object_classification["failure_reasons"]
    ] + [f"fixture:{reason}" for reason in fixture_classification["failure_reasons"]]
    if not contact_observed:
        failure_reasons.append("object_fixture_contact_not_observed")
    classification = {
        **object_classification,
        "stable": stable,
        "fallen": bool(
            object_classification["fallen"] or fixture_classification["fallen"]
        ),
        "target_pose_retained": bool(
            object_classification["target_pose_retained"]
            and fixture_classification["target_pose_retained"]
        ),
        "settled": bool(
            object_classification["settled"] and fixture_classification["settled"]
        ),
        "failure_reasons": failure_reasons,
    }
    fixture_assembly = {
        "separate_bodies": True,
        "object_fixture_contact_observed": bool(contact_observed),
        "combined_mass_properties_audited": bool(mass_audit["valid"]),
        "object_body_id_distinct_from_fixture": object_body != fixture_body,
        "mass_properties": mass_audit,
        "object_classification": object_classification,
        "fixture_classification": fixture_classification,
    }
    result = {
        "schema": "make-it-stand-pybullet-case-v1",
        "case": case.name,
        "provenance": case.provenance,
        "collision_mesh": str(case.collision_mesh),
        "visual_mesh": str(case.visual_mesh),
        "fixture_mesh": str(case.fixture_mesh),
        "expected_static_audit": {
            "stable": case.expected_stable,
            "signed_margin_cells": case.expected_margin_cells,
        },
        "object_vhacd": object_vhacd_metadata,
        "fixture_vhacd": fixture_vhacd_metadata,
        "initial_position_m": object_initial,
        "fixture_initial_position_m": fixture_initial,
        "support_scene": {
            "type": config.scene,
            "support_height_m": scene.height_m,
            "support_body": "tabletop" if config.scene == "table" else "ground",
        },
        "classification": classification,
        "fixture_assembly": fixture_assembly,
        "trajectory": object_trajectory,
        "fixture_trajectory": fixture_trajectory,
    }
    _json_dump(case_root / "result.json", result)
    return result


def _write_summary(
    output: Path, results: list[dict[str, Any]], config: SimulationConfig
) -> None:
    agreement_values = []
    for result in results:
        expected = result["expected_static_audit"]["stable"]
        if expected is not None:
            agreement_values.append(expected == result["classification"]["stable"])
    summary = {
        "schema": "make-it-stand-pybullet-summary-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "configuration": asdict(config),
        "cases_total": len(results),
        "stable_cases": sum(result["classification"]["stable"] for result in results),
        "fallen_cases": sum(result["classification"]["fallen"] for result in results),
        "static_dynamic_agreement": (
            float(np.mean(agreement_values)) if agreement_values else None
        ),
        "cases": [
            {
                "name": result["case"],
                "expected_stable": result["expected_static_audit"]["stable"],
                "expected_margin_cells": result["expected_static_audit"][
                    "signed_margin_cells"
                ],
                **result["classification"],
                **(
                    {"fixture_assembly": result["fixture_assembly"]}
                    if "fixture_assembly" in result
                    else {}
                ),
            }
            for result in results
        ],
    }
    _json_dump(output / "summary.json", summary)
    fieldnames = [
        "name",
        "expected_stable",
        "expected_margin_cells",
        "stable",
        "fallen",
        "target_pose_retained",
        "final_tilt_deg",
        "maximum_tilt_deg",
        "maximum_horizontal_drift_m",
        "settled",
        "failure_reasons",
    ]
    with (output / "summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for case in summary["cases"]:
            row = {key: case.get(key) for key in fieldnames}
            row["failure_reasons"] = ";".join(case["failure_reasons"])
            writer.writerow(row)
    _make_contact_sheet(output, results)


def _make_sequential_video(output: Path, results: list[dict[str, Any]]) -> None:
    videos = [output / result["case"] / "simulation.mp4" for result in results]
    videos = [path for path in videos if path.is_file()]
    if not videos:
        return
    concat = output / "video_playlist.txt"
    lines = []
    for path in videos:
        escaped = str(path.resolve()).replace("'", "'\\''")
        lines.append(f"file '{escaped}'")
    concat.write_text("\n".join(lines) + "\n", encoding="utf-8")
    command = [
        "ffmpeg",
        "-hide_banner",
        "-loglevel",
        "error",
        "-y",
        "-f",
        "concat",
        "-safe",
        "0",
        "-i",
        str(concat),
        "-c",
        "copy",
        str(output / "all_cases_sequential.mp4"),
    ]
    try:
        subprocess.run(command, check=True, timeout=120.0)
    except (
        FileNotFoundError,
        subprocess.CalledProcessError,
        subprocess.TimeoutExpired,
    ):
        print(
            "warning: per-case videos were created, but ffmpeg concatenation failed",
            file=sys.stderr,
        )


def _make_contact_sheet(output: Path, results: list[dict[str, Any]]) -> None:
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        return
    tiles: list[tuple[str, bool, Any]] = []
    for result in results:
        path = output / result["case"] / "final.png"
        if path.is_file():
            tiles.append(
                (result["case"], result["classification"]["stable"], Image.open(path))
            )
    if not tiles:
        return
    tile_width, tile_height = 480, 400
    columns = 2
    rows = math.ceil(len(tiles) / columns)
    sheet = Image.new("RGB", (columns * tile_width, rows * tile_height), "white")
    draw = ImageDraw.Draw(sheet)
    font = _display_font(22)
    for index, (name, stable, image) in enumerate(tiles):
        column, row = index % columns, index // columns
        resized = image.resize((tile_width, 360))
        x, y = column * tile_width, row * tile_height
        sheet.paste(resized, (x, y))
        status = "STABLE" if stable else "FAILED"
        color = (20, 125, 55) if stable else (180, 35, 35)
        label = name.replace("_", " ")
        draw.text((x + 12, y + 369), f"{label}  {status}", fill=color, font=font)
    sheet.save(output / "final_contact_sheet.png")
    for _, _, image in tiles:
        image.close()


def _select_cases(cases: Iterable[Case], requested: Sequence[str]) -> list[Case]:
    cases = list(cases)
    if not requested or "all" in requested:
        return cases
    lookup = {case.name: case for case in cases}
    missing = sorted(set(requested) - set(lookup))
    if missing:
        raise ValueError(f"Unknown cases: {', '.join(missing)}")
    return [lookup[name] for name in requested]


def _config_from_arguments(arguments: argparse.Namespace) -> SimulationConfig:
    preset = VHACD_PRESETS[arguments.preset]
    target_extent_m = (
        arguments.model_size_cm / 100.0
        if arguments.model_size_cm is not None
        else arguments.target_extent_m
    )
    if target_extent_m is None:
        target_extent_m = 0.16 if arguments.scene == "table" else 0.10
    camera_distance_m = arguments.camera_distance_m
    if camera_distance_m is None:
        camera_distance_m = 0.50 if arguments.scene == "table" else 0.30
    return SimulationConfig(
        scene=arguments.scene,
        preset=arguments.preset,
        target_extent_m=target_extent_m,
        table_width_m=arguments.table_width_m,
        table_depth_m=arguments.table_depth_m,
        table_height_m=arguments.table_height_m,
        table_thickness_m=arguments.table_thickness_m,
        camera_distance_m=camera_distance_m,
        camera_yaw_deg=arguments.camera_yaw_deg,
        camera_pitch_deg=arguments.camera_pitch_deg,
        mass_kg=arguments.mass_kg,
        initial_tilt_deg=arguments.initial_tilt_deg,
        initial_tilt_azimuth_deg=arguments.initial_tilt_azimuth_deg,
        duration_s=arguments.duration,
        time_step_s=1.0 / arguments.physics_hz,
        sample_hz=arguments.sample_hz,
        video_fps=arguments.video_fps,
        object_friction=arguments.friction,
        ground_friction=arguments.friction,
        vhacd_resolution=(
            arguments.vhacd_resolution
            if arguments.vhacd_resolution is not None
            else int(preset["resolution"])
        ),
        vhacd_concavity=(
            arguments.vhacd_concavity
            if arguments.vhacd_concavity is not None
            else float(preset["concavity"])
        ),
        vhacd_timeout_s=arguments.vhacd_timeout,
    )


def _cases_from_arguments(arguments: argparse.Namespace) -> list[Case]:
    cases = (
        load_manifest(arguments.manifest.resolve())
        if arguments.manifest
        else discover_r11_cases(arguments.r11_root.resolve())
    )
    return _select_cases(cases, arguments.case)


def batch_main(arguments: argparse.Namespace) -> int:
    pybullet, _ = _require_modules()
    config = _config_from_arguments(arguments)
    cases = _cases_from_arguments(arguments)
    if arguments.dry_run:
        checks = []
        for case in cases:
            properties = _load_mesh_properties(
                case, config.target_extent_m, config.mass_kg
            )
            checks.append(
                {
                    "case": case.name,
                    "collision_mesh": str(case.collision_mesh),
                    "faces": properties["faces"],
                    "watertight": properties["watertight"],
                    "mass_properties_source": properties["mass_properties_source"],
                }
            )
        print(
            json.dumps(
                {
                    "status": "ready",
                    "python": sys.executable,
                    "cases": checks,
                    "configuration": asdict(config),
                },
                indent=2,
            )
        )
        return 0
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    client = pybullet.connect(pybullet.DIRECT)
    if client < 0:
        raise RuntimeError("Could not connect to PyBullet DIRECT mode")
    results: list[dict[str, Any]] = []
    try:
        for index, case in enumerate(cases, start=1):
            print(f"[{index:02d}/{len(cases):02d}] {case.name}", flush=True)
            try:
                vhacd, metadata = ensure_vhacd(
                    case,
                    output / "collision_cache",
                    config,
                    force=arguments.force_vhacd,
                )
                if case.fixture_mesh is None:
                    result = simulate_case(
                        pybullet,
                        client,
                        case,
                        vhacd,
                        metadata,
                        output,
                        config,
                        render=not arguments.no_render,
                        record_video=arguments.record_video,
                    )
                else:
                    fixture_case = Case(
                        name=f"{case.name}__fixture",
                        collision_mesh=case.fixture_mesh,
                        visual_mesh=case.fixture_mesh,
                        provenance=case.provenance,
                    )
                    fixture_vhacd, fixture_metadata = ensure_vhacd(
                        fixture_case,
                        output / "collision_cache",
                        config,
                        force=arguments.force_vhacd,
                    )
                    result = simulate_fixture_case(
                        pybullet,
                        client,
                        case,
                        vhacd,
                        metadata,
                        fixture_vhacd,
                        fixture_metadata,
                        output,
                        config,
                        render=not arguments.no_render,
                        record_video=arguments.record_video,
                    )
                results.append(result)
                status = "stable" if result["classification"]["stable"] else "failed"
                print(
                    f"    {status}: final tilt "
                    f"{result['classification']['final_tilt_deg']:.2f} deg",
                    flush=True,
                )
            except Exception as error:
                failure = {
                    "schema": "make-it-stand-pybullet-case-v1",
                    "case": case.name,
                    "provenance": case.provenance,
                    "expected_static_audit": {
                        "stable": case.expected_stable,
                        "signed_margin_cells": case.expected_margin_cells,
                    },
                    "classification": {
                        "stable": False,
                        "fallen": False,
                        "target_pose_retained": False,
                        "failure_reasons": [f"simulation_error: {error}"],
                    },
                }
                _json_dump(output / case.name / "result.json", failure)
                results.append(failure)
                print(f"    error: {error}", file=sys.stderr, flush=True)
    finally:
        pybullet.disconnect(physicsClientId=client)
    _write_summary(output, results, config)
    if arguments.record_video:
        _make_sequential_video(output, results)
    print(f"Summary: {output / 'summary.json'}")
    return (
        0
        if all(
            "simulation_error"
            not in "".join(result["classification"]["failure_reasons"])
            for result in results
        )
        else 2
    )


def replay_main(arguments: argparse.Namespace) -> int:
    pybullet, _ = _require_modules()
    config = _config_from_arguments(arguments)
    cases = _cases_from_arguments(arguments)
    output = arguments.output.resolve()
    client = pybullet.connect(pybullet.GUI)
    if client < 0:
        raise RuntimeError(
            "Could not open the PyBullet GUI (is an X display available?)"
        )
    pybullet.configureDebugVisualizer(
        pybullet.COV_ENABLE_GUI, 1, physicsClientId=client
    )
    pybullet.resetDebugVisualizerCamera(
        cameraDistance=config.camera_distance_m,
        cameraYaw=config.camera_yaw_deg,
        cameraPitch=config.camera_pitch_deg,
        cameraTargetPosition=[
            0.0,
            0.0,
            (
                config.table_height_m + 0.18 * config.target_extent_m
                if config.scene == "table"
                else 0.45 * config.target_extent_m
            ),
        ],
        physicsClientId=client,
    )
    try:
        repeat = True
        while repeat:
            for case in cases:
                print(f"Replaying {case.name}", flush=True)
                vhacd, metadata = ensure_vhacd(
                    case, output / "collision_cache", config, force=False
                )
                simulate_case(
                    pybullet,
                    client,
                    case,
                    vhacd,
                    metadata,
                    output,
                    config,
                    render=False,
                    record_video=False,
                    realtime=True,
                )
                time.sleep(arguments.pause)
            repeat = arguments.loop
        if not arguments.exit_on_complete:
            print(
                "Simulation complete. The final scene will remain open; "
                "close the PyBullet window or press Ctrl+C to exit.",
                flush=True,
            )
            try:
                while pybullet.isConnected(physicsClientId=client):
                    time.sleep(0.1)
            except KeyboardInterrupt:
                print("Closing PyBullet replay.", flush=True)
    finally:
        if pybullet.isConnected(physicsClientId=client):
            pybullet.disconnect(physicsClientId=client)
    return 0


def vhacd_main(arguments: argparse.Namespace) -> int:
    pybullet, _ = _require_modules()
    client = pybullet.connect(pybullet.DIRECT)
    if client < 0:
        return 3
    try:
        pybullet.vhacd(
            str(arguments.input),
            str(arguments.output),
            str(arguments.log),
            resolution=arguments.resolution,
            concavity=arguments.concavity,
            maxNumVerticesPerCH=arguments.max_vertices,
        )
    finally:
        pybullet.disconnect(physicsClientId=client)
    return 0 if arguments.output.is_file() and arguments.output.stat().st_size else 4


def _add_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--r11-root", type=Path, default=DEFAULT_R11_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument("--duration", type=float, default=8.0)
    parser.add_argument("--physics-hz", type=float, default=240.0)
    parser.add_argument("--sample-hz", type=float, default=30.0)
    parser.add_argument("--video-fps", type=float, default=12.0)
    parser.add_argument(
        "--target-extent-m",
        type=float,
        default=None,
        help="model maximum extent in metres (10 cm ground / 16 cm table default)",
    )
    parser.add_argument(
        "--model-size-cm",
        type=float,
        default=None,
        help="convenient override for model maximum extent in centimetres",
    )
    parser.add_argument("--table-width-m", type=float, default=0.50)
    parser.add_argument("--table-depth-m", type=float, default=0.38)
    parser.add_argument("--table-height-m", type=float, default=0.72)
    parser.add_argument("--table-thickness-m", type=float, default=0.035)
    parser.add_argument(
        "--camera-distance-m",
        type=float,
        default=None,
        help="camera distance; smaller values zoom in",
    )
    parser.add_argument("--camera-yaw-deg", type=float, default=38.0)
    parser.add_argument("--camera-pitch-deg", type=float, default=-12.0)
    parser.add_argument("--mass-kg", type=float, default=0.20)
    parser.add_argument("--initial-tilt-deg", type=float, default=0.0)
    parser.add_argument("--initial-tilt-azimuth-deg", type=float, default=0.0)
    parser.add_argument("--friction", type=float, default=0.80)
    parser.add_argument(
        "--scene",
        choices=("ground", "table"),
        default="ground",
        help="place the object on the legacy ground or on a physical table",
    )
    parser.add_argument(
        "--preset",
        choices=sorted(VHACD_PRESETS),
        default="default",
        help=(
            "V-HACD preset: default preserves the original 100k/0.0025 "
            "protocol; robust uses 1M/0.001 for a sensitivity check"
        ),
    )
    parser.add_argument(
        "--vhacd-resolution",
        type=int,
        default=None,
        help="override the selected preset's V-HACD resolution",
    )
    parser.add_argument(
        "--vhacd-concavity",
        type=float,
        default=None,
        help="override the selected preset's V-HACD concavity",
    )
    parser.add_argument("--vhacd-timeout", type=float, default=180.0)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    batch = commands.add_parser("batch", help="run deterministic headless simulations")
    _add_common_arguments(batch)
    batch.add_argument("--force-vhacd", action="store_true")
    batch.add_argument("--no-render", action="store_true")
    batch.add_argument(
        "--dry-run",
        action="store_true",
        help="validate dependencies, case discovery, and mesh mass properties only",
    )
    batch.add_argument(
        "--record-video",
        action="store_true",
        help="write one MP4 per case and an all-cases sequential MP4",
    )
    batch.set_defaults(function=batch_main)
    replay = commands.add_parser(
        "replay", help="show one/all simulations in PyBullet GUI"
    )
    _add_common_arguments(replay)
    replay.add_argument("--pause", type=float, default=1.5)
    replay.add_argument("--loop", action="store_true")
    replay.add_argument(
        "--exit-on-complete",
        action="store_true",
        help="close the GUI after playback instead of holding the final scene",
    )
    replay.set_defaults(function=replay_main)
    worker = commands.add_parser("vhacd", help=argparse.SUPPRESS)
    worker.add_argument("--input", type=Path, required=True)
    worker.add_argument("--output", type=Path, required=True)
    worker.add_argument("--log", type=Path, required=True)
    worker.add_argument("--resolution", type=int, required=True)
    worker.add_argument("--concavity", type=float, required=True)
    worker.add_argument("--max-vertices", type=int, required=True)
    worker.set_defaults(function=vhacd_main)
    return parser


def main() -> int:
    parser = build_parser()
    arguments = parser.parse_args()
    if getattr(arguments, "duration", 1.0) <= 0.0:
        parser.error("--duration must be positive")
    if getattr(arguments, "physics_hz", 1.0) <= 0.0:
        parser.error("--physics-hz must be positive")
    if (
        getattr(arguments, "model_size_cm", None) is not None
        and getattr(arguments, "target_extent_m", None) is not None
    ):
        parser.error("use only one of --model-size-cm and --target-extent-m")
    positive_dimensions = (
        "model_size_cm",
        "target_extent_m",
        "table_width_m",
        "table_depth_m",
        "table_height_m",
        "table_thickness_m",
        "camera_distance_m",
    )
    for name in positive_dimensions:
        value = getattr(arguments, name, None)
        if value is not None and value <= 0.0:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    if not math.isfinite(getattr(arguments, "initial_tilt_deg", 0.0)) or not (
        0.0 <= getattr(arguments, "initial_tilt_deg", 0.0) < 90.0
    ):
        parser.error("--initial-tilt-deg must lie in [0, 90)")
    if not math.isfinite(getattr(arguments, "initial_tilt_azimuth_deg", 0.0)):
        parser.error("--initial-tilt-azimuth-deg must be finite")
    return int(arguments.function(arguments))


if __name__ == "__main__":
    raise SystemExit(main())
