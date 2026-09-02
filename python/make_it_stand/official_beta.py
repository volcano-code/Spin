"""Explicit importer for the archived Make-It-Stand beta presets.

This module is intentionally separate from :class:`BalanceModel`.  Loading a
``.mis`` file here opts into the archived unit-cube occupancy and BBW cache
contract; ordinary mesh construction continues to use the modern voxelizer,
weight solver, and strict source gate.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
import os
from pathlib import Path
import re
from typing import Any

import numpy as np

from .model import (
    BalanceResult,
    ModelConfig,
    RunConfig,
    _load_native,
    _native_model_config,
    _native_run_config,
)
from .paper_compat import PaperNormalization, released_beta_run_config


_DEFAULT_GAMMA_RADIANS = math.radians(5.0)
_BETA_MINIMUM_STABILITY_INSET = 1.0e-4


@dataclass(frozen=True)
class ReleasedBetaMis:
    """Supported single-object subset of the archived ``.mis`` grammar."""

    mesh: str
    mesh_rotation_wxyz: tuple[float, float, float, float]
    flatten_threshold: float
    vertex_handles: tuple[tuple[float, float, float], ...]


def _finite_numbers(
    fields: list[str],
    count: int,
    path: Path,
    line_number: int,
    tag: str,
) -> tuple[float, ...]:
    if len(fields) != count + 1:
        raise ValueError(
            f"{path}:{line_number}: {tag} requires exactly {count} value(s)"
        )
    try:
        values = tuple(float(item) for item in fields[1:])
    except ValueError as error:
        raise ValueError(
            f"{path}:{line_number}: {tag} contains a non-number"
        ) from error
    if not all(math.isfinite(item) for item in values):
        raise ValueError(f"{path}:{line_number}: {tag} must be finite")
    return values


def parse_released_beta_mis(path: os.PathLike[str] | str) -> ReleasedBetaMis:
    """Parse the archived, single-standing-object preset grammar strictly."""

    source = Path(path)
    singles: dict[str, Any] = {}
    handles: list[tuple[float, float, float]] = []
    for line_number, raw_line in enumerate(
        source.read_text(encoding="utf-8", errors="strict").splitlines(),
        start=1,
    ):
        fields = raw_line.split("#", maxsplit=1)[0].split()
        if not fields:
            continue
        tag = fields[0]
        if tag == "mesh":
            if len(fields) != 2 or not fields[1]:
                raise ValueError(f"{source}:{line_number}: mesh requires one name")
            value: Any = fields[1]
        elif tag == "mesh_rotation":
            value = _finite_numbers(fields, 4, source, line_number, tag)
        elif tag == "flatten_threshold":
            value = _finite_numbers(fields, 1, source, line_number, tag)[0]
            if value < 0.0:
                raise ValueError(
                    f"{source}:{line_number}: flatten_threshold must be nonnegative"
                )
        elif tag == "vertex_handle":
            values = _finite_numbers(fields, 3, source, line_number, tag)
            handles.append(values)  # type: ignore[arg-type]
            continue
        else:
            raise ValueError(
                f"{source}:{line_number}: unsupported compatibility tag {tag!r}"
            )
        if tag in singles:
            raise ValueError(f"{source}:{line_number}: duplicate {tag} tag")
        singles[tag] = value

    missing = [
        tag for tag in ("mesh", "mesh_rotation", "flatten_threshold")
        if tag not in singles
    ]
    if missing:
        raise ValueError(f"{source}: missing tag(s): {', '.join(missing)}")
    mesh_name = str(singles["mesh"])
    mesh_path = Path(mesh_name)
    if mesh_path.is_absolute() or ".." in mesh_path.parts:
        raise ValueError(f"{source}: mesh name must stay within data/meshes")
    return ReleasedBetaMis(
        mesh=mesh_name,
        mesh_rotation_wxyz=singles["mesh_rotation"],
        flatten_threshold=float(singles["flatten_threshold"]),
        vertex_handles=tuple(handles),
    )


def _parse_obj(path: Path) -> tuple[np.ndarray, np.ndarray]:
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []
    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8", errors="strict").splitlines(),
        start=1,
    ):
        fields = raw_line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if fields[0] == "v":
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: incomplete vertex")
            try:
                point = tuple(float(item) for item in fields[1:4])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid vertex") from error
            if not all(math.isfinite(item) for item in point):
                raise ValueError(f"{path}:{line_number}: non-finite vertex")
            vertices.append(point)  # type: ignore[arg-type]
        elif fields[0] == "f":
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: incomplete face")
            indices: list[int] = []
            for field in fields[1:]:
                try:
                    raw_index = int(field.split("/", maxsplit=1)[0])
                except ValueError as error:
                    raise ValueError(
                        f"{path}:{line_number}: invalid face index"
                    ) from error
                if raw_index == 0:
                    raise ValueError(f"{path}:{line_number}: OBJ index zero is invalid")
                index = raw_index - 1 if raw_index > 0 else len(vertices) + raw_index
                if index < 0 or index >= len(vertices):
                    raise ValueError(f"{path}:{line_number}: face index out of range")
                indices.append(index)
            for corner in range(1, len(indices) - 1):
                faces.append((indices[0], indices[corner], indices[corner + 1]))
        elif fields[0] in {"vn", "vt", "o", "g", "s", "usemtl", "mtllib"}:
            continue

    if not vertices or not faces:
        raise ValueError(f"{path}: OBJ must contain vertices and faces")
    return (
        np.ascontiguousarray(vertices, dtype=np.float64),
        np.ascontiguousarray(faces, dtype=np.int32),
    )


def _infer_cache(config_directory: Path, stem: str) -> tuple[int, Path, Path]:
    pattern = re.compile(rf"^{re.escape(stem)}_(\d+)\.vox$")
    candidates: list[tuple[int, Path]] = []
    for path in config_directory.glob(f"{stem}_*.vox"):
        match = pattern.fullmatch(path.name)
        if match is not None:
            candidates.append((int(match.group(1)), path))
    if not candidates:
        raise ValueError(f"no released-beta .vox cache matches {stem!r}")
    resolution, vox_path = max(candidates, key=lambda item: item[0])
    bbw_path = vox_path.with_suffix(".bbw")
    if not bbw_path.is_file():
        raise ValueError(f"missing matching released-beta BBW cache: {bbw_path.name}")
    return resolution, vox_path, bbw_path


def _parse_vox(path: Path, resolution: int) -> np.ndarray:
    tokens = path.read_text(encoding="utf-8", errors="strict").split()
    expected = resolution ** 3
    if len(tokens) != expected:
        raise ValueError(f"{path}: expected {expected} cell IDs, got {len(tokens)}")
    try:
        values = np.asarray([int(token) for token in tokens], dtype=np.int64)
    except ValueError as error:
        raise ValueError(f"{path}: .vox contains a non-integer token") from error
    if np.any(values < -1):
        raise ValueError(f"{path}: .vox cell IDs must be -1 or nonnegative")
    occupied = values[values >= 0]
    if occupied.size == 0 or not np.array_equal(
        np.sort(occupied), np.arange(occupied.size, dtype=np.int64)
    ):
        raise ValueError(f"{path}: occupied IDs must be unique and contiguous")
    if occupied.size > np.iinfo(np.int32).max:
        raise ValueError(f"{path}: occupied IDs exceed int32 capacity")
    return np.ascontiguousarray(values, dtype=np.int32)


def _parse_bbw(path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    columns: int | None = None
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8", errors="strict").splitlines(),
        start=1,
    ):
        fields = line.split()
        if not fields:
            continue
        try:
            row = [float(item) for item in fields]
        except ValueError as error:
            raise ValueError(f"{path}:{line_number}: invalid BBW value") from error
        if columns is None:
            columns = len(row)
        if not row or len(row) != columns:
            raise ValueError(f"{path}:{line_number}: inconsistent BBW column count")
        if not all(math.isfinite(item) for item in row):
            raise ValueError(f"{path}:{line_number}: non-finite BBW value")
        rows.append(row)
    if not rows:
        raise ValueError(f"{path}: BBW cache is empty")
    return np.ascontiguousarray(rows, dtype=np.float64)


def _active_node_mask(cell_ids_z_fast: np.ndarray, resolution: int) -> np.ndarray:
    occupied = np.asarray(cell_ids_z_fast).reshape((resolution,) * 3) >= 0
    active = np.zeros((resolution + 1,) * 3, dtype=np.bool_)
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                active[
                    dx:dx + resolution,
                    dy:dy + resolution,
                    dz:dz + resolution,
                ] |= occupied
    return active


def normalize_and_map_beta_bbw(
    cell_ids_z_fast: Any,
    beta_weights_z_fast: Any,
    resolution: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (beta-row-normalized, current-coordinate-order) BBW arrays."""

    cells = np.asarray(cell_ids_z_fast)
    weights = np.asarray(beta_weights_z_fast, dtype=np.float64)
    if cells.shape != (resolution ** 3,):
        raise ValueError("cell_ids_z_fast must contain exactly N^3 entries")
    active = _active_node_mask(cells, resolution)
    active_coordinates = np.argwhere(active)
    if (
        weights.ndim != 2 or
        weights.shape[0] != active_coordinates.shape[0] or
        weights.shape[1] <= 0 or
        not np.isfinite(weights).all() or
        np.any(weights < 0.0)
    ):
        raise ValueError("BBW must be finite, nonnegative, and active-node by handle")
    row_sums = weights.sum(axis=1)
    if not np.isfinite(row_sums).all() or np.any(row_sums <= 0.0):
        raise ValueError("every BBW row must have a positive finite sum")
    normalized = np.ascontiguousarray(weights / row_sums[:, np.newaxis])

    # np.argwhere traverses (x,y,z) in C order, exactly the beta's z-fast
    # row order. Sort those coordinates by the current x-fast global corner
    # key. No serialized node or occupied-cell ID participates in the map.
    side = resolution + 1
    current_keys = (
        active_coordinates[:, 0]
        + side * (active_coordinates[:, 1] + side * active_coordinates[:, 2])
    )
    current_order = np.argsort(current_keys, kind="stable")
    mapped = np.ascontiguousarray(normalized[current_order])
    return normalized, mapped


def _gravity_from_raw_quaternion(
    quaternion_wxyz: tuple[float, float, float, float],
) -> tuple[float, float, float]:
    # Match Eigen::Quaternion::conjugate().toRotationMatrix() in the beta.
    w, x, y, z = quaternion_wxyz
    x, y, z = -x, -y, -z
    rotation = np.asarray(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ],
        dtype=np.float64,
    )
    gravity = rotation @ np.asarray((0.0, -1.0, 0.0), dtype=np.float64)
    if not np.isfinite(gravity).all() or float(np.linalg.norm(gravity)) <= 1.0e-15:
        raise ValueError("mesh_rotation does not define a finite gravity direction")
    return tuple(float(item) for item in gravity)


@dataclass(frozen=True)
class ReleasedBetaArchive:
    """Validated archived preset plus its exact occupancy/BBW payload."""

    preset_name: str
    specification: ReleasedBetaMis
    source_vertices: np.ndarray
    normalized_vertices: np.ndarray
    faces: np.ndarray
    voxel_resolution: int
    cell_ids_z_fast: np.ndarray
    beta_weights_z_fast: np.ndarray
    beta_normalized_weights: np.ndarray
    coordinate_mapped_weights: np.ndarray
    gravity: tuple[float, float, float]

    @property
    def occupied_cell_count(self) -> int:
        return int(np.count_nonzero(self.cell_ids_z_fast >= 0))

    @property
    def active_node_count(self) -> int:
        return int(self.beta_weights_z_fast.shape[0])

    @classmethod
    def load(
        cls,
        mis_path: os.PathLike[str] | str,
        *,
        beta_root: os.PathLike[str] | str | None = None,
    ) -> "ReleasedBetaArchive":
        source = Path(mis_path).resolve()
        root = (
            Path(beta_root).resolve()
            if beta_root is not None
            else source.parent.parent.parent
        )
        expected_configs = (root / "data" / "configs").resolve()
        if source.parent != expected_configs:
            raise ValueError(".mis file must be inside the selected data/configs directory")
        specification = parse_released_beta_mis(source)
        mesh_path = (root / "data" / "meshes" / specification.mesh).resolve()
        mesh_directory = (root / "data" / "meshes").resolve()
        if mesh_path.parent != mesh_directory or not mesh_path.is_file():
            raise ValueError("preset mesh is missing or escapes data/meshes")
        source_vertices, faces = _parse_obj(mesh_path)
        normalization = PaperNormalization.from_vertices(source_vertices)
        normalized_vertices = normalization.to_normalized(source_vertices)
        resolution, vox_path, bbw_path = _infer_cache(source.parent, source.stem)
        cell_ids = _parse_vox(vox_path, resolution)
        beta_weights = _parse_bbw(bbw_path)
        expected_columns = 1 + len(specification.vertex_handles)
        if beta_weights.shape[1] != expected_columns:
            raise ValueError(
                ".bbw column count must equal one support plus vertex_handle count"
            )
        beta_normalized, mapped = normalize_and_map_beta_bbw(
            cell_ids,
            beta_weights,
            resolution,
        )
        return cls(
            preset_name=source.stem,
            specification=specification,
            source_vertices=source_vertices,
            normalized_vertices=normalized_vertices,
            faces=faces,
            voxel_resolution=resolution,
            cell_ids_z_fast=cell_ids,
            beta_weights_z_fast=beta_weights,
            beta_normalized_weights=beta_normalized,
            coordinate_mapped_weights=mapped,
            gravity=_gravity_from_raw_quaternion(
                specification.mesh_rotation_wxyz
            ),
        )

    def build_model(
        self,
        *,
        hull_depth: int = 2,
        native_module: Any | None = None,
    ) -> "ReleasedBetaArchiveModel":
        if isinstance(hull_depth, bool) or not isinstance(hull_depth, int) or hull_depth < 0:
            raise ValueError("hull_depth must be a nonnegative integer")
        native = native_module if native_module is not None else _load_native()
        settings = ModelConfig(
            path="<released-beta-archive>",
            voxel_resolution=self.voxel_resolution,
            minimum_wall_thickness=0.0,
            fixed_unit_cube_voxel_grid=True,
        )
        native_model = native.BalanceModel.from_released_beta_archive(
            self.normalized_vertices,
            self.faces,
            self.voxel_resolution,
            self.cell_ids_z_fast,
            self.beta_weights_z_fast,
            hull_depth,
            _native_model_config(native, settings),
        )
        native_model.set_gravity(self.gravity)
        original_normalized = PaperNormalization.from_vertices(
            self.source_vertices
        ).to_normalized(self.source_vertices)
        if (
            self.normalized_vertices.shape == original_normalized.shape and
            np.array_equal(self.normalized_vertices, original_normalized)
        ):
            native_model.set_support_from_mesh(
                self.specification.flatten_threshold,
                _DEFAULT_GAMMA_RADIANS,
                _BETA_MINIMUM_STABILITY_INSET,
                True,
            )
        else:
            native_model.set_support_from_reference_mesh(
                original_normalized,
                self.specification.flatten_threshold,
                _DEFAULT_GAMMA_RADIANS,
                _BETA_MINIMUM_STABILITY_INSET,
                True,
            )
        for handle in self.specification.vertex_handles:
            native_model.add_handle(handle)
        return ReleasedBetaArchiveModel(self, native_model, native)

    def with_reconstructed_normalized_surface(
        self,
        path: os.PathLike[str] | str,
    ) -> "ReleasedBetaArchive":
        """Replace only the normalized outer surface for labelled reconstruction.

        Occupancy, BBW, handles, gravity, and the untouched source vertices
        remain bound to the original archive.  The caller is responsible for
        recording the replacement path/hash and must not classify the result
        as an untouched official input.
        """

        vertices, faces = _parse_obj(Path(path).expanduser().resolve())
        tolerance = 4096.0 * np.finfo(np.float64).eps
        if (
            np.any(vertices < -tolerance) or
            np.any(vertices > 1.0 + tolerance)
        ):
            raise ValueError(
                "reconstructed normalized surface must remain in the unit cube"
            )
        return replace(
            self,
            normalized_vertices=np.ascontiguousarray(vertices),
            faces=np.ascontiguousarray(faces),
        )


class ReleasedBetaArchiveModel:
    """Owning model built only through :meth:`ReleasedBetaArchive.build_model`."""

    def __init__(self, archive: ReleasedBetaArchive, native_model: Any, native: Any):
        self._archive = archive
        self._native = native_model
        self._native_module = native

    @property
    def archive(self) -> ReleasedBetaArchive:
        return self._archive

    @property
    def native(self) -> Any:
        return self._native

    def optimize(self, config: RunConfig | None = None) -> BalanceResult:
        settings = released_beta_run_config() if config is None else config
        return BalanceResult(
            self._native.optimize(
                _native_run_config(self._native_module, settings)
            )
        )

    def optimize_released_beta(self, **overrides: Any) -> BalanceResult:
        return self.optimize(released_beta_run_config(**overrides))


def load_released_beta_archive(
    mis_path: os.PathLike[str] | str,
    *,
    beta_root: os.PathLike[str] | str | None = None,
) -> ReleasedBetaArchive:
    """Load and validate one archived single-object beta preset."""

    return ReleasedBetaArchive.load(mis_path, beta_root=beta_root)


__all__ = [
    "ReleasedBetaArchive",
    "ReleasedBetaArchiveModel",
    "ReleasedBetaMis",
    "load_released_beta_archive",
    "normalize_and_map_beta_bbw",
    "parse_released_beta_mis",
]
