"""Paper-compatible normalization and support-flattening production wrapper.

The released beta implementation normalizes every input into a unit box
before voxelization and projects the vertices selected by a positive flatten
threshold onto the support plane when it exports the result.  The audited
native facade intentionally treats a positive threshold as hypothetical,
because its sealed manufacturing report describes the unprojected geometry.

This module bridges those two contracts without weakening either one:

* optimization runs in the beta coordinate convention (longest extent 0.95,
  bounding-box centre at ``(0.5, 0.5, 0.5)``);
* pivots and exported geometry are mapped between caller and normalized
  coordinates by one explicit, invertible similarity transform;
* support-band vertices are really projected for diagnostic OBJ/STL export;
* a post-projection boundary audit is reported, but the projected result is
  never promoted to sealed manufacturing-ready status without the remaining
  native manufacturing checks being rerun on that exact geometry.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
import os
from pathlib import Path
import struct
from typing import Any, Sequence

import numpy as np

from .io import write_json_atomic, write_text_atomic
from .model import BalanceModel, BalanceResult, ModelConfig, RunConfig


_TARGET_CENTER = (0.5, 0.5, 0.5)
_TARGET_LONGEST_EXTENT = 0.95


PAPER_RUN_CONFIG = RunConfig(
    # The paper's released implementation initializes these user-facing
    # weights to 0.75 and 20.  The paper profile differs from the beta below
    # in its printed 1% fixed-alpha block-coordinate policy, not these weights.
    initial_mu=0.75,
    shape_scale_lambda=20.0,
    required_relative_decrease=0.01,
    freeze_targets_after_initial_refresh=False,
    released_beta_semantics=False,
    released_beta_outer_loop=False,
)


def paper_run_config(**overrides: Any) -> RunConfig:
    """Return the paper block-coordinate profile with explicit overrides."""

    return replace(PAPER_RUN_CONFIG, **overrides)


RELEASED_BETA_RUN_CONFIG = RunConfig(
    initial_mu=0.75,
    shape_scale_lambda=20.0,
    initial_step=1.0,
    reduction_ratio=0.8,
    minimum_step=0.4,
    required_relative_decrease=0.03,
    mu_decrement=0.05,
    freeze_targets_after_initial_refresh=True,
    released_beta_semantics=True,
    released_beta_outer_loop=True,
)


def released_beta_run_config(**overrides: Any) -> RunConfig:
    """Return the audited released-beta parameter/semantics profile.

    The beta GUI has no counterparts for this implementation's finite block
    and iteration caps, gradient tolerance, or numerical mass guards. Those
    bounded-execution settings retain their modern defaults and may be
    overridden explicitly here.
    """

    return replace(RELEASED_BETA_RUN_CONFIG, **overrides)


def _vertices(
    value: Any,
    name: str = "vertices",
    *,
    allow_empty: bool = False,
) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if (
        result.ndim != 2 or
        result.shape[1] != 3 or
        (result.shape[0] == 0 and not allow_empty)
    ):
        qualifier = "N >= 0" if allow_empty else "N > 0"
        raise ValueError(f"{name} must have shape (N, 3) with {qualifier}")
    if not np.isfinite(result).all():
        raise ValueError(f"{name} must contain only finite values")
    return np.ascontiguousarray(result)


def _faces(value: Any) -> np.ndarray:
    original = np.asarray(value)
    if original.ndim != 2 or original.shape[1] != 3 or original.shape[0] == 0:
        raise ValueError("faces must have shape (M, 3) with M > 0")
    if original.dtype.kind not in {"i", "u"}:
        raise TypeError("faces must use an integer dtype")
    if original.size and (
        int(original.min()) < 0 or int(original.max()) > np.iinfo(np.int32).max
    ):
        raise ValueError("face indices must be non-negative int32 values")
    return np.ascontiguousarray(original, dtype=np.int32)


def _unit_vector(value: Sequence[float], name: str) -> np.ndarray:
    vector = np.asarray(value, dtype=np.float64)
    if vector.shape != (3,) or not np.isfinite(vector).all():
        raise ValueError(f"{name} must contain exactly three finite values")
    norm = float(np.linalg.norm(vector))
    if norm <= 1.0e-15:
        raise ValueError(f"{name} must be non-zero")
    return vector / norm


@dataclass(frozen=True)
class PaperNormalization:
    """Invertible beta-coordinate similarity transform."""

    source_center: tuple[float, float, float]
    scale: float
    target_center: tuple[float, float, float] = _TARGET_CENTER
    target_longest_extent: float = _TARGET_LONGEST_EXTENT

    @classmethod
    def from_vertices(cls, vertices: Any) -> "PaperNormalization":
        points = _vertices(vertices)
        minimum = points.min(axis=0)
        maximum = points.max(axis=0)
        extent = maximum - minimum
        longest = float(extent.max())
        if not math.isfinite(longest) or longest <= 0.0:
            raise ValueError("source mesh must have a positive bounding-box extent")
        center = 0.5 * (minimum + maximum)
        return cls(
            source_center=tuple(float(item) for item in center),
            scale=_TARGET_LONGEST_EXTENT / longest,
        )

    def to_normalized(self, points: Any) -> np.ndarray:
        values = _vertices(points, "points", allow_empty=True)
        return (
            self.scale
            * (values - np.asarray(self.source_center, dtype=np.float64))
            + np.asarray(self.target_center, dtype=np.float64)
        )

    def to_original(self, points: Any) -> np.ndarray:
        values = _vertices(points, "points", allow_empty=True)
        return (
            (values - np.asarray(self.target_center, dtype=np.float64))
            / self.scale
            + np.asarray(self.source_center, dtype=np.float64)
        )

    def point_to_normalized(self, point: Sequence[float]) -> tuple[float, float, float]:
        values = np.asarray(point, dtype=np.float64)
        if values.shape != (3,) or not np.isfinite(values).all():
            raise ValueError("point must contain exactly three finite values")
        result = (
            self.scale
            * (values - np.asarray(self.source_center, dtype=np.float64))
            + np.asarray(self.target_center, dtype=np.float64)
        )
        return tuple(float(item) for item in result)

    def point_to_original(self, point: Sequence[float]) -> tuple[float, float, float]:
        values = np.asarray(point, dtype=np.float64)
        if values.shape != (3,) or not np.isfinite(values).all():
            raise ValueError("point must contain exactly three finite values")
        result = (
            (values - np.asarray(self.target_center, dtype=np.float64))
            / self.scale
            + np.asarray(self.source_center, dtype=np.float64)
        )
        return tuple(float(item) for item in result)

    def to_dict(self) -> dict[str, Any]:
        return {
            "source_center": list(self.source_center),
            "scale": self.scale,
            "target_center": list(self.target_center),
            "target_longest_extent": self.target_longest_extent,
        }


@dataclass(frozen=True)
class SupportFlattening:
    """Frozen support band selected from the normalized source mesh."""

    gravity: tuple[float, float, float]
    contact_threshold: float
    minimum_stability_inset: float
    support_height: float
    vertex_indices: tuple[int, ...]

    @classmethod
    def from_vertices(
        cls,
        vertices: Any,
        gravity: Sequence[float],
        contact_threshold: float,
        minimum_stability_inset: float = 0.0,
    ) -> "SupportFlattening":
        points = _vertices(vertices)
        direction = _unit_vector(gravity, "gravity")
        threshold = float(contact_threshold)
        if not math.isfinite(threshold) or threshold < 0.0:
            raise ValueError("contact_threshold must be finite and non-negative")
        inset = float(minimum_stability_inset)
        if not math.isfinite(inset) or inset < 0.0:
            raise ValueError(
                "minimum_stability_inset must be finite and non-negative"
            )
        heights = -points @ direction
        support_height = float(heights.min()) + threshold
        # This strict comparison matches Support::updateStandingZone in the
        # released beta. A zero threshold therefore records no export-time
        # projection, as expected for an actual-contact support.
        indices = np.flatnonzero(heights < support_height).astype(np.int64)
        return cls(
            gravity=tuple(float(item) for item in direction),
            contact_threshold=threshold,
            minimum_stability_inset=inset,
            support_height=support_height,
            vertex_indices=tuple(int(item) for item in indices),
        )

    @property
    def applied(self) -> bool:
        return self.contact_threshold > 0.0 and bool(self.vertex_indices)

    def project(self, vertices: Any) -> np.ndarray:
        result = _vertices(vertices).copy()
        if not self.applied:
            return result
        indices = np.asarray(self.vertex_indices, dtype=np.int64)
        if int(indices.max()) >= result.shape[0]:
            raise ValueError(
                "support vertex indices do not address the returned outer mesh"
            )
        gravity = np.asarray(self.gravity, dtype=np.float64)
        heights = -result[indices] @ gravity
        result[indices] += (
            (heights - self.support_height)[:, np.newaxis] * gravity
        )
        return result

    def to_dict(self) -> dict[str, Any]:
        return {
            "gravity": list(self.gravity),
            "contact_threshold_normalized": self.contact_threshold,
            "minimum_stability_inset_normalized": (
                self.minimum_stability_inset
            ),
            "support_height_normalized": self.support_height,
            "support_vertex_count": len(self.vertex_indices),
            "support_vertex_indices": list(self.vertex_indices),
            "projection_applied": self.applied,
        }


def _mesh_text(vertices: np.ndarray, faces: np.ndarray, suffix: str) -> str:
    extension = suffix.lower()
    if extension == ".obj":
        lines = [
            f"v {point[0]:.17g} {point[1]:.17g} {point[2]:.17g}"
            for point in vertices
        ]
        lines.extend(
            f"f {face[0] + 1} {face[1] + 1} {face[2] + 1}"
            for face in faces
        )
        return "\n".join(lines) + "\n"
    if extension == ".stl":
        lines = ["solid make_it_stand_paper_compat"]
        for face in faces:
            triangle = vertices[face]
            normal = np.cross(triangle[1] - triangle[0], triangle[2] - triangle[0])
            length = float(np.linalg.norm(normal))
            if length > 0.0:
                normal /= length
            else:
                normal[:] = 0.0
            lines.append(
                f"  facet normal {normal[0]:.17g} {normal[1]:.17g} {normal[2]:.17g}"
            )
            lines.append("    outer loop")
            for point in triangle:
                lines.append(
                    f"      vertex {point[0]:.17g} {point[1]:.17g} {point[2]:.17g}"
                )
            lines.extend(("    endloop", "  endfacet"))
        lines.append("endsolid make_it_stand_paper_compat")
        return "\n".join(lines) + "\n"
    raise ValueError("diagnostic mesh output must end in .obj or .stl")


def _load_stl_indexed(path: os.PathLike[str] | str) -> tuple[np.ndarray, np.ndarray]:
    """Load ASCII/binary STL and exactly weld its triangle-soup vertices.

    STL has no indexed-vertex representation. libigl therefore exposes some
    valid STL files as disconnected per-facet triangle soup, which the strict
    source-boundary gate correctly rejects. Exact welding here is a format
    normalization, not a geometric repair: only bit-identical coordinates are
    merged, and the reconstructed indexed mesh still passes through the full
    native boundary audit.
    """

    source = Path(path)
    payload = source.read_bytes()
    triangles: list[tuple[tuple[float, float, float], ...]] = []
    is_binary = False
    if len(payload) >= 84:
        triangle_count = struct.unpack_from("<I", payload, 80)[0]
        is_binary = len(payload) == 84 + 50 * triangle_count
    if is_binary:
        triangle_count = struct.unpack_from("<I", payload, 80)[0]
        for index in range(triangle_count):
            values = struct.unpack_from("<12fH", payload, 84 + 50 * index)
            triangles.append(
                (
                    tuple(float(item) for item in values[3:6]),
                    tuple(float(item) for item in values[6:9]),
                    tuple(float(item) for item in values[9:12]),
                )
            )
    else:
        try:
            text = payload.decode("utf-8")
        except UnicodeDecodeError as error:
            raise ValueError("STL is neither canonical binary nor UTF-8 ASCII") from error
        vertices_in_facets: list[tuple[float, float, float]] = []
        for line_number, line in enumerate(text.splitlines(), start=1):
            fields = line.split()
            if not fields or fields[0].lower() != "vertex":
                continue
            if len(fields) != 4:
                raise ValueError(f"invalid STL vertex at line {line_number}")
            try:
                point = tuple(float(item) for item in fields[1:4])
            except ValueError as error:
                raise ValueError(
                    f"invalid STL vertex number at line {line_number}"
                ) from error
            vertices_in_facets.append(point)  # type: ignore[arg-type]
        if len(vertices_in_facets) % 3 != 0:
            raise ValueError("ASCII STL does not contain complete triangle facets")
        triangles = [
            tuple(vertices_in_facets[index:index + 3])  # type: ignore[arg-type]
            for index in range(0, len(vertices_in_facets), 3)
        ]

    if not triangles:
        raise ValueError("STL contains no triangles")

    vertices: list[tuple[float, float, float]] = []
    lookup: dict[tuple[float, float, float], int] = {}
    faces: list[tuple[int, int, int]] = []
    for triangle in triangles:
        face: list[int] = []
        for raw_point in triangle:
            # Canonicalize signed zero while preserving every other exact
            # coordinate value.
            point = tuple(0.0 if item == 0.0 else item for item in raw_point)
            if not all(math.isfinite(item) for item in point):
                raise ValueError("STL contains a non-finite vertex")
            vertex = lookup.get(point)
            if vertex is None:
                vertex = len(vertices)
                lookup[point] = vertex
                vertices.append(point)
            face.append(vertex)
        faces.append(tuple(face))  # type: ignore[arg-type]
    return _vertices(vertices), _faces(faces)


class PaperCompatibilityResult:
    """Owning result whose mesh coordinates have been restored for the caller."""

    def __init__(
        self,
        native_result: BalanceResult,
        normalization: PaperNormalization,
        flattenings: Sequence[SupportFlattening],
        native_module: Any,
    ):
        self._native_result = native_result
        self._normalization = normalization
        self._flattenings = tuple(flattenings)
        normalized_outer = _vertices(native_result.outer_vertices)
        for flattening in self._flattenings:
            normalized_outer = flattening.project(normalized_outer)
        self._normalized_outer = normalized_outer
        self._outer = normalization.to_original(normalized_outer)
        self._outer_faces = _faces(native_result.outer_faces)
        self._inner = normalization.to_original(native_result.inner_vertices)
        self._inner_faces = np.asarray(native_result.inner_faces, dtype=np.int32)
        if self._inner_faces.ndim != 2 or self._inner_faces.shape[1] != 3:
            raise ValueError("native inner faces must have shape (M, 3)")
        self._post_flatten_outer_boundary_audit = dict(
            native_module.audit_triangle_boundary(
                self._normalized_outer,
                self._outer_faces,
            )
        )

    @property
    def status(self) -> str:
        return self._native_result.status

    @property
    def succeeded(self) -> bool:
        return self._native_result.succeeded

    @property
    def native_result(self) -> BalanceResult:
        return self._native_result

    @property
    def normalization(self) -> PaperNormalization:
        return self._normalization

    @property
    def flattenings(self) -> tuple[SupportFlattening, ...]:
        return self._flattenings

    @property
    def flatten_applied(self) -> bool:
        return any(item.applied for item in self._flattenings)

    @property
    def outer_vertices(self) -> np.ndarray:
        return self._outer.copy()

    @property
    def outer_faces(self) -> np.ndarray:
        return self._outer_faces.copy()

    @property
    def inner_vertices(self) -> np.ndarray:
        return self._inner.copy()

    @property
    def inner_faces(self) -> np.ndarray:
        return self._inner_faces.copy()

    @property
    def normalized_outer_vertices(self) -> np.ndarray:
        return self._normalized_outer.copy()

    @property
    def post_flatten_outer_boundary_audit(self) -> dict[str, Any]:
        return dict(self._post_flatten_outer_boundary_audit)

    @property
    def native_manufacturing_ready(self) -> bool:
        return self._native_result.manufacturing_ready

    @property
    def native_manufacturing_failures(self) -> tuple[str, ...]:
        return self._native_result.manufacturing_failures

    @property
    def manufacturing_ready(self) -> bool:
        # The native integrity seal binds the exact normalized coordinates.
        # Even an invertible similarity (and especially support projection)
        # creates a distinct snapshot that must be audited and sealed again.
        return False

    @property
    def manufacturing_failures(self) -> tuple[str, ...]:
        failures = list(self._native_result.manufacturing_failures)
        failures.append("coordinate_restoration_requires_manufacturing_reaudit")
        if self.flatten_applied:
            failures = [
                item for item in failures
                if item != "source_hypothetical_flattening_not_applied"
            ]
            failures.append("post_export_flattening_requires_manufacturing_reaudit")
            if not bool(
                self._post_flatten_outer_boundary_audit.get(
                    "strict_boundary_valid",
                    False,
                )
            ):
                failures.append("post_export_flattened_outer_boundary_invalid")
        return tuple(dict.fromkeys(failures))

    @property
    def diagnostics(self) -> dict[str, Any]:
        return {
            "schema": "make-it-stand-paper-compat-result-v1",
            "schema_version": 1,
            "coordinate_contract": {
                "optimization": "paper_normalized",
                "returned_meshes": "source_coordinates",
            },
            "normalization": self._normalization.to_dict(),
            "support_flattenings": [item.to_dict() for item in self._flattenings],
            "post_flatten_outer_boundary_audit": (
                self.post_flatten_outer_boundary_audit
            ),
            "manufacturing_ready": self.manufacturing_ready,
            "manufacturing_failures": list(self.manufacturing_failures),
            "native_manufacturing_ready": self.native_manufacturing_ready,
            "native_manufacturing_failures": list(
                self.native_manufacturing_failures
            ),
            "native_diagnostics": self._native_result.diagnostics,
        }

    def save_outer_mesh(self, path: os.PathLike[str] | str) -> None:
        destination = Path(path)
        write_text_atomic(
            destination,
            _mesh_text(self._outer, self._outer_faces, destination.suffix),
        )

    def save_inner_mesh(self, path: os.PathLike[str] | str) -> None:
        destination = Path(path)
        write_text_atomic(
            destination,
            _mesh_text(self._inner, self._inner_faces, destination.suffix),
        )

    def save_coherent_mesh(self, path: os.PathLike[str] | str) -> None:
        """Save restored outer plus inward cavity faces as diagnostic geometry."""

        vertices = np.vstack((self._outer, self._inner))
        faces = self._outer_faces
        if self._inner_faces.shape[0] != 0:
            faces = np.vstack(
                (
                    faces,
                    self._inner_faces + self._outer.shape[0],
                )
            )
        destination = Path(path)
        write_text_atomic(
            destination,
            _mesh_text(vertices, faces, destination.suffix),
        )

    def save_diagnostics(self, path: os.PathLike[str] | str) -> None:
        write_json_atomic(path, self.diagnostics)

    def save_printable_mesh(self, path: os.PathLike[str] | str) -> None:
        raise RuntimeError(
            "paper-compatible result is not sealed manufacturing-ready; "
            + ", ".join(self.manufacturing_failures)
        )


class PaperCompatibilityModel:
    """Production wrapper for beta coordinates, support, and user handles."""

    def __init__(
        self,
        mesh: os.PathLike[str] | str | Any,
        faces: Any | None = None,
        *,
        config: ModelConfig | None = None,
        native_module: Any | None = None,
    ):
        settings = config
        if settings is None:
            settings = ModelConfig(
                path=os.fspath(mesh)
                if faces is None and isinstance(mesh, (str, os.PathLike))
                else "<memory>"
            )

        if faces is None:
            if not isinstance(mesh, (str, os.PathLike)):
                raise TypeError("faces are required for an in-memory mesh")
            if Path(mesh).suffix.lower() == ".stl":
                source_vertices, source_faces = _load_stl_indexed(mesh)
            else:
                # The native loader is the single source of truth for all
                # indexed formats. This probe exposes the exact mesh so it can
                # be normalized before the real voxel/deformation model.
                probe = BalanceModel(
                    mesh,
                    config=settings,
                    native_module=native_module,
                )
                source_vertices = _vertices(probe.native.vertices)
                source_faces = _faces(probe.native.faces)
                native_module = probe._native_module
        else:
            source_vertices = _vertices(mesh)
            source_faces = _faces(faces)

        if source_faces.size and int(source_faces.max()) >= source_vertices.shape[0]:
            raise ValueError("face indices are outside the source vertex array")

        self._normalization = PaperNormalization.from_vertices(source_vertices)
        self._source_vertices = source_vertices.copy()
        self._source_faces = source_faces.copy()
        self._normalized_vertices = self._normalization.to_normalized(
            source_vertices
        )
        normalized_settings = replace(
            settings,
            path="<paper-normalized-memory>",
            fixed_unit_cube_voxel_grid=True,
            minimum_wall_thickness=(
                settings.minimum_wall_thickness * self._normalization.scale
            ),
        )
        self._model = BalanceModel(
            self._normalized_vertices,
            self._source_faces,
            config=normalized_settings,
            native_module=native_module,
        )
        self._native_module = self._model._native_module
        self._gravity = np.asarray((0.0, -1.0, 0.0), dtype=np.float64)
        self._flattenings: list[SupportFlattening] = []

    @property
    def native(self) -> Any:
        return self._model.native

    @property
    def normalization(self) -> PaperNormalization:
        return self._normalization

    @property
    def source_vertices(self) -> np.ndarray:
        return self._source_vertices.copy()

    @property
    def normalized_vertices(self) -> np.ndarray:
        return self._normalized_vertices.copy()

    @property
    def faces(self) -> np.ndarray:
        return self._source_faces.copy()

    def set_gravity(self, gravity: Sequence[float]) -> "PaperCompatibilityModel":
        direction = _unit_vector(gravity, "gravity")
        self._model.set_gravity(direction)
        self._gravity = direction
        return self

    def set_support_from_mesh(
        self,
        contact_threshold: float = 0.02,
        gamma_radians: float = 0.0,
        *,
        beta_minimum_inset: bool = False,
    ) -> "PaperCompatibilityModel":
        if not isinstance(beta_minimum_inset, bool):
            raise TypeError("beta_minimum_inset must be a boolean")
        minimum_inset = 1.0e-4 if beta_minimum_inset else 0.0
        flattening = SupportFlattening.from_vertices(
            self._normalized_vertices,
            self._gravity,
            contact_threshold,
            minimum_inset,
        )
        self._model.set_support_from_mesh(
            contact_threshold,
            gamma_radians,
            minimum_inset,
            True,
        )
        self._flattenings = [flattening]
        return self

    def set_support_polygon(
        self,
        polygon: Sequence[Sequence[float]],
        gamma_radians: float = 0.0,
    ) -> "PaperCompatibilityModel":
        normalized = self._normalization.to_normalized(polygon)
        self._model.set_support_polygon(normalized, gamma_radians)
        self._flattenings = []
        return self

    def add_standing_constraint_from_mesh(
        self,
        gravity_direction: Sequence[float],
        contact_threshold: float = 0.02,
        gamma_radians: float = 0.0,
        *,
        beta_minimum_inset: bool = False,
    ) -> "PaperCompatibilityModel":
        direction = _unit_vector(gravity_direction, "gravity_direction")
        if not isinstance(beta_minimum_inset, bool):
            raise TypeError("beta_minimum_inset must be a boolean")
        minimum_inset = 1.0e-4 if beta_minimum_inset else 0.0
        flattening = SupportFlattening.from_vertices(
            self._normalized_vertices,
            direction,
            contact_threshold,
            minimum_inset,
        )
        self._model.add_standing_constraint_from_mesh(
            direction,
            contact_threshold,
            gamma_radians,
            minimum_inset,
            True,
        )
        self._flattenings.append(flattening)
        return self

    def add_standing_constraint_polygon(
        self,
        polygon: Sequence[Sequence[float]],
        gravity_direction: Sequence[float],
        gamma_radians: float = 0.0,
    ) -> "PaperCompatibilityModel":
        normalized = self._normalization.to_normalized(polygon)
        self._model.add_standing_constraint_polygon(
            normalized,
            gravity_direction,
            gamma_radians,
        )
        return self

    def add_handle(
        self,
        rest_pivot: Sequence[float],
        fixed_rotation: Sequence[Sequence[float]] | None = None,
    ) -> "PaperCompatibilityModel":
        self._model.add_handle(
            self._normalization.point_to_normalized(rest_pivot),
            fixed_rotation,
        )
        return self

    def add_normalized_handle(
        self,
        rest_pivot: Sequence[float],
        fixed_rotation: Sequence[Sequence[float]] | None = None,
    ) -> "PaperCompatibilityModel":
        """Add an official ``.mis`` pivot that is already beta-normalized."""

        self._model.add_handle(rest_pivot, fixed_rotation)
        return self

    def add_suspension_constraint(
        self,
        attachment_point: Sequence[float],
        allowed_deviation_radians: float,
    ) -> "PaperCompatibilityModel":
        self._model.add_suspension_constraint(
            self._normalization.point_to_normalized(attachment_point),
            allowed_deviation_radians,
        )
        return self

    def clear_handles(self) -> "PaperCompatibilityModel":
        self._model.clear_handles()
        return self

    def clear_constraints(self) -> "PaperCompatibilityModel":
        self._model.clear_constraints()
        self._flattenings = []
        return self

    def optimize(self, config: RunConfig | None = None) -> PaperCompatibilityResult:
        return PaperCompatibilityResult(
            self._model.optimize(config),
            self._normalization,
            self._flattenings,
            self._native_module,
        )

    @staticmethod
    def beta_run_config(**overrides: Any) -> RunConfig:
        """Get a fresh released-beta profile with optional bounded overrides."""

        return released_beta_run_config(**overrides)

    def optimize_released_beta(
        self,
        **overrides: Any,
    ) -> PaperCompatibilityResult:
        """Optimize with the explicit released-beta compatibility profile."""

        return self.optimize(released_beta_run_config(**overrides))


__all__ = [
    "PAPER_RUN_CONFIG",
    "PaperCompatibilityModel",
    "PaperCompatibilityResult",
    "PaperNormalization",
    "RELEASED_BETA_RUN_CONFIG",
    "SupportFlattening",
    "paper_run_config",
    "released_beta_run_config",
]
