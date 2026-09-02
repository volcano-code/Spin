"""Immutable Python facade for the Spin-It Faster formula reproduction."""

from __future__ import annotations

from dataclasses import dataclass, field
import importlib
import json
import math
import os
from pathlib import Path
from typing import Any, Sequence

import numpy as np


SPIN_ARTIFACT_SCHEMA = "spin-it-faster-artifact-v1"
SPIN_MODES = frozenset({"STATIC_BALANCE", "FIDGET_SPINNER", "YO_YO", "SPINNING_TOP"})
SPIN_EVIDENCE = frozenset({"PAPER_FORMULA_REPRODUCTION", "RECONSTRUCTED_ASSET"})
STATIC_OBJECTIVES = frozenset({"MAXIMUM_MASS", "LOWER_CENTER_OF_MASS"})


def _native_module() -> Any:
    errors: list[BaseException] = []
    for name in ("make_it_stand._make_it_stand", "_make_it_stand"):
        try:
            return importlib.import_module(name)
        except (ImportError, ModuleNotFoundError) as error:
            errors.append(error)
    raise RuntimeError(
        "the Spin-It Faster native extension is not installed; configure "
        "with MAKE_IT_STAND_BUILD_PYTHON=ON and "
        "MAKE_IT_STAND_BUILD_SPIN_IT_FASTER=ON"
    ) from errors[-1]


def _finite_tuple(values: Sequence[float], length: int, name: str) -> tuple[float, ...]:
    if len(values) != length:
        raise ValueError(f"{name} must contain exactly {length} values")
    result = tuple(float(value) for value in values)
    if not all(math.isfinite(value) for value in result):
        raise ValueError(f"{name} must contain only finite values")
    return result


def _vertices(values: Sequence[Sequence[float]], name: str) -> tuple[tuple[float, float, float], ...]:
    result = tuple(_finite_tuple(row, 3, name) for row in values)
    if len(result) < 4:
        raise ValueError(f"{name} must contain at least four vertices")
    return result  # type: ignore[return-value]


def _faces(values: Sequence[Sequence[int]], name: str) -> tuple[tuple[int, int, int], ...]:
    result: list[tuple[int, int, int]] = []
    for row in values:
        if len(row) != 3 or any(isinstance(value, bool) or not isinstance(value, int) for value in row):
            raise ValueError(f"{name} must contain integer triangles")
        result.append((int(row[0]), int(row[1]), int(row[2])))
    if len(result) < 4:
        raise ValueError(f"{name} must contain at least four faces")
    return tuple(result)


@dataclass(frozen=True)
class SpinSolverConfig:
    maximum_iterations: int = 60
    maximum_restarts: int = 8
    maximum_function_evaluations_per_restart: int = 60
    equation_tolerance: float = 1e-12
    physical_relative_tolerance: float = 1e-6
    deterministic_seed: int = 0

    def __post_init__(self) -> None:
        for name in (
            "maximum_iterations",
            "maximum_restarts",
            "maximum_function_evaluations_per_restart",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ValueError(f"{name} must be a positive integer")
        if isinstance(self.deterministic_seed, bool) or not isinstance(self.deterministic_seed, int):
            raise ValueError("deterministic_seed must be an integer")
        for name in ("equation_tolerance", "physical_relative_tolerance"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class SpinItFasterRequest:
    outer_vertices: tuple[tuple[float, float, float], ...]
    outer_faces: tuple[tuple[int, int, int], ...]
    design_vertices: tuple[tuple[float, float, float], ...]
    design_faces: tuple[tuple[int, int, int], ...]
    fixed_shell_moments: tuple[float, ...] = (0.0,) * 10
    target_origin: tuple[float, float, float] = (0.0, 0.0, 0.0)
    target_axes: tuple[tuple[float, float, float], ...] = (
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    )
    mode: str = "STATIC_BALANCE"
    static_objective: str = "MAXIMUM_MASS"
    evidence: str = "PAPER_FORMULA_REPRODUCTION"
    material_density: float = 1.0
    original_units_per_normalized_unit: float = 1.0
    minimum_wall_thickness: float = 0.0
    center_of_mass_weight: float = 0.0
    settings: SpinSolverConfig = field(default_factory=SpinSolverConfig)

    def __post_init__(self) -> None:
        object.__setattr__(self, "outer_vertices", _vertices(self.outer_vertices, "outer_vertices"))
        object.__setattr__(self, "outer_faces", _faces(self.outer_faces, "outer_faces"))
        object.__setattr__(self, "design_vertices", _vertices(self.design_vertices, "design_vertices"))
        object.__setattr__(self, "design_faces", _faces(self.design_faces, "design_faces"))
        object.__setattr__(self, "fixed_shell_moments", _finite_tuple(self.fixed_shell_moments, 10, "fixed_shell_moments"))
        object.__setattr__(self, "target_origin", _finite_tuple(self.target_origin, 3, "target_origin"))
        axes = tuple(_finite_tuple(row, 3, "target_axes") for row in self.target_axes)
        if len(axes) != 3:
            raise ValueError("target_axes must contain exactly three rows")
        object.__setattr__(self, "target_axes", axes)
        if self.mode not in SPIN_MODES:
            raise ValueError(f"unsupported Spin-It Faster mode: {self.mode!r}")
        if self.static_objective not in STATIC_OBJECTIVES:
            raise ValueError(f"unsupported static objective: {self.static_objective!r}")
        if self.evidence not in SPIN_EVIDENCE:
            raise ValueError(f"unsupported evidence level: {self.evidence!r}")
        for name in (
            "material_density",
            "original_units_per_normalized_unit",
            "minimum_wall_thickness",
            "center_of_mass_weight",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class SpinQuadricEvaluation:
    variable_raw_moments: tuple[float, ...]
    total_raw_moments: tuple[float, ...]
    geometry_hash: str
    topology_hash: str
    interface_type: str


@dataclass(frozen=True)
class SpinIteration:
    restart: int
    iteration: int
    function_evaluations: int
    normalized_equation_residual: float
    physical_relative_residual: float
    objective: float
    multiplier_mu: float
    accepted_geodesic_step: float
    quadric_hash: str
    material_geometry_hash: str
    material_topology_hash: str


class SpinItFasterResult:
    """Read-only Python view over an immutable native result snapshot."""

    __slots__ = ("_native", "_diagnostics", "_trace")

    def __init__(self, native: Any):
        diagnostics = json.loads(str(native.diagnostics_json))
        if diagnostics.get("schema") != SPIN_ARTIFACT_SCHEMA or diagnostics.get("schema_version") != 1:
            raise ValueError("native Spin-It Faster diagnostics use an unsupported schema")
        object.__setattr__(self, "_native", native)
        object.__setattr__(self, "_diagnostics", diagnostics)
        object.__setattr__(self, "_trace", tuple(SpinIteration(**entry) for entry in diagnostics["trace"]))

    def __setattr__(self, name: str, value: Any) -> None:
        raise AttributeError("SpinItFasterResult is immutable")

    @property
    def status(self) -> str:
        return str(self._native.status)

    @property
    def evidence(self) -> str:
        return str(self._native.evidence)

    @property
    def succeeded(self) -> bool:
        return bool(self._native.succeeded)

    @property
    def geometry_sealed(self) -> bool:
        return bool(self._native.geometry_sealed)

    @property
    def diagnostics(self) -> dict[str, Any]:
        return json.loads(json.dumps(self._diagnostics))

    @property
    def trace(self) -> tuple[SpinIteration, ...]:
        return self._trace

    @property
    def outer_vertices(self) -> np.ndarray:
        value = np.asarray(self._native.outer_vertices).copy()
        value.flags.writeable = False
        return value

    @property
    def outer_faces(self) -> np.ndarray:
        value = np.asarray(self._native.outer_faces).copy()
        value.flags.writeable = False
        return value

    @property
    def material_vertices(self) -> np.ndarray:
        value = np.asarray(self._native.material_vertices).copy()
        value.flags.writeable = False
        return value

    @property
    def material_faces(self) -> np.ndarray:
        value = np.asarray(self._native.material_faces).copy()
        value.flags.writeable = False
        return value

    @property
    def cavity_vertices(self) -> np.ndarray:
        value = np.asarray(self._native.cavity_vertices).copy()
        value.flags.writeable = False
        return value

    @property
    def cavity_faces(self) -> np.ndarray:
        value = np.asarray(self._native.cavity_faces).copy()
        value.flags.writeable = False
        return value

    def save_outer_mesh(self, path: os.PathLike[str] | str) -> None:
        self._native.save_outer_mesh(os.fspath(path))

    def save_material_mesh(self, path: os.PathLike[str] | str) -> None:
        self._native.save_material_mesh(os.fspath(path))

    def save_cavity_mesh(self, path: os.PathLike[str] | str) -> None:
        self._native.save_cavity_mesh(os.fspath(path))

    def save_diagnostics(self, path: os.PathLike[str] | str) -> None:
        self._native.save_diagnostics(os.fspath(path))


def _native_request(request: SpinItFasterRequest) -> tuple[Any, Any]:
    native = _native_module()
    result = native.SpinDesignRequest()
    outer = native.SpinTriangleMesh()
    outer.vertices = np.asarray(request.outer_vertices, dtype=np.float64)
    outer.faces = np.asarray(request.outer_faces, dtype=np.int32)
    design = native.SpinTriangleMesh()
    design.vertices = np.asarray(request.design_vertices, dtype=np.float64)
    design.faces = np.asarray(request.design_faces, dtype=np.int32)
    frame = native.SpinTargetFrame()
    frame.origin = np.asarray(request.target_origin, dtype=np.float64)
    frame.axes = np.asarray(request.target_axes, dtype=np.float64)
    settings = native.SpinSolverSettings()
    for name in SpinSolverConfig.__dataclass_fields__:
        setattr(settings, name, getattr(request.settings, name))
    result.outer_surface = outer
    result.design_domain = design
    result.fixed_shell_moments = native.RawMassMoments10.from_values(
        np.asarray(request.fixed_shell_moments, dtype=np.float64)
    )
    result.target_frame = frame
    result.mode = getattr(native.SpinMode, request.mode)
    result.static_objective = getattr(native.StaticObjective, request.static_objective)
    result.material_density = request.material_density
    result.original_units_per_normalized_unit = request.original_units_per_normalized_unit
    result.minimum_wall_thickness = request.minimum_wall_thickness
    result.center_of_mass_weight = request.center_of_mass_weight
    result.settings = settings
    return native, result


def evaluate_spin_it_faster_quadric(
    request: SpinItFasterRequest,
    coefficients: Sequence[float],
) -> SpinQuadricEvaluation:
    native, native_request = _native_request(request)
    values = _finite_tuple(coefficients, 10, "quadric coefficients")
    evaluated = native.evaluate_spin_it_faster_quadric(
        native_request,
        np.asarray(values, dtype=np.float64),
    )
    return SpinQuadricEvaluation(
        variable_raw_moments=tuple(float(value) for value in evaluated["variable_raw_moments"]),
        total_raw_moments=tuple(float(value) for value in evaluated["total_raw_moments"]),
        geometry_hash=str(evaluated["geometry_hash"]),
        topology_hash=str(evaluated["topology_hash"]),
        interface_type=str(evaluated["interface_type"]),
    )


def integrate_raw_mass_moments(
    vertices: Sequence[Sequence[float]],
    faces: Sequence[Sequence[int]],
    density: float = 1.0,
) -> tuple[float, ...]:
    native = _native_module()
    normalized_vertices = _vertices(vertices, "vertices")
    normalized_faces = _faces(faces, "faces")
    values = native.integrate_raw_mass_moments(
        np.asarray(normalized_vertices, dtype=np.float64),
        np.asarray(normalized_faces, dtype=np.int32),
        float(density),
    )
    return tuple(float(value) for value in values)


def solve_spin_it_faster(
    request: SpinItFasterRequest,
    *,
    initial_quadrics: Sequence[Sequence[float]] = (),
) -> SpinItFasterResult:
    native, native_request = _native_request(request)
    initial = np.asarray(
        tuple(_finite_tuple(row, 10, "initial quadric") for row in initial_quadrics),
        dtype=np.float64,
    ).reshape((-1, 10))
    artifact = native.solve_spin_it_faster(
        native_request,
        getattr(native.SpinItFasterEvidence, request.evidence),
        initial,
    )
    return SpinItFasterResult(artifact)
