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
SPIN_MODES = frozenset({
    "STATIC_BALANCE",
    "BUOYANT_STABILITY",
    "FIDGET_SPINNER",
    "YO_YO",
    "SPINNING_TOP",
})
SPIN_EVIDENCE = frozenset({"PAPER_FORMULA_REPRODUCTION", "RECONSTRUCTED_ASSET"})
STATIC_OBJECTIVES = frozenset({
    "MAXIMUM_MASS", "LOWER_CENTER_OF_MASS", "PARETO_TRADEOFF",
})


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
    equation_profile: str = "LEGACY_NATIVE_NORMALIZED_V1"

    def __post_init__(self) -> None:
        if self.equation_profile not in {"LEGACY_NATIVE_NORMALIZED_V1", "DUAL_MONOMIAL_RAW_V1"}:
            raise ValueError("unsupported equation_profile")
        if self.equation_profile == "DUAL_MONOMIAL_RAW_V1" and self.maximum_function_evaluations_per_restart > 60:
            raise ValueError("strict equation profile permits at most 60 evaluations per restart")
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
class SpinEvaluatorConfig:
    """Explicit native geometry controls; omitting this keeps legacy defaults."""
    derivative_method: str = "LEGACY_CONTINUOUS_IMPLICIT_SHAPE"
    icosphere_subdivisions: int = 3
    paraboloid_grid_resolution: int = 24
    cylinder_angular_resolution: int = 96
    truncation_radius: float = 0.0
    maximum_radius_doublings: int = 8
    retain_fixed_normal_provenance: bool = False
    minimum_interface_gradient: float = 1e-12

    def __post_init__(self):
        if self.derivative_method not in {"LEGACY_CONTINUOUS_IMPLICIT_SHAPE", "DISCRETE_ELLIPSOID"}:
            raise ValueError("unsupported derivative_method")
        for name, minimum in (("icosphere_subdivisions", 0), ("paraboloid_grid_resolution", 1),
                              ("cylinder_angular_resolution", 8), ("maximum_radius_doublings", 0)):
            value = getattr(self, name)
            if type(value) is not int or value < minimum:
                raise ValueError(name + " must be an integer in the native supported range")
        if self.icosphere_subdivisions > 7:
            raise ValueError("icosphere_subdivisions exceeds native maximum 7")
        if type(self.retain_fixed_normal_provenance) is not bool:
            raise ValueError("retain_fixed_normal_provenance must be bool")
        for name in ("truncation_radius", "minimum_interface_gradient"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(name + " must be finite")
            if value < 0 or (name == "minimum_interface_gradient" and value == 0):
                raise ValueError(name + " is outside native supported range")
            object.__setattr__(self, name, float(value))


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
    static_pareto_weight: float = 0.0
    evidence: str = "PAPER_FORMULA_REPRODUCTION"
    material_density: float = 1.0
    original_units_per_normalized_unit: float = 1.0
    minimum_wall_thickness: float = 0.0
    center_of_mass_weight: float = 0.0
    liquid_density: float = 0.0
    displaced_liquid_volume: float = 0.0
    center_of_buoyancy_target: tuple[float, float, float] = (0.0, 0.0, 0.0)
    waterline_area: float = 0.0
    waterline_central_second_moment: tuple[tuple[float, float], ...] = (
        (0.0, 0.0),
        (0.0, 0.0),
    )
    fixed_shell_volume: float = 0.0
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
        object.__setattr__(
            self,
            "center_of_buoyancy_target",
            _finite_tuple(
                self.center_of_buoyancy_target,
                3,
                "center_of_buoyancy_target",
            ),
        )
        waterline = tuple(
            _finite_tuple(row, 2, "waterline_central_second_moment")
            for row in self.waterline_central_second_moment
        )
        if len(waterline) != 2:
            raise ValueError(
                "waterline_central_second_moment must contain two rows"
            )
        if abs(waterline[0][1] - waterline[1][0]) > 1e-12:
            raise ValueError(
                "waterline_central_second_moment must be symmetric"
            )
        object.__setattr__(self, "waterline_central_second_moment", waterline)
        for name in (
            "material_density",
            "original_units_per_normalized_unit",
            "minimum_wall_thickness",
            "center_of_mass_weight",
            "static_pareto_weight",
            "liquid_density",
            "displaced_liquid_volume",
            "waterline_area",
            "fixed_shell_volume",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)
        if self.static_objective == "PARETO_TRADEOFF":
            if not 0.0 <= self.static_pareto_weight <= 1.0:
                raise ValueError("static_pareto_weight must lie in [0,1]")
        elif self.static_pareto_weight != 0.0:
            raise ValueError(
                "static_pareto_weight is supported only for PARETO_TRADEOFF"
            )
        if self.mode == "BUOYANT_STABILITY":
            if self.liquid_density <= 0.0:
                raise ValueError("liquid_density must be positive")
            if self.displaced_liquid_volume <= 0.0:
                raise ValueError("displaced_liquid_volume must be positive")
            if self.waterline_area <= 0.0:
                raise ValueError("waterline_area must be positive")
            if self.fixed_shell_volume < 0.0:
                raise ValueError("fixed_shell_volume must be nonnegative")
            if min(np.linalg.eigvalsh(np.asarray(waterline))) <= 0.0:
                raise ValueError(
                    "waterline_central_second_moment must be positive definite"
                )
        elif any((
            self.liquid_density,
            self.displaced_liquid_volume,
            *self.center_of_buoyancy_target,
            self.waterline_area,
            *waterline[0],
            *waterline[1],
            self.fixed_shell_volume,
        )):
            raise ValueError(
                "buoyancy inputs require mode='BUOYANT_STABILITY'"
            )


@dataclass(frozen=True)
class SpinQuadricEvaluation:
    variable_raw_moments: tuple[float, ...]
    total_raw_moments: tuple[float, ...]
    geometry_hash: str
    topology_hash: str
    interface_type: str
    interface_subtype: str
    base_request_hash: str = ""
    request_hash: str = ""
    derivative_method: str = ""
    derivative_defined: bool = False
    derivative_diagnostic: str = ""
    derivative_scientific_hash: str = ""
    provenance_topology_hash: str = ""
    vertex_provenance: tuple[str, ...] = ()
    material_vertices: tuple[tuple[float, ...], ...] = ()
    material_faces: tuple[tuple[int, ...], ...] = ()
    outer_vertices: tuple[tuple[float, ...], ...] = ()
    outer_faces: tuple[tuple[int, ...], ...] = ()


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
    equation_profile: str = "LEGACY_NATIVE_NORMALIZED_V1"
    solver_coordinate_identity: str = "native_matrix_entry_unit_sphere_v1"
    raw_equation_residual: float | None = None
    quadric_coefficients: tuple[float, ...] | None = None

    def __post_init__(self) -> None:
        if self.quadric_coefficients is not None:
            object.__setattr__(self, "quadric_coefficients", _finite_tuple(
                self.quadric_coefficients, 10, "trace quadric coefficients"))


@dataclass(frozen=True)
class SpinSolverEvent:
    type: str
    restart: int
    iteration: int
    evaluation_index: int
    restart_evaluation_index: int
    step: float | None
    outcome: str
    reason: str
    detail: str
    quadric_coefficients: tuple[float, ...] | None
    quadric_hash: str
    normalized_equation_residual: float | None
    physical_relative_residual: float | None
    objective: float | None
    multiplier_mu: float | None
    material_geometry_hash: str | None
    material_topology_hash: str | None
    equation_profile: str = "LEGACY_NATIVE_NORMALIZED_V1"
    solver_coordinate_identity: str = "native_matrix_entry_unit_sphere_v1"
    raw_equation_residual: float | None = None


def _immutable_spin_solver_event(entry: dict[str, Any]) -> SpinSolverEvent:
    coefficients = entry["quadric_coefficients"]
    immutable_coefficients: tuple[float, ...] | None = None
    if coefficients is not None:
        immutable_coefficients = tuple(float(value) for value in coefficients)
        if (
            len(immutable_coefficients) != 10
            or not all(math.isfinite(value) for value in immutable_coefficients)
        ):
            raise ValueError(
                "native Spin-It Faster event coefficients must be 10 finite values"
            )
    return SpinSolverEvent(
        **{**entry, "quadric_coefficients": immutable_coefficients}
    )


class SpinItFasterResult:
    """Read-only Python view over an immutable native result snapshot."""

    __slots__ = ("_native", "_diagnostics", "_trace", "_events")

    def __init__(self, native: Any):
        diagnostics = json.loads(str(native.diagnostics_json))
        if diagnostics.get("schema") != SPIN_ARTIFACT_SCHEMA or diagnostics.get("schema_version") != 1:
            raise ValueError("native Spin-It Faster diagnostics use an unsupported schema")
        object.__setattr__(self, "_native", native)
        object.__setattr__(self, "_diagnostics", diagnostics)
        object.__setattr__(self, "_trace", tuple(SpinIteration(**entry) for entry in diagnostics["trace"]))
        object.__setattr__(
            self,
            "_events",
            tuple(
                _immutable_spin_solver_event(entry)
                for entry in diagnostics["events"]
            ),
        )

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
    def residual_evaluation_seconds(self) -> float:
        return float(self._native.residual_evaluation_seconds)

    @property
    def derivative_evaluation_seconds(self) -> float:
        return float(self._native.derivative_evaluation_seconds)

    @property
    def solve_evaluator_wall_seconds(self) -> float:
        return float(self._native.solve_evaluator_wall_seconds)

    @property
    def solve_function_evaluations(self) -> int:
        return int(self._native.solve_function_evaluations)

    @property
    def post_solve_material_evaluation_count(self) -> int:
        return int(self._native.post_solve_material_evaluation_count)

    @property
    def post_solve_material_evaluation_wall_seconds(self) -> float:
        return float(self._native.post_solve_material_evaluation_wall_seconds)

    @property
    def event_hash(self) -> str:
        return str(self._native.event_hash)

    @property
    def equation_profile(self) -> str:
        return self._diagnostics.get("equation_profile", "LEGACY_NATIVE_NORMALIZED_V1")

    @property
    def solver_coordinate_identity(self) -> str:
        return self._diagnostics.get("solver_coordinate_identity", "native_matrix_entry_unit_sphere_v1")

    @property
    def raw_equation_residual(self) -> float | None:
        return self._diagnostics.get("raw_equation_residual")

    @property
    def evaluation_timings(self) -> tuple[dict[str, Any], ...]:
        return tuple(dict(entry) for entry in self._native.evaluation_timings)

    @property
    def diagnostics(self) -> dict[str, Any]:
        return json.loads(json.dumps(self._diagnostics))

    @property
    def derivative_metadata(self) -> dict[str, str]:
        return {"method": str(self._native.derivative_method),
                "diagnostic": str(self._native.derivative_diagnostic),
                "scientific_hash": str(self._native.derivative_scientific_hash)}

    @property
    def final_raw_moments(self) -> tuple[float, ...] | None:
        """Actual converged export in native target coordinates, or no final claim."""
        return tuple(self._diagnostics['final_raw_moments']) if self.succeeded else None

    @property
    def solution_coefficients(self) -> tuple[float, ...] | None:
        return tuple(self._diagnostics['quadric_coefficients']) if self.succeeded else None

    @property
    def trace(self) -> tuple[SpinIteration, ...]:
        return self._trace

    @property
    def events(self) -> tuple[SpinSolverEvent, ...]:
        return self._events

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
    def design_vertices(self) -> np.ndarray:
        value = np.asarray(self._native.design_vertices).copy()
        value.flags.writeable = False
        return value

    @property
    def design_faces(self) -> np.ndarray:
        value = np.asarray(self._native.design_faces).copy()
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

    def save_design_mesh(self, path: os.PathLike[str] | str) -> None:
        self._native.save_design_mesh(os.fspath(path))

    def save_material_mesh(self, path: os.PathLike[str] | str) -> None:
        self._native.save_material_mesh(os.fspath(path))

    def save_cavity_mesh(self, path: os.PathLike[str] | str) -> None:
        self._native.save_cavity_mesh(os.fspath(path))

    def save_diagnostics(self, path: os.PathLike[str] | str) -> None:
        self._native.save_diagnostics(os.fspath(path))

    def save_viewer_bundle(self, directory: os.PathLike[str] | str) -> tuple[str, ...]:
        root = os.fspath(directory)
        os.makedirs(root, exist_ok=True)
        paths = tuple(os.path.join(root, name) for name in (
            "outer.obj", "design.obj", "material.obj", "cavity.obj", "diagnostics.json"
        ))
        self.save_outer_mesh(paths[0])
        self.save_design_mesh(paths[1])
        self.save_material_mesh(paths[2])
        self.save_cavity_mesh(paths[3])
        self.save_diagnostics(paths[4])
        return paths


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
        value = getattr(request.settings, name)
        if name == "equation_profile":
            value = getattr(native.SpinEquationProfile, value)
        setattr(settings, name, value)
    result.outer_surface = outer
    result.design_domain = design
    result.fixed_shell_moments = native.RawMassMoments10.from_values(
        np.asarray(request.fixed_shell_moments, dtype=np.float64)
    )
    result.target_frame = frame
    result.mode = getattr(native.SpinMode, request.mode)
    result.static_objective = getattr(native.StaticObjective, request.static_objective)
    result.static_pareto_weight = request.static_pareto_weight
    result.material_density = request.material_density
    result.original_units_per_normalized_unit = request.original_units_per_normalized_unit
    result.minimum_wall_thickness = request.minimum_wall_thickness
    result.center_of_mass_weight = request.center_of_mass_weight
    buoyancy = native.BuoyantStabilityData()
    buoyancy.liquid_density = request.liquid_density
    buoyancy.displaced_liquid_volume = request.displaced_liquid_volume
    buoyancy.center_of_buoyancy_target = np.asarray(
        request.center_of_buoyancy_target,
        dtype=np.float64,
    )
    buoyancy.waterline_area = request.waterline_area
    buoyancy.waterline_central_second_moment = np.asarray(
        request.waterline_central_second_moment,
        dtype=np.float64,
    )
    buoyancy.fixed_shell_volume = request.fixed_shell_volume
    result.buoyancy = buoyancy
    result.settings = settings
    return native, result


def _native_evaluator_settings(native, config):
    if not isinstance(config, SpinEvaluatorConfig):
        raise ValueError("expected immutable SpinEvaluatorConfig")
    if not hasattr(native, "SpinEvaluatorSettings"):
        raise RuntimeError("native Faster evaluator settings support is unavailable")
    settings = native.SpinEvaluatorSettings()
    triangulation = native.SpinTriangulationSettings()
    for name in SpinEvaluatorConfig.__dataclass_fields__:
        value = getattr(config, name)
        if name == "derivative_method":
            settings.derivative_method = getattr(native.SpinDerivativeMethod, value)
        elif name == "minimum_interface_gradient":
            settings.minimum_interface_gradient = value
        else:
            setattr(triangulation, name, value)
    settings.triangulation = triangulation
    return settings


def evaluate_spin_it_faster_quadric(
    request: SpinItFasterRequest,
    coefficients: Sequence[float],
    *, evaluator_settings: SpinEvaluatorConfig | None = None,
) -> SpinQuadricEvaluation:
    native, native_request = _native_request(request)
    values = _finite_tuple(coefficients, 10, "quadric coefficients")
    extra = () if evaluator_settings is None else (_native_evaluator_settings(native, evaluator_settings),)
    evaluated = native.evaluate_spin_it_faster_quadric(
        native_request,
        np.asarray(values, dtype=np.float64),
        *extra,
    )
    return SpinQuadricEvaluation(
        variable_raw_moments=tuple(float(value) for value in evaluated["variable_raw_moments"]),
        total_raw_moments=tuple(float(value) for value in evaluated["total_raw_moments"]),
        geometry_hash=str(evaluated["geometry_hash"]),
        topology_hash=str(evaluated["topology_hash"]),
        interface_type=str(evaluated["interface_type"]),
        interface_subtype=str(evaluated["interface_subtype"]),
        **({} if evaluator_settings is None else {
            **{name: evaluated[name] for name in ("base_request_hash", "request_hash", "derivative_method",
                "derivative_defined", "derivative_diagnostic", "derivative_scientific_hash", "provenance_topology_hash")},
            "vertex_provenance": tuple(evaluated["vertex_provenance"]),
            **{name: tuple(tuple(row) for row in evaluated[name].tolist()) for name in
               ("material_vertices", "material_faces", "outer_vertices", "outer_faces")},
        }),
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
    evaluator_settings: SpinEvaluatorConfig | None = None,
) -> SpinItFasterResult:
    native, native_request = _native_request(request)
    initial = np.asarray(
        tuple(_finite_tuple(row, 10, "initial quadric") for row in initial_quadrics),
        dtype=np.float64,
    ).reshape((-1, 10))
    extra = () if evaluator_settings is None else (_native_evaluator_settings(native, evaluator_settings),)
    artifact = native.solve_spin_it_faster(
        native_request,
        getattr(native.SpinItFasterEvidence, request.evidence),
        initial,
        *extra,
    )
    return SpinItFasterResult(artifact)
