"""Immutable Python facade for the Spin-It Faster formula reproduction."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Any, Callable, Sequence

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


@dataclass(frozen=True)
class Spin2014QuadricObjectiveConfig:
    """Original 2014 formulas in declared units; moments already carry density."""
    objective: str = "TOP_MASS_AND_LOWERING"
    inertia_weight: float = 1.0
    center_weight: float = 1.0
    declared_length_unit: str = "model_length"
    declared_mass_unit: str = "model_mass"
    solver_mass_unit_in_declared_units: float = 1.0
    solver_length_unit_in_declared_units: float = 1.0

    def __post_init__(self) -> None:
        if self.objective not in {"YOYO", "TOP_INERTIA_ONLY", "TOP_LOWERING_ONLY",
                                  "TOP_MASS_AND_LOWERING", "STATIC_LOWERING"}:
            raise ValueError("unknown Spin-It 2014 objective")
        for name in ("declared_length_unit", "declared_mass_unit"):
            if not isinstance(getattr(self, name), str) or not getattr(self, name):
                raise ValueError(name + " must be a nonempty unit label")
        for name in ("inertia_weight", "center_weight",
                     "solver_mass_unit_in_declared_units", "solver_length_unit_in_declared_units"):
            x = getattr(self, name)
            if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x):
                raise ValueError(name + " must be finite")
            if x < 0 or (name.startswith("solver_") and x == 0):
                raise ValueError(name + " is outside the supported range")
            object.__setattr__(self, name, float(x))


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
    capture_evaluation_payload: bool = False

    def __post_init__(self) -> None:
        if self.equation_profile not in {"LEGACY_NATIVE_NORMALIZED_V1", "DUAL_MONOMIAL_RAW_V1"}:
            raise ValueError("unsupported equation_profile")
        if self.equation_profile == "DUAL_MONOMIAL_RAW_V1" and self.maximum_function_evaluations_per_restart > 60:
            raise ValueError("strict equation profile permits at most 60 evaluations per restart")
        if type(self.capture_evaluation_payload) is not bool:
            raise ValueError("capture_evaluation_payload must be bool")
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
    maximum_icosphere_subdivisions: int = -1
    paraboloid_grid_resolution: int = 24
    cylinder_angular_resolution: int = 96
    hyperbolic_cylinder_grid_resolution: int = 24
    parabolic_cylinder_grid_resolution: int = 24
    truncation_radius: float = 0.0
    maximum_radius_doublings: int = 8
    retain_fixed_normal_provenance: bool = False
    minimum_interface_gradient: float = 1e-12
    cap_radius_policy: str = "REJECT_UNCERTIFIED"
    paper_projective_hyperboloid: bool = False
    algebraic_finite_map: bool = False
    fixed_normal_paraboloid: bool = False
    exact_binary64_self_intersection: bool = False
    exact_binary64_source_boundary: bool = False
    exact_binary64_material_boundary: bool = False
    rational_source_coordinates: bool = False
    compact_rational_source_coordinates: bool = False

    def __post_init__(self):
        if type(self.compact_rational_source_coordinates) is not bool:
            raise ValueError("compact_rational_source_coordinates must be bool")
        if self.compact_rational_source_coordinates and not self.rational_source_coordinates:
            raise ValueError("compact rational source requires rational_source_coordinates")
        if type(self.rational_source_coordinates) is not bool:
            raise ValueError("rational_source_coordinates must be bool")
        if self.rational_source_coordinates and (
                self.derivative_method != 'DISCRETE_FULL_RANK_QUADRIC' or
                self.algebraic_finite_map or self.paper_projective_hyperboloid or
                self.fixed_normal_paraboloid):
            raise ValueError('rational source requires an unmixed discrete full-rank profile')
        if type(self.exact_binary64_self_intersection) is not bool:
            raise ValueError("exact_binary64_self_intersection must be bool")
        if type(self.exact_binary64_source_boundary) is not bool:
            raise ValueError("exact_binary64_source_boundary must be bool")
        if type(self.exact_binary64_material_boundary) is not bool:
            raise ValueError("exact_binary64_material_boundary must be bool")
        if self.derivative_method not in {"LEGACY_CONTINUOUS_IMPLICIT_SHAPE", "DISCRETE_FULL_RANK_QUADRIC",
                                          "DISCRETE_REDUCED_STRATUM", "DISCRETE_ELLIPSOID"}:
            raise ValueError("unsupported derivative_method")
        if self.cap_radius_policy not in {"REJECT_UNCERTIFIED", "GROW_TO_CERTIFIED_V1"}:
            raise ValueError("unsupported cap_radius_policy")
        for name, minimum in (("icosphere_subdivisions", 0),
                              ("maximum_icosphere_subdivisions", -1),
                              ("paraboloid_grid_resolution", 1),
                              ("cylinder_angular_resolution", 8),
                              ("hyperbolic_cylinder_grid_resolution", 2),
                              ("parabolic_cylinder_grid_resolution", 1),
                              ("maximum_radius_doublings", 0)):
            value = getattr(self, name)
            if type(value) is not int or value < minimum:
                raise ValueError(name + " must be an integer in the native supported range")
        if self.icosphere_subdivisions > 7:
            raise ValueError("icosphere_subdivisions exceeds native maximum 7")
        if self.maximum_icosphere_subdivisions > 7:
            raise ValueError("maximum_icosphere_subdivisions exceeds native maximum 7")
        if (self.maximum_icosphere_subdivisions != -1 and
                self.maximum_icosphere_subdivisions < self.icosphere_subdivisions):
            raise ValueError(
                "maximum_icosphere_subdivisions must be -1 or at least icosphere_subdivisions"
            )
        if type(self.retain_fixed_normal_provenance) is not bool:
            raise ValueError("retain_fixed_normal_provenance must be bool")
        if type(self.paper_projective_hyperboloid) is not bool:
            raise ValueError("paper_projective_hyperboloid must be bool")
        if type(self.algebraic_finite_map) is not bool:
            raise ValueError("algebraic_finite_map must be bool")
        if type(self.fixed_normal_paraboloid) is not bool:
            raise ValueError("fixed_normal_paraboloid must be bool")
        if self.fixed_normal_paraboloid and (
                self.derivative_method != 'DISCRETE_REDUCED_STRATUM' or
                self.paper_projective_hyperboloid or self.algebraic_finite_map or
                self.maximum_icosphere_subdivisions not in (-1,self.icosphere_subdivisions)):
            raise ValueError('fixed-normal paraboloid requires an unmixed reduced fixed-normal profile')
        for name in ("truncation_radius", "minimum_interface_gradient"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(name + " must be finite")
            if value < 0 or (name == "minimum_interface_gradient" and value == 0):
                raise ValueError(name + " is outside native supported range")
            object.__setattr__(self, name, float(value))

    def identity_dict(self):
        """Keep saved legacy payloads unchanged unless exact boundary policy is opted in.

        Each true option is bound to new evidence. Source admission and solver
        material-region proof are separate; neither changes a mesh or confers
        manufacturing approval.
        """
        result = asdict(self)
        if not self.exact_binary64_source_boundary:
            result.pop("exact_binary64_source_boundary")
        if not self.exact_binary64_material_boundary:
            result.pop("exact_binary64_material_boundary")
        if not self.rational_source_coordinates:
            result.pop("rational_source_coordinates")
        if not self.compact_rational_source_coordinates:
            result.pop("compact_rational_source_coordinates")
        return result


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
    fixed_shell_density: float = 0.0
    settings: SpinSolverConfig = field(default_factory=SpinSolverConfig)
    static_objective_measure: str = "DENSITY_WEIGHTED_MASS_V1"
    spinning_top_objective: str = "INERTIA_AND_LOWERING"

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
        if self.spinning_top_objective not in ("INERTIA_AND_LOWERING","LOWER_CENTER_OF_MASS"):
            raise ValueError("unsupported spinning-top objective")
        if self.spinning_top_objective == "LOWER_CENTER_OF_MASS" and (
                self.mode != "SPINNING_TOP" or self.center_of_mass_weight != 0.0):
            raise ValueError("pure-height objective requires spinning-top mode and zero extra weight")
        if self.static_objective not in STATIC_OBJECTIVES:
            raise ValueError(f"unsupported static objective: {self.static_objective!r}")
        if self.static_objective_measure not in ("DENSITY_WEIGHTED_MASS_V1","UNIFORM_DENSITY_VOLUME_V1"):
            raise ValueError("unsupported static objective measure")
        if self.static_objective_measure != "DENSITY_WEIGHTED_MASS_V1" and self.mode != "STATIC_BALANCE":
            raise ValueError("volume objective measure is supported only for static balance")
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
            "fixed_shell_density",
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
            if self.fixed_shell_density < 0.0:
                raise ValueError("fixed_shell_density must be nonnegative")
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
            self.fixed_shell_density,
        )):
            raise ValueError(
                "buoyancy inputs require mode='BUOYANT_STABILITY'"
            )

    @property
    def request_digest(self) -> str:
        """Canonical digest of this immutable request, including solver settings."""
        from .verification.spin_common import digest
        return digest(asdict(self))


def with_measured_fixed_shell_geometry(
    request: SpinItFasterRequest,
    vertices: Sequence[Sequence[float]],
    faces: Sequence[Sequence[int]],
    *,
    density: float | None = None,
) -> SpinItFasterRequest:
    """Return a request whose fixed-shell moments are integrated from a mesh.

    The request API stores fixed-shell mass moments rather than a second mesh.
    This adapter closes that input path from actual geometry. If the fixed
    shell uses the same material as the design region, omitting ``density``
    retains the historical single-density behavior.
    """
    if not isinstance(request, SpinItFasterRequest):
        raise TypeError("request must be a SpinItFasterRequest")
    material_density = request.material_density if density is None else density
    if isinstance(material_density, bool) or not isinstance(material_density, (int, float)) or \
            not math.isfinite(material_density) or material_density <= 0.0:
        raise ValueError("fixed-shell density must be finite and positive")
    moments = integrate_raw_mass_moments(vertices, faces, density=float(material_density))
    return replace(request, fixed_shell_moments=moments)


def with_measured_buoyancy_geometry(
    request: SpinItFasterRequest,
    displacement_vertices: Sequence[Sequence[float]],
    displacement_faces: Sequence[Sequence[int]],
    fixed_shell_vertices: Sequence[Sequence[float]],
    fixed_shell_faces: Sequence[Sequence[int]],
    *,
    waterline_origin: Sequence[float],
    waterline_normal: Sequence[float],
    liquid_density: float,
    fixed_shell_density: float | None = None,
) -> SpinItFasterRequest:
    """Build prescribed-waterline inputs from actual closed geometry.

    ``displacement_*`` must be the caller-designated closed boundary for the
    displacement model (including the intended sealed/floodable semantics).
    The waterline normal must point along the request target frame's positive
    vertical axis because the Faster planar family is parameterized as
    ``a + bx*x + by*y - z``.  Fixed-shell moments and volume are independently
    integrated from ``fixed_shell_*`` at ``fixed_shell_density`` (defaulting
    to the design material density for legacy calls) in model coordinates.
    Native evaluation transforms these fixed moments to the target frame once;
    hydrostatic quantities below are already expressed in the target frame.

    This constructs measured solver inputs; it does not infer sealed cavities,
    solve the plane family, prove physical stability, or attest geometry
    provenance.  Those remain separate evidence requirements.
    """
    if not isinstance(request, SpinItFasterRequest):
        raise TypeError("request must be a SpinItFasterRequest")
    if isinstance(liquid_density, bool) or not isinstance(liquid_density, (int, float)) or \
            not math.isfinite(liquid_density) or liquid_density <= 0.0:
        raise ValueError("liquid_density must be finite and positive")
    shell_density = request.material_density if fixed_shell_density is None else fixed_shell_density
    if isinstance(shell_density, bool) or not isinstance(shell_density, (int, float)) or \
            not math.isfinite(shell_density) or shell_density <= 0.0:
        raise ValueError("fixed_shell_density must be finite and positive")
    axes = np.asarray(request.target_axes, dtype=np.float64)
    origin = np.asarray(request.target_origin, dtype=np.float64)
    if not np.allclose(axes.T @ axes, np.eye(3), rtol=0.0, atol=1e-10) or \
            not math.isclose(float(np.linalg.det(axes)), 1.0, rel_tol=0.0, abs_tol=1e-10):
        raise ValueError("target_axes must be a right-handed orthonormal frame")

    world_origin = np.asarray(_finite_tuple(waterline_origin, 3, "waterline_origin"))
    world_normal = np.asarray(_finite_tuple(waterline_normal, 3, "waterline_normal"))
    normal_length = float(np.linalg.norm(world_normal))
    if not math.isfinite(normal_length) or normal_length <= 0.0:
        raise ValueError("waterline_normal must be nonzero and finite")
    target_normal = axes.T @ (world_normal / normal_length)
    if not np.allclose(target_normal, (0.0, 0.0, 1.0), rtol=0.0, atol=1e-10):
        raise ValueError("waterline normal must align with target-frame vertical")

    displacement = _vertices(displacement_vertices, "displacement_vertices")
    faces = _faces(displacement_faces, "displacement_faces")
    target_displacement = (np.asarray(displacement) - origin) @ axes
    target_waterline_origin = (world_origin - origin) @ axes
    hydrostatic = compute_hydrostatic_geometry(
        target_displacement,
        faces,
        waterline_origin=target_waterline_origin,
        waterline_normal=(0.0, 0.0, 1.0),
    )

    fixed_vertices = _vertices(fixed_shell_vertices, "fixed_shell_vertices")
    fixed_faces = _faces(fixed_shell_faces, "fixed_shell_faces")
    fixed_moments = integrate_raw_mass_moments(
        fixed_vertices, fixed_faces, density=float(shell_density)
    )
    fixed_volume = float(integrate_raw_mass_moments(
        fixed_vertices, fixed_faces, density=1.0
    )[0])
    return replace(
        request,
        fixed_shell_moments=fixed_moments,
        mode="BUOYANT_STABILITY",
        liquid_density=float(liquid_density),
        displaced_liquid_volume=hydrostatic.displaced_volume,
        center_of_buoyancy_target=hydrostatic.center_of_buoyancy,
        waterline_area=hydrostatic.waterline_area,
        waterline_central_second_moment=hydrostatic.waterline_central_second_moment,
        fixed_shell_volume=fixed_volume,
        fixed_shell_density=float(shell_density),
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
    stratum_derivative_defined: bool = False
    stratum_derivative_diagnostic: str = ""
    stratum_moment_jacobian: tuple[tuple[float, ...], ...] = ()
    stratum_coefficient_basis: tuple[tuple[float, ...], ...] = ()
    derivative_scientific_hash: str = ""
    provenance_topology_hash: str = ""
    boolean_backend_identity: str = ""
    quadric_solid_is_positive_region: bool = False
    vertex_provenance: tuple[str, ...] = ()
    interface_vertices: tuple[tuple[float, ...], ...] = ()
    interface_faces: tuple[tuple[int, ...], ...] = ()
    material_vertices: tuple[tuple[float, ...], ...] = ()
    material_faces: tuple[tuple[int, ...], ...] = ()
    outer_vertices: tuple[tuple[float, ...], ...] = ()
    outer_faces: tuple[tuple[int, ...], ...] = ()


@dataclass(frozen=True)
class SpinQuadricTeacherDiagnostic:
    """Same-call moments and ambient derivative from the default evaluator.

    This diagnostic export preserves an undefined derivative as ``None``.
    It is not a certificate that a defined derivative matches finite
    differences; callers must retain and separately validate that distinction.
    """

    coefficients: tuple[float, ...]
    variable_raw_moments: tuple[float, ...]
    total_raw_moments: tuple[float, ...]
    fixed_raw_moments: tuple[float, ...]
    moment_jacobian: tuple[tuple[float, ...], ...] | None
    derivative_defined: bool
    derivative_method: str
    derivative_diagnostic: str
    equation_profile: str
    geometry_hash: str
    topology_hash: str
    base_request_hash: str
    evaluator_request_hash: str
    derivative_scientific_hash: str
    provenance_topology_hash: str
    residual_evaluation_seconds: float
    derivative_evaluation_seconds: float


@dataclass(frozen=True)
class SpinFullRankQuadricEvaluation:
    """Discrete Boolean-to-moment derivative on a supported full-rank stratum."""
    variable_raw_moments: tuple[float, ...]
    total_raw_moments: tuple[float, ...]
    geometry_hash: str
    boolean_backend_identity: str
    source_scientific_hash: str
    scientific_hash: str
    derivative_defined: bool
    derivative_diagnostic: str
    coefficient_jacobian: tuple[tuple[float, ...], ...] | None
    material_vertices: tuple[tuple[float, ...], ...]
    material_faces: tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class SpinFullRankQuadricSolveResult:
    """Opt-in result from the canonical discrete-Jacobian Faster solver."""
    status: str
    diagnostic: str
    final_geometry_verified: bool
    verification_diagnostic: str
    input_hash: str
    output_hash: str
    function_evaluations: int
    restart_count: int
    normalized_equation_residual: float
    physical_relative_residual: float
    independent_raw_residual: float | None
    multiplier_mu: float
    multipliers: tuple[float, ...]
    quadric_coefficients: tuple[float, ...] | None
    equation_profile: str
    event_hash: str
    trace_hash: str
    trace: tuple[SpinIteration, ...]
    events: tuple[SpinSolverEvent, ...]
    restart_function_evaluations: tuple[int, ...]
    evaluation_timings: tuple[SpinFullRankEvaluationTiming, ...]
    final_geometry_hash: str | None
    final_derivative_defined: bool
    final_derivative_diagnostic: str
    final_variable_raw_moments: tuple[float, ...] | None
    final_total_raw_moments: tuple[float, ...] | None
    independent_total_raw_moments: tuple[float, ...] | None
    material_vertices: tuple[tuple[float, ...], ...] | None
    material_faces: tuple[tuple[int, ...], ...] | None


class _RationalAuthorityView:
    @property
    def rational_authority(self):
        """Fresh view of the same-call attachment; not a geometry/seal verdict."""
        return json.loads(self.rational_authority_json)


@dataclass(frozen=True)
class _RationalSpinQuadricEvaluation(_RationalAuthorityView, SpinQuadricEvaluation):
    rational_authority_json: str = ""


@dataclass(frozen=True)
class _RationalSpinFullRankEvaluation(_RationalAuthorityView, SpinFullRankQuadricEvaluation):
    rational_authority_json: str = ""
    normalized_coefficients: tuple[float, ...] = ()
    branch_identity: str = ""


class _MaterialProofView:
    @property
    def material_proof(self):
        """Independent view of the same-call tape, not a validation verdict."""
        return json.loads(self.material_proof_json)


@dataclass(frozen=True)
class _ProofSpinFullRankEvaluation(_MaterialProofView, SpinFullRankQuadricEvaluation):
    material_proof_json: str = ""


@dataclass(frozen=True)
class _ProofRationalSpinFullRankEvaluation(_MaterialProofView, _RationalSpinFullRankEvaluation):
    material_proof_json: str = ""


def _same_call_material_proof_json(value, design_face_count):
    """Validate attachment coverage/identity only; no second geometry operation."""
    from fractions import Fraction
    proof = value.get('material_proof_attachment')
    if not isinstance(proof, dict) or proof.get('schema') != 'same-call-cgal-material-proof-v1' or \
            not value.get('geometry_hash') or proof.get('geometry_hash') != value['geometry_hash']:
        raise RuntimeError('same-call material proof identity missing or mismatched')
    points = proof.get('exact_material_vertices', ())
    birth = np.asarray(proof.get('material_birth_faces', ()))
    source = np.asarray(proof.get('source_faces', ()))
    vertices = np.asarray(value['material_vertices'])
    faces = np.asarray(value['material_faces'])
    if len(points) != len(vertices) or birth.shape != (len(faces),) or \
            source.ndim != 2 or source.shape[1] != 3 or not len(source) or \
            not np.issubdtype(birth.dtype, np.integer) or not np.issubdtype(source.dtype, np.integer) or \
            np.any(source < 0) or np.any(birth < 0) or np.any(birth >= design_face_count + len(source)):
        raise RuntimeError('same-call material proof coverage invalid')
    try:
        canonical = [[str(Fraction(x)) for x in row] for row in points]
        rounded = np.asarray([[float(Fraction(x)) for x in row] for row in canonical])
    except (ValueError, TypeError, ZeroDivisionError, OverflowError) as exc:
        raise RuntimeError('same-call material proof coordinates invalid') from exc
    if rounded.shape != vertices.shape or not np.array_equal(rounded, vertices):
        raise RuntimeError('same-call material proof coordinate view mismatch')
    topology = proof.get('provenance_topology_hash')
    if not isinstance(topology, str) or not topology:
        raise RuntimeError('same-call material proof provenance missing')
    return json.dumps(dict(schema=proof['schema'], geometry_hash=proof['geometry_hash'],
        provenance_topology_hash=topology, exact_material_vertices=canonical,
        material_birth_faces=birth.tolist(), source_faces=source.tolist()),
        sort_keys=True, separators=(',', ':'), allow_nan=False)


@dataclass(frozen=True)
class _RationalSpinFullRankSolveResult(_RationalAuthorityView, SpinFullRankQuadricSolveResult):
    rational_authority_json: str = ""


def _rational_result_type(legacy, extended, value, *, required=False):
    """Extend only opt-in results; asdict of legacy results remains unchanged.

    The attachment stays canonical JSON rather than a caller-mutable dict.
    Actual material checks are still done by the independent seal/checker.
    """
    authority = value.get('rational_authority')
    if authority is None and 'rational_authority_json' in value:
        authority = json.loads(value['rational_authority_json'])
    if authority is None:
        if required:
            raise RuntimeError('native rational source authority attachment missing')
        return legacy, {}
    expected = value.get('geometry_hash', value.get('final_geometry_hash'))
    if not isinstance(authority, dict) or authority.get('schema') != \
            'rational-quadric-production-authority-v1' or \
            authority.get('coordinate_space') != 'TARGET_FRAME' or \
            authority.get('analytic_coordinates_exact') is not False or \
            not expected or authority.get('geometry_hash') != expected:
        raise RuntimeError('native rational source authority attachment identity mismatch')
    return extended, {'rational_authority_json': json.dumps(authority,
        sort_keys=True, separators=(',', ':'), allow_nan=False)}


@dataclass(frozen=True)
class HydrostaticGeometry:
    """Actual submerged closed-mesh and waterline measurements.

    This is a geometric observation, not a claim that a buoyant Faster input
    has an author-published waterline, geometry, or stable physical motion.
    """
    displaced_volume: float
    center_of_buoyancy: tuple[float, float, float]
    waterline_area: float
    waterline_centroid: tuple[float, float, float]
    waterline_central_second_moment: tuple[tuple[float, float], tuple[float, float]]
    submerged_vertices: tuple[tuple[float, float, float], ...]
    submerged_faces: tuple[tuple[int, int, int], ...]


@dataclass(frozen=True)
class HydrostaticEquilibrium:
    """A scalar volume-matched waterline for an explicitly supplied boundary.

    It is deliberately not a buoyancy optimization result: the caller still
    owns the displacement-boundary/floodability model and must independently
    assess stability at this fixed orientation.
    """
    geometry: HydrostaticGeometry
    waterline_signed_offset: float
    target_displaced_volume: float
    volume_residual: float
    iterations: int


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
    evaluation_payload: SpinSolverEvaluationPayload | None = None
    newton_system_payload: SpinNewtonSystemPayload | None = None


@dataclass(frozen=True)
class SpinNewtonSystemPayload:
    """Saved linear-system diagnostics, including unsuccessful rank checks."""
    schema: str
    matrix: tuple[tuple[float, ...], ...]
    right_hand_side: tuple[float, ...]
    computed_rank: int
    equilibrated_rank: int
    multiplier_mu: float | None
    constraint_multipliers: tuple[float, ...]
    row_scales: tuple[float, ...]
    column_scales: tuple[float, ...]
    linear_residual_infinity_norm: float
    linear_backward_error: float


@dataclass(frozen=True)
class SpinSolverEvaluationPayload:
    schema: str
    solver_coordinates: tuple[float, ...]
    solver_raw_moments: tuple[float, ...]
    variable_raw_moments: tuple[float, ...] | None
    solver_moment_jacobian: tuple[tuple[float, ...], ...] | None
    derivative_defined: bool
    jacobian_coordinate_identity: str
    evaluation_map_identity: str
    boolean_backend_identity: str
    derivative_identity: str
    derivative_diagnostic: str
    geometry_hash: str
    topology_hash: str
    evaluator_wall_seconds: float
    residual_evaluation_seconds: float
    derivative_evaluation_seconds: float


@dataclass(frozen=True)
class SpinFullRankEvaluationTiming:
    evaluation_index: int
    restart: int
    restart_evaluation_index: int
    evaluator_returned: bool
    wall_seconds: float
    reported_residual_seconds: float | None
    reported_derivative_seconds: float | None


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
    payload_entry = entry.get("evaluation_payload")
    payload = None
    if payload_entry is not None:
        coordinates = tuple(float(value) for value in payload_entry["solver_coordinates"])
        solver_moments = tuple(
            float(value) for value in payload_entry["solver_raw_moments"]
        )
        raw_variable_moments = payload_entry["variable_raw_moments"]
        variable_moments = None if raw_variable_moments is None else tuple(
            float(value) for value in raw_variable_moments
        )
        raw_jacobian = payload_entry["solver_moment_jacobian"]
        jacobian = None if raw_jacobian is None else tuple(
            tuple(float(value) for value in row) for row in np.asarray(raw_jacobian).tolist()
        )
        derivative_defined = bool(payload_entry["derivative_defined"])
        if (len(coordinates) != 10 or len(solver_moments) != 10 or
                (variable_moments is not None and len(variable_moments) != 10)):
            raise ValueError("captured evaluation payload must contain ten coordinates and moment values")
        if not all(math.isfinite(value) for value in coordinates + solver_moments):
            raise ValueError("captured evaluation payload coordinates/moments must be finite")
        if variable_moments is not None and not all(
            math.isfinite(value) for value in variable_moments
        ):
            raise ValueError("captured variable moments must be finite")
        if derivative_defined != (jacobian is not None):
            raise ValueError("captured Jacobian availability disagrees with derivative_defined")
        if jacobian is not None and (
            len(jacobian) != 10 or any(len(row) != 10 for row in jacobian) or
            not all(math.isfinite(value) for row in jacobian for value in row)
        ):
            raise ValueError("captured solver Jacobian must be a finite 10-by-10 matrix")
        times = tuple(float(payload_entry[name]) for name in (
            "evaluator_wall_seconds", "residual_evaluation_seconds",
            "derivative_evaluation_seconds",
        ))
        if not all(math.isfinite(value) and value >= 0.0 for value in times):
            raise ValueError("captured evaluation payload timings must be finite and nonnegative")
        payload = SpinSolverEvaluationPayload(
            schema=str(payload_entry["schema"]),
            solver_coordinates=coordinates,
            solver_raw_moments=solver_moments,
            variable_raw_moments=variable_moments,
            solver_moment_jacobian=jacobian,
            derivative_defined=derivative_defined,
            jacobian_coordinate_identity=str(payload_entry["jacobian_coordinate_identity"]),
            evaluation_map_identity=str(payload_entry["evaluation_map_identity"]),
            boolean_backend_identity=str(payload_entry["boolean_backend_identity"]),
            derivative_identity=str(payload_entry["derivative_identity"]),
            derivative_diagnostic=str(payload_entry["derivative_diagnostic"]),
            geometry_hash=str(payload_entry["geometry_hash"]),
            topology_hash=str(payload_entry["topology_hash"]),
            evaluator_wall_seconds=times[0],
            residual_evaluation_seconds=times[1],
            derivative_evaluation_seconds=times[2],
        )
        if payload.schema != "spin-evaluation-payload-v1":
            raise ValueError("unsupported captured evaluation payload schema")
    system_entry = entry.get("newton_system_payload")
    system = None
    if system_entry is not None:
        if system_entry["schema"] != "spin-newton-system-payload-v2":
            raise ValueError("unsupported captured Newton system schema")
        # Diagnostic infinities in rejected systems must remain diagnostic
        # failures, not be discarded or promoted to valid scientific outputs.
        system = SpinNewtonSystemPayload(
            schema=str(system_entry["schema"]),
            matrix=tuple(tuple(float(x) for x in row) for row in system_entry["matrix"]),
            right_hand_side=tuple(float(x) for x in system_entry["right_hand_side"]),
            computed_rank=int(system_entry["computed_rank"]),
            equilibrated_rank=int(system_entry["equilibrated_rank"]),
            multiplier_mu=(None if system_entry["multiplier_mu"] is None
                           else float(system_entry["multiplier_mu"])),
            constraint_multipliers=tuple(float(x) for x in system_entry["constraint_multipliers"]),
            row_scales=tuple(float(x) for x in system_entry["row_scales"]),
            column_scales=tuple(float(x) for x in system_entry["column_scales"]),
            linear_residual_infinity_norm=float(system_entry["linear_residual_infinity_norm"]),
            linear_backward_error=float(system_entry["linear_backward_error"]))
    return SpinSolverEvent(
        **{**entry, "quadric_coefficients": immutable_coefficients,
           "evaluation_payload": payload, "newton_system_payload": system}
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
    result.static_objective_measure = getattr(native.StaticObjectiveMeasure,request.static_objective_measure)
    result.material_density = request.material_density
    result.original_units_per_normalized_unit = request.original_units_per_normalized_unit
    result.minimum_wall_thickness = request.minimum_wall_thickness
    result.center_of_mass_weight = request.center_of_mass_weight
    result.spinning_top_objective = getattr(native.SpinningTopObjective,request.spinning_top_objective)
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
    buoyancy.fixed_shell_density = request.fixed_shell_density
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
        elif name == "cap_radius_policy":
            settings.cap_radius_policy = getattr(native.SpinCapRadiusPolicy, value)
        elif name == "minimum_interface_gradient":
            settings.minimum_interface_gradient = value
        elif name == "exact_binary64_self_intersection":
            if hasattr(settings,name):
                setattr(settings,name,value)
            elif value:
                raise RuntimeError("native exact binary64 self-intersection policy unavailable")
        elif name == "exact_binary64_source_boundary":
            if hasattr(settings, name):
                setattr(settings, name, value)
            elif value:
                raise RuntimeError("native exact binary64 source boundary policy unavailable")
        elif name == "exact_binary64_material_boundary":
            if hasattr(settings, name):
                setattr(settings, name, value)
            elif value:
                raise RuntimeError("native exact binary64 material boundary policy unavailable")
        elif name == "paper_projective_hyperboloid" and not hasattr(triangulation, name):
            if value:
                raise RuntimeError("native paper-projective hyperboloid support is unavailable")
            # Old native modules still implement the requested legacy default.
        elif name == "algebraic_finite_map" and not hasattr(triangulation, name):
            if value:
                raise RuntimeError("native algebraic finite-map support is unavailable")
        elif name == "fixed_normal_paraboloid" and not hasattr(triangulation, name):
            if value:
                raise RuntimeError("native fixed-normal paraboloid support is unavailable")
        elif name == "rational_source_coordinates" and not hasattr(triangulation, name):
            if value:
                raise RuntimeError("native rational source coordinate support is unavailable")
        elif name == "compact_rational_source_coordinates" and not hasattr(triangulation, name):
            if value:
                raise RuntimeError("native compact rational source support is unavailable")
        else:
            setattr(triangulation, name, value)
    settings.triangulation = triangulation
    return settings


def solve_spin_it_faster_buoyant_plane(
    request: SpinItFasterRequest,
    *,
    initial_parameters: Sequence[float] = (0.0, 0.0, 0.0),
    maximum_iterations: int = 60,
    maximum_line_search_steps: int = 24,
    constraint_tolerance: float = 1e-10,
) -> dict[str, Any]:
    """Run Faster's prescribed-waterline planar-family solver.

    The request should normally come from ``with_measured_buoyancy_geometry``.
    The solver operates on the actual plane/Boolean material evaluation in the
    request design domain. ``equilibrium_constraints_verified`` and the
    separate ``stable_at_prescribed_waterline`` result are intentionally not
    conflated. This is not a free-waterline or fluid-dynamics simulation.
    """
    if not isinstance(request, SpinItFasterRequest):
        raise TypeError("request must be a SpinItFasterRequest")
    if request.mode != "BUOYANT_STABILITY":
        raise ValueError("request mode must be BUOYANT_STABILITY")
    if type(maximum_iterations) is not int or maximum_iterations <= 0:
        raise ValueError("maximum_iterations must be a positive integer")
    if type(maximum_line_search_steps) is not int or maximum_line_search_steps <= 0:
        raise ValueError("maximum_line_search_steps must be a positive integer")
    if isinstance(constraint_tolerance, bool) or not isinstance(constraint_tolerance, (int, float)) or \
            not math.isfinite(constraint_tolerance) or constraint_tolerance <= 0.0:
        raise ValueError("constraint_tolerance must be finite and positive")
    initial = np.asarray(_finite_tuple(initial_parameters, 3, "initial_parameters"), dtype=np.float64)

    native, native_request = _native_request(request)
    solver = getattr(native, "solve_buoyant_stability_by_plane", None)
    settings_type = getattr(native, "BuoyantPlaneBalanceSettings", None)
    if solver is None or settings_type is None:
        raise RuntimeError("native prescribed-waterline plane solver is unavailable")
    settings = settings_type()
    settings.maximum_iterations = maximum_iterations
    settings.maximum_line_search_steps = maximum_line_search_steps
    settings.constraint_tolerance = float(constraint_tolerance)
    value = solver(native_request, initial, settings)
    result = dict(value)
    for key in ("plane_coefficients", "total_raw_moments", "variable_raw_moments",
                "constraints", "normalized_constraints", "variable_material_vertices",
                "variable_material_faces"):
        result[key] = np.asarray(result[key]).copy()
    result["trace"] = tuple(dict(event) for event in result["trace"])
    return result


def evaluate_spin_it_faster_plane_halfspaces(
    request: SpinItFasterRequest,
    *,
    coefficients: Sequence[float],
) -> dict[str, Any]:
    """Return both oriented half-spaces of a plane as an exact design partition.

    Each side is independently evaluated by the native plane/Boolean path.
    The returned boundary meshes are in target coordinates; exact partition
    checking is performed after transforming them back into request
    coordinates. This is useful for sealing the retained and removed regions
    of a solved buoyancy plane. It does not seal an artifact or infer fluid
    semantics.
    """
    if not isinstance(request, SpinItFasterRequest):
        raise TypeError("request must be a SpinItFasterRequest")
    coeff = np.asarray(_finite_tuple(coefficients, 4, "coefficients"), dtype=np.float64)
    if float(np.linalg.norm(coeff)) <= 1e-15:
        raise ValueError("plane coefficients must be nonzero")
    native, native_request = _native_request(request)
    evaluate = getattr(native, "evaluate_spin_it_faster_plane_halfspace", None)
    exact_region = getattr(native, "make_exact_material_region", None)
    partition_check = getattr(native, "check_exact_material_partition", None)
    if evaluate is None or exact_region is None or partition_check is None:
        raise RuntimeError("native exact plane-partition support is unavailable")

    positive = dict(evaluate(native_request, coeff))
    negative = dict(evaluate(native_request, -coeff))
    axes = np.asarray(request.target_axes, dtype=np.float64)
    origin = np.asarray(request.target_origin, dtype=np.float64)

    # Independent Boolean calls for c and -c can choose different diagonals
    # for the same rounded cut polygon. The polygon vertices can agree while
    # the resulting piecewise-planar surfaces differ by a tiny exact sliver.
    # Reuse one side's interface triangulation (with opposite winding) on the
    # other side. No vertex is moved or welded approximately; exact partition
    # verification below remains the deciding gate.
    def stitch_shared_interface(first, second):
        first_vertices = np.asarray(first["vertices"], dtype=np.float64)
        second_vertices = np.asarray(second["vertices"], dtype=np.float64)
        first_faces = np.asarray(first["faces"], dtype=np.int32)
        second_faces = np.asarray(second["faces"], dtype=np.int32)
        plane = np.asarray(coeff, dtype=np.float64)
        normal_norm = float(np.linalg.norm(plane[1:]))
        if not math.isfinite(normal_norm) or normal_norm <= 1e-15:
            raise ValueError("plane halfspace coefficients have an invalid normal")
        distance_tolerance = 1e-10 * max(
            1.0,
            float(np.max(np.abs(first_vertices))) if first_vertices.size else 1.0,
            float(np.max(np.abs(second_vertices))) if second_vertices.size else 1.0,
        )

        def interface_face_ids(vertices, faces):
            signed = (plane[0] + vertices @ plane[1:]) / normal_norm
            return [index for index, face in enumerate(faces)
                    if float(np.max(np.abs(signed[face]))) <= distance_tolerance]

        first_ids = interface_face_ids(first_vertices, first_faces)
        second_ids = interface_face_ids(second_vertices, second_faces)
        if not first_ids or not second_ids:
            raise ValueError("both complementary meshes must expose a plane interface")

        # Exact coordinate identity is required for every shared interface
        # vertex. Approximate nearest-vertex matching could silently alter the
        # material geometry and is deliberately not used.
        second_vertex_ids = {
            tuple(float(value) for value in point): index
            for index, point in enumerate(second_vertices)
        }
        mapped_faces = []
        for face_id in first_ids:
            mapped = []
            for point in first_vertices[first_faces[face_id]]:
                key = tuple(float(value) for value in point)
                if key not in second_vertex_ids:
                    raise ValueError("complementary plane interfaces do not share exact vertex coordinates")
                mapped.append(second_vertex_ids[key])
            mapped_faces.append(mapped)

        def patch_boundary_edges(vertices, faces, face_ids):
            edge_count = {}
            for face_id in face_ids:
                face = faces[face_id]
                for corner in range(3):
                    a = tuple(float(value) for value in vertices[int(face[corner])])
                    b = tuple(float(value) for value in vertices[int(face[(corner + 1) % 3])])
                    edge = tuple(sorted((a, b)))
                    edge_count[edge] = edge_count.get(edge, 0) + 1
            if any(count not in (1, 2) for count in edge_count.values()):
                raise ValueError("plane interface patch has a non-manifold triangulation")
            return sorted(edge for edge, count in edge_count.items() if count == 1)

        if patch_boundary_edges(first_vertices, first_faces, first_ids) != \
                patch_boundary_edges(second_vertices, second_faces, second_ids):
            raise ValueError("complementary plane interfaces have different exact boundary edges")

        for face in mapped_faces:
            tri = second_vertices[face]
            if float(np.linalg.norm(np.cross(tri[1] - tri[0], tri[2] - tri[0]))) == 0.0:
                raise ValueError("degenerate plane interface triangle")
        # The two material interiors lie on opposite sides of the same plane,
        # so the outward normal of their common interface must be reversed.
        mapped_faces = [face[::-1] for face in mapped_faces]

        removed = set(second_ids)
        second["faces"] = np.asarray(
            [face.tolist() for index, face in enumerate(second_faces) if index not in removed] +
            mapped_faces,
            dtype=np.int32,
        )
        return {"canonical_interface_triangles": len(mapped_faces),
                "removed_interface_triangles": len(second_ids),
                "vertex_coordinates_moved": False,
                "exact_coordinate_matching": True}

    interface_stitch = stitch_shared_interface(positive, negative)

    def model_region(result):
        target_vertices = np.asarray(result["vertices"], dtype=np.float64)
        model_vertices = target_vertices @ axes.T + origin
        return exact_region(model_vertices, np.asarray(result["faces"], dtype=np.int32))

    design_region = exact_region(
        np.asarray(request.design_vertices, dtype=np.float64),
        np.asarray(request.design_faces, dtype=np.int32),
    )
    positive_region, negative_region = model_region(positive), model_region(negative)
    audit = partition_check(design_region, [positive_region, negative_region])
    if not audit["completed"] or not audit["complete_partition"] or \
            not audit["no_positive_volume_overlap"]:
        raise ValueError("plane half-spaces do not form an exact complementary partition: " +
                         str(audit))

    design_target = (np.asarray(request.design_vertices, dtype=np.float64) - origin) @ axes
    expected = np.asarray(integrate_raw_mass_moments(
        design_target, request.design_faces, request.material_density), dtype=np.float64)
    actual_positive = np.asarray(positive["raw_moments"], dtype=np.float64)
    actual_negative = np.asarray(negative["raw_moments"], dtype=np.float64)
    span = float(np.ptp(design_target, axis=0).max())
    density = float(request.material_density)
    scales = np.asarray((density * span**3,) + (density * span**4,) * 3 +
                        (density * span**5,) * 6, dtype=np.float64)
    normalized_error = float(np.max(np.abs(actual_positive + actual_negative - expected) / scales))
    if not math.isfinite(normalized_error) or normalized_error > 1e-10:
        raise ValueError("plane partition moment sum differs from the design domain")

    # The returned meshes, rather than just evaluator-side analytical values,
    # are the objects that will be sealed. Independently reintegrate each and
    # require agreement with its native halfspace evaluation.
    positive_mesh_moments = np.asarray(integrate_raw_mass_moments(
        positive["vertices"], positive["faces"].tolist(), density), dtype=np.float64)
    negative_mesh_moments = np.asarray(integrate_raw_mass_moments(
        negative["vertices"], negative["faces"].tolist(), density), dtype=np.float64)
    positive_mesh_error = float(np.max(np.abs(positive_mesh_moments - actual_positive) / scales))
    negative_mesh_error = float(np.max(np.abs(negative_mesh_moments - actual_negative) / scales))
    if max(positive_mesh_error, negative_mesh_error) > 1e-10:
        raise ValueError("stitched plane meshes disagree with native halfspace mass moments")
    stitched_sum_error = float(np.max(np.abs(
        positive_mesh_moments + negative_mesh_moments - expected) / scales))
    if stitched_sum_error > 1e-10:
        raise ValueError("stitched plane mesh moments do not sum to the design-domain moments")

    def expose(value):
        return {
            "vertices": np.asarray(value["vertices"], dtype=np.float64).copy(),
            "faces": np.asarray(value["faces"], dtype=np.int32).copy(),
            "raw_moments": tuple(float(x) for x in value["raw_moments"]),
            "geometry_hash": str(value["geometry_hash"]),
            "cap_separation_certificate": str(value["cap_separation_certificate"]),
        }

    return {
        "status": "EXACT_COMPLEMENTARY_PARTITION",
        "positive": expose(positive),
        "negative": expose(negative),
        "interface_stitch": interface_stitch,
        "independent_mesh_moment_error": max(positive_mesh_error, negative_mesh_error),
        "independent_partition_moment_error": stitched_sum_error,
        "exact_partition_audit": dict(audit),
        "normalized_moment_partition_error": normalized_error,
        "design_raw_moments": tuple(float(x) for x in expected),
        "coefficient_convention": "BETA_GT_ZERO_IN_TARGET_FRAME_V1",
    }


def solve_spin_it_faster_static_plane_family(
    request: SpinItFasterRequest,
    *,
    initial_coefficients: Sequence[float],
    maximum_iterations: int = 60,
    maximum_line_search_steps: int = 24,
    constraint_tolerance: float = 1e-10,
) -> dict[str, Any]:
    """Solve the Faster static plane-family problem for its selected objective.

    The request's static objective dispatches to maximum mass, lower center of
    mass, or the weighted Pareto KKT system. The returned variable boundary is
    the actual native Boolean result; callers should independently integrate
    it before treating the result as a material artifact. The reported KKT
    flag is numerical verification for this plane family, not a global
    optimality certificate.

    Legacy calls weight mass. ``UNIFORM_DENSITY_VOLUME_V1`` explicitly weights
    geometric volume as in Figure 11 and validates the actual fixed shell's
    uniform-density moments; arbitrary supplied fixed mass is not sufficient.
    Strict full-domain optima carry a separate ``boundary_optimality`` witness;
    the actual accepted plane and its interior-equation residual are preserved.
    """
    if not isinstance(request, SpinItFasterRequest):
        raise TypeError("request must be a SpinItFasterRequest")
    if request.mode != "STATIC_BALANCE":
        raise ValueError("request mode must be STATIC_BALANCE")
    if type(maximum_iterations) is not int or maximum_iterations < 0:
        raise ValueError("maximum_iterations must be a nonnegative integer")
    if type(maximum_line_search_steps) is not int or maximum_line_search_steps <= 0:
        raise ValueError("maximum_line_search_steps must be a positive integer")
    if isinstance(constraint_tolerance, bool) or not isinstance(constraint_tolerance, (int, float)) or \
            not math.isfinite(constraint_tolerance) or constraint_tolerance <= 0.0:
        raise ValueError("constraint_tolerance must be finite and positive")
    initial = np.asarray(
        _finite_tuple(initial_coefficients, 4, "initial_coefficients"),
        dtype=np.float64,
    )

    native, native_request = _native_request(request)
    solver = getattr(native, "solve_static_balance_by_plane_family", None)
    if solver is None:
        raise RuntimeError("native static plane-family solver is unavailable")
    value = solver(
        native_request,
        initial,
        maximum_iterations,
        maximum_line_search_steps,
        float(constraint_tolerance),
    )
    result = dict(value)
    for key in (
        "plane_coefficients", "lambda", "constraints", "normalized_constraints",
        "variable_raw_moments", "total_raw_moments", "variable_material_vertices",
        "variable_material_faces",
    ):
        result[key] = np.asarray(result[key]).copy()
    result["trace"] = tuple(dict(event) for event in result["trace"])
    return result


def solve_spin_it_faster_planar_fidget(
    *,
    design_region_loops: Sequence[Sequence[Sequence[float]]],
    fixed_material_loops: Sequence[Sequence[Sequence[float]]] = (),
    fixed_areal_density: float | None = None,
    design_areal_density: float,
    initial_center: Sequence[float],
    maximum_iterations: int = 100,
    maximum_function_evaluations: int = 60,
    residual_tolerance: float = 1e-10,
    initial_trust_radius: float = 0.25,
    maximum_trust_radius: float = 1.0,
    extrusion_interval: Sequence[float] | None = None,
    extrusion_volume_density: float | None = None,
    fixed_extrusion_layers: Sequence[dict[str, Any]] = (),
) -> dict[str, Any]:
    """Solve Faster's 2D disk-hollowing balance problem from actual loops.

    Boundary loops follow the native convention: outer boundaries are
    counter-clockwise and holes clockwise, with material on the left. Fixed
    material and design-region moments are integrated separately from their
    supplied polygon geometry; callers do not hand-fill moment vectors.
    The exact solved material boundary is returned as oriented line/circular-
    arc cycles. If ``extrusion_interval`` and ``extrusion_volume_density`` are
    supplied, exact volume moments for its analytic extrusion are returned.
    Fixed layers are exact polygonal prisms and can directly define the fixed
    moments used by the optimizer. If legacy fixed loops are also supplied,
    their moments are checked against the layer projection.
    """
    def normalize_loops(values: Sequence[Sequence[Sequence[float]]], name: str):
        loops = []
        for loop_index, loop in enumerate(values):
            if len(loop) < 3:
                raise ValueError(f"{name}[{loop_index}] must have at least three points")
            points = []
            for point_index, point in enumerate(loop):
                if len(point) != 2 or any(isinstance(value, (bool, np.bool_)) for value in point):
                    raise ValueError(f"{name}[{loop_index}][{point_index}] must be a finite 2D point")
                points.append(_finite_tuple(point, 2, f"{name}[{loop_index}][{point_index}]"))
            loops.append(tuple(points))
        if not loops:
            raise ValueError(f"{name} must contain at least one boundary loop")
        return tuple(loops)

    for name, value in (("design_areal_density", design_areal_density),):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or \
                not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    if fixed_areal_density is not None and (
            isinstance(fixed_areal_density, bool) or
            not isinstance(fixed_areal_density, (int, float)) or
            not math.isfinite(fixed_areal_density) or fixed_areal_density <= 0.0):
        raise ValueError("fixed_areal_density must be finite and positive when supplied")
    for name, value in (("maximum_iterations", maximum_iterations),
                        ("maximum_function_evaluations", maximum_function_evaluations)):
        if type(value) is not int or value <= 0:
            raise ValueError(f"{name} must be a positive integer")
    for name, value in (("residual_tolerance", residual_tolerance),
                        ("initial_trust_radius", initial_trust_radius),
                        ("maximum_trust_radius", maximum_trust_radius)):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or \
                not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    if maximum_trust_radius < initial_trust_radius:
        raise ValueError("maximum_trust_radius must be at least initial_trust_radius")

    design_loops = normalize_loops(design_region_loops, "design_region_loops")
    if any(isinstance(value, (bool, np.bool_)) for value in initial_center):
        raise ValueError("initial_center must contain finite real coordinates, not booleans")
    center = _finite_tuple(initial_center, 2, "initial_center")
    if (extrusion_interval is None) != (extrusion_volume_density is None):
        raise ValueError("extrusion_interval and extrusion_volume_density must be supplied together")
    native_extrusion = None
    if extrusion_interval is not None:
        native_extrusion = _finite_tuple(extrusion_interval, 2, "extrusion_interval")
        if native_extrusion[1] <= native_extrusion[0]:
            raise ValueError("extrusion_interval must be increasing")
        if isinstance(extrusion_volume_density, bool) or not isinstance(extrusion_volume_density, (int, float)) or \
                not math.isfinite(extrusion_volume_density) or extrusion_volume_density <= 0.0:
            raise ValueError("extrusion_volume_density must be finite and positive")
        if abs(extrusion_volume_density*(native_extrusion[1]-native_extrusion[0])-
               design_areal_density) > 1e-10*max(1.0, abs(design_areal_density)):
            raise ValueError(
                "design_areal_density must equal extrusion_volume_density times extrusion height")
    if not isinstance(fixed_extrusion_layers, Sequence) or isinstance(
            fixed_extrusion_layers, (str, bytes)):
        raise ValueError("fixed_extrusion_layers must be a sequence of layer dictionaries")
    native_layers = []
    for layer_index, layer in enumerate(fixed_extrusion_layers):
        if not isinstance(layer, dict) or set(layer) != {"loops", "z_interval", "volume_density"}:
            raise ValueError(
                f"fixed_extrusion_layers[{layer_index}] requires exactly loops, z_interval, and volume_density")
        layer_loops = normalize_loops(
            layer["loops"], f"fixed_extrusion_layers[{layer_index}].loops")
        layer_z = _finite_tuple(
            layer["z_interval"], 2, f"fixed_extrusion_layers[{layer_index}].z_interval")
        if layer_z[1] <= layer_z[0]:
            raise ValueError(f"fixed_extrusion_layers[{layer_index}].z_interval must be increasing")
        layer_density = layer["volume_density"]
        if isinstance(layer_density, bool) or not isinstance(layer_density, (int, float)) or \
                not math.isfinite(layer_density) or layer_density <= 0.0:
            raise ValueError(
                f"fixed_extrusion_layers[{layer_index}].volume_density must be finite and positive")
        native_layers.append({"loops": layer_loops, "z_interval": layer_z,
                              "volume_density": float(layer_density)})
    if fixed_material_loops:
        if fixed_areal_density is None:
            raise ValueError("fixed_areal_density is required with fixed_material_loops")
        fixed_loops = normalize_loops(fixed_material_loops, "fixed_material_loops")
    else:
        fixed_loops = ()
        if not native_layers:
            raise ValueError("fixed material must be supplied as loops or extrusion layers")
        if fixed_areal_density is not None:
            raise ValueError("omit fixed_areal_density when fixed layers define the fixed material")
    if native_layers and native_extrusion is None:
        raise ValueError("fixed_extrusion_layers require an extrusion_interval for the design material")
    native = _native_module()
    value = native.solve_planar_fidget_spinner(
        design_loops, fixed_loops,
        None if fixed_areal_density is None else float(fixed_areal_density),
        float(design_areal_density), np.asarray(center, dtype=np.float64),
        maximum_iterations, maximum_function_evaluations,
        float(residual_tolerance), float(initial_trust_radius),
        float(maximum_trust_radius), native_extrusion,
        None if extrusion_volume_density is None else float(extrusion_volume_density),
        native_layers,
    )
    result = dict(value)
    result["boundary_cycles"] = tuple(
        tuple(dict(piece) for piece in cycle)
        for cycle in result["boundary_cycles"]
    )
    for key in ("center", "fixed_material_moments", "initial_material_moments",
                "removed_material_moments", "remaining_material_moments",
                "normalized_first_moment_residual"):
        result[key] = np.asarray(result[key], dtype=np.float64).copy()
    if "extruded_design_material_raw_moments" in result:
        result["extruded_design_material_raw_moments"] = np.asarray(
            result["extruded_design_material_raw_moments"], dtype=np.float64
        ).copy()
    if "extruded_fixed_layer_raw_moments" in result:
        result["extruded_fixed_layer_raw_moments"] = tuple(
            np.asarray(layer, dtype=np.float64).copy()
            for layer in result["extruded_fixed_layer_raw_moments"]
        )
    if "fixed_extrusion_layers" in result:
        result["fixed_extrusion_layers"] = tuple(
            {
                "loops": tuple(tuple(tuple(float(value) for value in point)
                                      for point in loop)
                               for loop in layer["loops"]),
                "z_interval": np.asarray(layer["z_interval"], dtype=np.float64).copy(),
                "volume_density": float(layer["volume_density"]),
            }
            for layer in result["fixed_extrusion_layers"]
        )
    if "extruded_complete_material_raw_moments" in result:
        result["extruded_complete_material_raw_moments"] = np.asarray(
            result["extruded_complete_material_raw_moments"], dtype=np.float64
        ).copy()
    if "extruded_design_material_raw_moments" in result:
        entities = [{
            "role": "editable_design",
            "boundary_cycles": result["boundary_cycles"],
            "z_interval": np.asarray(result["extrusion_interval"], dtype=np.float64).copy(),
            "volume_density": float(result["extrusion_volume_density"]),
            "raw_mass_moments": result["extruded_design_material_raw_moments"].copy(),
            "geometry_representation": "exact_line_arc_profile_swept_between_planes",
        }]
        layers = result.get("fixed_extrusion_layers", ())
        layer_moments = result.get("extruded_fixed_layer_raw_moments", ())
        if len(layers) != len(layer_moments):
            raise RuntimeError("fixed extrusion entities and moment records do not match")
        for layer, moments in zip(layers, layer_moments):
            entities.append({
                "role": "fixed_material",
                "loops": layer["loops"],
                "z_interval": layer["z_interval"].copy(),
                "volume_density": layer["volume_density"],
                "raw_mass_moments": moments.copy(),
                "geometry_representation": "exact_polygonal_prism",
            })
        result["extruded_material_entities"] = tuple(entities)
        if "extruded_complete_material_raw_moments" not in result:
            result["extruded_complete_material_raw_moments"] = sum(
                (entity["raw_mass_moments"] for entity in entities), np.zeros(10)
            )
    return result


def export_spin_it_faster_planar_material_step(
    solver_result: dict[str, Any], *, occt_executable: str | os.PathLike[str],
    brep_path: str | os.PathLike[str], step_path: str | os.PathLike[str],
) -> dict[str, Any]:
    """Export one converged exact planar result through the pinned OCCT bridge.

    The bridge consumes the returned line/arc cycles and fixed polygonal
    prisms directly. It does not triangulate the circular boundary. Per-entity
    densities remain explicit in the returned sidecar evidence; STEP stores
    the distinct solids, while this function independently compares reloaded
    geometric integrals against the solver's density-weighted ten moments.
    """
    if not isinstance(solver_result, dict) or solver_result.get("status") != "CONVERGED":
        raise ValueError("CAD export requires a converged planar solver result")
    if solver_result.get("geometry_scope") != "exact_2d_boundary_and_analytic_layered_extrusion" or \
            not solver_result.get("boundary_defined"):
        raise ValueError("CAD export requires the exact solved planar boundary")
    entities = solver_result.get("extruded_material_entities")
    if not isinstance(entities, (tuple, list)) or not entities or \
            entities[0].get("role") != "editable_design":
        raise ValueError("CAD export requires the ordered editable-plus-fixed material entities")
    expected_complete = np.asarray(
        solver_result.get("extruded_complete_material_raw_moments", ()), dtype=np.float64)
    if expected_complete.shape != (10,) or not np.all(np.isfinite(expected_complete)):
        raise ValueError("CAD export requires ten finite solver mass moments")

    def number(value: Any, label: str) -> str:
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.number)):
            raise ValueError(f"{label} must be a finite real number")
        numeric = float(value)
        if not math.isfinite(numeric):
            raise ValueError(f"{label} must be a finite real number")
        return format(numeric, ".17g")

    lines = ["PLANAR_LAYERED_EXTRUSION_V2", str(len(entities))]
    roles: list[str] = []
    densities: list[float] = []
    entity_ids: list[str] = []
    source_entity_moments: list[np.ndarray] = []
    for entity_index, entity in enumerate(entities):
        if not isinstance(entity, dict):
            raise ValueError(f"material entity {entity_index} must be a mapping")
        role = entity.get("role")
        expected_role = "editable_design" if entity_index == 0 else "fixed_material"
        if role != expected_role:
            raise ValueError("material entities must contain one design region followed by fixed layers")
        density = entity.get("volume_density")
        density_token = number(density, f"material entity {entity_index} density")
        if float(density_token) <= 0.0:
            raise ValueError("material entity densities must be positive")
        z_interval = _finite_tuple(entity.get("z_interval", ()), 2,
                                   f"material entity {entity_index} z_interval")
        if z_interval[1] <= z_interval[0]:
            raise ValueError("material entity z intervals must be increasing")
        raw = np.asarray(entity.get("raw_mass_moments", ()), dtype=np.float64)
        if raw.shape != (10,) or not np.all(np.isfinite(raw)):
            raise ValueError(f"material entity {entity_index} lacks ten finite mass moments")
        source_entity_moments.append(raw)
        roles.append(role)
        densities.append(float(density_token))

        if role == "editable_design":
            cycles = entity.get("boundary_cycles")
            if not isinstance(cycles, (tuple, list)) or not cycles or \
                    entity.get("geometry_representation") != "exact_line_arc_profile_swept_between_planes":
                raise ValueError("editable material must be represented by exact line/arc boundary cycles")
            serialized_cycles = []
            for cycle_index, cycle in enumerate(cycles):
                serialized = []
                for piece_index, piece in enumerate(cycle):
                    if piece.get("type") == "LINE_SEGMENT":
                        start = _finite_tuple(piece.get("start", ()), 2, "line start")
                        end = _finite_tuple(piece.get("end", ()), 2, "line end")
                        serialized.append("L " + " ".join((
                            number(start[0], "line start x"), number(start[1], "line start y"),
                            number(end[0], "line end x"), number(end[1], "line end y"))))
                    elif piece.get("type") == "CIRCULAR_ARC":
                        center = _finite_tuple(piece.get("center", ()), 2, "arc center")
                        radius, start_angle, sweep = (piece.get("radius"),
                                                      piece.get("start_angle"), piece.get("sweep"))
                        serialized.append("A " + " ".join((
                            number(center[0], "arc center x"), number(center[1], "arc center y"),
                            number(radius, "arc radius"), number(start_angle, "arc start angle"),
                            number(sweep, "arc sweep"))))
                    else:
                        raise ValueError(
                            f"unsupported exact boundary piece at cycle {cycle_index}, piece {piece_index}")
                if not serialized:
                    raise ValueError("editable material contains an empty boundary cycle")
                serialized_cycles.append(serialized)
        else:
            if entity.get("geometry_representation") != "exact_polygonal_prism":
                raise ValueError("fixed material must be represented by exact polygonal prisms")
            loops = entity.get("loops")
            if not isinstance(loops, (tuple, list)) or not loops:
                raise ValueError("fixed material prism requires at least one boundary loop")
            serialized_cycles = []
            for loop in loops:
                points = tuple(_finite_tuple(point, 2, "fixed polygon point") for point in loop)
                if len(points) < 3:
                    raise ValueError("fixed polygon loops require at least three points")
                serialized_cycles.append([
                    "L " + " ".join((number(points[index][0], "fixed x0"),
                                      number(points[index][1], "fixed y0"),
                                      number(points[(index + 1) % len(points)][0], "fixed x1"),
                                      number(points[(index + 1) % len(points)][1], "fixed y1")))
                    for index in range(len(points))
                ])
        z_min_token = number(z_interval[0], "z minimum")
        z_max_token = number(z_interval[1], "z maximum")
        identity_payload = ["ENTITY_GEOMETRY_DENSITY_V1", role, density_token,
                            z_min_token, z_max_token, str(len(serialized_cycles))]
        for cycle in serialized_cycles:
            identity_payload.append(str(len(cycle)))
            identity_payload.extend(cycle)
        entity_id = hashlib.sha256(
            ("\n".join(identity_payload) + "\n").encode("ascii")).hexdigest()
        entity_ids.append(entity_id)
        lines.append("ENTITY " + " ".join((role, entity_id, density_token,
                    z_min_token, z_max_token, str(len(serialized_cycles)))))
        for cycle in serialized_cycles:
            lines.append(str(len(cycle)))
            lines.extend(cycle)
    serialized_input = "\n".join(lines) + "\n"
    source_digest = hashlib.sha256(serialized_input.encode("ascii")).hexdigest()

    executable = Path(occt_executable).resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError("occt_executable must be an executable pinned OCCT bridge")
    brep = Path(brep_path).absolute()
    step = Path(step_path).absolute()
    identity_manifest = step.with_suffix(step.suffix + ".identity.json")
    if brep == step:
        raise ValueError("B-rep and STEP output paths must differ")
    if brep.exists() or step.exists() or identity_manifest.exists():
        raise FileExistsError("CAD output paths must be new; refusing to overwrite existing evidence")
    if not brep.parent.is_dir() or not step.parent.is_dir():
        raise ValueError("CAD output parent directories must already exist")
    with tempfile.TemporaryDirectory(prefix="spin-planar-cad-") as temporary:
        input_path = Path(temporary) / "planar_material.txt"
        input_path.write_text(serialized_input, encoding="ascii")
        completed = subprocess.run(
            [str(executable), "--planar-layered-extrusion", str(input_path), str(brep), str(step)],
            check=False, capture_output=True, text=True, timeout=60)
    if completed.returncode != 0:
        raise RuntimeError("OCCT planar material export failed: " + completed.stderr.strip())
    try:
        cad_result = json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise RuntimeError("OCCT planar material export returned invalid evidence JSON") from error
    if cad_result.get("status") != "PLANAR_LAYERED_STEP_BREP_ROUNDTRIP_VALIDATED":
        raise RuntimeError("OCCT planar material export did not report a validated round trip")
    if cad_result.get("material_entity_ids") != entity_ids or \
            cad_result.get("entity_roles") != roles or \
            not np.array_equal(np.asarray(cad_result.get("entity_densities", ())), densities):
        raise RuntimeError("OCCT output material identities/roles/densities do not match the solver result")
    cad_entity_moments = np.asarray(cad_result.get("entity_mass_moments", ()), dtype=np.float64)
    source_entity_array = np.asarray(source_entity_moments, dtype=np.float64)
    if cad_entity_moments.shape != source_entity_array.shape or \
            not np.all(np.isfinite(cad_entity_moments)):
        raise RuntimeError("OCCT did not return a finite mass-moment record for every entity")
    entity_scale = np.maximum(1.0, np.max(np.abs(source_entity_array), axis=1))
    entity_error = float(np.max(np.abs(cad_entity_moments-source_entity_array) /
                                entity_scale[:, None]))
    if entity_error > 1e-6:
        raise RuntimeError(f"OCCT entity mass moments differ from solver result by {entity_error:.6g}")
    reloaded = np.asarray(cad_result.get("raw_mass_moments", {}).get("step", ()), dtype=np.float64)
    if reloaded.shape != (10,) or not np.all(np.isfinite(reloaded)):
        raise RuntimeError("OCCT STEP reload lacks ten finite density-weighted moments")
    normalized_error = float(np.max(np.abs(reloaded-expected_complete)) /
                             max(1.0, float(np.max(np.abs(expected_complete)))))
    if normalized_error > 1e-6:
        raise RuntimeError(f"OCCT STEP mass moments differ from solver result by {normalized_error:.6g}")
    if not brep.is_file() or not step.is_file() or brep.stat().st_size == 0 or step.stat().st_size == 0:
        raise RuntimeError("OCCT reported success without producing nonempty B-rep and STEP files")
    cad_result["solver_geometry_moment_error"] = normalized_error
    cad_result["entity_geometry_moment_error"] = entity_error
    cad_result["serialized_geometry_sha256"] = source_digest
    solid_counts = cad_result.get("entity_solid_counts")
    if not isinstance(solid_counts, list) or len(solid_counts) != len(entities) or any(
            type(count) is not int or count <= 0 for count in solid_counts):
        raise RuntimeError("OCCT output lacks a positive solid count for each material identity")
    cad_result["material_identity_sidecar"] = tuple(
        {"entity_id": entity_id, "role": role, "volume_density": density,
         "solid_count": solid_count,
         "raw_mass_moments": tuple(float(v) for v in entity_moments)}
        for entity_id, role, density, solid_count, entity_moments in zip(
            entity_ids, roles, densities, solid_counts, cad_entity_moments))
    cad_result["step_material_identity_scope"] = "geometry-density-bound-ordered-sidecar-v1"
    def file_record(path: Path) -> dict[str, Any]:
        digest = hashlib.sha256()
        with path.open("rb") as source:
            for block in iter(lambda: source.read(1024 * 1024), b""):
                digest.update(block)
        return {"path": str(path), "bytes": path.stat().st_size,
                "sha256": digest.hexdigest()}

    manifest = {
        "schema": "PLANAR_STEP_MATERIAL_IDENTITY_V1",
        "scope": cad_result["step_material_identity_scope"],
        "serialized_geometry_sha256": source_digest,
        "files": {"brep": file_record(brep), "step": file_record(step)},
        "entities": [
            {"entity_id": entity_id, "role": role,
             "volume_density": density, "solid_count": solid_count,
             "density_weighted_raw_moments": [float(v) for v in moments]}
            for entity_id, role, density, solid_count, moments in zip(
                entity_ids, roles, densities, solid_counts, cad_entity_moments)
        ],
        "solver_density_weighted_raw_moments": [float(v) for v in expected_complete],
        "step_reloaded_density_weighted_raw_moments": [float(v) for v in reloaded],
    }
    manifest_bytes = (json.dumps(manifest, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")
    descriptor, temporary_manifest = tempfile.mkstemp(
        prefix=identity_manifest.name + ".", suffix=".tmp", dir=str(identity_manifest.parent))
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(manifest_bytes)
            output.flush()
            os.fsync(output.fileno())
        os.link(temporary_manifest, identity_manifest)
    except BaseException:
        try:
            os.unlink(temporary_manifest)
        except FileNotFoundError:
            pass
        raise
    else:
        os.unlink(temporary_manifest)
    cad_result["material_identity_manifest_path"] = str(identity_manifest)
    cad_result["material_identity_manifest_sha256"] = hashlib.sha256(manifest_bytes).hexdigest()
    return cad_result


def evaluate_spin_it_faster_quadric(
    request: SpinItFasterRequest,
    coefficients: Sequence[float],
    *, evaluator_settings: SpinEvaluatorConfig | None = None,
    coefficients_are_normalized: bool = False,
) -> SpinQuadricEvaluation:
    """Evaluate raw coefficients, or explicitly replay a validated saved unit vector.

    The default retains historical raw-input normalization. Unit-vector replay
    preserves the represented endpoint and does not supply optimizer initial guesses.
    """
    if type(coefficients_are_normalized) is not bool:
        raise ValueError("coefficients_are_normalized must be a boolean")
    values = _finite_tuple(coefficients, 10, "quadric coefficients")
    native, native_request = _native_request(request)
    extra = () if evaluator_settings is None else (_native_evaluator_settings(native, evaluator_settings),)
    if coefficients_are_normalized:
        extra = (None if evaluator_settings is None else extra[0], True)
    evaluated = native.evaluate_spin_it_faster_quadric(
        native_request,
        np.asarray(values, dtype=np.float64),
        *extra,
    )
    result_type, authority = _rational_result_type(SpinQuadricEvaluation,
        _RationalSpinQuadricEvaluation, evaluated,
        required=bool(evaluator_settings and evaluator_settings.rational_source_coordinates))
    return result_type(
        **authority,
        variable_raw_moments=tuple(float(value) for value in evaluated["variable_raw_moments"]),
        total_raw_moments=tuple(float(value) for value in evaluated["total_raw_moments"]),
        geometry_hash=str(evaluated["geometry_hash"]),
        topology_hash=str(evaluated["topology_hash"]),
        interface_type=str(evaluated["interface_type"]),
        interface_subtype=str(evaluated["interface_subtype"]),
        **({} if evaluator_settings is None else {
            **{name: evaluated[name] for name in ("base_request_hash", "request_hash", "derivative_method",
                "derivative_defined", "derivative_diagnostic", "stratum_derivative_defined",
                "stratum_derivative_diagnostic", "derivative_scientific_hash", "provenance_topology_hash",
                "boolean_backend_identity")},
            "stratum_moment_jacobian": tuple(tuple(float(value) for value in row)
                for row in evaluated["stratum_moment_jacobian"].tolist()),
            "stratum_coefficient_basis": tuple(tuple(float(value) for value in row)
                for row in evaluated["stratum_coefficient_basis"].tolist()),
            "quadric_solid_is_positive_region": bool(evaluated["quadric_solid_is_positive_region"]),
            "vertex_provenance": tuple(evaluated["vertex_provenance"]),
            **{name: tuple(tuple(row) for row in evaluated[name].tolist()) for name in
               ("interface_vertices", "interface_faces", "material_vertices", "material_faces",
                "outer_vertices", "outer_faces")},
        }),
    )


def evaluate_spin_it_faster_quadric_diagnostic(
    request: SpinItFasterRequest,
    coefficients: Sequence[float],
) -> SpinQuadricTeacherDiagnostic:
    """Export moments and the default evaluator's same-call ambient Jacobian.

    This is a separate, opt-in diagnostic endpoint. It intentionally does not
    accept evaluator overrides: the data contract is the same default legacy
    evaluator used by the ordinary Spin-It Faster solver. A derivative marked
    undefined is returned as ``None`` and must not be supervised as zero.
    """
    values = _finite_tuple(coefficients, 10, "quadric coefficients")
    native, native_request = _native_request(request)
    if not hasattr(native, "evaluate_spin_it_faster_quadric_diagnostic"):
        raise RuntimeError("native same-call quadric diagnostic export is unavailable")
    evaluated = native.evaluate_spin_it_faster_quadric_diagnostic(
        native_request,
        np.asarray(values, dtype=np.float64),
    )
    return _quadric_teacher_diagnostic_from_native(evaluated, request, values)


def _quadric_teacher_diagnostic_from_native(
    evaluated: dict[str, Any],
    request: SpinItFasterRequest,
    values: tuple[float, ...],
) -> SpinQuadricTeacherDiagnostic:
    variable = tuple(float(value) for value in evaluated["variable_raw_moments"])
    total = tuple(float(value) for value in evaluated["total_raw_moments"])
    if len(variable) != 10 or len(total) != 10 or not np.isfinite(variable + total).all():
        raise RuntimeError("native diagnostic returned invalid ten-component moments")
    jacobian_value = evaluated["moment_jacobian"]
    jacobian = None if jacobian_value is None else tuple(
        tuple(float(value) for value in row) for row in jacobian_value.tolist()
    )
    derivative_defined = bool(evaluated["derivative_defined"])
    if derivative_defined:
        if jacobian is None or len(jacobian) != 10 or any(len(row) != 10 for row in jacobian):
            raise RuntimeError("native diagnostic marked a non-10x10 Jacobian as defined")
        if not np.isfinite(np.asarray(jacobian, dtype=np.float64)).all():
            raise RuntimeError("native diagnostic returned a non-finite defined Jacobian")
    elif jacobian is not None:
        raise RuntimeError("native diagnostic must not serialize an undefined Jacobian")
    residual_seconds = float(evaluated["residual_evaluation_seconds"])
    derivative_seconds = float(evaluated["derivative_evaluation_seconds"])
    if (not math.isfinite(residual_seconds) or residual_seconds < 0.0 or
            not math.isfinite(derivative_seconds) or derivative_seconds < 0.0):
        raise RuntimeError("native diagnostic returned invalid timing values")
    return SpinQuadricTeacherDiagnostic(
        coefficients=values,
        variable_raw_moments=variable,
        total_raw_moments=total,
        fixed_raw_moments=tuple(total[i] - variable[i] for i in range(10)),
        moment_jacobian=jacobian,
        derivative_defined=derivative_defined,
        derivative_method=str(evaluated["derivative_method"]),
        derivative_diagnostic=str(evaluated["derivative_diagnostic"]),
        equation_profile=request.settings.equation_profile,
        geometry_hash=str(evaluated["geometry_hash"]),
        topology_hash=str(evaluated["topology_hash"]),
        base_request_hash=str(evaluated["base_request_hash"]),
        evaluator_request_hash=str(evaluated["evaluator_request_hash"]),
        derivative_scientific_hash=str(evaluated["derivative_scientific_hash"]),
        provenance_topology_hash=str(evaluated["provenance_topology_hash"]),
        residual_evaluation_seconds=residual_seconds,
        derivative_evaluation_seconds=derivative_seconds,
    )


def evaluate_spin_it_faster_fullrank_quadric(
    request: SpinItFasterRequest,
    coefficients: Sequence[float],
    *, evaluator_settings: SpinEvaluatorConfig | None = None,
    diagnostic_stage_observer: Callable[[str], None] | None = None,
    capture_material_proof: bool = False,
) -> SpinFullRankQuadricEvaluation:
    """Evaluate the opt-in source-motion/actual-polyhedron derivative path.

    This does not replace the legacy public solver path; unsupported strata are
    represented by ``derivative_defined=False`` rather than a zero Jacobian.
    ``evaluator_settings`` controls the triangulation used to build the
    discrete surface; the canonical path selects its discrete derivative
    internally, so legacy derivative-method and interface-gradient knobs do
    not change this evaluator's derivative semantics.
    The optional observer receives diagnostic-only production stage names;
    it is not hashed or included in the scientific result. Sink exceptions
    abort the call rather than granting a result from incomplete recording.
    Opt-in proof capture carries the existing exact-material Boolean tape;
    it does not perform another evaluation or grant material/branch admission.
    """
    if diagnostic_stage_observer is not None and not callable(diagnostic_stage_observer):
        raise ValueError('diagnostic stage observer must be callable or None')
    if not isinstance(capture_material_proof, bool):
        raise ValueError('capture_material_proof must be boolean')
    if capture_material_proof and not (evaluator_settings and
            evaluator_settings.exact_binary64_source_boundary and
            evaluator_settings.exact_binary64_material_boundary):
        raise ValueError('proof capture requires explicit exact source and material boundary policies')
    native, native_request = _native_request(request)
    values = _finite_tuple(coefficients, 10, "quadric coefficients")
    extra = () if evaluator_settings is None else (_native_evaluator_settings(native, evaluator_settings),)
    if diagnostic_stage_observer is not None:
        settings = extra[0] if extra else native.SpinEvaluatorSettings()
        settings.diagnostic_stage_observer = diagnostic_stage_observer
        extra = (settings,)
    evaluated = native.evaluate_spin_it_faster_fullrank_quadric(
        native_request,
        np.asarray(values, dtype=np.float64),
        *extra,
        **({'capture_material_proof': True} if capture_material_proof else {}),
    )
    jacobian = evaluated["coefficient_jacobian"]
    result_type, authority = _rational_result_type(SpinFullRankQuadricEvaluation,
        _RationalSpinFullRankEvaluation, evaluated,
        required=bool(evaluator_settings and evaluator_settings.rational_source_coordinates))
    if result_type is _RationalSpinFullRankEvaluation:
        authority['normalized_coefficients'] = _finite_tuple(
            evaluated['normalized_coefficients'], 10, 'actual native normalized coefficients')
        branch = str(evaluated['branch_identity'])
        if len(branch)!=64 or any(x not in '0123456789abcdef' for x in branch):
            raise RuntimeError('native rational production branch identity missing')
        authority['branch_identity'] = branch
    if capture_material_proof:
        authority['material_proof_json'] = _same_call_material_proof_json(evaluated, len(request.design_faces))
        result_type = (_ProofRationalSpinFullRankEvaluation
            if result_type is _RationalSpinFullRankEvaluation else _ProofSpinFullRankEvaluation)
    return result_type(
        **authority,
        variable_raw_moments=tuple(float(value) for value in evaluated["variable_raw_moments"]),
        total_raw_moments=tuple(float(value) for value in evaluated["total_raw_moments"]),
        geometry_hash=str(evaluated["geometry_hash"]),
        boolean_backend_identity=str(evaluated["boolean_backend_identity"]),
        source_scientific_hash=str(evaluated["source_scientific_hash"]),
        scientific_hash=str(evaluated["scientific_hash"]),
        derivative_defined=bool(evaluated["derivative_defined"]),
        derivative_diagnostic=str(evaluated["derivative_diagnostic"]),
        coefficient_jacobian=None if jacobian is None else tuple(
            tuple(float(value) for value in row) for row in jacobian.tolist()
        ),
        material_vertices=tuple(
            tuple(float(value) for value in row)
            for row in evaluated["material_vertices"].tolist()
        ),
        material_faces=tuple(
            tuple(int(value) for value in row)
            for row in evaluated["material_faces"].tolist()
        ),
    )


def solve_spin_it_faster_fullrank_quadric(
    request: SpinItFasterRequest,
    *,
    initial_quadrics: Sequence[Sequence[float]] = (),
    evaluator_settings: SpinEvaluatorConfig | None = None,
) -> SpinFullRankQuadricSolveResult:
    """Run the opt-in canonical solver with the discrete Boolean-moment Jacobian.

    This does not replace ``solve_spin_it_faster``. The canonical solver
    requires the separately identified dual-monomial/raw equation profile and
    exposes a final geometry only after its independent moment and KKT audits
    pass.
    """
    if request.settings.equation_profile != "DUAL_MONOMIAL_RAW_V1":
        raise ValueError(
            "the canonical full-rank solver requires equation_profile="
            "'DUAL_MONOMIAL_RAW_V1'"
        )
    native, native_request = _native_request(request)
    initial = np.asarray(
        tuple(_finite_tuple(row, 10, "initial quadric") for row in initial_quadrics),
        dtype=np.float64,
    ).reshape((-1, 10))
    extra = () if evaluator_settings is None else (
        _native_evaluator_settings(native, evaluator_settings),
    )
    value = native.solve_spin_it_faster_fullrank_quadric(
        native_request, initial, *extra
    )

    return _fullrank_solve_result_from_native(value,
        rational_source_coordinates=bool(evaluator_settings and evaluator_settings.rational_source_coordinates))


def _fullrank_solve_result_from_native(value: dict[str, Any], *,
    rational_source_coordinates: bool = False) -> SpinFullRankQuadricSolveResult:
    trace = tuple(
        SpinIteration(
            restart=int(item["restart"]),
            iteration=int(item["iteration"]),
            function_evaluations=int(item["function_evaluations"]),
            normalized_equation_residual=float(item["normalized_equation_residual"]),
            physical_relative_residual=float(item["physical_relative_residual"]),
            objective=float(item["objective"]),
            multiplier_mu=float(item["multiplier_mu"]),
            accepted_geodesic_step=float(item["accepted_geodesic_step"]),
            quadric_hash=str(item["quadric_hash"]),
            material_geometry_hash=str(item["material_geometry_hash"]),
            material_topology_hash=str(item["material_topology_hash"]),
            equation_profile=str(value["equation_profile"]),
            solver_coordinate_identity="literal_monomial_unit_sphere_v1",
            raw_equation_residual=(
                None if item["raw_equation_residual"] is None
                else float(item["raw_equation_residual"])
            ),
            quadric_coefficients=(
                None if item["quadric_coefficients"] is None
                else tuple(float(x) for x in item["quadric_coefficients"])
            ),
        )
        for item in value["trace"]
    )
    events = tuple(_immutable_spin_solver_event(item) for item in value["events"])
    evaluation_timings = tuple(
        SpinFullRankEvaluationTiming(
            evaluation_index=int(item["evaluation_index"]),
            restart=int(item["restart"]),
            restart_evaluation_index=int(item["restart_evaluation_index"]),
            evaluator_returned=bool(item["evaluator_returned"]),
            wall_seconds=float(item["wall_seconds"]),
            reported_residual_seconds=(
                None if item["reported_residual_seconds"] is None
                else float(item["reported_residual_seconds"])
            ),
            reported_derivative_seconds=(
                None if item["reported_derivative_seconds"] is None
                else float(item["reported_derivative_seconds"])
            ),
        )
        for item in value["evaluation_timings"]
    )

    def optional_vector(key: str) -> tuple[float, ...] | None:
        values = value[key]
        return None if values is None else tuple(float(x) for x in values)

    material_vertices = value["material_vertices"]
    material_faces = value["material_faces"]
    result_type, authority = _rational_result_type(SpinFullRankQuadricSolveResult,
        _RationalSpinFullRankSolveResult, value,
        required=rational_source_coordinates and bool(value['final_geometry_verified']))
    return result_type(
        **authority,
        status=str(value["status"]),
        diagnostic=str(value["diagnostic"]),
        final_geometry_verified=bool(value["final_geometry_verified"]),
        verification_diagnostic=str(value["verification_diagnostic"]),
        input_hash=str(value["input_hash"]),
        output_hash=str(value["output_hash"]),
        function_evaluations=int(value["function_evaluations"]),
        restart_count=int(value["restart_count"]),
        normalized_equation_residual=float(value["normalized_equation_residual"]),
        physical_relative_residual=float(value["physical_relative_residual"]),
        independent_raw_residual=(
            None if value["independent_raw_residual"] is None
            else float(value["independent_raw_residual"])
        ),
        multiplier_mu=float(value["multiplier_mu"]),
        multipliers=tuple(float(x) for x in value["multipliers"]),
        quadric_coefficients=(
            None if value["quadric_coefficients"] is None
            else tuple(float(x) for x in value["quadric_coefficients"])
        ),
        equation_profile=str(value["equation_profile"]),
        event_hash=str(value["event_hash"]),
        trace_hash=str(value["trace_hash"]),
        trace=trace,
        events=events,
        restart_function_evaluations=tuple(
            int(x) for x in value["restart_function_evaluations"]
        ),
        evaluation_timings=evaluation_timings,
        final_geometry_hash=(
            None if value["final_geometry_hash"] is None
            else str(value["final_geometry_hash"])
        ),
        final_derivative_defined=bool(value["final_derivative_defined"]),
        final_derivative_diagnostic=str(value["final_derivative_diagnostic"]),
        final_variable_raw_moments=optional_vector("final_variable_raw_moments"),
        final_total_raw_moments=optional_vector("final_total_raw_moments"),
        independent_total_raw_moments=optional_vector("independent_total_raw_moments"),
        material_vertices=(
            None if material_vertices is None else tuple(
                tuple(float(x) for x in row) for row in material_vertices
            )
        ),
        material_faces=(
            None if material_faces is None else tuple(
                tuple(int(x) for x in row) for row in material_faces
            )
        ),
    )


@dataclass(frozen=True)
class Spin2014QuadricMaterialRun:
    """Attachment to the ordinary solver result, not a manufacturing seal."""
    objective_profile: str
    objective_identity: str
    objective_settings: Spin2014QuadricObjectiveConfig
    measured_request: SpinItFasterRequest
    fixed_material_density: float
    solve: SpinFullRankQuadricSolveResult
    material_export_complete: bool
    material_export_diagnostic: str
    cavity_vertices: tuple[tuple[float, ...], ...] | None
    cavity_faces: tuple[tuple[int, ...], ...] | None
    cavity_raw_moments: tuple[float, ...] | None


@dataclass(frozen=True)
class Spin2014StaticQuadricMaterialRun:
    """Plane candidate and freshly evaluated square objective, not a seal."""
    objective_profile: str
    objective_identity: str
    objective_settings: Spin2014QuadricObjectiveConfig
    measured_request: SpinItFasterRequest
    fixed_material_density: float
    target_length: float
    initialization: str
    initial_coefficients: tuple[float, ...]
    candidate: dict[str, Any]
    squared_height: dict[str, Any]
    material_verified: bool = False


def solve_spin2014_static_quadric_hollowing(
    request: SpinItFasterRequest, objective_settings: Spin2014QuadricObjectiveConfig,
    *, fixed_material_density: float,
    initial_coefficients: Sequence[float] | None = None,
    maximum_iterations: int = 60, maximum_line_search_steps: int = 24,
    constraint_tolerance: float = 1e-10,
) -> Spin2014StaticQuadricMaterialRun:
    """Obtain an h-minimizing candidate; independently evaluate 2014 gamma*h².

    Default initial plane comes only from the actual input full-material height.
    Explicit planes are diagnostic inputs, not frozen-batch initialization.
    Returned local conditions do not establish support or material acceptance.
    """
    if not isinstance(request, SpinItFasterRequest):
        raise TypeError("request must be a SpinItFasterRequest")
    if request.mode != "STATIC_BALANCE" or objective_settings.objective != "STATIC_LOWERING":
        raise ValueError("static plane adapter requires STATIC_BALANCE and STATIC_LOWERING")
    if type(maximum_iterations) is not int or maximum_iterations < 0:
        raise ValueError("maximum_iterations must be a nonnegative integer")
    if type(maximum_line_search_steps) is not int or maximum_line_search_steps <= 0:
        raise ValueError("maximum_line_search_steps must be a positive integer")
    if isinstance(constraint_tolerance, bool) or not isinstance(constraint_tolerance, (int,float)) or \
            not math.isfinite(constraint_tolerance) or not 0 < constraint_tolerance <= 1e-10:
        raise ValueError("static constraint tolerance must be finite in (0,1e-10]")
    native, native_request = _native_request(request)
    settings = native.Spin2014QuadricObjectiveSettings()
    for name,value in asdict(objective_settings).items():
        setattr(settings,name,getattr(native.Spin2014Objective,value) if name=="objective" else value)
    initial = None if initial_coefficients is None else np.asarray(
        _finite_tuple(initial_coefficients,4,"initial_coefficients"),dtype=np.float64)
    value = native.solve_spin2014_static_quadric_hollowing(native_request,settings,
        fixed_material_density,initial,maximum_iterations,maximum_line_search_steps,
        float(constraint_tolerance))
    # Preserve the ordinary candidate's trace, flags and geometry separately.
    # Only the new squared_height attachment describes the new objective gate.
    return Spin2014StaticQuadricMaterialRun(
        objective_profile=str(value["objective_profile"]),
        objective_identity=str(value["objective_identity"]),objective_settings=objective_settings,
        measured_request=replace(request,fixed_shell_moments=tuple(float(x) for x in value["fixed_shell_moments"])),
        fixed_material_density=float(value["fixed_material_density"]),
        target_length=float(value["target_length"]),initialization=str(value["initialization"]),
        initial_coefficients=tuple(float(x) for x in value["initial_coefficients"]),
        candidate=dict(value["candidate"]),squared_height=dict(value["squared_height"]))


def solve_spin2014_quadric_hollowing(
    request: SpinItFasterRequest,
    objective_settings: Spin2014QuadricObjectiveConfig,
    *, fixed_material_density: float,
    evaluator_settings: SpinEvaluatorConfig | None = None,
    initial_quadrics: Sequence[Sequence[float]] = (),
    event_sink: Any = None,
    coefficient_policy: str | None = None,
    authority_evaluation_budget_per_attempt: int = 60,
) -> Spin2014QuadricMaterialRun | dict[str, Any]:
    """2014 objectives through actual CGAL discrete materials and strict Newton.

    Fixed moments are always derived from geometry by native admission. Nonzero
    supplied moments are consistency checks only. Explicit starting quadrics
    are for diagnostics; functional batches must use request/seed-only guesses.
    Availability of both terminal exports is not exact partition or float mesh
    acceptance. Static squared-height uses its separate plane entry.

    Explicit MPFR1024 returns a versioned authoritative numerical candidate,
    not the legacy float material-run type or a final material seal. Authority
    observers receive each event before continuing, including before initial
    material construction. Omitted policy preserves the historical entry.
    """
    if coefficient_policy not in (None, "COEFFICIENT_MPFR1024_V1"):
        raise ValueError("unsupported coefficient_policy")
    if coefficient_policy is None and authority_evaluation_budget_per_attempt != 60:
        raise ValueError("authority budget requires explicit coefficient_policy")
    if event_sink is not None and not callable(event_sink):
        raise TypeError("event_sink must be callable or None")
    if coefficient_policy == "COEFFICIENT_MPFR1024_V1":
        if len(initial_quadrics):
            raise ValueError("authority production initialization is seed-only")
    native, native_request = _native_request(request)
    settings = native.Spin2014QuadricObjectiveSettings()
    for name, value in asdict(objective_settings).items():
        setattr(settings, name, getattr(native.Spin2014Objective, value)
                if name == "objective" else value)
    evaluator = evaluator_settings or SpinEvaluatorConfig(
        derivative_method="DISCRETE_FULL_RANK_QUADRIC", icosphere_subdivisions=2)
    if coefficient_policy == "COEFFICIENT_MPFR1024_V1":
        if evaluator_settings is None:
            evaluator = replace(evaluator, rational_source_coordinates=True,
                compact_rational_source_coordinates=True,
                exact_binary64_source_boundary=True,
                exact_binary64_material_boundary=True)
        value = native.solve_spin2014_authority_quadric_hollowing(
            native_request, settings, fixed_material_density,
            _native_evaluator_settings(native, evaluator),
            evaluation_budget_per_attempt=authority_evaluation_budget_per_attempt,
            absolute_linear_residual_limit="1/1000000000000", event_sink=event_sink)
        return {**dict(value), "coefficient_policy": coefficient_policy,
                "precision_bits": 1024, "objective_profile": "SPIN2014_OBJECTIVES_QUADRIC_V1",
                "objective_settings": asdict(objective_settings),
                "fixed_material_density": float(fixed_material_density)}
    initial = np.asarray(tuple(_finite_tuple(row, 10, "initial quadric")
                               for row in initial_quadrics), dtype=np.float64).reshape((-1, 10))
    value = native.solve_spin2014_quadric_hollowing(
        native_request, settings, fixed_material_density,
        _native_evaluator_settings(native, evaluator), initial, event_sink)

    def matrix(key, scalar):
        x = value[key]
        return None if x is None else tuple(tuple(scalar(v) for v in row) for row in x)

    moments = value["cavity_raw_moments"]
    return Spin2014QuadricMaterialRun(
        objective_profile=str(value["objective_profile"]),
        objective_identity=str(value["objective_identity"]),
        objective_settings=objective_settings,
        measured_request=replace(request, fixed_shell_moments=tuple(
            float(x) for x in value["fixed_shell_moments"])),
        fixed_material_density=float(value["fixed_material_density"]),
        solve=_fullrank_solve_result_from_native(value["solve"],
            rational_source_coordinates=evaluator.rational_source_coordinates),
        material_export_complete=bool(value["material_export_complete"]),
        material_export_diagnostic=str(value["material_export_diagnostic"]),
        cavity_vertices=matrix("cavity_vertices", float),
        cavity_faces=matrix("cavity_faces", int),
        cavity_raw_moments=None if moments is None else tuple(float(x) for x in moments))


def solve_spin_it_faster_with_external_moment_evaluator(
    request: SpinItFasterRequest,
    evaluator: Any,
    *,
    initial_quadrics: Sequence[Sequence[float]] = (),
    event_sink: Any = None,
    support_basis: Any = None,
    objective_settings: Spin2014QuadricObjectiveConfig | None = None,
) -> dict[str, Any]:
    """Run strict Newton against a caller-owned actual geometry/moment map.

    ``evaluator(q)`` must return target-frame raw moments and their matching
    native-coefficient Jacobian from one concrete geometry evaluation. The
    result is numerical Newton evidence only; it never seals or independently
    validates the final material geometry.

    With ``support_basis`` (orthonormal dual coefficient columns), callbacks
    instead provide ``native_tangent_basis`` and
    ``native_directional_jacobian``. Explicit initial quadrics are required;
    this entry never projects full-space guesses or invents ambient derivatives.
    A versioned ``supported_space`` attachment retains all directional calls.
    """
    if objective_settings is not None:
        if type(objective_settings) is not Spin2014QuadricObjectiveConfig:
            raise TypeError('expected Spin2014QuadricObjectiveConfig')
        mode = {'YOYO': 'YO_YO', 'STATIC_LOWERING': 'STATIC_BALANCE'}.get(
            objective_settings.objective, 'SPINNING_TOP')
        if request.mode != mode:
            raise ValueError('original objective and request mode disagree')
        if mode == 'STATIC_BALANCE' or support_basis is not None:
            raise ValueError('original objective requires rotational full-space entry; static uses its plane solver')
    if not callable(evaluator):
        raise TypeError("evaluator must be callable")
    if event_sink is not None and not callable(event_sink):
        raise TypeError("event_sink must be callable or None")
    if request.settings.equation_profile != "DUAL_MONOMIAL_RAW_V1":
        raise ValueError(
            "external strict evaluator requires equation_profile="
            "'DUAL_MONOMIAL_RAW_V1'"
        )
    native, native_request = _native_request(request)
    initial = np.asarray(
        tuple(_finite_tuple(row, 10, "initial quadric") for row in initial_quadrics),
        dtype=np.float64,
    ).reshape((-1, 10))
    args = (native_request, evaluator, initial, event_sink)
    if objective_settings is not None:
        if not hasattr(native, 'Spin2014QuadricObjectiveSettings'):
            raise RuntimeError('native build lacks Spin-It 2014 objective support')
        settings = native.Spin2014QuadricObjectiveSettings()
        for name, value in asdict(objective_settings).items():
            setattr(settings, name, getattr(native.Spin2014Objective, value)
                    if name == 'objective' else value)
        return dict(native.solve_spin_it_faster_external_moment_evaluator(
            *args, objective_settings=settings))
    if support_basis is None:
        return dict(native.solve_spin_it_faster_external_moment_evaluator(*args))
    return dict(native.solve_spin_it_faster_external_moment_evaluator(
        *args, np.asarray(support_basis, dtype=np.float64)
    ))


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


def compute_hydrostatic_geometry(
    vertices: Sequence[Sequence[float]],
    faces: Sequence[Sequence[int]],
    *,
    waterline_origin: Sequence[float],
    waterline_normal: Sequence[float],
) -> HydrostaticGeometry:
    """Measure the actual negative-side waterline clipping of a closed mesh.

    The waterline has the equation ``dot(x-origin, normal) == 0`` and the
    displaced region is its negative side.  Coplanar source faces and a
    non-intersecting waterline are rejected rather than given an ambiguous
    numerical value.
    """
    native = _native_module()
    value = native.compute_hydrostatic_geometry(
        np.asarray(_vertices(vertices, "vertices"), dtype=np.float64),
        np.asarray(_faces(faces, "faces"), dtype=np.int32),
        np.asarray(_finite_tuple(waterline_origin, 3, "waterline_origin"), dtype=np.float64),
        np.asarray(_finite_tuple(waterline_normal, 3, "waterline_normal"), dtype=np.float64),
    )
    mesh = value.submerged_boundary
    return HydrostaticGeometry(
        displaced_volume=float(value.displaced_volume),
        center_of_buoyancy=tuple(float(x) for x in value.center_of_buoyancy),
        waterline_area=float(value.waterline_area),
        waterline_centroid=tuple(float(x) for x in value.waterline_centroid),
        waterline_central_second_moment=tuple(
            tuple(float(x) for x in row) for row in value.waterline_central_second_moment
        ),
        submerged_vertices=tuple(tuple(float(x) for x in row) for row in mesh.vertices.tolist()),
        submerged_faces=tuple(tuple(int(x) for x in row) for row in mesh.faces.tolist()),
    )


def find_hydrostatic_equilibrium(
    vertices: Sequence[Sequence[float]],
    faces: Sequence[Sequence[int]],
    *,
    reference_origin: Sequence[float],
    waterline_normal: Sequence[float],
    target_displaced_volume: float,
    relative_volume_tolerance: float = 1e-10,
    maximum_iterations: int = 80,
) -> HydrostaticEquilibrium:
    """Solve only the displaced-volume equation for a closed boundary.

    This does not infer sealed cavities, floodability, or dynamic stability,
    and therefore cannot turn a reconstructed mesh into a verified Faster
    floating-object experiment.
    """
    native = _native_module()
    value = native.find_hydrostatic_equilibrium(
        np.asarray(_vertices(vertices, "vertices"), dtype=np.float64),
        np.asarray(_faces(faces, "faces"), dtype=np.int32),
        np.asarray(_finite_tuple(reference_origin, 3, "reference_origin"), dtype=np.float64),
        np.asarray(_finite_tuple(waterline_normal, 3, "waterline_normal"), dtype=np.float64),
        float(target_displaced_volume), float(relative_volume_tolerance), int(maximum_iterations),
    )
    mesh = value.geometry.submerged_boundary
    geometry = HydrostaticGeometry(
        displaced_volume=float(value.geometry.displaced_volume),
        center_of_buoyancy=tuple(float(x) for x in value.geometry.center_of_buoyancy),
        waterline_area=float(value.geometry.waterline_area),
        waterline_centroid=tuple(float(x) for x in value.geometry.waterline_centroid),
        waterline_central_second_moment=tuple(
            tuple(float(x) for x in row) for row in value.geometry.waterline_central_second_moment
        ),
        submerged_vertices=tuple(tuple(float(x) for x in row) for row in mesh.vertices.tolist()),
        submerged_faces=tuple(tuple(int(x) for x in row) for row in mesh.faces.tolist()),
    )
    return HydrostaticEquilibrium(
        geometry=geometry,
        waterline_signed_offset=float(value.waterline_signed_offset),
        target_displaced_volume=float(value.target_displaced_volume),
        volume_residual=float(value.volume_residual),
        iterations=int(value.iterations),
    )


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
