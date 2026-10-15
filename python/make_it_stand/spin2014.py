"""Python facade for the adaptive-voxel Spin-It 2014 reproduction."""

from __future__ import annotations

from dataclasses import dataclass
import importlib
import json
import math
import os
from typing import Any, Sequence

import numpy as np


SPIN_2014_ARTIFACT_SCHEMA = "spin-it-2014-artifact-v1"
SPIN_2014_MODES = frozenset({"STATIC_BALANCE", "YO_YO", "SPINNING_TOP"})
SPIN_2014_OBJECTIVES = frozenset({
    "YOYO",
    "TOP_INERTIA_ONLY",
    "TOP_LOWERING_ONLY",
    "TOP_MASS_AND_LOWERING",
    "STATIC_LOWERING",
})
SPIN_2014_EVIDENCE = frozenset({"PAPER_FORMULA_REPRODUCTION", "RECONSTRUCTED_ASSET"})

_SPIN_2014_OBJECTIVE_MODES = {
    "YOYO": "YO_YO",
    "TOP_INERTIA_ONLY": "SPINNING_TOP",
    "TOP_LOWERING_ONLY": "SPINNING_TOP",
    "TOP_MASS_AND_LOWERING": "SPINNING_TOP",
    "STATIC_LOWERING": "STATIC_BALANCE",
}
_SPIN_2014_MODE_DEFAULT_OBJECTIVES = {
    "STATIC_BALANCE": "STATIC_LOWERING",
    "YO_YO": "YOYO",
    "SPINNING_TOP": "TOP_MASS_AND_LOWERING",
}


def _native_module() -> Any:
    errors: list[BaseException] = []
    for name in ("make_it_stand._make_it_stand", "_make_it_stand"):
        try:
            return importlib.import_module(name)
        except (ImportError, ModuleNotFoundError) as error:
            errors.append(error)
    raise RuntimeError(
        "the Spin-It 2014 native extension is not installed; configure with "
        "MAKE_IT_STAND_BUILD_PYTHON=ON and MAKE_IT_STAND_BUILD_SPIN_IT_2014=ON"
    ) from errors[-1]


@dataclass(frozen=True)
class SpinIt2014Request:
    dimensions: tuple[int, int, int]
    minimum: tuple[float, float, float]
    maximum: tuple[float, float, float]
    material_mask: tuple[int, ...]
    carveable_mask: tuple[int, ...]
    mode: str = "SPINNING_TOP"
    evidence: str = "PAPER_FORMULA_REPRODUCTION"
    initial_level: int = 3
    maximum_level: int = 7
    maximum_adaptive_passes: int = 8
    maximum_leaf_count: int = 1000000
    initial_fill: float = 0.5
    fractional_epsilon: float = 0.05
    laplacian_weight: float = 1e-3
    center_of_mass_weight: float = 1.0
    boundary_density: float = 1.0
    interior_density: float = 1.0
    binary_constraint_tolerance: float = 1e-6
    slqp_maximum_iterations: int = 120
    slqp_projected_gradient_tolerance: float = 1e-7
    objective: str | None = None
    inertia_weight: float = 1.0
    formulation: str = "LEGACY_ANGLE_ELIMINATED"
    initial_phi: float = 0.0
    slqp_algorithm: str = "LEGACY_V3"
    slqp_restoration_trust_bound: float = 1.0
    slqp_lp_iteration_limit: int = 10000
    slqp_lp_time_limit: float = 10.0
    slqp_lp_feasibility_tolerance: float = 1e-9
    target_origin: tuple[float, float, float] = (0.0, 0.0, 0.0)
    target_axes: tuple[tuple[float, float, float], ...] = (
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    )

    def __post_init__(self) -> None:
        origin = np.asarray(self.target_origin, dtype=np.float64)
        axes = np.asarray(self.target_axes, dtype=np.float64)
        if origin.shape != (3,) or not np.all(np.isfinite(origin)):
            raise ValueError("target_origin must contain three finite values")
        if axes.shape != (3, 3) or not np.all(np.isfinite(axes)):
            raise ValueError("target_axes must be a finite 3x3 matrix")
        if not np.allclose(axes.T @ axes, np.eye(3), rtol=0.0, atol=1e-12) or not math.isclose(
            float(np.linalg.det(axes)), 1.0, rel_tol=0.0, abs_tol=1e-12
        ):
            raise ValueError("target_axes must be a right-handed orthonormal frame")
        object.__setattr__(self, "target_origin", tuple(float(value) for value in origin))
        object.__setattr__(
            self, "target_axes", tuple(tuple(float(value) for value in row) for row in axes)
        )
        if self.slqp_algorithm not in {"LEGACY_V3", "AUDITED_KKT_RESTORATION_V1"}:
            raise ValueError("unsupported SLQP algorithm")
        if self.slqp_algorithm == "AUDITED_KKT_RESTORATION_V1":
            if self.formulation != "PAPER_GIVENS":
                raise ValueError("audited SLQP requires PAPER_GIVENS")
            if any(not math.isfinite(v) or v <= 0 for v in
                   (self.slqp_restoration_trust_bound, self.slqp_lp_time_limit)):
                raise ValueError("audit trust bound and LP time limit must be positive and finite")
            if isinstance(self.slqp_lp_iteration_limit, bool) or not isinstance(self.slqp_lp_iteration_limit, int) or self.slqp_lp_iteration_limit <= 0:
                raise ValueError("LP iteration limit must be a positive integer")
            if not math.isfinite(self.slqp_lp_feasibility_tolerance) or not 1e-10 <= self.slqp_lp_feasibility_tolerance <= 1e-4:
                raise ValueError("LP feasibility tolerance must be in [1e-10,1e-4]")
        if self.formulation not in {"LEGACY_ANGLE_ELIMINATED", "PAPER_GIVENS"}:
            raise ValueError("unsupported Spin-It 2014 formulation")
        if not math.isfinite(self.initial_phi) or (
            (self.formulation == "LEGACY_ANGLE_ELIMINATED" or self.mode == "STATIC_BALANCE")
            and self.initial_phi != 0.0
        ):
            raise ValueError("phi must be finite; legacy and static require initial_phi=0")
        if len(self.dimensions) != 3 or any(
            isinstance(value, bool) or not isinstance(value, int) or value <= 0
            for value in self.dimensions
        ):
            raise ValueError("dimensions must contain three positive integers")
        expected = math.prod(self.dimensions)
        for name in ("material_mask", "carveable_mask"):
            values = tuple(int(value) for value in getattr(self, name))
            if len(values) != expected or any(value not in (0, 1) for value in values):
                raise ValueError(f"{name} must contain dimensions-product binary values")
            object.__setattr__(self, name, values)
        if any(c and not m for c, m in zip(self.carveable_mask, self.material_mask)):
            raise ValueError("carveable_mask must be a subset of material_mask")
        if self.mode not in SPIN_2014_MODES:
            raise ValueError("mode is not part of the Spin-It 2014 objective set")
        if self.objective is not None:
            if self.objective not in SPIN_2014_OBJECTIVES:
                raise ValueError("objective is not part of the Spin-It 2014 objective set")
            if _SPIN_2014_OBJECTIVE_MODES[self.objective] != self.mode:
                raise ValueError("Spin-It 2014 objective does not match mode")
        if self.evidence not in SPIN_2014_EVIDENCE:
            raise ValueError("unsupported Spin-It 2014 evidence class")
        if len(self.minimum) != 3 or len(self.maximum) != 3:
            raise ValueError("minimum and maximum must contain three values")
        if any(
            not math.isfinite(float(lo)) or not math.isfinite(float(hi)) or hi <= lo
            for lo, hi in zip(self.minimum, self.maximum)
        ):
            raise ValueError("grid bounds must be finite and strictly increasing")
        if (
            isinstance(self.maximum_leaf_count, bool)
            or not isinstance(self.maximum_leaf_count, int)
            or self.maximum_leaf_count <= 0
        ):
            raise ValueError("maximum_leaf_count must be a positive integer")
        if (
            not math.isfinite(float(self.binary_constraint_tolerance))
            or self.binary_constraint_tolerance <= 0.0
        ):
            raise ValueError("binary_constraint_tolerance must be positive and finite")
        if (
            isinstance(self.slqp_maximum_iterations, bool)
            or not isinstance(self.slqp_maximum_iterations, int)
            or self.slqp_maximum_iterations <= 0
        ):
            raise ValueError("slqp_maximum_iterations must be a positive integer")
        if (
            not math.isfinite(float(self.slqp_projected_gradient_tolerance))
            or self.slqp_projected_gradient_tolerance <= 0.0
        ):
            raise ValueError("slqp_projected_gradient_tolerance must be positive and finite")
        for name in ("inertia_weight", "center_of_mass_weight"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and nonnegative")


class SpinIt2014Result:
    __slots__ = ("_native", "_diagnostics")

    def __init__(self, native: Any):
        diagnostics = json.loads(str(native.diagnostics_json))
        if diagnostics.get("schema") != SPIN_2014_ARTIFACT_SCHEMA:
            raise ValueError("native Spin-It 2014 diagnostics use an unsupported schema")
        object.__setattr__(self, "_native", native)
        object.__setattr__(self, "_diagnostics", diagnostics)

    def __setattr__(self, name: str, value: Any) -> None:
        raise AttributeError("SpinIt2014Result is immutable")

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
    def formulation(self) -> str:
        return self._diagnostics.get("formulation", "LEGACY_ANGLE_ELIMINATED")

    @property
    def initial_phi(self) -> float:
        return float(self._diagnostics.get("initial_phi", 0.0))

    @property
    def final_phi(self) -> float:
        return float(self._diagnostics.get("final_phi", 0.0))

    @property
    def diagnostics(self) -> dict[str, Any]:
        return json.loads(json.dumps(self._diagnostics))

    def _array(self, name: str) -> np.ndarray:
        value = np.asarray(getattr(self._native, name)).copy()
        value.flags.writeable = False
        return value

    @property
    def outer_vertices(self) -> np.ndarray: return self._array("outer_vertices")
    @property
    def outer_faces(self) -> np.ndarray: return self._array("outer_faces")
    @property
    def design_vertices(self) -> np.ndarray: return self._array("design_vertices")
    @property
    def design_faces(self) -> np.ndarray: return self._array("design_faces")
    @property
    def material_vertices(self) -> np.ndarray: return self._array("material_vertices")
    @property
    def material_faces(self) -> np.ndarray: return self._array("material_faces")
    @property
    def cavity_vertices(self) -> np.ndarray: return self._array("cavity_vertices")
    @property
    def cavity_faces(self) -> np.ndarray: return self._array("cavity_faces")

    def save_viewer_bundle(self, directory: os.PathLike[str] | str) -> tuple[str, ...]:
        root = os.fspath(directory)
        os.makedirs(root, exist_ok=True)
        paths = tuple(os.path.join(root, name) for name in (
            "outer.obj", "design.obj", "material.obj", "cavity.obj", "diagnostics.json"
        ))
        self._native.save_outer_mesh(paths[0])
        self._native.save_design_mesh(paths[1])
        self._native.save_material_mesh(paths[2])
        self._native.save_cavity_mesh(paths[3])
        self._native.save_diagnostics(paths[4])
        return paths


def solve_spin_it_2014(request: SpinIt2014Request) -> SpinIt2014Result:
    native = _native_module()
    artifact = native.solve_spin_it_2014(
        np.asarray(request.dimensions, dtype=np.int32),
        np.asarray(request.minimum, dtype=np.float64),
        np.asarray(request.maximum, dtype=np.float64),
        list(request.material_mask),
        list(request.carveable_mask),
        getattr(native.SpinMode, request.mode),
        getattr(native.SpinIt2014Evidence, request.evidence),
        request.initial_level,
        request.maximum_level,
        request.maximum_adaptive_passes,
        request.initial_fill,
        request.fractional_epsilon,
        request.laplacian_weight,
        request.center_of_mass_weight,
        request.boundary_density,
        request.interior_density,
        request.maximum_leaf_count,
        request.binary_constraint_tolerance,
        request.slqp_maximum_iterations,
        request.slqp_projected_gradient_tolerance,
        (
            getattr(native.Spin2014Objective, request.objective)
            if request.objective is not None
            else None
        ),
        request.inertia_weight,
        request.formulation,
        request.initial_phi,
        request.slqp_algorithm,
        request.slqp_restoration_trust_bound,
        request.slqp_lp_iteration_limit,
        request.slqp_lp_time_limit,
        request.slqp_lp_feasibility_tolerance,
        np.asarray(request.target_origin, dtype=np.float64),
        np.asarray(request.target_axes, dtype=np.float64),
    )
    return SpinIt2014Result(artifact)


def solve_spin_it_2014_cage(
    hollowing: SpinIt2014Result,
    request: SpinIt2014Request,
    *,
    maximum_handle_count: int = 12,
) -> dict[str, Any]:
    """Run the paper's joint cage/fill stage from a hollowing artifact."""
    if not isinstance(hollowing, SpinIt2014Result):
        raise TypeError("hollowing must be a SpinIt2014Result")
    if request.formulation == "PAPER_GIVENS" or hollowing.formulation == "PAPER_GIVENS":
        raise ValueError("paper Givens formulation is unsupported by cage deformation")
    if (
        isinstance(maximum_handle_count, bool)
        or not isinstance(maximum_handle_count, int)
        or maximum_handle_count < 2
    ):
        raise ValueError("maximum_handle_count must be an integer of at least two")
    native = _native_module()
    return dict(native.solve_spin_it_2014_cage(
        hollowing._native,
        np.asarray(request.dimensions, dtype=np.int32),
        np.asarray(request.minimum, dtype=np.float64),
        np.asarray(request.maximum, dtype=np.float64),
        list(request.material_mask),
        list(request.carveable_mask),
        request.center_of_mass_weight,
        maximum_handle_count,
        getattr(
            native.Spin2014Objective,
            request.objective or _SPIN_2014_MODE_DEFAULT_OBJECTIVES[request.mode],
        ),
        request.inertia_weight,
    ))
