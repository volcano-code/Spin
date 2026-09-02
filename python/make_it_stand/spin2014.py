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
SPIN_2014_EVIDENCE = frozenset({"PAPER_FORMULA_REPRODUCTION", "RECONSTRUCTED_ASSET"})


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
    initial_fill: float = 0.5
    fractional_epsilon: float = 0.05
    laplacian_weight: float = 1e-3
    center_of_mass_weight: float = 1.0
    boundary_density: float = 1.0
    interior_density: float = 1.0

    def __post_init__(self) -> None:
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
        if self.evidence not in SPIN_2014_EVIDENCE:
            raise ValueError("unsupported Spin-It 2014 evidence class")
        if len(self.minimum) != 3 or len(self.maximum) != 3:
            raise ValueError("minimum and maximum must contain three values")
        if any(
            not math.isfinite(float(lo)) or not math.isfinite(float(hi)) or hi <= lo
            for lo, hi in zip(self.minimum, self.maximum)
        ):
            raise ValueError("grid bounds must be finite and strictly increasing")


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
    )
    return SpinIt2014Result(artifact)
