"""Strict schema-v1 configuration and friendly owning-model wrapper."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
import importlib
import math
from numbers import Real
import os
from pathlib import Path
from typing import Any, Mapping, Sequence

from .io import DataFormatError, loads_json_strict, prepare_output_path, read_json


CONFIG_SCHEMA = "make-it-stand-balance-config-v1"
CONFIG_SCHEMA_VERSION = 1
DIAGNOSTICS_SCHEMA = "make-it-stand-r9-balance-result-v1"


class ConfigurationError(ValueError):
    """A stable, path-qualified configuration validation error."""


def _fail(path: str, message: str) -> None:
    raise ConfigurationError(f"{path}: {message}")


def _mapping(value: Any, path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        _fail(path, "expected an object")
    if not all(isinstance(key, str) for key in value):
        _fail(path, "object keys must be strings")
    return value


def _keys(
    value: Mapping[str, Any],
    path: str,
    *,
    allowed: set[str],
    required: set[str] = frozenset(),
) -> None:
    unknown = sorted(set(value) - allowed)
    missing = sorted(required - set(value))
    if unknown:
        _fail(path, f"unknown field(s): {', '.join(unknown)}")
    if missing:
        _fail(path, f"missing field(s): {', '.join(missing)}")


def _string(value: Any, path: str, *, nonempty: bool = True) -> str:
    if not isinstance(value, str):
        _fail(path, "expected a string")
    if nonempty and not value:
        _fail(path, "must not be empty")
    return value


def _bool(value: Any, path: str) -> bool:
    if not isinstance(value, bool):
        _fail(path, "expected a boolean")
    return value


def _number(
    value: Any,
    path: str,
    *,
    minimum: float | None = None,
    maximum: float | None = None,
    strictly_positive: bool = False,
) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        _fail(path, "expected a finite number")
    result = float(value)
    if not math.isfinite(result):
        _fail(path, "expected a finite number")
    if strictly_positive and result <= 0.0:
        _fail(path, "must be greater than zero")
    if minimum is not None and result < minimum:
        _fail(path, f"must be at least {minimum:g}")
    if maximum is not None and result > maximum:
        _fail(path, f"must be at most {maximum:g}")
    return result


def _integer(
    value: Any,
    path: str,
    *,
    minimum: int | None = None,
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        _fail(path, "expected an integer")
    if minimum is not None and value < minimum:
        _fail(path, f"must be at least {minimum}")
    return value


def _vector3(value: Any, path: str) -> tuple[float, float, float]:
    if isinstance(value, (str, bytes, bytearray, Mapping)):
        _fail(path, "expected an array of exactly three numbers")
    try:
        if len(value) != 3:
            _fail(path, "expected an array of exactly three numbers")
        items = tuple(value[index] for index in range(3))
    except (TypeError, IndexError, KeyError):
        _fail(path, "expected an array of exactly three numbers")
    return tuple(  # type: ignore[return-value]
        _number(item, f"{path}[{index}]")
        for index, item in enumerate(items)
    )


def _matrix3(value: Any, path: str) -> tuple[tuple[float, float, float], ...]:
    if isinstance(value, (str, bytes, bytearray, Mapping)):
        _fail(path, "expected a 3 by 3 array")
    try:
        if len(value) != 3:
            _fail(path, "expected a 3 by 3 array")
        rows = tuple(value[index] for index in range(3))
    except (TypeError, IndexError, KeyError):
        _fail(path, "expected a 3 by 3 array")
    return tuple(
        _vector3(row, f"{path}[{index}]")
        for index, row in enumerate(rows)
    )


def _nonzero(vector: Sequence[float], path: str) -> None:
    if math.sqrt(sum(component * component for component in vector)) <= 1.0e-15:
        _fail(path, "must be non-zero")


@dataclass(frozen=True)
class ModelConfig:
    path: str
    voxel_resolution: int = 50
    minimum_wall_thickness: float = 0.0
    bbw_max_iterations: int = 200
    bbw_tolerance: float = 1.0e-8


@dataclass(frozen=True)
class SupportConfig:
    kind: str
    gamma_radians: float
    polygon: tuple[tuple[float, float, float], ...] = ()
    contact_threshold: float | None = None


@dataclass(frozen=True)
class StandingConstraintConfig:
    kind: str
    gravity: tuple[float, float, float]
    gamma_radians: float
    polygon: tuple[tuple[float, float, float], ...] = ()
    contact_threshold: float | None = None


@dataclass(frozen=True)
class SuspensionConfig:
    attachment_point: tuple[float, float, float]
    allowed_deviation_radians: float


@dataclass(frozen=True)
class HandleConfig:
    rest_pivot: tuple[float, float, float]
    fixed_rotation: tuple[tuple[float, float, float], ...] = (
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    )


@dataclass(frozen=True)
class RunConfig:
    initial_mu: float = 0.2
    maximum_carving_blocks: int = 20
    maximum_deformation_iterations_per_block: int = 200
    gradient_tolerance: float = 1.0e-8
    shape_scale_lambda: float = 1.0
    density: float = 1.0
    minimum_positive_mass: float = 0.0
    minimum_positive_cell_mass: float = 0.0
    minimum_carving_removals_to_continue: int = 1
    initial_step: float = 1.0
    minimum_step: float = 0.4
    reduction_ratio: float = 0.8
    required_relative_decrease: float = 0.01
    mu_decrement: float = 0.05
    enable_carving: bool = field(default=True, kw_only=True)
    enable_deformation: bool = field(default=True, kw_only=True)
    enable_laplacian: bool = field(default=True, kw_only=True)


@dataclass(frozen=True)
class OutputConfig:
    diagnostics: str = "diagnostics.json"
    outer_mesh: str | None = None
    inner_mesh: str | None = None
    printable_mesh: str | None = None


@dataclass(frozen=True)
class BalanceConfig:
    model: ModelConfig
    gravity: tuple[float, float, float]
    support: SupportConfig | None
    handles: tuple[HandleConfig, ...]
    suspensions: tuple[SuspensionConfig, ...]
    run: RunConfig = field(default_factory=RunConfig)
    outputs: OutputConfig = field(default_factory=OutputConfig)
    schema: str = CONFIG_SCHEMA
    schema_version: int = CONFIG_SCHEMA_VERSION
    # R10 is keyword-only so every pre-R10 positional constructor keeps its
    # exact field mapping.
    standing_constraints: tuple[StandingConstraintConfig, ...] = field(
        default_factory=tuple,
        kw_only=True,
    )


_MODEL_FIELDS = {
    "path", "voxel_resolution", "minimum_wall_thickness",
    "bbw_max_iterations", "bbw_tolerance",
}
_RUN_FIELDS = set(RunConfig.__dataclass_fields__)
_OUTPUT_FIELDS = {"diagnostics", "outer_mesh", "inner_mesh", "printable_mesh"}


def validate_config(value: Any) -> BalanceConfig:
    """Validate and normalize one exact R9 JSON configuration object."""

    root = _mapping(value, "config")
    _keys(
        root,
        "config",
        allowed={
            "schema", "schema_version", "model", "gravity", "support",
            "standing_constraints", "handles", "suspensions", "run", "outputs",
        },
        required={"schema", "schema_version", "model", "gravity", "outputs"},
    )
    schema = _string(root["schema"], "config.schema")
    if schema != CONFIG_SCHEMA:
        _fail("config.schema", f"expected {CONFIG_SCHEMA!r}")
    version = _integer(root["schema_version"], "config.schema_version")
    if version != CONFIG_SCHEMA_VERSION:
        _fail("config.schema_version", f"expected {CONFIG_SCHEMA_VERSION}")

    raw_model = _mapping(root["model"], "config.model")
    _keys(raw_model, "config.model", allowed=_MODEL_FIELDS, required={"path"})
    model = ModelConfig(
        path=_string(raw_model["path"], "config.model.path"),
        voxel_resolution=_integer(
            raw_model.get("voxel_resolution", 50),
            "config.model.voxel_resolution",
            minimum=3,
        ),
        minimum_wall_thickness=_number(
            raw_model.get("minimum_wall_thickness", 0.0),
            "config.model.minimum_wall_thickness",
            minimum=0.0,
        ),
        bbw_max_iterations=_integer(
            raw_model.get("bbw_max_iterations", 200),
            "config.model.bbw_max_iterations",
            minimum=1,
        ),
        bbw_tolerance=_number(
            raw_model.get("bbw_tolerance", 1.0e-8),
            "config.model.bbw_tolerance",
            strictly_positive=True,
        ),
    )

    gravity = _vector3(root["gravity"], "config.gravity")
    _nonzero(gravity, "config.gravity")

    support: SupportConfig | None = None
    if root.get("support") is not None:
        raw_support = _mapping(root["support"], "config.support")
        _keys(
            raw_support,
            "config.support",
            allowed={"kind", "polygon", "contact_threshold", "gamma_radians"},
            required={"kind", "gamma_radians"},
        )
        kind = _string(raw_support["kind"], "config.support.kind")
        gamma = _number(
            raw_support["gamma_radians"],
            "config.support.gamma_radians",
            minimum=0.0,
            maximum=None,
        )
        if gamma >= math.pi / 2.0:
            _fail("config.support.gamma_radians", "must be less than pi/2")
        if kind == "polygon":
            if "contact_threshold" in raw_support:
                _fail("config.support", "contact_threshold is valid only for kind 'mesh'")
            raw_polygon = raw_support.get("polygon")
            if not isinstance(raw_polygon, (list, tuple)) or len(raw_polygon) < 3:
                _fail("config.support.polygon", "expected at least three 3D points")
            polygon = tuple(
                _vector3(point, f"config.support.polygon[{index}]")
                for index, point in enumerate(raw_polygon)
            )
            support = SupportConfig(kind=kind, gamma_radians=gamma, polygon=polygon)
        elif kind == "mesh":
            if "polygon" in raw_support:
                _fail("config.support", "polygon is valid only for kind 'polygon'")
            threshold = _number(
                raw_support.get("contact_threshold", 0.0),
                "config.support.contact_threshold",
                minimum=0.0,
            )
            support = SupportConfig(
                kind=kind,
                gamma_radians=gamma,
                contact_threshold=threshold,
            )
        else:
            _fail("config.support.kind", "expected 'polygon' or 'mesh'")

    standing_constraints: list[StandingConstraintConfig] = []
    raw_standings = root.get("standing_constraints", [])
    if not isinstance(raw_standings, list):
        _fail("config.standing_constraints", "expected an array")
    for index, item in enumerate(raw_standings):
        path = f"config.standing_constraints[{index}]"
        raw_standing = _mapping(item, path)
        _keys(
            raw_standing,
            path,
            allowed={
                "kind", "gravity", "polygon", "contact_threshold",
                "gamma_radians",
            },
            required={"kind", "gravity", "gamma_radians"},
        )
        kind = _string(raw_standing["kind"], f"{path}.kind")
        constraint_gravity = _vector3(
            raw_standing["gravity"], f"{path}.gravity"
        )
        _nonzero(constraint_gravity, f"{path}.gravity")
        gamma = _number(
            raw_standing["gamma_radians"],
            f"{path}.gamma_radians",
            minimum=0.0,
        )
        if gamma >= math.pi / 2.0:
            _fail(f"{path}.gamma_radians", "must be less than pi/2")
        if kind == "polygon":
            if "contact_threshold" in raw_standing:
                _fail(path, "contact_threshold is valid only for kind 'mesh'")
            raw_polygon = raw_standing.get("polygon")
            if not isinstance(raw_polygon, (list, tuple)) or len(raw_polygon) < 3:
                _fail(f"{path}.polygon", "expected at least three 3D points")
            polygon = tuple(
                _vector3(point, f"{path}.polygon[{point_index}]")
                for point_index, point in enumerate(raw_polygon)
            )
            standing_constraints.append(
                StandingConstraintConfig(
                    kind=kind,
                    gravity=constraint_gravity,
                    gamma_radians=gamma,
                    polygon=polygon,
                )
            )
        elif kind == "mesh":
            if "polygon" in raw_standing:
                _fail(path, "polygon is valid only for kind 'polygon'")
            threshold = _number(
                raw_standing.get("contact_threshold", 0.0),
                f"{path}.contact_threshold",
                minimum=0.0,
            )
            standing_constraints.append(
                StandingConstraintConfig(
                    kind=kind,
                    gravity=constraint_gravity,
                    gamma_radians=gamma,
                    contact_threshold=threshold,
                )
            )
        else:
            _fail(f"{path}.kind", "expected 'polygon' or 'mesh'")

    handles: list[HandleConfig] = []
    raw_handles = root.get("handles", [])
    if not isinstance(raw_handles, list):
        _fail("config.handles", "expected an array")
    for index, item in enumerate(raw_handles):
        path = f"config.handles[{index}]"
        raw_handle = _mapping(item, path)
        _keys(
            raw_handle,
            path,
            allowed={"rest_pivot", "fixed_rotation"},
            required={"rest_pivot"},
        )
        handles.append(
            HandleConfig(
                rest_pivot=_vector3(raw_handle["rest_pivot"], f"{path}.rest_pivot"),
                fixed_rotation=_matrix3(
                    raw_handle.get(
                        "fixed_rotation",
                        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                    ),
                    f"{path}.fixed_rotation",
                ),
            )
        )

    suspensions: list[SuspensionConfig] = []
    raw_suspensions = root.get("suspensions", [])
    if not isinstance(raw_suspensions, list):
        _fail("config.suspensions", "expected an array")
    for index, item in enumerate(raw_suspensions):
        path = f"config.suspensions[{index}]"
        raw_suspension = _mapping(item, path)
        _keys(
            raw_suspension,
            path,
            allowed={"attachment_point", "allowed_deviation_radians"},
            required={"attachment_point", "allowed_deviation_radians"},
        )
        suspensions.append(
            SuspensionConfig(
                attachment_point=_vector3(
                    raw_suspension["attachment_point"],
                    f"{path}.attachment_point",
                ),
                allowed_deviation_radians=_number(
                    raw_suspension["allowed_deviation_radians"],
                    f"{path}.allowed_deviation_radians",
                    minimum=0.0,
                    maximum=None,
                ),
            )
        )
        if suspensions[-1].allowed_deviation_radians >= math.pi / 2.0:
            _fail(
                f"{path}.allowed_deviation_radians",
                "must be less than pi/2",
            )
    if support is None and not standing_constraints and not suspensions:
        _fail(
            "config",
            "at least one support, standing, or suspension constraint is required",
        )

    raw_run = _mapping(root.get("run", {}), "config.run")
    _keys(raw_run, "config.run", allowed=_RUN_FIELDS)
    run_defaults = RunConfig()
    run_values: dict[str, Any] = {}
    integer_run_fields = {
        "maximum_carving_blocks",
        "maximum_deformation_iterations_per_block",
        "minimum_carving_removals_to_continue",
    }
    boolean_run_fields = {
        "enable_carving",
        "enable_deformation",
        "enable_laplacian",
    }
    for name in _RUN_FIELDS:
        raw = raw_run.get(name, getattr(run_defaults, name))
        if name in boolean_run_fields:
            run_values[name] = _bool(raw, f"config.run.{name}")
        elif name in integer_run_fields:
            minimum = 0 if name == "minimum_carving_removals_to_continue" else 1
            run_values[name] = _integer(raw, f"config.run.{name}", minimum=minimum)
        else:
            run_values[name] = _number(raw, f"config.run.{name}")
    positive = {
        "density", "initial_step", "minimum_step", "mu_decrement",
    }
    for name in positive - {"bbw_tolerance"}:
        if run_values[name] <= 0.0:
            _fail(f"config.run.{name}", "must be greater than zero")
    for name in (
        "initial_mu", "gradient_tolerance", "shape_scale_lambda", "minimum_positive_mass",
        "minimum_positive_cell_mass", "required_relative_decrease",
    ):
        if run_values[name] < 0.0:
            _fail(f"config.run.{name}", "must be non-negative")
    if not 0.0 < run_values["reduction_ratio"] < 1.0:
        _fail("config.run.reduction_ratio", "must be strictly between zero and one")
    if run_values["minimum_step"] > run_values["initial_step"]:
        _fail("config.run.minimum_step", "must not exceed initial_step")
    if run_values["initial_mu"] > 1.0:
        _fail("config.run.initial_mu", "must not exceed one")
    if run_values["required_relative_decrease"] >= 1.0:
        _fail("config.run.required_relative_decrease", "must be less than one")
    if run_values["mu_decrement"] > 1.0:
        _fail("config.run.mu_decrement", "must not exceed one")
    if not run_values["enable_carving"] and not run_values["enable_deformation"]:
        _fail(
            "config.run",
            "enable_carving and enable_deformation must not both be false",
        )
    run = RunConfig(**run_values)

    raw_outputs = _mapping(root["outputs"], "config.outputs")
    _keys(raw_outputs, "config.outputs", allowed=_OUTPUT_FIELDS)
    output_values: dict[str, str | None] = {}
    for name in _OUTPUT_FIELDS:
        raw = raw_outputs.get(name, "diagnostics.json" if name == "diagnostics" else None)
        if raw is not None:
            raw = _string(raw, f"config.outputs.{name}")
            candidate = Path(raw)
            if candidate.is_absolute() or ".." in candidate.parts:
                _fail(
                    f"config.outputs.{name}",
                    "must be a relative path without parent traversal",
                )
        output_values[name] = raw
    outputs = OutputConfig(**output_values)  # type: ignore[arg-type]
    if outputs.diagnostics is None:
        _fail("config.outputs.diagnostics", "must be a non-empty relative path")
    selected_outputs = [
        Path(path).as_posix()
        for path in (
            outputs.diagnostics,
            outputs.outer_mesh,
            outputs.inner_mesh,
            outputs.printable_mesh,
        )
        if path is not None
    ]
    if any(path in {"", "."} for path in selected_outputs):
        _fail("config.outputs", "output paths must name files")
    if len(set(selected_outputs)) != len(selected_outputs):
        _fail("config.outputs", "output paths must be distinct")
    selected_parts = [Path(path).parts for path in selected_outputs]
    for index, left in enumerate(selected_parts):
        for right in selected_parts[index + 1:]:
            shared = min(len(left), len(right))
            if left[:shared] == right[:shared]:
                _fail(
                    "config.outputs",
                    "an output file path must not contain another output path",
                )

    return BalanceConfig(
        model=model,
        gravity=gravity,
        support=support,
        standing_constraints=tuple(standing_constraints),
        handles=tuple(handles),
        suspensions=tuple(suspensions),
        run=run,
        outputs=outputs,
    )


def _schema_tree(value: Any) -> Any:
    """Convert dataclass output to the JSON container forms accepted by v1."""

    if isinstance(value, Mapping):
        return {key: _schema_tree(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_schema_tree(item) for item in value]
    if isinstance(value, list):
        return [_schema_tree(item) for item in value]
    return value


def validate_balance_config(value: Any) -> BalanceConfig:
    """Revalidate a public ``BalanceConfig`` through the canonical v1 schema.

    Frozen dataclasses prevent direct assignment but ``dataclasses.replace``
    can still construct an invalid value.  Public execution entry points call
    this function so the parsed-JSON and programmatic APIs have one contract.
    """

    if not isinstance(value, BalanceConfig):
        _fail("config", "expected a BalanceConfig instance")
    try:
        raw = _schema_tree(asdict(value))
    except (TypeError, ValueError) as error:
        raise ConfigurationError(f"config: invalid dataclass tree: {error}") from error
    support = raw.get("support")
    if isinstance(support, dict):
        if support.get("kind") == "polygon":
            if support.get("contact_threshold") is not None:
                _fail(
                    "config.support",
                    "contact_threshold is valid only for kind 'mesh'",
                )
            support.pop("contact_threshold", None)
        elif support.get("kind") == "mesh":
            if support.get("polygon"):
                _fail(
                    "config.support",
                    "polygon is valid only for kind 'polygon'",
                )
            support.pop("polygon", None)
    raw_standings = raw.get("standing_constraints", [])
    if not isinstance(raw_standings, (list, tuple)):
        _fail("config.standing_constraints", "expected an array")
    for index, standing in enumerate(raw_standings):
        if not isinstance(standing, dict):
            _fail(
                f"config.standing_constraints[{index}]",
                "expected an object",
            )
        if standing.get("kind") == "polygon":
            if standing.get("contact_threshold") is not None:
                _fail(
                    "config.standing_constraints",
                    "contact_threshold is valid only for kind 'mesh'",
                )
            standing.pop("contact_threshold", None)
        elif standing.get("kind") == "mesh":
            if standing.get("polygon"):
                _fail(
                    "config.standing_constraints",
                    "polygon is valid only for kind 'polygon'",
                )
            standing.pop("polygon", None)
    return validate_config(raw)


def _validation_document(
    *,
    model: Mapping[str, Any] | None = None,
    run: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "schema": CONFIG_SCHEMA,
        "schema_version": CONFIG_SCHEMA_VERSION,
        "model": model if model is not None else {"path": "<validation>"},
        "gravity": [0.0, -1.0, 0.0],
        "support": {"kind": "mesh", "gamma_radians": 0.0},
        "handles": [],
        "suspensions": [],
        "run": run if run is not None else {},
        "outputs": {},
    }


def _validate_model_config(value: Any) -> ModelConfig:
    if not isinstance(value, ModelConfig):
        _fail("config.model", "expected a ModelConfig instance")
    try:
        raw = _schema_tree(asdict(value))
    except (TypeError, ValueError) as error:
        raise ConfigurationError(f"config.model: invalid dataclass: {error}") from error
    return validate_config(_validation_document(model=raw)).model


def _validate_run_config(value: Any) -> RunConfig:
    if not isinstance(value, RunConfig):
        _fail("config.run", "expected a RunConfig instance")
    try:
        raw = _schema_tree(asdict(value))
    except (TypeError, ValueError) as error:
        raise ConfigurationError(f"config.run: invalid dataclass: {error}") from error
    return validate_config(_validation_document(run=raw)).run


def load_config(path: os.PathLike[str] | str) -> BalanceConfig:
    source = Path(path)
    config = validate_config(read_json(source))
    mesh_path = Path(config.model.path)
    if not mesh_path.is_absolute():
        mesh_path = (source.resolve().parent / mesh_path).resolve()
        config = replace(
            config,
            model=replace(config.model, path=os.fspath(mesh_path)),
        )
    return config


def _load_native() -> Any:
    errors: list[BaseException] = []
    for module_name in ("make_it_stand._make_it_stand", "_make_it_stand"):
        try:
            return importlib.import_module(module_name)
        except (ImportError, ModuleNotFoundError) as error:
            errors.append(error)
    message = (
        "the optional _make_it_stand extension is not installed; build the "
        "R9 Python-enabled target or install the wheel"
    )
    raise RuntimeError(message) from errors[-1]


def _native_model_config(native: Any, config: ModelConfig) -> Any:
    result = native.BalanceModelConfig()
    for name in (
        "voxel_resolution", "minimum_wall_thickness",
        "bbw_max_iterations", "bbw_tolerance",
    ):
        setattr(result, name, getattr(config, name))
    return result


def _native_run_config(native: Any, config: RunConfig) -> Any:
    result = native.BalanceRunConfig()
    for name in _RUN_FIELDS:
        setattr(result, name, getattr(config, name))
    return result


class BalanceResult:
    """Read-only Python wrapper around an owning native result snapshot."""

    def __init__(self, native_result: Any):
        self._native = native_result

    @property
    def status(self) -> str:
        return str(self._native.status)

    @property
    def succeeded(self) -> bool:
        return bool(self._native.succeeded)

    @property
    def manufacturing_ready(self) -> bool:
        return bool(self._native.manufacturing_ready)

    @property
    def manufacturing_failures(self) -> tuple[str, ...]:
        return tuple(str(item) for item in self._native.manufacturing_failures)

    @property
    def outer_vertices(self) -> Any:
        return self._native.outer_vertices

    @property
    def outer_faces(self) -> Any:
        return self._native.outer_faces

    @property
    def inner_vertices(self) -> Any:
        return self._native.inner_vertices

    @property
    def inner_faces(self) -> Any:
        return self._native.inner_faces

    @property
    def manufacturing_report(self) -> dict[str, Any]:
        value = loads_json_strict(str(self._native.manufacturing_report_json))
        if not isinstance(value, dict):
            raise DataFormatError("native manufacturing_report_json must contain an object")
        return value

    @property
    def diagnostics(self) -> dict[str, Any]:
        value = loads_json_strict(str(self._native.diagnostics_json))
        if not isinstance(value, dict):
            raise DataFormatError("native diagnostics_json must contain an object")
        if value.get("schema") != DIAGNOSTICS_SCHEMA:
            raise DataFormatError(
                f"native diagnostics_json has unsupported schema; "
                f"expected {DIAGNOSTICS_SCHEMA!r}"
            )
        version = value.get("schema_version")
        if type(version) is not int or version != 1:
            raise DataFormatError(
                "native diagnostics_json has unsupported schema_version; expected 1"
            )
        if not isinstance(value.get("status"), str) or not value["status"]:
            raise DataFormatError(
                "native diagnostics_json must contain a non-empty status"
            )
        if not isinstance(value.get("succeeded"), bool):
            raise DataFormatError(
                "native diagnostics_json must contain a boolean succeeded field"
            )
        return value

    @property
    def performance_report(self) -> dict[str, Any]:
        value = loads_json_strict(str(self._native.performance_report_json))
        if not isinstance(value, dict):
            raise DataFormatError(
                "native performance_report_json must contain an object"
            )
        if set(value) != {
            "schema", "schema_version", "clock", "scope", "phases",
            "total_seconds",
        }:
            raise DataFormatError(
                "native performance_report_json has an invalid field set"
            )
        if value.get("schema") != "make-it-stand-r10-performance-report-v1":
            raise DataFormatError(
                "native performance_report_json has unsupported schema"
            )
        if value.get("schema_version") != 1 or value.get("clock") != "steady_clock":
            raise DataFormatError(
                "native performance_report_json has unsupported version or clock"
            )
        if not isinstance(value.get("scope"), str) or not value["scope"]:
            raise DataFormatError(
                "native performance_report_json has invalid scope"
            )
        phases = value.get("phases")
        required_phases = {
            "input_snapshot_seconds",
            "voxelization_seconds",
            "domain_and_handles_seconds",
            "deformation_system_and_bbw_seconds",
            "target_and_optimizer_setup_seconds",
            "optimization_seconds",
            "manufacturing_audit_seconds",
            "diagnostics_payload_setup_seconds",
            "diagnostics_serialization_seconds",
        }
        if not isinstance(phases, dict) or set(phases) != required_phases:
            raise DataFormatError(
                "native performance_report_json has invalid phases"
            )
        numbers = tuple(phases.values()) + (value.get("total_seconds"),)
        if any(
            isinstance(item, bool)
            or not isinstance(item, (int, float))
            or not math.isfinite(float(item))
            or float(item) < 0.0
            for item in numbers
        ):
            raise DataFormatError(
                "native performance_report_json contains invalid durations"
            )
        total = float(value["total_seconds"])
        tolerance = max(1.0e-12, total * 1.0e-9)
        if sum(float(item) for item in phases.values()) > total + tolerance:
            raise DataFormatError(
                "native performance_report_json phases exceed total duration"
            )
        return value

    @property
    def native(self) -> Any:
        return self._native

    def save_outer_mesh(self, path: os.PathLike[str] | str) -> None:
        self._native.save_outer_mesh(os.fspath(path))

    def save_inner_mesh(self, path: os.PathLike[str] | str) -> None:
        self._native.save_inner_mesh(os.fspath(path))

    def save_printable_mesh(self, path: os.PathLike[str] | str) -> None:
        self._native.save_printable_mesh(os.fspath(path))

    def save_diagnostics(self, path: os.PathLike[str] | str) -> None:
        self._native.save_diagnostics(os.fspath(path))


class BalanceModel:
    """Friendly builder that keeps the native binding surface intentionally thin."""

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
            if faces is None and isinstance(mesh, (str, os.PathLike)):
                settings = ModelConfig(path=os.fspath(mesh))
            else:
                settings = ModelConfig(path="<memory>")
        settings = _validate_model_config(settings)
        native = native_module if native_module is not None else _load_native()
        self._native_module = native
        native_config = _native_model_config(native, settings)
        if faces is None:
            self._native = native.BalanceModel(os.fspath(mesh), native_config)
        else:
            self._native = native.BalanceModel(mesh, faces, native_config)

    @property
    def native(self) -> Any:
        return self._native

    def set_gravity(self, gravity: Sequence[float]) -> "BalanceModel":
        vector = _vector3(gravity, "gravity")
        _nonzero(vector, "gravity")
        self._native.set_gravity(vector)
        return self

    def set_support_polygon(
        self,
        polygon: Sequence[Sequence[float]],
        gamma_radians: float,
    ) -> "BalanceModel":
        points = tuple(
            _vector3(point, f"support_polygon[{index}]")
            for index, point in enumerate(polygon)
        )
        if len(points) < 3:
            _fail("support_polygon", "expected at least three points")
        gamma = _number(gamma_radians, "gamma_radians", minimum=0.0)
        if gamma >= math.pi / 2.0:
            _fail("gamma_radians", "must be less than pi/2")
        self._native.set_support_polygon(points, gamma)
        return self

    def set_support_from_mesh(
        self,
        contact_threshold: float,
        gamma_radians: float,
    ) -> "BalanceModel":
        threshold = _number(contact_threshold, "contact_threshold", minimum=0.0)
        gamma = _number(gamma_radians, "gamma_radians", minimum=0.0)
        if gamma >= math.pi / 2.0:
            _fail("gamma_radians", "must be less than pi/2")
        self._native.set_support_from_mesh(threshold, gamma)
        return self

    def add_standing_constraint_polygon(
        self,
        polygon: Sequence[Sequence[float]],
        gravity_direction: Sequence[float],
        gamma_radians: float,
    ) -> "BalanceModel":
        points = tuple(
            _vector3(point, f"standing_polygon[{index}]")
            for index, point in enumerate(polygon)
        )
        if len(points) < 3:
            _fail("standing_polygon", "expected at least three points")
        gravity = _vector3(gravity_direction, "gravity_direction")
        _nonzero(gravity, "gravity_direction")
        gamma = _number(gamma_radians, "gamma_radians", minimum=0.0)
        if gamma >= math.pi / 2.0:
            _fail("gamma_radians", "must be less than pi/2")
        self._native.add_standing_constraint_polygon(points, gravity, gamma)
        return self

    def add_standing_constraint_from_mesh(
        self,
        gravity_direction: Sequence[float],
        contact_threshold: float,
        gamma_radians: float,
    ) -> "BalanceModel":
        gravity = _vector3(gravity_direction, "gravity_direction")
        _nonzero(gravity, "gravity_direction")
        threshold = _number(contact_threshold, "contact_threshold", minimum=0.0)
        gamma = _number(gamma_radians, "gamma_radians", minimum=0.0)
        if gamma >= math.pi / 2.0:
            _fail("gamma_radians", "must be less than pi/2")
        self._native.add_standing_constraint_from_mesh(
            gravity,
            threshold,
            gamma,
        )
        return self

    def add_handle(
        self,
        rest_pivot: Sequence[float],
        fixed_rotation: Sequence[Sequence[float]] | None = None,
    ) -> "BalanceModel":
        pivot = _vector3(rest_pivot, "rest_pivot")
        if fixed_rotation is None:
            self._native.add_handle(pivot)
        else:
            self._native.add_handle(pivot, _matrix3(fixed_rotation, "fixed_rotation"))
        return self

    def add_suspension_constraint(
        self,
        attachment_point: Sequence[float],
        allowed_deviation_radians: float,
    ) -> "BalanceModel":
        point = _vector3(attachment_point, "attachment_point")
        angle = _number(
            allowed_deviation_radians,
            "allowed_deviation_radians",
            minimum=0.0,
        )
        if angle >= math.pi / 2.0:
            _fail("allowed_deviation_radians", "must be less than pi/2")
        self._native.add_suspension_constraint(point, angle)
        return self

    def clear_handles(self) -> "BalanceModel":
        self._native.clear_handles()
        return self

    def clear_constraints(self) -> "BalanceModel":
        self._native.clear_constraints()
        return self

    def optimize(self, config: RunConfig | None = None) -> BalanceResult:
        settings = _validate_run_config(config if config is not None else RunConfig())
        return BalanceResult(
            self._native.optimize(_native_run_config(self._native_module, settings))
        )


def _build_model(config: BalanceConfig, *, native_module: Any | None = None) -> BalanceModel:
    model = BalanceModel(
        config.model.path,
        config=config.model,
        native_module=native_module,
    )
    model.set_gravity(config.gravity)
    if config.support is not None:
        if config.support.kind == "polygon":
            model.set_support_polygon(config.support.polygon, config.support.gamma_radians)
        else:
            assert config.support.contact_threshold is not None
            model.set_support_from_mesh(
                config.support.contact_threshold,
                config.support.gamma_radians,
            )
    for standing in config.standing_constraints:
        if standing.kind == "polygon":
            model.add_standing_constraint_polygon(
                standing.polygon,
                standing.gravity,
                standing.gamma_radians,
            )
        else:
            assert standing.contact_threshold is not None
            model.add_standing_constraint_from_mesh(
                standing.gravity,
                standing.contact_threshold,
                standing.gamma_radians,
            )
    for suspension in config.suspensions:
        model.add_suspension_constraint(
            suspension.attachment_point,
            suspension.allowed_deviation_radians,
        )
    for handle in config.handles:
        model.add_handle(handle.rest_pivot, handle.fixed_rotation)
    return model


def build_model(config: BalanceConfig, *, native_module: Any | None = None) -> BalanceModel:
    """Build a native model from a fully revalidated public configuration."""

    return _build_model(validate_balance_config(config), native_module=native_module)


def run_config(
    config: BalanceConfig,
    *,
    output_directory: os.PathLike[str] | str = ".",
    native_module: Any | None = None,
) -> BalanceResult:
    """Run one validated configuration and write only explicitly named outputs."""

    config = validate_balance_config(config)
    outputs = config.outputs
    destinations = {
        "diagnostics": prepare_output_path(output_directory, outputs.diagnostics),
        "outer_mesh": (
            prepare_output_path(output_directory, outputs.outer_mesh)
            if outputs.outer_mesh is not None else None
        ),
        "inner_mesh": (
            prepare_output_path(output_directory, outputs.inner_mesh)
            if outputs.inner_mesh is not None else None
        ),
        "printable_mesh": (
            prepare_output_path(output_directory, outputs.printable_mesh)
            if outputs.printable_mesh is not None else None
        ),
    }
    result = _build_model(config, native_module=native_module).optimize(config.run)

    # Publish the run record before optional mesh exports. Each native save is
    # atomic for its own destination, but the set of files is intentionally not
    # presented as one cross-file transaction. In particular, a guarded
    # printable export may reject a non-ready result without erasing the run's
    # diagnostics.
    result.save_diagnostics(
        prepare_output_path(output_directory, outputs.diagnostics)
    )
    if destinations["outer_mesh"] is not None:
        assert outputs.outer_mesh is not None
        result.save_outer_mesh(prepare_output_path(output_directory, outputs.outer_mesh))
    if destinations["inner_mesh"] is not None:
        assert outputs.inner_mesh is not None
        result.save_inner_mesh(prepare_output_path(output_directory, outputs.inner_mesh))
    if destinations["printable_mesh"] is not None:
        assert outputs.printable_mesh is not None
        result.save_printable_mesh(
            prepare_output_path(output_directory, outputs.printable_mesh)
        )
    return result


BalanceModelConfig = ModelConfig
BalanceRunConfig = RunConfig


__all__ = [
    "BalanceConfig", "BalanceModel", "BalanceModelConfig", "BalanceResult",
    "BalanceRunConfig", "CONFIG_SCHEMA",
    "CONFIG_SCHEMA_VERSION", "ConfigurationError", "DIAGNOSTICS_SCHEMA",
    "HandleConfig", "ModelConfig", "OutputConfig", "RunConfig",
    "StandingConstraintConfig", "SupportConfig", "SuspensionConfig",
    "build_model", "load_config",
    "run_config", "validate_balance_config", "validate_config",
]
