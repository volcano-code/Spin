"""Auditable experiment records for the R10 paper-alignment phase.

This module deliberately separates values copied from the 2013 paper from
measurements produced by this repository.  It records *paper-governed* or
*configuration-anchored* experiments without promoting them to exact
reproductions.  Documents are serialized with one deterministic JSON
encoding and bind manifests, inputs, outputs, and run records by SHA-256.

The module has no dependency on the optional native extension.  A caller can
therefore use :class:`RunRecorder` around either the C++ executable or the R9
Python API while retaining the same evidence contract.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from contextlib import contextmanager
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import platform
import re
import stat
import time
from typing import Any, Iterable, Iterator, Mapping, Sequence

from .io import DataFormatError, loads_json_strict


MANIFEST_SCHEMA = "make-it-stand-r10-experiment-manifest-v1"
RUN_RECORD_SCHEMA = "make-it-stand-r10-run-record-v1"
AGGREGATE_SCHEMA = "make-it-stand-r10-aggregate-v1"
SCHEMA_VERSION = 1
CANONICALIZATION = "utf8-json-sorted-keys-compact-v1"
OUTER_TIMING_SCOPE = "model_setup_start-through_run_record_finalization"

CLAIM_LEVELS = frozenset(
    {"paper_governed", "configuration_anchored", "local_robustness", "smoke"}
)
RUN_STATUSES = frozenset({"succeeded", "failed", "skipped", "invalid"})
METRIC_SOURCES = frozenset({"paper_reference", "local_measurement"})

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_GIT_OBJECT_ID = re.compile(r"^(?:[0-9a-f]{40}|[0-9a-f]{64})$")
_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_WINDOWS_ABSOLUTE = re.compile(r"^[A-Za-z]:[\\/]")
_EMBEDDED_HOST_PATH = re.compile(
    r"(?:^|[\s=]|-I|-L)(?:/(?!/)[A-Za-z0-9._-]+(?:/[A-Za-z0-9._-]+)+|"
    r"[A-Za-z]:[\\/]|\\\\[^\\\s]+[\\/])"
)
_PROHIBITED_HOST_KEYS = frozenset(
    {"absolute_path", "cwd", "host", "host_name", "hostname", "working_directory"}
)


class R10RecordError(ValueError):
    """Raised when an R10 evidence document violates its stable contract."""


def _fail(path: str, message: str) -> None:
    raise R10RecordError(f"{path}: {message}")


def _mapping(value: Any, path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        _fail(path, "expected an object")
    if any(not isinstance(key, str) for key in value):
        _fail(path, "object keys must be strings")
    return value


def _sequence(value: Any, path: str) -> Sequence[Any]:
    if isinstance(value, (str, bytes, bytearray, Mapping)):
        _fail(path, "expected an array")
    try:
        return tuple(value)
    except TypeError:
        _fail(path, "expected an array")


def _keys(
    value: Mapping[str, Any],
    path: str,
    *,
    required: set[str],
    optional: set[str] = frozenset(),
) -> None:
    missing = sorted(required - set(value))
    unknown = sorted(set(value) - required - optional)
    if missing:
        _fail(path, f"missing field(s): {', '.join(missing)}")
    if unknown:
        _fail(path, f"unknown field(s): {', '.join(unknown)}")


def _string(value: Any, path: str, *, identifier: bool = False) -> str:
    if not isinstance(value, str) or not value:
        _fail(path, "expected a non-empty string")
    if identifier and not _IDENTIFIER.fullmatch(value):
        _fail(path, "expected a portable identifier")
    return value


def _bool(value: Any, path: str) -> bool:
    if not isinstance(value, bool):
        _fail(path, "expected a boolean")
    return value


def _integer(value: Any, path: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        _fail(path, "expected an integer")
    if value < minimum:
        _fail(path, f"must be at least {minimum}")
    return value


def _number(value: Any, path: str, *, minimum: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        _fail(path, "expected a finite number")
    result = float(value)
    if not math.isfinite(result):
        _fail(path, "expected a finite number")
    if minimum is not None and result < minimum:
        _fail(path, f"must be at least {minimum:g}")
    return result


def _sha256(value: Any, path: str) -> str:
    text = _string(value, path)
    if not _SHA256.fullmatch(text):
        _fail(path, "expected 64 lowercase hexadecimal SHA-256 characters")
    return text


def _git_object_id(value: Any, path: str) -> str:
    text = _string(value, path)
    if not _GIT_OBJECT_ID.fullmatch(text):
        _fail(path, "expected a 40- or 64-character lowercase Git object ID")
    return text


def _reject_non_json(value: Any, path: str = "document") -> None:
    """Reject values for which Python's permissive encoder is ambiguous."""

    if value is None or isinstance(value, (str, bool, int)):
        return
    if isinstance(value, float):
        if not math.isfinite(value):
            _fail(path, "non-finite JSON numbers are forbidden")
        return
    if isinstance(value, Mapping):
        for key, item in value.items():
            if not isinstance(key, str):
                _fail(path, "object keys must be strings")
            _reject_non_json(item, f"{path}.{key}")
        return
    if isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _reject_non_json(item, f"{path}[{index}]")
        return
    _fail(path, f"unsupported JSON value type: {type(value).__name__}")


def _reject_host_paths(value: Any, path: str = "document") -> None:
    """Prevent local machine paths and identity from leaking into records."""

    if isinstance(value, Mapping):
        for key, item in value.items():
            if key in _PROHIBITED_HOST_KEYS:
                _fail(f"{path}.{key}", "host-specific metadata is forbidden")
            _reject_host_paths(item, f"{path}.{key}")
        return
    if isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _reject_host_paths(item, f"{path}[{index}]")
        return
    if isinstance(value, str):
        if (
            value.startswith(("/", "\\", "file://"))
            or _WINDOWS_ABSOLUTE.match(value)
        ):
            _fail(path, "absolute host paths are forbidden")
        if _EMBEDDED_HOST_PATH.search(value):
            _fail(path, "absolute host paths are forbidden")


def canonical_json_bytes(value: Any) -> bytes:
    """Return the sole canonical byte representation used by R10.

    The representation is compact UTF-8 JSON with sorted keys, no NaN or
    infinity, and no trailing newline.  It is intentionally a project-local
    scheme (named by :data:`CANONICALIZATION`), not a claim of RFC 8785
    conformance.
    """

    _reject_non_json(value)
    _reject_host_paths(value)
    try:
        text = json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        )
    except (TypeError, ValueError) as error:
        raise R10RecordError(f"document: cannot encode canonical JSON: {error}") from error
    return text.encode("utf-8")


def document_sha256(value: Any) -> str:
    """Hash the canonical representation of one evidence document."""

    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def read_canonical_json(path: os.PathLike[str] | str) -> Any:
    """Read a canonical R10 document, rejecting alternate byte encodings."""

    source = Path(path)
    raw = source.read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise DataFormatError(f"JSON is not valid UTF-8: {source}") from error
    value = loads_json_strict(text)
    if raw != canonical_json_bytes(value):
        raise DataFormatError(
            f"R10 JSON is not canonical ({CANONICALIZATION}): {source}"
        )
    return value


def write_canonical_json(path: os.PathLike[str] | str, value: Any) -> None:
    """Atomically write a canonical R10 document without a trailing newline."""

    destination = Path(path)
    payload = canonical_json_bytes(value)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, destination)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def validate_relative_path(value: Any, path: str) -> str:
    """Validate one portable POSIX relative artifact path."""

    text = _string(value, path)
    if text.startswith(("/", "\\")) or _WINDOWS_ABSOLUTE.match(text):
        _fail(path, "must be a relative path")
    if "\\" in text:
        _fail(path, "must use POSIX '/' separators")
    pure = PurePosixPath(text)
    if text in {".", ".."} or any(part in {"", ".", ".."} for part in pure.parts):
        _fail(path, "must not contain empty, current, or parent components")
    if pure.as_posix() != text:
        _fail(path, "must be in canonical POSIX form")
    return text


def _validate_binding(value: Any, path: str, *, require_bytes: bool) -> None:
    binding = _mapping(value, path)
    required = {"artifact_id", "relative_path", "sha256", "role"}
    if require_bytes:
        required.add("byte_count")
    _keys(
        binding,
        path,
        required=required,
        optional={"byte_count", "media_type", "description"} - required,
    )
    _string(binding["artifact_id"], f"{path}.artifact_id", identifier=True)
    validate_relative_path(binding["relative_path"], f"{path}.relative_path")
    _sha256(binding["sha256"], f"{path}.sha256")
    _string(binding["role"], f"{path}.role", identifier=True)
    if "byte_count" in binding:
        _integer(binding["byte_count"], f"{path}.byte_count")
    if "media_type" in binding:
        _string(binding["media_type"], f"{path}.media_type")
    if "description" in binding:
        _string(binding["description"], f"{path}.description")


def _validate_metric(value: Any, path: str, *, source: str) -> None:
    metric = _mapping(value, path)
    _keys(
        metric,
        path,
        required={"name", "value", "unit", "source", "semantics"},
        optional={"paper_location", "comparability"},
    )
    _string(metric["name"], f"{path}.name", identifier=True)
    _number(metric["value"], f"{path}.value")
    _string(metric["unit"], f"{path}.unit")
    observed_source = _string(metric["source"], f"{path}.source")
    if observed_source != source:
        _fail(f"{path}.source", f"expected {source!r}")
    _string(metric["semantics"], f"{path}.semantics")
    if source == "paper_reference":
        if "paper_location" not in metric or "comparability" not in metric:
            _fail(path, "paper metrics require paper_location and comparability")
        _string(metric["paper_location"], f"{path}.paper_location")
        if metric["comparability"] != "reference_only":
            _fail(f"{path}.comparability", "must be 'reference_only'")
    elif "paper_location" in metric or "comparability" in metric:
        _fail(path, "local metrics must not masquerade as paper reference values")


def _validate_claim(claim: Any, path: str) -> str:
    item = _mapping(claim, path)
    _keys(
        item,
        path,
        required={"level", "statement", "exact_reproduction_claim", "limitations"},
    )
    level = _string(item["level"], f"{path}.level")
    if level not in CLAIM_LEVELS:
        _fail(f"{path}.level", f"expected one of {sorted(CLAIM_LEVELS)}")
    _string(item["statement"], f"{path}.statement")
    if _bool(item["exact_reproduction_claim"], f"{path}.exact_reproduction_claim"):
        _fail(
            f"{path}.exact_reproduction_claim",
            "R10 evidence cannot claim an exact paper reproduction",
        )
    limitations = _sequence(item["limitations"], f"{path}.limitations")
    if not limitations:
        _fail(f"{path}.limitations", "at least one limitation is required")
    for index, limitation in enumerate(limitations):
        _string(limitation, f"{path}.limitations[{index}]")
    return level


def validate_manifest(value: Any) -> Mapping[str, Any]:
    """Validate an exact R10 experiment-manifest-v1 document."""

    root = _mapping(value, "manifest")
    _keys(
        root,
        "manifest",
        required={
            "schema", "schema_version", "canonicalization", "experiment_id",
            "title", "claim", "paper_reference", "protocol", "provenance",
            "inputs",
        },
    )
    if root["schema"] != MANIFEST_SCHEMA:
        _fail("manifest.schema", f"expected {MANIFEST_SCHEMA!r}")
    if _integer(root["schema_version"], "manifest.schema_version", minimum=1) != 1:
        _fail("manifest.schema_version", "expected 1")
    if root["canonicalization"] != CANONICALIZATION:
        _fail("manifest.canonicalization", f"expected {CANONICALIZATION!r}")
    _string(root["experiment_id"], "manifest.experiment_id", identifier=True)
    _string(root["title"], "manifest.title")
    _validate_claim(root["claim"], "manifest.claim")

    paper = _mapping(root["paper_reference"], "manifest.paper_reference")
    _keys(
        paper,
        "manifest.paper_reference",
        required={"citation_key", "scope", "reported_metrics"},
    )
    _string(paper["citation_key"], "manifest.paper_reference.citation_key")
    _string(paper["scope"], "manifest.paper_reference.scope")
    metrics = _sequence(
        paper["reported_metrics"], "manifest.paper_reference.reported_metrics"
    )
    for index, metric in enumerate(metrics):
        _validate_metric(
            metric,
            f"manifest.paper_reference.reported_metrics[{index}]",
            source="paper_reference",
        )

    protocol = _mapping(root["protocol"], "manifest.protocol")
    _keys(
        protocol,
        "manifest.protocol",
        required={
            "variant", "repeats", "seed_policy", "preprocessing",
            "algorithm", "metric_definitions",
        },
    )
    _string(protocol["variant"], "manifest.protocol.variant", identifier=True)
    _integer(protocol["repeats"], "manifest.protocol.repeats", minimum=1)
    _mapping(protocol["seed_policy"], "manifest.protocol.seed_policy")
    _mapping(protocol["preprocessing"], "manifest.protocol.preprocessing")
    algorithm = _mapping(protocol["algorithm"], "manifest.protocol.algorithm")
    if "warmup_runs_per_measured_repeat" not in algorithm:
        _fail(
            "manifest.protocol.algorithm.warmup_runs_per_measured_repeat",
            "is required",
        )
    _integer(
        algorithm["warmup_runs_per_measured_repeat"],
        "manifest.protocol.algorithm.warmup_runs_per_measured_repeat",
        minimum=0,
    )
    if "native_extension_artifact_id" not in algorithm:
        _fail(
            "manifest.protocol.algorithm.native_extension_artifact_id",
            "is required",
        )
    native_extension_artifact_id = _string(
        algorithm["native_extension_artifact_id"],
        "manifest.protocol.algorithm.native_extension_artifact_id",
        identifier=True,
    )
    definitions = _sequence(
        protocol["metric_definitions"], "manifest.protocol.metric_definitions"
    )
    if not definitions:
        _fail("manifest.protocol.metric_definitions", "must not be empty")
    metric_names: set[str] = set()
    for index, raw_definition in enumerate(definitions):
        path = f"manifest.protocol.metric_definitions[{index}]"
        definition = _mapping(raw_definition, path)
        _keys(
            definition,
            path,
            required={"name", "unit", "direction", "definition", "source"},
        )
        name = _string(definition["name"], f"{path}.name", identifier=True)
        if name in metric_names:
            _fail(f"{path}.name", "duplicate metric name")
        metric_names.add(name)
        _string(definition["unit"], f"{path}.unit")
        if definition["direction"] not in {"lower_is_better", "higher_is_better", "descriptive"}:
            _fail(f"{path}.direction", "unsupported direction")
        _string(definition["definition"], f"{path}.definition")
        if definition["source"] != "local_measurement":
            _fail(f"{path}.source", "must be 'local_measurement'")

    provenance = _mapping(root["provenance"], "manifest.provenance")
    _keys(
        provenance,
        "manifest.provenance",
        required={
            "repository_commit_sha", "dirty_worktree", "implementation_profile",
            "resources", "deviations",
        },
    )
    _git_object_id(
        provenance["repository_commit_sha"],
        "manifest.provenance.repository_commit_sha",
    )
    _bool(provenance["dirty_worktree"], "manifest.provenance.dirty_worktree")
    _string(provenance["implementation_profile"], "manifest.provenance.implementation_profile")
    resources = _sequence(provenance["resources"], "manifest.provenance.resources")
    for index, resource in enumerate(resources):
        _validate_binding(resource, f"manifest.provenance.resources[{index}]", require_bytes=False)
    native_extensions = tuple(
        resource
        for resource in resources
        if resource.get("artifact_id") == native_extension_artifact_id
    )
    if len(native_extensions) != 1:
        _fail(
            "manifest.provenance.resources",
            "must contain exactly one declared native extension artifact",
        )
    native_extension = native_extensions[0]
    _validate_binding(
        native_extension,
        "manifest.provenance.resources[native_extension]",
        require_bytes=True,
    )
    if native_extension.get("role") != "native_extension":
        _fail(
            "manifest.provenance.resources[native_extension].role",
            "must be 'native_extension'",
        )
    if any(
        resource.get("role") == "native_extension" and
        resource.get("artifact_id") != native_extension_artifact_id
        for resource in resources
    ):
        _fail(
            "manifest.provenance.resources",
            "contains an undeclared native extension artifact",
        )
    deviations = _sequence(provenance["deviations"], "manifest.provenance.deviations")
    if not deviations:
        _fail("manifest.provenance.deviations", "at least one known deviation is required")
    for index, deviation in enumerate(deviations):
        _string(deviation, f"manifest.provenance.deviations[{index}]")

    inputs = _sequence(root["inputs"], "manifest.inputs")
    if not inputs:
        _fail("manifest.inputs", "at least one SHA-256-bound input is required")
    artifact_ids: set[str] = set()
    for index, binding in enumerate(inputs):
        _validate_binding(binding, f"manifest.inputs[{index}]", require_bytes=True)
        artifact_id = binding["artifact_id"]
        if artifact_id in artifact_ids:
            _fail(f"manifest.inputs[{index}].artifact_id", "duplicate artifact identifier")
        artifact_ids.add(artifact_id)

    canonical_json_bytes(root)
    return root


def _validate_hardware(value: Any, path: str) -> None:
    hardware = _mapping(value, path)
    _keys(
        hardware,
        path,
        required={"cpu_model", "logical_cpu_count", "memory_bytes", "os", "machine"},
        optional={"accelerator"},
    )
    _string(hardware["cpu_model"], f"{path}.cpu_model")
    _integer(hardware["logical_cpu_count"], f"{path}.logical_cpu_count", minimum=1)
    _integer(hardware["memory_bytes"], f"{path}.memory_bytes", minimum=1)
    _string(hardware["os"], f"{path}.os")
    _string(hardware["machine"], f"{path}.machine")
    if "accelerator" in hardware:
        _string(hardware["accelerator"], f"{path}.accelerator")


def _validate_build(value: Any, path: str) -> None:
    build = _mapping(value, path)
    _keys(
        build,
        path,
        required={"build_type", "compiler", "compiler_version", "commit_sha", "dirty_worktree", "flags"},
    )
    _string(build["build_type"], f"{path}.build_type")
    _string(build["compiler"], f"{path}.compiler")
    _string(build["compiler_version"], f"{path}.compiler_version")
    _git_object_id(build["commit_sha"], f"{path}.commit_sha")
    _bool(build["dirty_worktree"], f"{path}.dirty_worktree")
    flags = _sequence(build["flags"], f"{path}.flags")
    for index, flag in enumerate(flags):
        _string(flag, f"{path}.flags[{index}]")


def _validate_timing(value: Any, path: str) -> None:
    timing = _mapping(value, path)
    _keys(
        timing,
        path,
        required={"clock", "scope", "warmup_runs", "phases", "total_seconds"},
    )
    if timing["clock"] != "monotonic_performance_counter":
        _fail(f"{path}.clock", "expected 'monotonic_performance_counter'")
    if timing["scope"] != OUTER_TIMING_SCOPE:
        _fail(f"{path}.scope", f"expected {OUTER_TIMING_SCOPE!r}")
    _integer(timing["warmup_runs"], f"{path}.warmup_runs")
    phases = _sequence(timing["phases"], f"{path}.phases")
    names: set[str] = set()
    phase_total = 0.0
    for index, raw_phase in enumerate(phases):
        phase_path = f"{path}.phases[{index}]"
        phase = _mapping(raw_phase, phase_path)
        _keys(phase, phase_path, required={"name", "seconds"})
        name = _string(phase["name"], f"{phase_path}.name", identifier=True)
        if name in names:
            _fail(f"{phase_path}.name", "duplicate timing phase")
        names.add(name)
        phase_total += _number(
            phase["seconds"],
            f"{phase_path}.seconds",
            minimum=0.0,
        )
    total = _number(
        timing["total_seconds"],
        f"{path}.total_seconds",
        minimum=0.0,
    )
    tolerance = max(1.0e-12, total * 1.0e-9)
    if phase_total > total + tolerance:
        _fail(f"{path}.phases", "phase durations exceed total duration")


def validate_run_record(
    value: Any,
    *,
    manifest: Mapping[str, Any] | None = None,
) -> Mapping[str, Any]:
    """Validate an exact R10 run-record-v1 and optionally bind its manifest."""

    root = _mapping(value, "run_record")
    _keys(
        root,
        "run_record",
        required={
            "schema", "schema_version", "canonicalization", "experiment_id",
            "run_id", "repeat_index", "seed", "manifest_sha256", "claim_level",
            "status", "execution", "paper_reference_metrics", "local_metrics",
            "artifacts", "notes",
        },
    )
    if root["schema"] != RUN_RECORD_SCHEMA:
        _fail("run_record.schema", f"expected {RUN_RECORD_SCHEMA!r}")
    if _integer(root["schema_version"], "run_record.schema_version", minimum=1) != 1:
        _fail("run_record.schema_version", "expected 1")
    if root["canonicalization"] != CANONICALIZATION:
        _fail("run_record.canonicalization", f"expected {CANONICALIZATION!r}")
    _string(root["experiment_id"], "run_record.experiment_id", identifier=True)
    _string(root["run_id"], "run_record.run_id", identifier=True)
    _integer(root["repeat_index"], "run_record.repeat_index")
    _integer(root["seed"], "run_record.seed")
    _sha256(root["manifest_sha256"], "run_record.manifest_sha256")
    if root["claim_level"] not in CLAIM_LEVELS:
        _fail("run_record.claim_level", f"expected one of {sorted(CLAIM_LEVELS)}")
    if root["status"] not in RUN_STATUSES:
        _fail("run_record.status", f"expected one of {sorted(RUN_STATUSES)}")

    execution = _mapping(root["execution"], "run_record.execution")
    _keys(execution, "run_record.execution", required={"hardware", "build", "timing"})
    _validate_hardware(execution["hardware"], "run_record.execution.hardware")
    _validate_build(execution["build"], "run_record.execution.build")
    _validate_timing(execution["timing"], "run_record.execution.timing")

    paper_metrics = _sequence(root["paper_reference_metrics"], "run_record.paper_reference_metrics")
    for index, metric in enumerate(paper_metrics):
        _validate_metric(metric, f"run_record.paper_reference_metrics[{index}]", source="paper_reference")
    local_metrics = _sequence(root["local_metrics"], "run_record.local_metrics")
    names: set[str] = set()
    for index, metric in enumerate(local_metrics):
        path = f"run_record.local_metrics[{index}]"
        _validate_metric(metric, path, source="local_measurement")
        if metric["name"] in names:
            _fail(f"{path}.name", "duplicate local metric name")
        names.add(metric["name"])

    artifacts = _sequence(root["artifacts"], "run_record.artifacts")
    for index, artifact in enumerate(artifacts):
        _validate_binding(artifact, f"run_record.artifacts[{index}]", require_bytes=True)
    notes = _sequence(root["notes"], "run_record.notes")
    for index, note in enumerate(notes):
        _string(note, f"run_record.notes[{index}]")

    if manifest is not None:
        checked_manifest = validate_manifest(manifest)
        expected_hash = document_sha256(checked_manifest)
        if root["manifest_sha256"] != expected_hash:
            _fail("run_record.manifest_sha256", "does not bind the supplied manifest")
        if root["experiment_id"] != checked_manifest["experiment_id"]:
            _fail("run_record.experiment_id", "does not match the supplied manifest")
        if root["claim_level"] != checked_manifest["claim"]["level"]:
            _fail("run_record.claim_level", "does not match the supplied manifest")
        build = root["execution"]["build"]
        provenance = checked_manifest["provenance"]
        if build["commit_sha"] != provenance["repository_commit_sha"]:
            _fail(
                "run_record.execution.build.commit_sha",
                "does not match manifest provenance",
            )
        if build["dirty_worktree"] != provenance["dirty_worktree"]:
            _fail(
                "run_record.execution.build.dirty_worktree",
                "does not match manifest provenance",
            )
        algorithm = checked_manifest["protocol"]["algorithm"]
        declared_warmups = algorithm.get("warmup_runs_per_measured_repeat")
        if declared_warmups is None:
            _fail(
                "manifest.protocol.algorithm.warmup_runs_per_measured_repeat",
                "is required to bind run-record timing",
            )
        _integer(
            declared_warmups,
            "manifest.protocol.algorithm.warmup_runs_per_measured_repeat",
            minimum=0,
        )
        if root["execution"]["timing"]["warmup_runs"] != declared_warmups:
            _fail(
                "run_record.execution.timing.warmup_runs",
                "does not match manifest protocol",
            )
        if list(root["paper_reference_metrics"]) != list(
            checked_manifest["paper_reference"]["reported_metrics"]
        ):
            _fail(
                "run_record.paper_reference_metrics",
                "must copy the manifest paper-reference metrics exactly",
            )
        definitions = {
            item["name"]: item
            for item in checked_manifest["protocol"]["metric_definitions"]
        }
        for index, metric in enumerate(root["local_metrics"]):
            definition = definitions.get(metric["name"])
            if definition is None:
                _fail(
                    f"run_record.local_metrics[{index}].name",
                    "metric is not declared by the manifest",
                )
            if metric["unit"] != definition["unit"]:
                _fail(
                    f"run_record.local_metrics[{index}].unit",
                    "does not match the manifest metric definition",
                )
            if metric["semantics"] != definition["definition"]:
                _fail(
                    f"run_record.local_metrics[{index}].semantics",
                    "does not match the manifest metric definition",
                )
        native_id = algorithm["native_extension_artifact_id"]
        manifest_native = next(
            resource
            for resource in provenance["resources"]
            if resource["artifact_id"] == native_id
        )
        record_native = tuple(
            artifact
            for artifact in root["artifacts"]
            if artifact.get("artifact_id") == native_id
        )
        if len(record_native) != 1 or record_native[0] != manifest_native:
            _fail(
                "run_record.artifacts",
                "must copy the manifest native extension binding exactly",
            )
        if any(
            artifact.get("role") == "native_extension" and
            artifact.get("artifact_id") != native_id
            for artifact in root["artifacts"]
        ):
            _fail(
                "run_record.artifacts",
                "contains an undeclared native extension artifact",
            )

    canonical_json_bytes(root)
    return root


def _memory_bytes() -> int:
    try:
        return int(os.sysconf("SC_PAGE_SIZE")) * int(os.sysconf("SC_PHYS_PAGES"))
    except (AttributeError, OSError, ValueError):
        return 1


def collect_hardware_metadata(*, accelerator: str | None = None) -> dict[str, Any]:
    """Collect reproducibility metadata without hostname, username, or paths."""

    cpu_model = platform.processor().strip()
    if not cpu_model:
        try:
            for line in Path("/proc/cpuinfo").read_text(encoding="utf-8").splitlines():
                if line.lower().startswith("model name") and ":" in line:
                    cpu_model = line.split(":", 1)[1].strip()
                    break
        except OSError:
            pass
    result: dict[str, Any] = {
        "cpu_model": cpu_model or "unreported",
        "logical_cpu_count": os.cpu_count() or 1,
        "memory_bytes": max(1, _memory_bytes()),
        "os": f"{platform.system()} {platform.release()}".strip(),
        "machine": platform.machine() or "unreported",
    }
    if accelerator:
        result["accelerator"] = accelerator
    _validate_hardware(result, "hardware")
    return result


def bind_artifact(
    root: os.PathLike[str] | str,
    relative_path: str,
    *,
    artifact_id: str,
    role: str,
    media_type: str | None = None,
    description: str | None = None,
) -> dict[str, Any]:
    """Hash one ordinary file below *root* without serializing *root* itself."""

    relative = validate_relative_path(relative_path, "relative_path")
    _string(artifact_id, "artifact_id", identifier=True)
    _string(role, "role", identifier=True)
    anchor = Path(root).resolve(strict=True)
    candidate = anchor
    for component in PurePosixPath(relative).parts:
        candidate /= component
        metadata = os.lstat(candidate)
        if stat.S_ISLNK(metadata.st_mode):
            _fail("relative_path", "symbolic-link artifacts are forbidden")
    resolved = candidate.resolve(strict=True)
    try:
        resolved.relative_to(anchor)
    except ValueError:
        _fail("relative_path", "artifact escapes the selected root")
    if not resolved.is_file():
        _fail("relative_path", "artifact must be an ordinary file")
    digest = hashlib.sha256()
    byte_count = 0
    with resolved.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
            byte_count += len(chunk)
    result: dict[str, Any] = {
        "artifact_id": artifact_id,
        "relative_path": relative,
        "sha256": digest.hexdigest(),
        "byte_count": byte_count,
        "role": role,
    }
    if media_type is not None:
        result["media_type"] = _string(media_type, "media_type")
    if description is not None:
        result["description"] = _string(description, "description")
    _validate_binding(result, "artifact", require_bytes=True)
    return result


@dataclass
class _PhaseMeasurement:
    name: str
    seconds: float


class RunRecorder:
    """Build a SHA-bound run record while measuring named phases.

    Parameters which describe the build are explicit because guessing compiler
    flags from a Python process would be misleading.  Hardware collection is
    automatic but can be overridden with a validated mapping for controlled
    or deterministic tests.
    """

    def __init__(
        self,
        manifest: Mapping[str, Any],
        *,
        run_id: str,
        repeat_index: int,
        seed: int,
        build: Mapping[str, Any],
        hardware: Mapping[str, Any] | None = None,
        timing_scope: str = OUTER_TIMING_SCOPE,
        warmup_runs: int = 0,
    ) -> None:
        self._manifest = dict(validate_manifest(manifest))
        self._run_id = _string(run_id, "run_id", identifier=True)
        self._repeat_index = _integer(repeat_index, "repeat_index")
        self._seed = _integer(seed, "seed")
        self._build = dict(_mapping(build, "build"))
        _validate_build(self._build, "build")
        self._hardware = dict(hardware or collect_hardware_metadata())
        _validate_hardware(self._hardware, "hardware")
        self._timing_scope = _string(timing_scope, "timing_scope")
        self._warmup_runs = _integer(warmup_runs, "warmup_runs")
        self._phases: list[_PhaseMeasurement] = []
        self._active_phase: str | None = None
        self._metrics: list[dict[str, Any]] = []
        native_id = self._manifest["protocol"]["algorithm"][
            "native_extension_artifact_id"
        ]
        native_binding = next(
            resource
            for resource in self._manifest["provenance"]["resources"]
            if resource["artifact_id"] == native_id
        )
        self._artifacts: list[dict[str, Any]] = [dict(native_binding)]
        self._notes: list[str] = []
        self._start = time.perf_counter()
        self._finished = False

    def _ensure_open(self) -> None:
        if self._finished:
            raise R10RecordError("run: recorder has already been finalized")

    @contextmanager
    def phase(self, name: str) -> Iterator[None]:
        """Measure one non-nested, uniquely named phase."""

        self._ensure_open()
        checked = _string(name, "phase.name", identifier=True)
        if self._active_phase is not None:
            raise R10RecordError("phase: nested timing phases are forbidden")
        if any(item.name == checked for item in self._phases):
            raise R10RecordError(f"phase.name: duplicate phase {checked!r}")
        self._active_phase = checked
        start = time.perf_counter()
        try:
            yield
        finally:
            elapsed = max(0.0, time.perf_counter() - start)
            self._phases.append(_PhaseMeasurement(checked, elapsed))
            self._active_phase = None

    def add_local_metric(self, name: str, value: float, *, semantics: str) -> None:
        self._ensure_open()
        definitions = {
            item["name"]: item
            for item in self._manifest["protocol"]["metric_definitions"]
        }
        checked_name = _string(name, "metric.name", identifier=True)
        if checked_name not in definitions:
            _fail("metric.name", "metric is not declared by the manifest")
        if any(item["name"] == checked_name for item in self._metrics):
            _fail("metric.name", "duplicate local metric")
        self._metrics.append(
            {
                "name": checked_name,
                "value": _number(value, "metric.value"),
                "unit": definitions[checked_name]["unit"],
                "source": "local_measurement",
                "semantics": _string(semantics, "metric.semantics"),
            }
        )

    def add_artifact(self, binding: Mapping[str, Any]) -> None:
        self._ensure_open()
        copied = dict(binding)
        _validate_binding(copied, "artifact", require_bytes=True)
        if any(item["artifact_id"] == copied["artifact_id"] for item in self._artifacts):
            _fail("artifact.artifact_id", "duplicate output artifact")
        self._artifacts.append(copied)

    def add_note(self, note: str) -> None:
        self._ensure_open()
        self._notes.append(_string(note, "note"))

    def finish(self, *, status: str) -> dict[str, Any]:
        self._ensure_open()
        if self._active_phase is not None:
            raise R10RecordError("run: cannot finish while a timing phase is active")
        if status not in RUN_STATUSES:
            _fail("status", f"expected one of {sorted(RUN_STATUSES)}")
        record: dict[str, Any] = {
            "schema": RUN_RECORD_SCHEMA,
            "schema_version": SCHEMA_VERSION,
            "canonicalization": CANONICALIZATION,
            "experiment_id": self._manifest["experiment_id"],
            "run_id": self._run_id,
            "repeat_index": self._repeat_index,
            "seed": self._seed,
            "manifest_sha256": document_sha256(self._manifest),
            "claim_level": self._manifest["claim"]["level"],
            "status": status,
            "execution": {
                "hardware": self._hardware,
                "build": self._build,
                "timing": {
                    "clock": "monotonic_performance_counter",
                    "scope": self._timing_scope,
                    "warmup_runs": self._warmup_runs,
                    "phases": [
                        {"name": phase.name, "seconds": phase.seconds}
                        for phase in self._phases
                    ],
                    # Finalized below, after all other record fields have been
                    # assembled, so the declared outer scope reaches the
                    # run-record finalization boundary.
                    "total_seconds": 0.0,
                },
            },
            "paper_reference_metrics": list(
                self._manifest["paper_reference"]["reported_metrics"]
            ),
            "local_metrics": list(self._metrics),
            "artifacts": list(self._artifacts),
            "notes": list(self._notes),
        }
        record["execution"]["timing"]["total_seconds"] = max(
            0.0, time.perf_counter() - self._start
        )
        validate_run_record(record, manifest=self._manifest)
        self._finished = True
        return record


def _quantile(sorted_values: Sequence[float], probability: float) -> float:
    """Linearly interpolated quantile (Hyndman-Fan type 7)."""

    if not sorted_values:
        raise R10RecordError("aggregate: cannot summarize an empty sample")
    position = (len(sorted_values) - 1) * probability
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return float(sorted_values[lower])
    weight = position - lower
    return float(sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight)


def _summary(values: Sequence[float]) -> dict[str, Any]:
    ordered = sorted(float(value) for value in values)
    q1 = _quantile(ordered, 0.25)
    median = _quantile(ordered, 0.5)
    q3 = _quantile(ordered, 0.75)
    return {
        "sample_count": len(ordered),
        "median": median,
        "q1": q1,
        "q3": q3,
        "iqr": q3 - q1,
        "minimum": ordered[0],
        "maximum": ordered[-1],
    }


def aggregate_run_records(
    manifest: Mapping[str, Any],
    run_records: Iterable[Mapping[str, Any]],
) -> dict[str, Any]:
    """Aggregate comparable local values; paper values remain reference-only."""

    checked_manifest = validate_manifest(manifest)
    records = [dict(validate_run_record(record, manifest=checked_manifest)) for record in run_records]
    if not records:
        _fail("run_records", "at least one run record is required")
    run_ids = [record["run_id"] for record in records]
    if len(set(run_ids)) != len(run_ids):
        _fail("run_records", "run_id values must be unique")
    repeat_indices = [record["repeat_index"] for record in records]
    if len(set(repeat_indices)) != len(repeat_indices):
        _fail("run_records", "repeat_index values must be unique")
    expected_repeats = checked_manifest["protocol"]["repeats"]
    if len(records) != expected_repeats:
        _fail(
            "run_records",
            f"expected exactly {expected_repeats} records declared by the manifest",
        )
    if sorted(repeat_indices) != list(range(expected_repeats)):
        _fail(
            "run_records",
            "repeat_index values must cover the manifest range from zero",
        )

    metric_values: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    timing_values: dict[str, list[float]] = defaultdict(list)
    for record in records:
        if record["status"] == "succeeded":
            for metric in record["local_metrics"]:
                key = (metric["name"], metric["unit"], metric["semantics"])
                metric_values[key].append(float(metric["value"]))
            timing = record["execution"]["timing"]
            timing_values["total"].append(float(timing["total_seconds"]))
            for phase in timing["phases"]:
                timing_values[phase["name"]].append(float(phase["seconds"]))

    aggregate = {
        "schema": AGGREGATE_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "canonicalization": CANONICALIZATION,
        "experiment_id": checked_manifest["experiment_id"],
        "manifest_sha256": document_sha256(checked_manifest),
        "claim_level": checked_manifest["claim"]["level"],
        "exact_reproduction_claim": False,
        "comparison_semantics": {
            "paper_values": "reference_only_not_pooled",
            "local_values": "same-manifest_repeated_measurements",
            "quantile_method": "linear_type_7",
        },
        "repeat_count": len(records),
        "successful_repeat_count": sum(record["status"] == "succeeded" for record in records),
        "status_counts": dict(sorted(Counter(record["status"] for record in records).items())),
        "run_record_sha256s": sorted(document_sha256(record) for record in records),
        "paper_reference_metrics": list(
            checked_manifest["paper_reference"]["reported_metrics"]
        ),
        "local_metric_aggregates": [
            {
                "name": name,
                "unit": unit,
                "source": "local_measurement",
                "semantics": semantics,
                **_summary(values),
            }
            for (name, unit, semantics), values in sorted(metric_values.items())
        ],
        "timing_aggregates": [
            {"name": name, "unit": "seconds", **_summary(values)}
            for name, values in sorted(timing_values.items())
        ],
    }
    validate_aggregate(aggregate, manifest=checked_manifest)
    return aggregate


def _validate_summary(value: Any, path: str) -> None:
    summary = _mapping(value, path)
    required = {"sample_count", "median", "q1", "q3", "iqr", "minimum", "maximum"}
    if path.endswith("local_metric_aggregates[]"):
        required |= {"name", "unit", "source", "semantics"}
    elif path.endswith("timing_aggregates[]"):
        required |= {"name", "unit"}
    _keys(summary, path, required=required)
    count = _integer(summary["sample_count"], f"{path}.sample_count", minimum=1)
    del count
    minimum = _number(summary["minimum"], f"{path}.minimum")
    q1 = _number(summary["q1"], f"{path}.q1")
    median = _number(summary["median"], f"{path}.median")
    q3 = _number(summary["q3"], f"{path}.q3")
    maximum = _number(summary["maximum"], f"{path}.maximum")
    iqr = _number(summary["iqr"], f"{path}.iqr", minimum=0.0)
    if not (minimum <= q1 <= median <= q3 <= maximum):
        _fail(path, "summary order must be minimum <= q1 <= median <= q3 <= maximum")
    if not math.isclose(iqr, q3 - q1, rel_tol=1.0e-12, abs_tol=1.0e-15):
        _fail(f"{path}.iqr", "must equal q3 - q1")


def validate_aggregate(
    value: Any,
    *,
    manifest: Mapping[str, Any] | None = None,
) -> Mapping[str, Any]:
    """Validate an exact R10 aggregate-v1 document."""

    root = _mapping(value, "aggregate")
    _keys(
        root,
        "aggregate",
        required={
            "schema", "schema_version", "canonicalization", "experiment_id",
            "manifest_sha256", "claim_level", "exact_reproduction_claim",
            "comparison_semantics", "repeat_count", "successful_repeat_count",
            "status_counts", "run_record_sha256s", "paper_reference_metrics",
            "local_metric_aggregates", "timing_aggregates",
        },
    )
    if root["schema"] != AGGREGATE_SCHEMA:
        _fail("aggregate.schema", f"expected {AGGREGATE_SCHEMA!r}")
    if _integer(root["schema_version"], "aggregate.schema_version", minimum=1) != 1:
        _fail("aggregate.schema_version", "expected 1")
    if root["canonicalization"] != CANONICALIZATION:
        _fail("aggregate.canonicalization", f"expected {CANONICALIZATION!r}")
    _string(root["experiment_id"], "aggregate.experiment_id", identifier=True)
    _sha256(root["manifest_sha256"], "aggregate.manifest_sha256")
    if root["claim_level"] not in CLAIM_LEVELS:
        _fail("aggregate.claim_level", f"expected one of {sorted(CLAIM_LEVELS)}")
    if _bool(root["exact_reproduction_claim"], "aggregate.exact_reproduction_claim"):
        _fail("aggregate.exact_reproduction_claim", "must be false")
    semantics = _mapping(root["comparison_semantics"], "aggregate.comparison_semantics")
    _keys(
        semantics,
        "aggregate.comparison_semantics",
        required={"paper_values", "local_values", "quantile_method"},
    )
    if semantics != {
        "paper_values": "reference_only_not_pooled",
        "local_values": "same-manifest_repeated_measurements",
        "quantile_method": "linear_type_7",
    }:
        _fail("aggregate.comparison_semantics", "unsupported comparison semantics")
    repeat_count = _integer(root["repeat_count"], "aggregate.repeat_count", minimum=1)
    successful_count = _integer(
        root["successful_repeat_count"], "aggregate.successful_repeat_count"
    )
    if successful_count > repeat_count:
        _fail("aggregate.successful_repeat_count", "cannot exceed repeat_count")
    status_counts = _mapping(root["status_counts"], "aggregate.status_counts")
    if any(status not in RUN_STATUSES for status in status_counts):
        _fail("aggregate.status_counts", "contains an unsupported status")
    for status, count in status_counts.items():
        _integer(count, f"aggregate.status_counts.{status}")
    if sum(status_counts.values()) != repeat_count:
        _fail("aggregate.status_counts", "counts must sum to repeat_count")
    hashes = _sequence(root["run_record_sha256s"], "aggregate.run_record_sha256s")
    if len(hashes) != repeat_count:
        _fail("aggregate.run_record_sha256s", "length must equal repeat_count")
    if list(hashes) != sorted(hashes) or len(set(hashes)) != len(hashes):
        _fail("aggregate.run_record_sha256s", "must be sorted and unique")
    for index, digest in enumerate(hashes):
        _sha256(digest, f"aggregate.run_record_sha256s[{index}]")
    for index, metric in enumerate(
        _sequence(root["paper_reference_metrics"], "aggregate.paper_reference_metrics")
    ):
        _validate_metric(metric, f"aggregate.paper_reference_metrics[{index}]", source="paper_reference")
    for index, summary in enumerate(
        _sequence(root["local_metric_aggregates"], "aggregate.local_metric_aggregates")
    ):
        _validate_summary(summary, "aggregate.local_metric_aggregates[]")
        if summary["source"] != "local_measurement":
            _fail(f"aggregate.local_metric_aggregates[{index}].source", "must be local_measurement")
    for summary in _sequence(root["timing_aggregates"], "aggregate.timing_aggregates"):
        _validate_summary(summary, "aggregate.timing_aggregates[]")
        if summary["unit"] != "seconds":
            _fail("aggregate.timing_aggregates[].unit", "must be seconds")

    if manifest is not None:
        checked = validate_manifest(manifest)
        if root["manifest_sha256"] != document_sha256(checked):
            _fail("aggregate.manifest_sha256", "does not bind the supplied manifest")
        if root["experiment_id"] != checked["experiment_id"]:
            _fail("aggregate.experiment_id", "does not match the supplied manifest")
        if root["claim_level"] != checked["claim"]["level"]:
            _fail("aggregate.claim_level", "does not match the supplied manifest")
        if list(root["paper_reference_metrics"]) != list(
            checked["paper_reference"]["reported_metrics"]
        ):
            _fail("aggregate.paper_reference_metrics", "must copy manifest values exactly")

    canonical_json_bytes(root)
    return root


__all__ = [
    "AGGREGATE_SCHEMA",
    "CANONICALIZATION",
    "CLAIM_LEVELS",
    "MANIFEST_SCHEMA",
    "METRIC_SOURCES",
    "OUTER_TIMING_SCOPE",
    "R10RecordError",
    "RUN_RECORD_SCHEMA",
    "RUN_STATUSES",
    "RunRecorder",
    "aggregate_run_records",
    "bind_artifact",
    "canonical_json_bytes",
    "collect_hardware_metadata",
    "document_sha256",
    "read_canonical_json",
    "validate_aggregate",
    "validate_manifest",
    "validate_relative_path",
    "validate_run_record",
    "write_canonical_json",
]
