"""Fail-closed provenance contracts for paper-simulation experiments."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
from typing import Any, Iterable, Mapping


class EvidenceClass(str, Enum):
    """Maximum claim supported by an experiment's available inputs."""

    EXACT_CACHE = "EXACT_CACHE"
    RECONSTRUCTED_PARAMETER = "RECONSTRUCTED_PARAMETER"
    RECONSTRUCTED_ASSET = "RECONSTRUCTED_ASSET"
    BLOCKED_MISSING_ASSET = "BLOCKED_MISSING_ASSET"


class FreshnessStatus(str, Enum):
    CURRENT = "CURRENT"
    STALE_RESULT = "STALE_RESULT"


@dataclass(frozen=True)
class FreshnessAssessment:
    status: FreshnessStatus
    mismatches: tuple[str, ...]


def classify_evidence(
    *,
    complete_cache: bool,
    parameter_override: bool = False,
    reconstructed_surface: bool = False,
) -> EvidenceClass:
    """Return the strongest defensible claim without silently upgrading it."""

    if not complete_cache:
        if parameter_override or reconstructed_surface:
            raise ValueError(
                "missing assets cannot be upgraded by parameter or surface overrides"
            )
        return EvidenceClass.BLOCKED_MISSING_ASSET
    if reconstructed_surface:
        return EvidenceClass.RECONSTRUCTED_ASSET
    if parameter_override:
        return EvidenceClass.RECONSTRUCTED_PARAMETER
    return EvidenceClass.EXACT_CACHE


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_json(value: Any) -> bytes:
    return json.dumps(
        value,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


def _relative_or_absolute(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return os.fspath(path)


def _file_records(
    repository_root: Path,
    paths: Iterable[os.PathLike[str] | str],
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    expanded: list[Path] = []

    def is_source_file(path: Path) -> bool:
        return (
            "__pycache__" not in path.parts
            and path.suffix not in {".pyc", ".pyo"}
        )

    for raw in paths:
        path = Path(raw).expanduser().resolve()
        if path.is_dir():
            expanded.extend(
                item
                for item in path.rglob("*")
                if item.is_file() and is_source_file(item)
            )
        else:
            expanded.append(path)
    for path in sorted(set(expanded), key=os.fspath):
        if not path.is_file():
            raise FileNotFoundError(path)
        records.append({
            "path": _relative_or_absolute(path, repository_root),
            "size": path.stat().st_size,
            "sha256": _sha256_file(path),
        })
    return records


def _git_value(root: Path, *arguments: str) -> str | None:
    try:
        completed = subprocess.run(
            ("git", "-C", os.fspath(root), *arguments),
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    value = completed.stdout.strip()
    return value or None


def capture_execution_fingerprint(
    repository_root: os.PathLike[str] | str,
    *,
    source_paths: Iterable[os.PathLike[str] | str],
    executable_paths: Iterable[os.PathLike[str] | str] = (),
    configuration: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Bind one result to exact source, executables, configuration and runtime."""

    root = Path(repository_root).expanduser().resolve()
    sources = _file_records(root, source_paths)
    executables = _file_records(root, executable_paths)
    config = dict(configuration or {})
    return {
        "schema": "make-it-stand-execution-fingerprint-v1",
        "schema_version": 1,
        "git_head": _git_value(root, "rev-parse", "HEAD"),
        "git_branch": _git_value(root, "branch", "--show-current"),
        "source_tree_sha256": _sha256_bytes(_canonical_json(sources)),
        "sources": sources,
        "executables_sha256": _sha256_bytes(_canonical_json(executables)),
        "executables": executables,
        "configuration_sha256": _sha256_bytes(_canonical_json(config)),
        "configuration": config,
        "runtime": {
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
    }


def assess_result_freshness(
    recorded: Mapping[str, Any],
    current: Mapping[str, Any],
) -> FreshnessAssessment:
    """Compare result provenance against the exact current execution context."""

    fields = (
        "git_head",
        "source_tree_sha256",
        "executables_sha256",
        "configuration_sha256",
    )
    mismatches = [field for field in fields if recorded.get(field) != current.get(field)]
    if recorded.get("executables") != current.get("executables"):
        mismatches.append("executables")
    return FreshnessAssessment(
        status=(
            FreshnessStatus.STALE_RESULT
            if mismatches
            else FreshnessStatus.CURRENT
        ),
        mismatches=tuple(dict.fromkeys(mismatches)),
    )


def _recorded_paths(
    records: Any,
    repository_root: Path,
) -> list[Path]:
    if not isinstance(records, list):
        raise ValueError("fingerprint file records must be a list")
    result: list[Path] = []
    for record in records:
        if not isinstance(record, Mapping) or not isinstance(record.get("path"), str):
            raise ValueError("fingerprint file record must contain a path")
        path = Path(record["path"])
        result.append(path if path.is_absolute() else repository_root / path)
    return result


def audit_manifest_freshness(
    manifest: Mapping[str, Any],
    repository_root: os.PathLike[str] | str,
) -> dict[str, Any]:
    """Return a manifest copy with a freshly evaluated fail-closed status."""

    result = dict(manifest)
    recorded = manifest.get("execution_fingerprint")
    if not isinstance(recorded, Mapping):
        result["freshness"] = {
            "status": FreshnessStatus.STALE_RESULT.value,
            "mismatches": ["missing_execution_fingerprint"],
        }
        return result
    root = Path(repository_root).expanduser().resolve()
    try:
        current = capture_execution_fingerprint(
            root,
            source_paths=_recorded_paths(recorded.get("sources"), root),
            executable_paths=_recorded_paths(recorded.get("executables"), root),
            configuration=recorded.get("configuration")
            if isinstance(recorded.get("configuration"), Mapping)
            else {},
        )
    except (OSError, ValueError) as error:
        result["freshness"] = {
            "status": FreshnessStatus.STALE_RESULT.value,
            "mismatches": ["fingerprinted_file_missing_or_invalid"],
            "error": str(error),
        }
        return result
    assessment = assess_result_freshness(recorded, current)
    result["freshness"] = {
        "status": assessment.status.value,
        "mismatches": list(assessment.mismatches),
    }
    return result


__all__ = [
    "EvidenceClass",
    "FreshnessAssessment",
    "FreshnessStatus",
    "assess_result_freshness",
    "audit_manifest_freshness",
    "capture_execution_fingerprint",
    "classify_evidence",
]
