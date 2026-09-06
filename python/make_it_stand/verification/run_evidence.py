"""Observed file/runtime provenance, deliberately *not* build attestation.

Use ``begin_run`` before executing a run, ``finish_run`` after it, and
``check_freshness`` on the resulting immutable PaperRunEvidenceV2. CURRENT
means the observed inputs and output files still match, not scientific success.
Failed attempts must be supplied in their original event order.

Source policy v1 covers this repository's code/build/test trees and experiment
scripts. It is a declared policy, not compiler-discovered dependency closure.
Dependency coverage is bounded to explicitly declared local roots. Git HEAD
and all three binary diffs are contextual observations: effective file content
defines freshness, so committing only a report cannot create a source hash cycle.

Runtime hashes measure files on disk, not already-loaded memory. Snapshots and
their integrity digests are unauthenticated and cannot prove author identity,
build-to-source correspondence, or changes reverted between observations.
"""

from __future__ import annotations

import base64
from dataclasses import dataclass
import hashlib
from importlib.machinery import ExtensionFileLoader
import json
import math
import os
from pathlib import Path
import platform
import stat
import subprocess
import sys
from types import ModuleType
from typing import Any


SOURCE_POLICY = "make-it-stand-effective-source-v1"
_TREES = {"src", "include", "bindings", "cmake", "python", "scripts", "tests"}
_EXCLUDED = {".git", "build", "__pycache__", ".pytest_cache", ".mypy_cache",
             "reports", "results", "output", "outputs"}
_ROOT_FILES = {"CMakeLists.txt", "CMakePresets.json", "CMakeUserPresets.json", "pyproject.toml",
               "setup.cfg", "setup.py", "Makefile", "requirements.txt"}
_CODE_SUFFIXES = {".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx", ".cmake", ".py", ".sh"}
_DISK_SCOPE = "Observed files on disk; not a measurement of loaded memory"
_BEGIN_SCHEMA = "make-it-stand-run-inputs-v2"
_PAPER_SCHEMA = "make-it-stand-paper-run-evidence-v2"


class EvidenceValidationError(ValueError):
    """An input or serialized observation cannot be accepted."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceValidationError(message)


def _json_value(value: Any) -> None:
    if value is None or type(value) in (str, bool, int):
        return
    if type(value) is float:
        _require(math.isfinite(value), "nonfinite JSON number")
        return
    if type(value) is list:
        for item in value:
            _json_value(item)
        return
    if type(value) is dict:
        for key, item in value.items():
            _require(type(key) is str, "JSON object keys must be strings")
            _json_value(item)
        return
    raise EvidenceValidationError("non-JSON value: " + type(value).__name__)


def _canonical(value: Any) -> str:
    _json_value(value)
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


def _digest(value: Any) -> str:
    return hashlib.sha256(_canonical(value).encode("utf-8")).hexdigest()


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        _require(key not in result, "duplicate JSON key: " + key)
        result[key] = value
    return result


def _parse(payload: str) -> dict:
    _require(type(payload) is str, "serialized evidence must be JSON text")
    try:
        result = json.loads(payload, object_pairs_hook=_pairs)
        _json_value(result)
    except (TypeError, RecursionError, json.JSONDecodeError) as exc:
        raise EvidenceValidationError("malformed JSON evidence") from exc
    _require(type(result) is dict, "evidence must be an object")
    return result


def _keys(value, expected, label):
    _require(type(value) is dict and set(value) == set(expected), "invalid " + label + " fields")


def _nonempty_object(value, label):
    _require(type(value) is dict and bool(value), label + " must be a nonempty JSON object")
    _json_value(value)


def _text(value, label):
    _require(type(value) is str and bool(value.strip()), label + " must be nonempty text")


def _path(raw, *, directory=False, resolve_runtime=False) -> Path:
    _require(isinstance(raw, (str, os.PathLike)) and bool(os.fspath(raw)), "missing path")
    path = Path(raw)
    if not path.is_absolute():
        path = Path.cwd() / path
    if resolve_runtime:
        path = path.resolve()
    for part in (path, *path.parents):
        _require(not part.is_symlink(), "symlink is unsupported: " + str(part))
    path = Path(os.path.abspath(path))
    _require(path.is_dir() if directory else path.is_file(), "missing or unsupported path: " + str(path))
    return path


def _record(path: Path, *, label=None) -> dict:
    path = _path(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        before = os.fstat(stream.fileno())
        _require(stat.S_ISREG(before.st_mode), "not a regular file: " + str(path))
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
        after = os.fstat(stream.fileno())
    now = path.stat()
    identity = lambda item: (item.st_dev, item.st_ino, item.st_size, item.st_mtime_ns, item.st_ctime_ns)
    _require(identity(before) == identity(after) == identity(now), "file changed while hashing: " + str(path))
    return {"path": str(path) if label is None else label, "size": after.st_size, "sha256": digest.hexdigest()}


def _records(paths, *, root=None, nonempty=True):
    _require(type(paths) in (list, tuple), "file inputs must be an explicit list")
    _require(bool(paths) or not nonempty, "file inputs cannot be empty")
    seen_paths, seen_inodes, result = set(), set(), []
    for raw in paths:
        path = _path(raw)
        inode = (path.stat().st_dev, path.stat().st_ino)
        _require(str(path) not in seen_paths and inode not in seen_inodes, "duplicate file input: " + str(path))
        seen_paths.add(str(path))
        seen_inodes.add(inode)
        result.append(_record(path, label=path.relative_to(root).as_posix() if root else None))
    return sorted(result, key=lambda record: record["path"])


def _selected(relative: str) -> bool:
    path = Path(relative)
    if set(path.parts) & _EXCLUDED or path.suffix in {".pyc", ".pyo"}:
        return False
    if len(path.parts) == 1:
        return path.name in _ROOT_FILES or path.suffix in _CODE_SUFFIXES or path.name.startswith("requirements")
    return path.parts[0] in _TREES or (path.parts[0] == "experiments" and
           (path.suffix in _CODE_SUFFIXES or path.name == "CMakeLists.txt"))


def _walk(root, *, source=False):
    result = []
    def onerror(exc):
        raise EvidenceValidationError("cannot scan input tree: " + str(exc)) from exc
    for directory, dirs, files in os.walk(root, followlinks=False, onerror=onerror):
        here = Path(directory)
        dirs[:] = sorted(name for name in dirs if name not in _EXCLUDED and
                         (not source or here != root or name in _TREES | {"experiments"}))
        for name in dirs:
            _require(not (here / name).is_symlink(), "symlink directory is unsupported: " + str(here / name))
        for name in sorted(files):
            path = here / name
            if path.suffix in {".pyc", ".pyo"}:
                continue
            if not source or _selected(path.relative_to(root).as_posix()):
                result.append(path)
    return result


def _git(root, *args) -> bytes:
    try:
        return subprocess.run(["git", "-C", str(root), *args], check=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
    except (OSError, subprocess.CalledProcessError) as exc:
        raise EvidenceValidationError("Git observation failed at " + str(root)) from exc


def _head(root):
    actual_root = _git(root, "rev-parse", "--show-toplevel").decode().strip()
    _require(Path(actual_root).resolve() == root, "input must be the actual Git repository root: " + str(root))
    return _git(root, "rev-parse", "--verify", "HEAD").decode().strip()


def _source(root):
    head = _head(root)
    files = _records(_walk(root, source=True), root=root)
    tracked = _git(root, "ls-tree", "-r", "--name-only", "-z", "HEAD").decode("utf-8").split("\0")
    indexed = _git(root, "ls-files", "--cached", "-z").decode("utf-8").split("\0")
    paths = sorted({record["path"] for record in files} | {path for path in tracked + indexed if path and _selected(path)})
    context = {"head": head, "diff_paths": paths}
    for name, args in (("head_to_worktree", ["HEAD"]), ("head_to_index", ["--cached", "HEAD"]),
                       ("index_to_worktree", [])):
        diff = _git(root, "--literal-pathspecs", "diff", "--binary", "--no-ext-diff", "--no-textconv", *args, "--", *paths)
        context[name + "_diff_base64"] = base64.b64encode(diff).decode("ascii")
    _require(head == _head(root), "Git HEAD changed while capturing source")
    return {"policy": SOURCE_POLICY, "selected_paths": [item["path"] for item in files],
            "files": files, "content_sha256": _digest(files), "git": context}


def _dependencies(roots):
    observed = {}
    for name, raw in sorted(roots.items()):
        root = _path(raw, directory=True)
        files = _records(_walk(root), root=root)
        revision = _head(root) if (root / ".git").exists() else None
        versions = {}
        for item in files:
            if Path(item["path"]).name.lower() in {"version", "version.txt", "cmakelists.txt", "pyproject.toml", "package.json"}:
                # Preserve bounded, actual metadata text, including build-declared versions.
                if item["size"] <= 1024 * 1024:
                    versions[item["path"]] = (root / item["path"]).read_text(encoding="utf-8")
        observed[name] = {"root": str(root), "files": files, "content_sha256": _digest(files),
                          "git_head": revision, "version_observations": versions,
                          "version_observation_status": "OBSERVED" if revision or versions else "UNAVAILABLE"}
    return {"coverage": "DECLARED_ROOTS_ONLY", "roots": observed}


def _module_names(modules):
    _require(type(modules) in (list, tuple), "native_modules must be an explicit list")
    names = []
    for requested in modules:
        if isinstance(requested, ModuleType):
            name = requested.__name__
            _require(sys.modules.get(name) is requested, "requested module object is not loaded: " + name)
        else:
            _text(requested, "native module name")
            name = requested
        _require(name not in names, "duplicate native module: " + name)
        names.append(name)
    return sorted(names)


def _runtime(names, pure_python):
    _require(type(pure_python) is bool, "pure_python must be a boolean")
    _require((not names) == pure_python, "explicit pure-Python run or nonempty native modules required")
    executable = _path(sys.executable, resolve_runtime=True)
    modules = {}
    for name in names:
        module = sys.modules.get(name)
        _require(isinstance(module, ModuleType), "requested native module is not loaded: " + name)
        spec = getattr(module, "__spec__", None)
        _require(spec is not None and isinstance(spec.loader, ExtensionFileLoader) and spec.name == name,
                 "requested module is not a native extension: " + name)
        path = _path(getattr(module, "__file__", ""))
        _require(path == _path(spec.origin), "module file does not match imported origin: " + name)
        modules[name] = _record(path)
    return {"executable": _record(executable), "native_modules": modules,
            "pure_python": pure_python, "version": sys.version, "platform": platform.platform(),
            "measurement_scope": _DISK_SCOPE}


_FILE_GROUPS = ("runner_files", "checker_files", "asset_files", "build_configuration_files")
_REQUEST_KEYS = {"repository_root", *_FILE_GROUPS, "parameters", "seed_protocol", "dependency_roots", "native_modules", "pure_python"}


def _validate_request(request):
    _keys(request, _REQUEST_KEYS, "capture request")
    _text(request["repository_root"], "repository root")
    _require(Path(request["repository_root"]).is_absolute(), "repository root must be absolute")
    _nonempty_object(request["parameters"], "parameters")
    _keys(request["seed_protocol"], {"seed", "protocol"}, "seed protocol")
    seed = request["seed_protocol"]["seed"]
    _require(type(seed) is int and seed >= 0, "seed must be a nonnegative integer, not a boolean")
    _text(request["seed_protocol"]["protocol"], "seed protocol description")
    for group in _FILE_GROUPS:
        paths = request[group]
        _require(type(paths) is list and bool(paths), group + " cannot be empty")
        _require(all(type(path) is str and Path(path).is_absolute() for path in paths), group + " needs absolute file paths")
        _require(paths == sorted(set(paths)), "duplicate or unordered " + group)
    roots = request["dependency_roots"]
    _nonempty_object(roots, "dependency roots")
    seen = []
    for name, raw in roots.items():
        _text(name, "dependency name")
        _text(raw, "dependency root")
        path = Path(raw)
        _require(path.is_absolute(), "dependency root must be absolute")
        _require(all(path != old and path not in old.parents and old not in path.parents for old in seen),
                 "duplicate or overlapping dependency roots")
        seen.append(path)
    names = request["native_modules"]
    _require(type(names) is list and all(type(name) is str and bool(name) for name in names), "invalid native module names")
    _require(names == sorted(set(names)), "duplicate native modules")
    _require(type(request["pure_python"]) is bool and ((not names) == request["pure_python"]), "invalid runtime mode")


def _observe(request):
    root = _path(request["repository_root"], directory=True)
    return {"source": _source(root), "dependencies": _dependencies(request["dependency_roots"]),
            "files": {group: _records(request[group]) for group in _FILE_GROUPS},
            "runtime": _runtime(request["native_modules"], request["pure_python"])}


def _seal(value):
    return {**value, "integrity_sha256": _digest(value)}


def _validate_integrity(value):
    _require(value.get("integrity_sha256") == _digest({key: item for key, item in value.items() if key != "integrity_sha256"}),
             "evidence integrity digest mismatch")


def _valid_revision(value):
    return (type(value) is str and len(value) in (40, 64) and
            all(char in "0123456789abcdef" for char in value))


def _validate_records(records, *, absolute=False):
    _require(type(records) is list and bool(records), "missing file content records")
    paths = []
    for item in records:
        _keys(item, {"path", "size", "sha256"}, "file record")
        _text(item["path"], "file path")
        path = Path(item["path"])
        _require(path.is_absolute() == absolute and ".." not in path.parts, "invalid file path identity")
        _require(type(item["size"]) is int and item["size"] >= 0, "invalid file size")
        _require(type(item["sha256"]) is str and len(item["sha256"]) == 64 and
                 all(char in "0123456789abcdef" for char in item["sha256"]), "invalid SHA256")
        paths.append(item["path"])
    _require(paths == sorted(set(paths)), "duplicate or unordered file records")


def _validate_begin(value):
    _keys(value, {"schema", "request", "observation", "integrity_sha256"}, "run inputs")
    _require(value["schema"] == _BEGIN_SCHEMA, "unsupported run input schema")
    _validate_integrity(value)
    _validate_request(value["request"])
    obs = value["observation"]
    _keys(obs, {"source", "dependencies", "files", "runtime"}, "observation")
    source = obs["source"]
    _keys(source, {"policy", "selected_paths", "files", "content_sha256", "git"}, "source")
    _require(source["policy"] == SOURCE_POLICY, "unsupported source selection policy")
    _validate_records(source["files"])
    _require(source["content_sha256"] == _digest(source["files"]), "source digest mismatch")
    _require(source["selected_paths"] == [item["path"] for item in source["files"]] and
             all(_selected(path) for path in source["selected_paths"]), "invalid source selection")
    context = source["git"]
    _keys(context, {"head", "diff_paths", "head_to_worktree_diff_base64", "head_to_index_diff_base64", "index_to_worktree_diff_base64"}, "Git context")
    _require(_valid_revision(context["head"]), "invalid observed Git HEAD")
    _require(type(context["diff_paths"]) is list and bool(context["diff_paths"]) and
             all(type(path) is str and _selected(path) for path in context["diff_paths"]), "invalid Git diff paths")
    _require(context["diff_paths"] == sorted(set(context["diff_paths"])), "duplicate or unordered Git diff paths")
    for key in context:
        if key.endswith("_base64"):
            try:
                base64.b64decode(context[key], validate=True)
            except (ValueError, TypeError) as exc:
                raise EvidenceValidationError("invalid binary Git diff") from exc
    _keys(obs["files"], _FILE_GROUPS, "input files")
    for group in _FILE_GROUPS:
        _validate_records(obs["files"][group], absolute=True)
        _require([item["path"] for item in obs["files"][group]] == value["request"][group], "file input coverage mismatch")
    deps = obs["dependencies"]
    _keys(deps, {"coverage", "roots"}, "dependencies")
    _require(deps["coverage"] == "DECLARED_ROOTS_ONLY", "invalid dependency coverage")
    _keys(deps["roots"], value["request"]["dependency_roots"], "dependency roots")
    for name, item in deps["roots"].items():
        _keys(item, {"root", "files", "content_sha256", "git_head", "version_observations", "version_observation_status"}, "dependency")
        _require(item["root"] == value["request"]["dependency_roots"][name], "dependency root mismatch")
        _validate_records(item["files"])
        _require(item["content_sha256"] == _digest(item["files"]), "dependency digest mismatch")
        _require(type(item["version_observations"]) is dict and item["version_observation_status"] in {"OBSERVED", "UNAVAILABLE"}, "invalid dependency version observations")
        _require(item["git_head"] is None or _valid_revision(item["git_head"]), "invalid dependency revision")
        _require(all(type(text) is str for text in item["version_observations"].values()), "invalid dependency version text")
        version_status = "OBSERVED" if item["git_head"] or item["version_observations"] else "UNAVAILABLE"
        _require(item["version_observation_status"] == version_status, "inconsistent version observation status")
    runtime = obs["runtime"]
    _keys(runtime, {"executable", "native_modules", "pure_python", "version", "platform", "measurement_scope"}, "runtime")
    _validate_records([runtime["executable"]], absolute=True)
    _keys(runtime["native_modules"], value["request"]["native_modules"], "native module observations")
    for record in runtime["native_modules"].values():
        _validate_records([record], absolute=True)
    _require(type(runtime["pure_python"]) is bool and runtime["pure_python"] == value["request"]["pure_python"] and runtime["measurement_scope"] == _DISK_SCOPE, "invalid runtime scope")
    _text(runtime["version"], "runtime version")
    _text(runtime["platform"], "runtime platform")


@dataclass(frozen=True)
class RunInputsV2:
    """Immutable begin/end observation. Not a completed scientific result."""

    payload_json: str

    def __post_init__(self):
        value = _parse(self.payload_json)
        _validate_begin(value)
        object.__setattr__(self, "payload_json", _canonical(value))

    def as_dict(self):
        return _parse(self.payload_json)

    def to_json(self):
        return self.payload_json


def begin_run(*, repository_root, runner_files, checker_files, asset_files,
              build_configuration_files, parameters, seed_protocol,
              dependency_roots, native_modules, pure_python) -> RunInputsV2:
    """Read actual inputs now; no caller-provided commit/hash/verification flags."""
    _require(type(dependency_roots) is dict and bool(dependency_roots), "dependency roots must be nonempty")
    request = {"repository_root": str(_path(repository_root, directory=True)),
               "runner_files": [item["path"] for item in _records(runner_files)],
               "checker_files": [item["path"] for item in _records(checker_files)],
               "asset_files": [item["path"] for item in _records(asset_files)],
               "build_configuration_files": [item["path"] for item in _records(build_configuration_files)],
               "parameters": parameters, "seed_protocol": seed_protocol,
               "dependency_roots": {name: str(_path(root, directory=True)) for name, root in dependency_roots.items()},
               "native_modules": _module_names(native_modules), "pure_python": pure_python}
    _validate_request(request)
    request = _parse(_canonical(request))
    return RunInputsV2(_canonical(_seal({"schema": _BEGIN_SCHEMA, "request": request, "observation": _observe(request)})))


def _validate_events(events):
    _require(type(events) is list and bool(events), "scientific events must be a nonempty ordered list")
    for event in events:
        _keys(event, {"iteration", "status", "data"}, "scientific event")
        _require(type(event["iteration"]) is int and event["iteration"] >= 0, "iteration must be a nonnegative integer")
        _require(event["status"] in ("ok", "failed"), "event status must preserve ok/failed attempts")
        _nonempty_object(event["data"], "scientific event data")


def _output_records(paths):
    records = _records(paths)
    for record in records:
        path = Path(record["path"])
        _require(record["size"] > 0, "empty scientific output: " + str(path))
        if path.suffix.lower() == ".json":
            try:
                with path.open(encoding="utf-8") as stream:
                    value = json.load(stream, object_pairs_hook=_pairs)
                _json_value(value)
            except ValueError as exc:
                raise EvidenceValidationError("invalid scientific JSON: " + str(path) + ": " + str(exc)) from exc
            _require(type(value) in (dict, list) and bool(value), "empty or malformed scientific JSON: " + str(path))
            _require(_record(path) == record, "scientific output changed while validating: " + str(path))
    return records


def _validate_paper(value):
    _keys(value, {"schema", "begin", "end", "scientific", "measurements", "build_source_attestation", "physical_validation", "integrity_sha256"}, "paper evidence")
    _require(value["schema"] == _PAPER_SCHEMA, "unsupported paper evidence schema")
    _validate_integrity(value)
    for key in ("begin", "end"):
        _validate_begin(value[key])
    _require(value["begin"]["request"] == value["end"]["request"], "begin/end input specification mismatch")
    _require(value["build_source_attestation"] == "NOT_VERIFIED" and value["physical_validation"] == "NOT_IN_SCOPE", "observation cannot grant attestation")
    scientific = value["scientific"]
    _keys(scientific, {"events", "outputs", "determinism_sha256"}, "scientific evidence")
    _validate_events(scientific["events"])
    _validate_records(scientific["outputs"], absolute=True)
    _require(all(item["size"] > 0 for item in scientific["outputs"]), "empty scientific output")
    _require(scientific["determinism_sha256"] == _digest({"events": scientific["events"], "outputs": scientific["outputs"]}), "scientific determinism digest mismatch")
    _nonempty_object(value["measurements"], "measurements")


@dataclass(frozen=True)
class PaperRunEvidenceV2:
    """Canonical JSON freezes every nested input and observation by value."""

    payload_json: str

    def __post_init__(self):
        value = _parse(self.payload_json)
        _validate_paper(value)
        object.__setattr__(self, "payload_json", _canonical(value))

    @classmethod
    def from_json(cls, payload):
        return cls(payload)

    def as_dict(self):
        return _parse(self.payload_json)

    def to_json(self):
        return self.payload_json


def finish_run(begin: RunInputsV2, *, events, output_files, measurements) -> PaperRunEvidenceV2:
    """Capture the end even if inputs changed; freshness retains that mismatch.

    Missing/unreadable end inputs raise EvidenceValidationError; they cannot
    create completed evidence. Measurements (timestamps, timings, resources)
    are kept separately from ordered scientific events and their digest.
    """
    _require(type(begin) is RunInputsV2, "finish_run requires a V2 begin observation")
    _validate_events(events)
    _nonempty_object(measurements, "measurements")
    initial = begin.as_dict()
    outputs = _output_records(output_files)
    scientific = {"events": events, "outputs": outputs}
    end = _seal({"schema": _BEGIN_SCHEMA, "request": initial["request"], "observation": _observe(initial["request"])})
    return PaperRunEvidenceV2(_canonical(_seal({"schema": _PAPER_SCHEMA, "begin": initial, "end": end,
        "scientific": {**scientific, "determinism_sha256": _digest(scientific)}, "measurements": measurements,
        "build_source_attestation": "NOT_VERIFIED", "physical_validation": "NOT_IN_SCOPE"})))


@dataclass(frozen=True)
class FreshnessAssessmentV2:
    status: str
    mismatches: tuple[str, ...] = ()
    contextual_changes: tuple[str, ...] = ()
    build_source_attestation: str = "NOT_VERIFIED"
    physical_validation: str = "NOT_IN_SCOPE"


def _differences(before, after, location="observation"):
    if before == after:
        return []
    if type(before) is dict and type(after) is dict:
        if "path" in before and "path" in after:
            location += "[" + str(before["path"]) + "]"
        return [item for key in sorted(set(before) | set(after))
                for item in _differences(before.get(key), after.get(key), location + "." + key)]
    if type(before) is list and type(after) is list and all(type(item) is dict and "path" in item for item in before + after):
        return _differences({item["path"]: item for item in before}, {item["path"]: item for item in after}, location)
    return [location + " changed"]


def _compare(before, after):
    # Copy first: preserve all contextual data in the actual immutable evidence.
    before, after = _parse(_canonical(before)), _parse(_canonical(after))
    contextual = _differences(before["source"].pop("git"), after["source"].pop("git"), "source.git")
    for name in before["dependencies"]["roots"]:
        old, new = before["dependencies"]["roots"][name], after["dependencies"]["roots"][name]
        contextual += _differences(old.pop("git_head"), new.pop("git_head"), "dependency." + name + ".git_head")
    return _differences(before, after), contextual


def check_freshness(evidence) -> FreshnessAssessmentV2:
    """Re-read effective inputs/output files. No `passed`/`verified` override."""
    try:
        if type(evidence) is PaperRunEvidenceV2:
            value = evidence.as_dict()
        elif type(evidence) is str:
            value = _parse(evidence)
        elif type(evidence) is dict:
            value = _parse(_canonical(evidence))
        else:
            raise EvidenceValidationError("missing or unsupported V2 evidence")
        _validate_paper(value)
    except (ValueError, TypeError, KeyError, RecursionError) as exc:
        return FreshnessAssessmentV2("INVALID_EVIDENCE", (str(exc),))
    mismatches, contextual = _compare(value["begin"]["observation"], value["end"]["observation"])
    try:
        fresh = _observe(value["begin"]["request"])
        new_mismatches, new_contextual = _compare(value["end"]["observation"], fresh)
        mismatches += new_mismatches
        contextual += new_contextual
        outputs = _output_records([item["path"] for item in value["scientific"]["outputs"]])
        mismatches += _differences(value["scientific"]["outputs"], outputs, "scientific.outputs")
    except (OSError, ValueError, TypeError) as exc:
        mismatches.append("recapture failed: " + str(exc))
    return FreshnessAssessmentV2("STALE_RESULT" if mismatches else "CURRENT",
                                 tuple(dict.fromkeys(mismatches)), tuple(dict.fromkeys(contextual)))
