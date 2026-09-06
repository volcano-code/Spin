"""Opt-in GNU/Linux compiler observations, NOT full build/source attestation.

Invoke as a CMake compiler launcher with --evidence-dir DIR --timeout SECONDS
-- COMPILER ARGS. Outputs must be absent: use a fresh build directory (or let
the build owner explicitly remove an old object before retrying). No shell,
response files, PCH, modules or unknown compiler flags are admitted.

The two dependency probes observe preprocessing closure, including system
headers. Before/after disk hashes cannot exclude transient reverted changes.
Inline assembly can read files outside this preprocessing closure; assembler
input coverage remains NOT_VERIFIED. Linking, build control, dynamically loaded
toolchain components, fresh runtime launch and scientific checks remain
NOT_VERIFIED. Receipts are unauthenticated.
The total invocation timeout is not aggregate memory accounting. Timing is
contextual only. This observer is not a sandbox for an untrusted compiler.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import time

if __package__:
    from .compiler_dependencies import ObservationError, parse_command, parse_dependencies
else:
    from compiler_dependencies import ObservationError, parse_command, parse_dependencies


ENVIRONMENT_KEYS = ("PATH", "CPATH", "C_INCLUDE_PATH", "CPLUS_INCLUDE_PATH", "OBJC_INCLUDE_PATH",
                    "LANG", "LC_ALL", "LC_CTYPE", "LC_MESSAGES", "SOURCE_DATE_EPOCH", "TMPDIR",
                    "GCC_COLORS", "LD_LIBRARY_PATH", "LIBRARY_PATH")
UNSUPPORTED_ENVIRONMENT = ("LD_PRELOAD", "LD_AUDIT", "GCC_EXEC_PREFIX", "COMPILER_PATH",
                           "DEPENDENCIES_OUTPUT", "SUNPRO_DEPENDENCIES", "GCC_COMPARE_DEBUG")
NOT_VERIFIED = ("build_source_attestation", "link_coverage", "build_control", "dynamic_compiler_components",
                "assembler_input_coverage", "fresh_runtime_launch", "scientific_checks")


def _write_json(path, data):
    with path.open("x", encoding="utf-8") as stream:
        json.dump(data, stream, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False)
        stream.write("\n")


def _path_links(path):
    """Observe each traversed symlink, respecting symlinks before '..'."""
    pending = list(path.parts[1:])
    current = Path(path.anchor)
    links = []
    while pending:
        part = pending.pop(0)
        if part == "..":
            current = current.parent
            continue
        candidate = current / part
        info = candidate.lstat()
        if stat.S_ISLNK(info.st_mode):
            target = os.readlink(candidate)
            links.append({"path": str(candidate), "target": target, "device": info.st_dev,
                          "inode": info.st_ino, "mtime_ns": info.st_mtime_ns, "ctime_ns": info.st_ctime_ns})
            if len(links) > 40:
                raise ObservationError("too many symlinks in input path")
            target_path = Path(target)
            if target_path.is_absolute():
                current = Path(target_path.anchor)
                pending = list(target_path.parts[1:]) + pending
            else:
                pending = list(target_path.parts) + pending
        else:
            current = candidate
    return current, links


def _snapshot(raw, cwd, deadline=None):
    path = Path(raw)
    if not path.is_absolute():
        path = cwd / path
    resolved, links = _path_links(path)
    digest = hashlib.sha256()
    if not stat.S_ISREG(path.stat().st_mode):
        raise ObservationError("input is not a regular file: " + str(path))
    with os.fdopen(os.open(path, os.O_RDONLY | os.O_NONBLOCK), "rb") as stream:
        before = os.fstat(stream.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise ObservationError("input is not a regular file: " + str(path))
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            if deadline is not None and time.monotonic() >= deadline:
                raise subprocess.TimeoutExpired(["hash", str(path)], 0)
            digest.update(block)
        after = os.fstat(stream.fileno())
    now = path.stat()
    identity = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
    if identity(before) != identity(after) or identity(after) != identity(now) or _path_links(path) != (resolved, links):
        raise ObservationError("file changed while hashing: " + str(path))
    return {"logical_path": str(raw), "resolved_path": str(resolved), "symlinks": links, "sha256": digest.hexdigest(),
            "size": after.st_size, "device": after.st_dev, "inode": after.st_ino,
            "mtime_ns": after.st_mtime_ns, "ctime_ns": after.st_ctime_ns}


def _stop_group(child):
    try:
        os.killpg(child.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        child.wait(timeout=0.2)
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(child.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    child.wait()


def _run(argv, *, cwd, directory, phase, deadline, receipt):
    """Execute real argv; all phase logs are exclusive, including failed phases."""
    operation = {"phase": phase, "argv": argv, "exit_status": None}
    receipt["operations"].append(operation)
    start = time.monotonic()
    with (directory / (phase + ".stdout")).open("xb") as out, (directory / (phase + ".stderr")).open("xb") as err:
        child = subprocess.Popen(argv, cwd=cwd, stdout=out, stderr=err, start_new_session=True)
        try:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise subprocess.TimeoutExpired(argv, 0)
            code = child.wait(timeout=remaining)
            operation["exit_status"] = code
        except BaseException:
            _stop_group(child)
            operation["exit_status"] = child.returncode
            operation["incomplete"] = True
            raise
        finally:
            operation["elapsed_seconds"] = time.monotonic() - start
    return code


def _compiler(executable, cwd, directory, phase, deadline, receipt):
    record = _snapshot(executable, cwd, deadline)
    with Path(executable).open("rb") as stream:
        if stream.read(4) != b"\x7fELF":
            raise ObservationError("direct GNU ELF compiler binary required; script wrappers are unsupported")
    code = _run([executable, "--version"], cwd=cwd, directory=directory, phase=phase, deadline=deadline, receipt=receipt)
    version = (directory / (phase + ".stdout")).read_text(errors="replace")
    if code or "Free Software Foundation" not in version:
        raise ObservationError("GNU compiler required; version probe failed or unsupported compiler")
    record["version"] = version
    return record


def _probe(command, cwd, directory, phase, deadline, receipt):
    code = _run(command.probe_argv, cwd=cwd, directory=directory, phase=phase, deadline=deadline, receipt=receipt)
    if code:
        raise ObservationError("dependency probe failed: " + phase)
    dependency_text = (directory / (phase + ".stdout")).read_text(errors="strict")
    # -fpch-preprocess exposes selected PCH even when only the .gch exists
    # in an earlier search directory and plain -M would report a later header.
    if "#pragma GCC pch_preprocess" in dependency_text:
        raise ObservationError("selected precompiled header is unsupported")
    paths = parse_dependencies(dependency_text)
    snapshots = [_snapshot(p, cwd, deadline) for p in paths]
    source = (cwd / command.source).resolve(strict=True)
    if str(source) not in {r["resolved_path"] for r in snapshots}:
        raise ObservationError("source absent from dependency closure")
    for record in snapshots:
        for name in (str(cwd / record["logical_path"]), record["resolved_path"]):
            if os.path.lexists(name + ".gch"):
                raise ObservationError("implicit precompiled header is unsupported: " + name + ".gch")
    return snapshots


def observe(argv, *, evidence_dir, timeout):
    """Observe an invocation and return its own execution status and receipt path.

    No caller-supplied hashes, receipt contents, or success booleans are accepted.
    """
    if type(timeout) not in (int, float) or not math.isfinite(timeout) or timeout <= 0:
        raise ObservationError("timeout must be a finite positive number of seconds")
    if not str(evidence_dir).strip():
        raise ObservationError("explicit evidence directory required")
    cwd = Path.cwd()
    root = Path(evidence_dir).absolute()
    root.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix="invocation-", dir=root))
    started = time.monotonic()
    deadline = started + timeout
    receipt = {"schema": "make-it-stand-compiler-observation-v1", "invocation_id": directory.name,
               "status": "IN_PROGRESS", "argv": list(argv), "cwd": str(cwd), "compiler_exit_status": None,
               "operations": [], "settings": {"timeout_seconds": timeout, "evidence_dir": str(root)},
               "environment": {k: os.environ.get(k) for k in ENVIRONMENT_KEYS},
               "limitations": ["before/after observations cannot exclude transient reverted changes",
                               "fresh output required; remove old outputs explicitly or use a fresh build directory",
                               "unauthenticated disk observations; no loaded-memory measurement",
                               "preprocessing closure only; inline assembly can introduce unobserved input files",
                               "elapsed time is contextual, not scientific determinism evidence",
                               "timeout is not aggregate 24 GiB memory accounting"],
               **{name: "NOT_VERIFIED" for name in NOT_VERIFIED}}
    _write_json(directory / "in-progress.json", receipt)
    result = 1
    try:
        if sys.platform != "linux":
            raise ObservationError("compiler observer requires GNU/Linux")
        if any(os.environ.get(k) for k in UNSUPPORTED_ENVIRONMENT):
            raise ObservationError("unsupported compiler-affecting environment: " + ", ".join(k for k in UNSUPPORTED_ENVIRONMENT if os.environ.get(k)))
        command = parse_command(argv)
        executable = shutil.which(argv[0])
        if executable is None:
            raise ObservationError("compiler executable not found")
        output = cwd / command.output
        if os.path.lexists(output):
            raise ObservationError("pre-existing output rejected; use a fresh build directory: " + str(output))
        protected = [output] + ([cwd / command.depfile] if command.depfile else [])
        for path in protected:
            if path.resolve().is_relative_to(root.resolve()):
                raise ObservationError("compiler output collides with evidence directory")
            if path.is_symlink():
                raise ObservationError("symlink dependency output unsupported")
        if len({str(p.resolve()) for p in protected}) != len(protected):
            raise ObservationError("object and dependency output collide")
        receipt["observer_code"] = [_snapshot(str(Path(__file__).resolve()), cwd, deadline),
                                    _snapshot(str(Path(__file__).with_name("compiler_dependencies.py").resolve()), cwd, deadline)]
        receipt["python"] = _snapshot(sys.executable, cwd, deadline)
        receipt["compiler_before"] = _compiler(executable, cwd, directory, "compiler-before", deadline, receipt)
        receipt["inputs_before"] = _probe(command, cwd, directory, "dependencies-before", deadline, receipt)
        inputs = receipt["inputs_before"] + receipt["observer_code"] + [receipt["python"], receipt["compiler_before"]]
        for path in protected:
            if str(path.resolve()) in {r["resolved_path"] for r in inputs} or (path.exists() and (path.stat().st_dev, path.stat().st_ino) in {(r["device"], r["inode"]) for r in inputs}):
                raise ObservationError("compiler output collides with observed input")
        if shutil.which(argv[0]) != executable:
            raise ObservationError("compiler executable search changed before compilation")
        # An exclusively created empty reservation prevents two observers even
        # with different evidence roots from accepting the same absent output.
        # Failed reservations/partial objects are retained for explicit cleanup.
        with output.open("xb"):
            pass
        receipt["output_reservation"] = _snapshot(command.output, cwd, deadline)
        code = _run(argv, cwd=cwd, directory=directory, phase="compile", deadline=deadline, receipt=receipt)
        receipt["compiler_exit_status"] = code
        result = code if code > 0 else 128 - code if code < 0 else 1
        post_errors = []
        if shutil.which(argv[0]) != executable:
            post_errors.append("compiler executable search changed")
        receipt["inputs_rehashed_after"] = []
        for record in receipt["inputs_before"]:
            try:
                observed = _snapshot(record["logical_path"], cwd, deadline)
            except (OSError, ValueError) as exc:
                observed = {"logical_path": record["logical_path"], "error": str(exc)}
            receipt["inputs_rehashed_after"].append(observed)
        if receipt["inputs_rehashed_after"] != receipt["inputs_before"]:
            post_errors.append("original dependency bytes/identity changed or missing")
        for name, action in (
            ("inputs_after", lambda: _probe(command, cwd, directory, "dependencies-after", deadline, receipt)),
            ("compiler_after", lambda: _compiler(executable, cwd, directory, "compiler-after", deadline, receipt)),
        ):
            try:
                receipt[name] = action()
            except (OSError, ValueError) as exc:
                post_errors.append(name + ": " + str(exc))
        if code:
            post_errors.append("compiler failed with exit status " + str(code))
        if receipt.get("inputs_after") != receipt["inputs_before"]:
            post_errors.append("dependency closure or input bytes/identity changed")
        if receipt.get("compiler_after") != receipt["compiler_before"]:
            post_errors.append("compiler identity/version changed")
        if [_snapshot(r["logical_path"], cwd, deadline) for r in receipt["observer_code"]] != receipt["observer_code"]:
            post_errors.append("observer code changed")
        if _snapshot(sys.executable, cwd, deadline) != receipt["python"]:
            post_errors.append("observer interpreter changed")
        if not code:
            receipt["object"] = _snapshot(command.output, cwd, deadline)
            if not receipt["object"]["size"]:
                post_errors.append("empty compiler output")
        if post_errors:
            raise ObservationError("; ".join(post_errors))
        receipt["status"] = "OBSERVED_SUCCESS"
        result = 0
    except BaseException as exc:
        receipt["status"] = "FAILED"
        receipt["error"] = type(exc).__name__ + ": " + str(exc)
        if isinstance(exc, subprocess.TimeoutExpired) and not receipt["compiler_exit_status"]:
            result = 124
        elif isinstance(exc, KeyboardInterrupt):
            result = 130
        if receipt["compiler_exit_status"] is None:
            for operation in receipt["operations"]:
                if operation["phase"] == "compile":
                    receipt["compiler_exit_status"] = operation["exit_status"]
    finally:
        receipt["elapsed_seconds"] = time.monotonic() - started
        _write_json(directory / "receipt.json", receipt)
        for operation in receipt["operations"]:
            for suffix, destination in (("stderr", sys.stderr), ("stdout", sys.stdout)):
                path = directory / (operation["phase"] + "." + suffix)
                if path.exists() and (suffix == "stderr" or operation["phase"] == "compile"):
                    destination.write(path.read_text(errors="replace"))
        if result:
            print("compiler observer: " + receipt.get("error", "incomplete") + "; receipt: " + str(directory / "receipt.json"), file=sys.stderr)
    return result, directory / "receipt.json"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-dir", required=True)
    parser.add_argument("--timeout", type=float, required=True)
    parser.add_argument("compiler_argv", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    compiler_argv = args.compiler_argv
    if not compiler_argv or compiler_argv[0] != "--":
        parser.error("separate compiler argv with --")
    def interrupted(signum, frame):
        raise KeyboardInterrupt("received signal " + str(signum))
    previous = signal.signal(signal.SIGTERM, interrupted)
    try:
        return observe(compiler_argv[1:], evidence_dir=args.evidence_dir, timeout=args.timeout)[0]
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    finally:
        signal.signal(signal.SIGTERM, previous)


if __name__ == "__main__":
    raise SystemExit(main())
