"""Fail-closed controller for the supplemental G7 post-SLQP diagnostic."""
import argparse
import base64
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time

SCOPE = "SUPPLEMENTAL_DIAGNOSTIC_ONLY"
WALL_SECONDS = 180
MEMORY_BYTES = 24 * 1024**3
MAX_DIAGNOSTIC_BYTES = 64 * 1024
MAX_PIPE_BYTES = 256 * 1024
MAX_CLEANUP_SECONDS = 5.0
CLEANUP_BUDGET_FRACTION = 0.25
G7_STAGE_INDEX = 19
G7_COMPLETED_PASSES = 20
G7_STAGE19_RECEIPT_BODY_SHA256 = "4db1c9163682612db66c4c03fd0b81eced08c205346fddff5027517cc3d7d9ce"
G7_STAGE19_CHECKPOINT_FILE_SHA256 = "cdafa89df9aace3ac24b58d558ad75ce54b2f6f9d0beb1db1d244f3a338a93bf"
G7_PASS20_TRACE_SHA256 = "7eb548bd94de5b5cdb622ac1734d80bf21f82a0bc6f44b9d46360d7ea80d0ca0"


def _sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _read_json(path):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return None


def _write_new(path, value):
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(value, stream, sort_keys=True, indent=2)
        stream.write("\n")


def _stale_sources(snapshot):
    recorded = _read_json(snapshot)
    if not isinstance(recorded, dict) or not recorded:
        return None
    stale = []
    for path, digest in recorded.items():
        if not isinstance(path, str) or not isinstance(digest, str) or len(digest) != 64:
            return None
        candidate = Path(path)
        if not candidate.is_file() or _sha256(candidate) != digest:
            stale.append(path)
    return sorted(stale)


def _verify_spec_identity(frozen):
    spec = _read_json(frozen / "spec.json")
    identity = spec.get("identity") if isinstance(spec, dict) else None
    if not isinstance(identity, dict) or not identity:
        return "UNAVAILABLE_SPEC_IDENTITY"
    for key, fact in identity.items():
        if not isinstance(key, str) or not key or not isinstance(fact, dict):
            return "UNAVAILABLE_SPEC_IDENTITY"
        path, digest, size = fact.get("path"), fact.get("sha256"), fact.get("size")
        if not isinstance(path, str) or not isinstance(digest, str) or len(digest) != 64 or not isinstance(size, int):
            return "UNAVAILABLE_SPEC_IDENTITY"
        candidate = Path(path)
        if not candidate.is_absolute() or str(candidate.resolve()) != path:
            return "UNAVAILABLE_SPEC_IDENTITY"
        if not candidate.is_file() or candidate.stat().st_size != size or _sha256(candidate) != digest:
            return "UNAVAILABLE_NATIVE_OR_REQUEST_IDENTITY"
    return None


def _g7_stage_chain_status(frozen):
    verification = _read_json(frozen / "verification.json")
    rows = verification.get("rows") if isinstance(verification, dict) else None
    if not isinstance(verification, dict) or verification.get("source_binding") != "VERIFIED":
        return "UNAVAILABLE_VERIFICATION"
    if not isinstance(rows, list) or not rows or not isinstance(rows[0], dict) or rows[0].get("checkpoint_chain") != "VERIFIED":
        return "UNAVAILABLE_CHECKPOINT_CHAIN"
    if rows[0].get("published_boundaries") != G7_COMPLETED_PASSES:
        return "UNAVAILABLE_G7_PASS_COUNT"
    stage = _read_json(frozen / "repeat-0/payload/attempts/attempt-00000019/payload/stage.json")
    if not isinstance(stage, dict) or stage.get("passes") != G7_COMPLETED_PASSES or not stage.get("checkpoint"):
        return "UNAVAILABLE_G7_STAGE19_BOUNDARY"
    checkpoint = frozen / "repeat-0/payload/checkpoints/stage-00000019.ckpt"
    receipt = frozen / "repeat-0/payload/checkpoints/commit-00000019.sha256"
    trace = frozen / "repeat-0/payload/attempts/attempt-00000020/payload/trace/pass-000020.slqptrace"
    if not checkpoint.is_file() or not receipt.is_file() or not trace.is_file():
        return "UNAVAILABLE_G7_CHAIN_ARTIFACT"
    try:
        if receipt.read_text(encoding="ascii").strip().split() != ["00000019", G7_STAGE19_RECEIPT_BODY_SHA256]:
            return "UNAVAILABLE_G7_CHECKPOINT_RECEIPT"
        if _sha256(checkpoint) != G7_STAGE19_CHECKPOINT_FILE_SHA256:
            return "UNAVAILABLE_G7_CHECKPOINT_IDENTITY"
        if _sha256(trace) != G7_PASS20_TRACE_SHA256:
            return "UNAVAILABLE_G7_TRACE_IDENTITY"
    except (OSError, UnicodeDecodeError):
        return "UNAVAILABLE_G7_CHAIN_ARTIFACT"
    return None


def _frozen_verifier_status(frozen):
    verifier = frozen / "verify_results.py"
    if not verifier.is_file():
        return "UNAVAILABLE_FROZEN_VERIFIER"
    try:
        completed = subprocess.run([sys.executable, str(verifier)], cwd=frozen, text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return "UNAVAILABLE_FROZEN_VERIFIER"
    return None if completed.returncode == 0 else "UNAVAILABLE_FROZEN_VERIFIER"


def _native_child(callback, sender):
    os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    resource.setrlimit(resource.RLIMIT_AS, (MEMORY_BYTES, MEMORY_BYTES))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    try:
        payload = callback()
        if not isinstance(payload, dict):
            raise TypeError("native diagnostic must return a compact dictionary")
        diagnostic = payload.pop("diagnostic_bytes", b"")
        if not isinstance(diagnostic, bytes) or len(diagnostic) > MAX_DIAGNOSTIC_BYTES:
            raise ValueError("diagnostic byte payload is invalid or too large")
        envelope = {"omp": os.environ.get("OMP_NUM_THREADS"), "openblas": os.environ.get("OPENBLAS_NUM_THREADS"),
                    "mkl": os.environ.get("MKL_NUM_THREADS"), "as": resource.getrlimit(resource.RLIMIT_AS)}
        if envelope != {"omp": "1", "openblas": "1", "mkl": "1", "as": (MEMORY_BYTES, MEMORY_BYTES)}:
            raise RuntimeError("resource envelope was not applied")
        encoded = json.dumps({"kind": "ok", "payload": payload,
                              "diagnostic_b64": base64.b64encode(diagnostic).decode("ascii")},
                             allow_nan=False, separators=(",", ":")).encode("utf-8")
        if len(encoded) > MAX_PIPE_BYTES:
            raise ValueError("compact diagnostic payload exceeds pipe cap")
    except BaseException as exc:
        encoded = json.dumps({"kind": "error", "detail": type(exc).__name__ + ": " + str(exc)}).encode("utf-8")
    try:
        sender.send_bytes(encoded)
    finally:
        sender.close()


class _ForkWorker:
    """Small POSIX child handle whose PID is captured directly at fork."""
    def __init__(self, pid):
        self.pid = pid
        self.exitcode = None

    def is_alive(self):
        if self.exitcode is not None:
            return False
        try:
            waited, status = os.waitpid(self.pid, os.WNOHANG)
        except ChildProcessError:
            self.exitcode = -1
            return False
        except InterruptedError:
            return True
        if waited == 0:
            return True
        self.exitcode = os.waitstatus_to_exitcode(status)
        return False

    def join(self, timeout):
        if timeout is None:
            raise ValueError("worker joins must have a finite timeout")
        deadline = time.monotonic() + max(0.0, timeout)
        while self.is_alive():
            remaining = deadline - time.monotonic()
            if remaining <= 0.0:
                return
            time.sleep(min(0.01, remaining))

    def terminate(self):
        os.kill(self.pid, signal.SIGTERM)

    def kill(self):
        os.kill(self.pid, signal.SIGKILL)


def _close_quietly(endpoint):
    if endpoint is not None:
        try:
            endpoint.close()
        except Exception:
            pass


class _WorkerCleanupFailed(RuntimeError):
    def __init__(self, pid):
        super().__init__("failed to close worker pipe and reap forked child")
        self.pid = pid


def _spawn_native_worker(callback, deadline):
    """Fork a worker and capture its PID before any fallible parent cleanup."""
    context = multiprocessing.get_context("fork")
    receiver, sender = context.Pipe(duplex=False)
    try:
        pid = os.fork()
    except BaseException:
        _close_quietly(receiver)
        _close_quietly(sender)
        raise
    if pid == 0:
        try:
            receiver.close()
            _native_child(callback, sender)
            os._exit(0)
        except BaseException:
            os._exit(1)

    process = _ForkWorker(pid)
    try:
        sender.close()
    except Exception as close_error:
        try:
            reaped = _cleanup_worker(process, deadline)
        except Exception:
            reaped = False
        _close_quietly(receiver)
        if not reaped:
            raise _WorkerCleanupFailed(process.pid) from close_error
        raise
    return receiver, process


def _join_worker_until(process, deadline):
    process.join(max(0.0, deadline - time.monotonic()))
    return not process.is_alive()


def _cleanup_worker(process, deadline):
    """Terminate and reap a worker, escalating to SIGKILL within deadline."""
    try:
        if not process.is_alive():
            return _join_worker_until(process, deadline)
    except Exception:
        return False
    try:
        process.terminate()
    except Exception:
        # Still try SIGKILL below; terminate can fail for an already-exiting
        # process or due to a transient OS-level error.
        pass
    terminate_deadline = time.monotonic() + max(0.0, deadline - time.monotonic()) / 2.0
    try:
        if _join_worker_until(process, terminate_deadline):
            return True
        try:
            process.kill()
        except ProcessLookupError:
            # The worker can exit between the TERM reap attempt and SIGKILL.
            # Reap the already-exited child instead of misreporting a leak.
            return _join_worker_until(process, deadline)
        return _join_worker_until(process, deadline)
    except Exception:
        return False


def _cleanup_failure(process):
    return "UNAVAILABLE_NATIVE_CLEANUP_FAILED", {"child_pid": process.pid}, None


def _bounded_native_call(callback, _wall_seconds=None):
    wall_seconds = WALL_SECONDS if _wall_seconds is None else _wall_seconds
    deadline = time.monotonic() + max(0.0, wall_seconds)
    cleanup_seconds = min(MAX_CLEANUP_SECONDS, max(0.0, wall_seconds) * CLEANUP_BUDGET_FRACTION)
    work_deadline = deadline - cleanup_seconds
    receiver = process = None
    try:
        receiver, process = _spawn_native_worker(callback, deadline)
        wire = None
        receive_failed = False
        while time.monotonic() < work_deadline:
            if receiver.poll(min(.1, max(0.0, work_deadline - time.monotonic()))):
                try:
                    wire = receiver.recv_bytes(MAX_PIPE_BYTES)
                except (EOFError, OSError):
                    receive_failed = True
                break
            if not process.is_alive():
                break
        if wire is None:
            wall_limited = process.is_alive() and not receive_failed
            if process.is_alive():
                if not _cleanup_worker(process, deadline):
                    return _cleanup_failure(process)
            elif not _join_worker_until(process, deadline):
                return _cleanup_failure(process)
            return ("UNAVAILABLE_RESOURCE_LIMIT_WALL" if wall_limited else "UNAVAILABLE_NATIVE_FAILURE"), None, None

        if not _join_worker_until(process, work_deadline):
            if not _cleanup_worker(process, deadline):
                return _cleanup_failure(process)
            return "UNAVAILABLE_RESOURCE_LIMIT_WALL", None, None
        if process.exitcode != 0:
            return "UNAVAILABLE_NATIVE_FAILURE", None, None
        result = json.loads(wire.decode("utf-8"))
        if result.get("kind") != "ok": return "UNAVAILABLE_NATIVE_FAILURE", result.get("detail"), None
        return "COMPLETED_SUPPLEMENTAL_DIAGNOSTIC", result["payload"], base64.b64decode(result["diagnostic_b64"])
    except _WorkerCleanupFailed as exc:
        return "UNAVAILABLE_NATIVE_CLEANUP_FAILED", {"child_pid": exc.pid}, None
    except Exception:
        if process is not None:
            try:
                if not _cleanup_worker(process, deadline):
                    return _cleanup_failure(process)
            except Exception:
                return _cleanup_failure(process)
        return "UNAVAILABLE_NATIVE_FAILURE", None, None
    finally:
        _close_quietly(receiver)


def run_controller(frozen_root, snapshot_path, batch_directory, native_invoker=None):
    """Fail closed; no command-line native adapter exists yet."""
    frozen, batch = Path(frozen_root).resolve(), Path(batch_directory).resolve()
    batch.mkdir(parents=True, exist_ok=False)
    result = {"scope": SCOPE, "resource_envelope": {"wall_seconds": WALL_SECONDS, "address_space_bytes": MEMORY_BYTES,
              "numerical_library_threads": 1}, "native_invoked": False, "binary_validation": "NOT_RUN",
              "scientific_acceptance": "NOT_GRANTED"}
    stale = _stale_sources(snapshot_path)
    if stale is None: result.update(status="UNAVAILABLE_SOURCE_SNAPSHOT", stale_sources=[])
    elif stale: result.update(status="UNAVAILABLE_STALE_SOURCE", stale_sources=stale)
    else:
        status = _verify_spec_identity(frozen) or _g7_stage_chain_status(frozen) or _frozen_verifier_status(frozen)
        if status: result.update(status=status, stale_sources=[])
        elif native_invoker is None: result.update(status="UNAVAILABLE_ADAPTER_NOT_CONFIGURED", stale_sources=[])
        else:
            result["native_invoked"] = True; status, payload, diagnostic = _bounded_native_call(native_invoker)
            result.update(status=status, payload=payload, stale_sources=[], resource_envelope_enforced=True)
            if status == "COMPLETED_SUPPLEMENTAL_DIAGNOSTIC" and diagnostic:
                began = time.monotonic_ns()
                with (batch / "postsolve-diagnostic-only.payload").open("xb") as stream: stream.write(diagnostic)
                result["diagnostic_file_write_ns"] = time.monotonic_ns() - began
    _write_new(batch / "result.json", result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--frozen-root", required=True); parser.add_argument("--source-snapshot", required=True); parser.add_argument("--output", required=True)
    args = parser.parse_args(); print(json.dumps(run_controller(args.frozen_root, args.source_snapshot, args.output), sort_keys=True))


if __name__ == "__main__": main()
