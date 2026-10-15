import hashlib
import json
from pathlib import Path
import sys
import pytest

sys.path.insert(0, str(Path(__file__).parent))
import postsolve_replay as replay
from postsolve_replay import (
    G7_STAGE19_CHECKPOINT_FILE_SHA256,
    G7_STAGE19_RECEIPT_BODY_SHA256,
    MAX_DIAGNOSTIC_BYTES,
    _bounded_native_call,
    _verify_spec_identity,
    run_controller,
)


class _FakeEndpoint:
    def __init__(self, wire=None, poll_error=None, close_error=None):
        self.wire = wire
        self.poll_error = poll_error
        self.close_error = close_error
        self.closed = False
        self.poll_timeouts = []

    def poll(self, timeout):
        self.poll_timeouts.append(timeout)
        if self.poll_error is not None:
            raise self.poll_error
        if self.wire is not None:
            return True
        replay.time.sleep(timeout)
        return False

    def recv_bytes(self, _maximum):
        return self.wire

    def close(self):
        self.closed = True
        if self.close_error is not None:
            raise self.close_error


class _FakeProcess:
    def __init__(self, *, exit_on_join=False, exit_after_kill=False):
        self.pid = 43210
        self.alive = True
        self.exit_on_join = exit_on_join
        self.exit_after_kill = exit_after_kill
        self.join_timeouts = []
        self.exitcode = None
        self.started = False
        self.terminated = False
        self.killed = False

    def is_alive(self):
        return self.alive

    def terminate(self):
        self.terminated = True

    def kill(self):
        self.killed = True

    def join(self, timeout):
        assert timeout is not None, "worker join must always be bounded"
        self.join_timeouts.append(timeout)
        if self.exit_on_join or (self.killed and self.exit_after_kill):
            self.alive = False
            self.exitcode = 0 if self.exit_on_join else -9


class _FakePipeContext:
    def __init__(self, receiver, sender):
        self.receiver = receiver
        self.sender = sender

    def Pipe(self, duplex):
        assert duplex is False
        return self.receiver, self.sender


def _install_fake_context(monkeypatch, process, receiver):
    sender = _FakeEndpoint()
    def spawn(_callback, _deadline):
        process.started = True
        sender.close()
        return receiver, process

    monkeypatch.setattr(replay, "_spawn_native_worker", spawn)
    return sender


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _minimal_snapshot(tmp_path):
    source = tmp_path / "source.cpp"
    source.write_text("frozen source\n")
    snapshot = tmp_path / "source-snapshot.json"
    snapshot.write_text(json.dumps({str(source): _sha(source)}))
    return source, snapshot


def test_stale_frozen_source_fails_closed_without_native_invocation(tmp_path):
    source, snapshot = _minimal_snapshot(tmp_path)
    frozen = tmp_path / "frozen"
    frozen.mkdir()
    calls = []
    source.write_text("current source differs\n")
    result = run_controller(frozen, snapshot, tmp_path / "new-batch", lambda: calls.append("native"))
    assert result["status"] == "UNAVAILABLE_STALE_SOURCE"
    assert not result["native_invoked"] and calls == []


def test_empty_or_malformed_snapshot_is_an_explicit_unavailable_result(tmp_path):
    frozen = tmp_path / "frozen"
    frozen.mkdir()
    snapshot = tmp_path / "source-snapshot.json"
    snapshot.write_text("[]")
    result = run_controller(frozen, snapshot, tmp_path / "new-batch")
    assert result["status"] == "UNAVAILABLE_SOURCE_SNAPSHOT"
    assert (tmp_path / "new-batch" / "result.json").is_file()


def test_compact_parent_drained_pipe_enforces_envelope_and_strips_bytes():
    status, payload, diagnostic = _bounded_native_call(lambda: {"native": "synthetic", "diagnostic_bytes": b"x"})
    assert status == "COMPLETED_SUPPLEMENTAL_DIAGNOSTIC"
    assert payload == {"native": "synthetic"}
    assert diagnostic == b"x"


def test_pipe_rejects_unbounded_diagnostic_payload():
    status, payload, diagnostic = _bounded_native_call(lambda: {"diagnostic_bytes": b"x" * (MAX_DIAGNOSTIC_BYTES + 1)})
    assert status == "UNAVAILABLE_NATIVE_FAILURE"
    assert payload.startswith("ValueError:") and diagnostic is None


def test_cleanup_escalates_from_terminate_to_kill_with_only_timed_joins():
    process = _FakeProcess(exit_after_kill=True)
    exited = replay._cleanup_worker(process, replay.time.monotonic() + 0.05)
    assert exited is True
    assert process.terminated is True
    assert process.killed is True
    assert process.join_timeouts
    assert all(timeout is not None and timeout >= 0 for timeout in process.join_timeouts)


def test_cleanup_reaps_worker_that_exits_during_kill_race():
    process = _FakeProcess()

    def exit_before_kill_signal():
        process.alive = False
        process.exitcode = 0
        raise ProcessLookupError("child exited before SIGKILL")

    process.kill = exit_before_kill_signal
    exited = replay._cleanup_worker(process, replay.time.monotonic() + 0.05)

    assert exited is True
    assert process.join_timeouts
    assert all(timeout is not None for timeout in process.join_timeouts)


def test_controller_timeout_reaps_sigterm_resistant_worker_and_closes_pipe(monkeypatch):
    process = _FakeProcess(exit_after_kill=True)
    receiver = _FakeEndpoint()
    sender = _install_fake_context(monkeypatch, process, receiver)
    monkeypatch.setattr(replay, "WALL_SECONDS", 0.03)

    status, payload, diagnostic = replay._bounded_native_call(lambda: {})

    assert (status, payload, diagnostic) == ("UNAVAILABLE_RESOURCE_LIMIT_WALL", None, None)
    assert process.terminated is True and process.killed is True
    assert process.join_timeouts and all(timeout is not None for timeout in process.join_timeouts)
    assert all(0 <= timeout <= 0.03 for timeout in process.join_timeouts)
    assert receiver.closed is True and sender.closed is True


def test_controller_reports_child_pid_when_worker_survives_sigkill(monkeypatch):
    process = _FakeProcess()
    receiver = _FakeEndpoint()
    sender = _install_fake_context(monkeypatch, process, receiver)
    monkeypatch.setattr(replay, "WALL_SECONDS", 0.03)

    status, detail, diagnostic = replay._bounded_native_call(lambda: {})

    assert status == "UNAVAILABLE_NATIVE_CLEANUP_FAILED"
    assert detail == {"child_pid": 43210}
    assert diagnostic is None
    assert process.terminated is True and process.killed is True
    assert all(0 <= timeout <= 0.03 for timeout in process.join_timeouts)
    assert receiver.closed is True and sender.closed is True


def test_valid_payload_is_accepted_only_after_bounded_worker_exit(monkeypatch):
    wire = json.dumps({"kind": "ok", "payload": {"native": "synthetic"},
                       "diagnostic_b64": "eA=="}).encode("utf-8")
    process = _FakeProcess(exit_on_join=True)
    receiver = _FakeEndpoint(wire)
    sender = _install_fake_context(monkeypatch, process, receiver)

    status, payload, diagnostic = replay._bounded_native_call(lambda: {})

    assert process.is_alive() is False
    assert process.join_timeouts and all(timeout is not None for timeout in process.join_timeouts)
    assert (status, payload, diagnostic) == (
        "COMPLETED_SUPPLEMENTAL_DIAGNOSTIC", {"native": "synthetic"}, b"x"
    )
    assert receiver.closed is True and sender.closed is True


def test_valid_payload_from_worker_reaped_after_work_deadline_is_wall_limited(monkeypatch):
    wire = json.dumps({"kind": "ok", "payload": {"native": "synthetic"},
                       "diagnostic_b64": "eA=="}).encode("utf-8")
    process = _FakeProcess(exit_after_kill=True)
    receiver = _FakeEndpoint(wire)
    sender = _install_fake_context(monkeypatch, process, receiver)
    monkeypatch.setattr(replay, "WALL_SECONDS", 0.04)

    status, payload, diagnostic = replay._bounded_native_call(lambda: {})

    assert status == "UNAVAILABLE_RESOURCE_LIMIT_WALL"
    assert payload is None and diagnostic is None
    assert process.terminated and process.killed
    assert all(0 <= timeout <= 0.04 for timeout in process.join_timeouts)
    assert receiver.closed and sender.closed


def test_process_start_failure_is_reported_and_closes_pipe(monkeypatch):
    receiver = _FakeEndpoint()
    sender = _FakeEndpoint()
    context = _FakePipeContext(receiver, sender)
    monkeypatch.setattr(replay.multiprocessing, "get_context", lambda method: context)

    def fail_fork():
        raise OSError("synthetic start failure")

    monkeypatch.setattr(replay.os, "fork", fail_fork)

    status, payload, diagnostic = replay._bounded_native_call(lambda: {})

    assert (status, payload, diagnostic) == ("UNAVAILABLE_NATIVE_FAILURE", None, None)
    assert receiver.closed and sender.closed


def test_controller_exception_reaps_started_worker_and_closes_pipe(monkeypatch):
    process = _FakeProcess(exit_after_kill=True)
    receiver = _FakeEndpoint(poll_error=RuntimeError("synthetic receiver failure"))
    sender = _install_fake_context(monkeypatch, process, receiver)
    monkeypatch.setattr(replay, "WALL_SECONDS", 0.04)

    status, payload, diagnostic = replay._bounded_native_call(lambda: {})

    assert (status, payload, diagnostic) == ("UNAVAILABLE_NATIVE_FAILURE", None, None)
    assert process.started and process.terminated and process.killed
    assert receiver.closed and sender.closed


def test_fork_worker_captures_pid_and_closes_parent_sender(monkeypatch):
    receiver, sender = _FakeEndpoint(), _FakeEndpoint()
    context = _FakePipeContext(receiver, sender)
    monkeypatch.setattr(replay.multiprocessing, "get_context", lambda method: context)
    monkeypatch.setattr(replay.os, "fork", lambda: 56789)

    child_receiver, process = replay._spawn_native_worker(lambda: {}, replay.time.monotonic() + 1.0)

    assert child_receiver is receiver
    assert process.pid == 56789
    assert sender.closed is True
    assert receiver.closed is False


def test_fork_worker_closes_pipe_when_fork_fails(monkeypatch):
    receiver, sender = _FakeEndpoint(), _FakeEndpoint()
    context = _FakePipeContext(receiver, sender)
    monkeypatch.setattr(replay.multiprocessing, "get_context", lambda method: context)

    def fail_fork():
        raise OSError("synthetic fork failure")

    monkeypatch.setattr(replay.os, "fork", fail_fork)
    with pytest.raises(OSError, match="synthetic fork failure"):
        replay._spawn_native_worker(lambda: {}, replay.time.monotonic() + 1.0)

    assert receiver.closed and sender.closed


def test_fork_worker_attempts_both_closes_and_preserves_fork_error(monkeypatch):
    receiver = _FakeEndpoint(close_error=RuntimeError("receiver close failure"))
    sender = _FakeEndpoint()
    context = _FakePipeContext(receiver, sender)
    monkeypatch.setattr(replay.multiprocessing, "get_context", lambda method: context)

    def fail_fork():
        raise OSError("original fork failure")

    monkeypatch.setattr(replay.os, "fork", fail_fork)
    with pytest.raises(OSError, match="original fork failure"):
        replay._spawn_native_worker(lambda: {}, replay.time.monotonic() + 1.0)

    assert receiver.closed and sender.closed


def test_parent_receiver_closes_when_sender_close_fails(monkeypatch):
    receiver = _FakeEndpoint()
    sender = _FakeEndpoint(close_error=OSError("sender close failure"))
    context = _FakePipeContext(receiver, sender)
    monkeypatch.setattr(replay.multiprocessing, "get_context", lambda method: context)
    monkeypatch.setattr(replay.os, "fork", lambda: 56789)
    monkeypatch.setattr(replay, "_cleanup_worker", lambda process, deadline: True)

    with pytest.raises(OSError, match="sender close failure"):
        replay._spawn_native_worker(lambda: {}, replay.time.monotonic() + 1.0)

    assert receiver.closed and sender.closed


def test_spawn_cleanup_failure_keeps_child_pid_and_explicit_status(monkeypatch):
    receiver = _FakeEndpoint()
    sender = _FakeEndpoint(close_error=OSError("sender close failure"))
    context = _FakePipeContext(receiver, sender)
    monkeypatch.setattr(replay.multiprocessing, "get_context", lambda method: context)
    monkeypatch.setattr(replay.os, "fork", lambda: 56789)
    monkeypatch.setattr(replay, "_cleanup_worker", lambda process, deadline: False)

    status, detail, diagnostic = replay._bounded_native_call(lambda: {}, _wall_seconds=0.1)

    assert status == "UNAVAILABLE_NATIVE_CLEANUP_FAILED"
    assert detail == {"child_pid": 56789}
    assert diagnostic is None
    assert receiver.closed and sender.closed


@pytest.mark.parametrize(
    ("wait_result", "wait_error", "alive", "exitcode"),
    [
        ((0, 0), None, True, None),
        ((43210, 0), None, False, 0),
        ((43210, 9), None, False, -9),
        (None, InterruptedError(), True, None),
        (None, ChildProcessError(), False, -1),
    ],
)
def test_fork_worker_waitpid_state_handling(monkeypatch, wait_result, wait_error, alive, exitcode):
    process = replay._ForkWorker(43210)

    def fake_waitpid(pid, options):
        assert pid == 43210
        assert options == replay.os.WNOHANG
        if wait_error is not None:
            raise wait_error
        return wait_result

    monkeypatch.setattr(replay.os, "waitpid", fake_waitpid)
    assert process.is_alive() is alive
    assert process.exitcode == exitcode


def test_labeled_spec_identity_uses_canonical_fact_path_not_mapping_key(tmp_path):
    source = (tmp_path / "native-module.so").resolve()
    source.write_bytes(b"frozen native module")
    frozen = tmp_path / "frozen"
    frozen.mkdir()
    (frozen / "spec.json").write_text(json.dumps({"identity": {
        "native_module": {"path": str(source), "sha256": _sha(source), "size": source.stat().st_size}
    }}))
    assert _verify_spec_identity(frozen) is None


def test_real_g7_receipt_body_and_checkpoint_file_hash_are_distinct():
    root = Path(__file__).resolve().parents[2] / "build/spin2014-g7-linear-candidate-20260919"
    receipt = root / "repeat-0/payload/checkpoints/commit-00000019.sha256"
    checkpoint = root / "repeat-0/payload/checkpoints/stage-00000019.ckpt"
    assert receipt.read_text(encoding="ascii").strip().split() == ["00000019", G7_STAGE19_RECEIPT_BODY_SHA256]
    assert _sha(checkpoint) == G7_STAGE19_CHECKPOINT_FILE_SHA256
    assert G7_STAGE19_RECEIPT_BODY_SHA256 != G7_STAGE19_CHECKPOINT_FILE_SHA256


def test_21833_element_compact_transport_is_accepted_while_bulk_is_rejected():
    count = 21833
    compact = {"terminal_state_dimension": count, "terminal_state_sha256": "a" * 64,
               "terminal_kkt": {"alpha_dimension": count, "alpha_sha256": "b" * 64,
                                "eta_dimension": count, "eta_sha256": "c" * 64}}
    status, payload, diagnostic = _bounded_native_call(lambda: compact)
    assert status == "COMPLETED_SUPPLEMENTAL_DIAGNOSTIC" and payload == compact and diagnostic == b""
    bulk = {"terminal_state": [1.0] * count,
            "terminal_kkt": {"lambda": [1.0] * count, "alpha": [1.0] * count, "eta": [1.0] * count}}
    status, payload, diagnostic = _bounded_native_call(lambda: bulk)
    assert status == "UNAVAILABLE_NATIVE_FAILURE" and payload.startswith("ValueError:") and diagnostic is None
