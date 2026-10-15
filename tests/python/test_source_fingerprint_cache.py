"""Metadata memoization must preserve source/artifact change detection."""
import importlib
import os
from pathlib import Path
from types import SimpleNamespace

import pytest


@pytest.fixture
def api():
    return importlib.import_module("make_it_stand.verification.source_fingerprint_cache")


def count_hashes(monkeypatch, api):
    actual = api.hashlib.file_digest
    calls = []

    def counted(stream, digest):
        calls.append(stream.name)
        return actual(stream, digest)

    monkeypatch.setattr(api.hashlib, "file_digest", counted)
    return calls


def test_unchanged_regular_file_reuses_verified_digest(api, tmp_path, monkeypatch):
    path = tmp_path / "trace.bin"
    path.write_bytes(b"sealed trace bytes")
    memo = api.FingerprintMemo()
    calls = count_hashes(monkeypatch, api)
    first = api.fingerprints({"trace": path}, memo)
    second = api.fingerprints({"trace": path}, memo)
    assert first == second
    assert len(calls) == 1
    assert memo.metrics == {"hashes": 1, "cache_hits": 1, "bytes_hashed": len(b"sealed trace bytes")}


def test_same_size_write_with_restored_mtime_is_rehashed(api, tmp_path, monkeypatch):
    path = tmp_path / "source.bin"
    path.write_bytes(b"alpha")
    memo = api.FingerprintMemo()
    calls = count_hashes(monkeypatch, api)
    before = api.fingerprints({"source": path}, memo)
    old = path.stat()
    path.write_bytes(b"bravo")
    os.utime(path, ns=(old.st_atime_ns, old.st_mtime_ns))
    after = api.fingerprints({"source": path}, memo)
    assert before["source"]["size"] == after["source"]["size"]
    assert before["source"]["sha256"] != after["source"]["sha256"]
    assert len(calls) == 2


def test_atomic_replacement_does_not_reuse_old_inode_digest(api, tmp_path, monkeypatch):
    path = tmp_path / "module.so"
    replacement = tmp_path / "replacement.so"
    path.write_bytes(b"same content")
    replacement.write_bytes(b"same content")
    memo = api.FingerprintMemo()
    calls = count_hashes(monkeypatch, api)
    first = api.fingerprints({"module": path}, memo)
    os.replace(replacement, path)
    second = api.fingerprints({"module": path}, memo)
    assert first == second
    assert len(calls) == 2


def test_cached_artifact_verification_rehashes_after_mutation(api, tmp_path, monkeypatch):
    path = tmp_path / "pass.slqptrace"
    path.write_bytes(b"original")
    memo = api.FingerprintMemo()
    calls = count_hashes(monkeypatch, api)
    recorded = api.fingerprints({"pass": path}, memo)
    api.verify_artifacts(recorded, memo)
    assert len(calls) == 1
    old = path.stat()
    path.write_bytes(b"mutated!")
    os.utime(path, ns=(old.st_atime_ns, old.st_mtime_ns))
    with pytest.raises(ValueError, match="STALE_CHECKPOINT_ARTIFACT"):
        api.verify_artifacts(recorded, memo)
    assert len(calls) == 2


def test_runtime_installer_patches_all_checkpoint_call_sites(api, tmp_path, monkeypatch):
    base = SimpleNamespace(_fingerprints=lambda files: None)
    checkpoint = SimpleNamespace(_fingerprints=lambda files: None,
                                 _verify_artifacts=lambda rows: None)
    supervisor = SimpleNamespace(_fingerprints=lambda files: None,
                                 _verify_artifacts=lambda rows: None)
    job = SimpleNamespace(_fingerprints=lambda files: None)
    original = [(obj, key, getattr(obj, key)) for obj, keys in (
        (base, ("_fingerprints",)),
        (checkpoint, ("_fingerprints", "_verify_artifacts")),
        (supervisor, ("_fingerprints", "_verify_artifacts")),
        (job, ("_fingerprints",))) for key in keys]
    for obj, key, value in original:
        monkeypatch.setattr(obj, key, value)

    memo = api.install_runtime_memo(base, checkpoint, supervisor, job_module=job)
    path = tmp_path / "input.bin"
    path.write_bytes(b"input")
    first = base._fingerprints({"input": path})
    assert checkpoint._fingerprints({"input": path}) == first
    checkpoint._verify_artifacts(first)
    assert supervisor._fingerprints({"input": path}) == first
    assert job._fingerprints({"input": path}) == first
    assert memo.metrics["hashes"] == 1
    assert memo.metrics["cache_hits"] == 4
