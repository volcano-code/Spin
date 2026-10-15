"""Process-local memoization for repeatedly verified immutable file digests.

Cached SHA-256 values are reused only while the resolved regular file keeps the
same device, inode, size, mtime and ctime. Any observable replacement or write
forces a fresh content hash. This is intended for local POSIX filesystems with
reliable ctime semantics; it is an optimization, not an adversarial defense
against privileged filesystem mutation.
"""
from __future__ import annotations

from collections import OrderedDict
import hashlib
import json
import os
from pathlib import Path
import stat
import threading


def _metadata(info):
    return (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns)


class FingerprintMemo:
    """Bounded, process-local cache with measurable hits and content reads."""

    def __init__(self, maximum_entries=100_000):
        if type(maximum_entries) is not int or maximum_entries < 1:
            raise ValueError("INVALID_FINGERPRINT_CACHE_CAPACITY")
        self.maximum_entries = maximum_entries
        self._records = OrderedDict()
        self._lock = threading.RLock()
        self._hashes = 0
        self._cache_hits = 0
        self._bytes_hashed = 0

    @property
    def metrics(self):
        with self._lock:
            return {"hashes": self._hashes, "cache_hits": self._cache_hits,
                    "bytes_hashed": self._bytes_hashed}

    def digest(self, path, stream, before):
        key = str(path)
        observed = _metadata(before)
        with self._lock:
            entry = self._records.get(key)
            if entry is not None and entry[0] == observed:
                self._records.move_to_end(key)
                self._cache_hits += 1
                digest = entry[1]
                cached = True
            else:
                digest = None
                cached = False
        if not cached:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        after = os.fstat(stream.fileno())
        if _metadata(before) != _metadata(after):
            with self._lock:
                self._records.pop(key, None)
            raise ValueError("FILE_CHANGED_DURING_CHECKPOINT_HASH")
        with self._lock:
            self._records[key] = (_metadata(after), digest)
            self._records.move_to_end(key)
            if not cached:
                self._hashes += 1
                self._bytes_hashed += before.st_size
            while len(self._records) > self.maximum_entries:
                self._records.popitem(last=False)
        return digest


def fingerprints(files, memo):
    """Return the legacy fingerprint shape, reusing only stat-identical hashes."""
    if not isinstance(files, dict):
        raise ValueError("CHECKPOINT_FILE_MAP_REQUIRED")
    result = {}
    for label, path in files.items():
        if type(label) is not str or not label:
            raise ValueError("CHECKPOINT_FILE_LABEL_REQUIRED")
        anchored = Path(path).resolve(strict=True)
        fd = os.open(anchored, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC)
        try:
            before = os.fstat(fd)
            if not stat.S_ISREG(before.st_mode):
                raise ValueError("REGULAR_FILE_REQUIRED")
            with os.fdopen(fd, "rb", closefd=False) as stream:
                digest = memo.digest(anchored, stream, before)
            result[label] = {"path": str(anchored), "size": before.st_size, "sha256": digest}
        finally:
            os.close(fd)
    return result


def verify_artifacts(recorded, memo):
    """Verify checkpoint artifact records using the same fail-closed schema."""
    if not isinstance(recorded, dict):
        raise ValueError("INVALID_ARTIFACT_MAP")
    files = {}
    for label, value in recorded.items():
        if not isinstance(value, dict) or set(value) != {"path", "size", "sha256"}:
            raise ValueError("INVALID_ARTIFACT_RECORD")
        if type(value["path"]) is not str or not Path(value["path"]).is_absolute():
            raise ValueError("INVALID_ARTIFACT_PATH")
        files[label] = value["path"]
    actual = fingerprints(files, memo)
    if json.dumps(actual, sort_keys=True, separators=(",", ":"), allow_nan=False) != json.dumps(
            recorded, sort_keys=True, separators=(",", ":"), allow_nan=False):
        raise ValueError("STALE_CHECKPOINT_ARTIFACT")


def install_runtime_memo(base_module, checkpoint_module, supervisor_module, job_module=None, memo=None):
    """Patch existing source-bound modules only in an explicitly versioned runner."""
    memo = FingerprintMemo() if memo is None else memo

    def cached_fingerprints(files):
        return fingerprints(files, memo)

    def cached_artifact_verification(recorded):
        return verify_artifacts(recorded, memo)

    base_module._fingerprints = cached_fingerprints
    checkpoint_module._fingerprints = cached_fingerprints
    checkpoint_module._verify_artifacts = cached_artifact_verification
    supervisor_module._fingerprints = cached_fingerprints
    supervisor_module._verify_artifacts = cached_artifact_verification
    if job_module is not None:
        job_module._fingerprints = cached_fingerprints
    return memo
