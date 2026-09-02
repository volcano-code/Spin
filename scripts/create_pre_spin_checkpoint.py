#!/usr/bin/env python3
"""Create a deterministic evidence bundle for an uncommitted Git worktree."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Any, Iterable


def _git(root: Path, *arguments: str, text: bool = True) -> str | bytes:
    completed = subprocess.run(
        ["git", *arguments],
        cwd=root,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=text,
    )
    return completed.stdout.strip() if text else completed.stdout


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as destination:
            destination.write(payload)
            destination.flush()
            os.fsync(destination.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _is_below(path: Path, directory: Path) -> bool:
    try:
        path.relative_to(directory)
    except ValueError:
        return False
    return True


def _untracked_paths(root: Path, output: Path) -> Iterable[Path]:
    raw = _git(
        root,
        "ls-files",
        "--others",
        "--exclude-standard",
        "-z",
        text=False,
    )
    assert isinstance(raw, bytes)
    relative_paths = {
        Path(value.decode("utf-8"))
        for value in raw.split(b"\0")
        if value
    }
    tracked_raw = _git(root, "ls-files", "-z", text=False)
    assert isinstance(tracked_raw, bytes)
    tracked_paths = {
        Path(value.decode("utf-8"))
        for value in tracked_raw.split(b"\0")
        if value
    }
    legacy_experiments = root / "experiments" / "paper_sim_v1"
    if legacy_experiments.is_dir():
        for candidate in legacy_experiments.rglob("*"):
            relative = candidate.relative_to(root)
            if candidate.is_file() and relative not in tracked_paths:
                relative_paths.add(relative)
    for relative in sorted(relative_paths, key=lambda value: os.fsencode(value)):
        candidate = root / relative
        if _is_below(candidate.resolve(strict=False), output):
            continue
        if candidate.is_symlink():
            resolved = candidate.resolve(strict=True)
            try:
                resolved.relative_to(root)
            except ValueError as error:
                raise ValueError(
                    f"untracked symlink escapes repository: {relative}"
                ) from error
        if candidate.is_file():
            yield relative


def create_checkpoint(repository: Path, output: Path) -> dict[str, Any]:
    root = repository.resolve(strict=True)
    destination = output.resolve(strict=False)
    try:
        destination.relative_to(root)
    except ValueError as error:
        raise ValueError("checkpoint output must be inside the repository") from error
    if destination == root:
        raise ValueError("checkpoint output must not be the repository root")

    top_level = Path(str(_git(root, "rev-parse", "--show-toplevel"))).resolve()
    if top_level != root:
        raise ValueError("repository must identify the Git worktree root")

    patch = _git(root, "diff", "--binary", "--no-ext-diff", text=False)
    assert isinstance(patch, bytes)
    untracked = []
    for relative in _untracked_paths(root, destination):
        path = root / relative
        untracked.append({
            "path": relative.as_posix(),
            "size": path.stat().st_size,
            "sha256": _sha256_file(path),
        })

    manifest: dict[str, Any] = {
        "schema": "spin-it-pre-checkpoint-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git": {
            "branch": str(_git(root, "branch", "--show-current")),
            "head": str(_git(root, "rev-parse", "HEAD")),
            "user_name": str(_git(root, "config", "--local", "--get", "user.name")),
            "user_email": str(_git(root, "config", "--local", "--get", "user.email")),
        },
        "tracked_patch_sha256": hashlib.sha256(patch).hexdigest(),
        "untracked_files": untracked,
    }

    manifest_bytes = (
        json.dumps(
            manifest,
            ensure_ascii=False,
            allow_nan=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode("utf-8")
    patch_hash = hashlib.sha256(patch).hexdigest()
    manifest_hash = hashlib.sha256(manifest_bytes).hexdigest()
    sums = (
        f"{patch_hash}  tracked.patch\n"
        f"{manifest_hash}  manifest.json\n"
    ).encode("ascii")

    _atomic_write(destination / "tracked.patch", patch)
    _atomic_write(destination / "manifest.json", manifest_bytes)
    _atomic_write(destination / "SHA256SUMS", sums)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    manifest = create_checkpoint(arguments.repository, arguments.output)
    print(json.dumps(manifest, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
