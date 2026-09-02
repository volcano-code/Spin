#!/usr/bin/env python3
"""Run and seal the four deterministic Spin-It Faster analytic fixtures."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any


CASE_IDS = (
    "static-planar",
    "fidget-circular-void",
    "yoyo-symmetric",
    "top-symmetric",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def git_value(repository: Path, *arguments: str) -> str:
    return subprocess.run(
        ["git", *arguments],
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def run_once(executable: Path) -> str:
    completed = subprocess.run(
        [str(executable)],
        check=True,
        capture_output=True,
        text=True,
    )
    if completed.stderr:
        raise RuntimeError(f"analytic executable wrote stderr: {completed.stderr}")
    return completed.stdout


def parse_hashes(output: str) -> list[str]:
    prefix = "[PASS] scenario=spin-it-faster-end-to-end "
    if not output.startswith(prefix) or not output.endswith("\n"):
        raise RuntimeError("analytic executable did not emit its PASS contract")
    hashes = output[len(prefix):].strip().split()
    if len(hashes) != len(CASE_IDS):
        raise RuntimeError("analytic executable emitted the wrong case count")
    if any(len(value) != 64 or set(value) - set("0123456789abcdef") for value in hashes):
        raise RuntimeError("analytic executable emitted an invalid SHA-256")
    return hashes


def main() -> int:
    script_directory = Path(__file__).resolve().parent
    repository = script_directory.parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository", type=Path, default=repository)
    parser.add_argument(
        "--executable",
        type=Path,
        default=repository / "build/spin-faster/spin-it-faster-end-to-end",
    )
    parser.add_argument("--output", type=Path, default=script_directory / "analytic-results")
    arguments = parser.parse_args()

    repository = arguments.repository.expanduser().resolve()
    executable = arguments.executable.expanduser().resolve()
    output_directory = arguments.output.expanduser().resolve()
    manifest_path = script_directory / "analytic_cases.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    observed_ids = tuple(case["id"] for case in manifest["cases"])
    if observed_ids != CASE_IDS:
        raise RuntimeError("analytic manifest order differs from native fixture order")
    if not executable.is_file():
        raise FileNotFoundError(f"missing analytic executable: {executable}")

    first = run_once(executable)
    second = run_once(executable)
    if first != second:
        raise RuntimeError("analytic repeats are not byte-identical")
    hashes = parse_hashes(first)

    diff = git_value(repository, "status", "--porcelain=v1", "--untracked-files=all")
    summary = {
        "schema": "spin-it-faster-analytic-results-v1",
        "schema_version": 1,
        "claim_level": manifest["claim_level"],
        "git_commit": git_value(repository, "rev-parse", "HEAD"),
        "worktree_clean": not bool(diff),
        "worktree_status_sha256": hashlib.sha256(diff.encode("utf-8")).hexdigest(),
        "manifest_sha256": sha256(manifest_path),
        "executable_sha256": sha256(executable),
        "stdout_sha256": hashlib.sha256(first.encode("utf-8")).hexdigest(),
        "repeat_count": 2,
        "byte_identical_repeats": True,
        "all_cases_passed": True,
        "equation_residual_limit": manifest["paper_equation_residual_limit"],
        "physical_relative_residual_limit": manifest["physical_relative_residual_limit"],
        "cases": [
            {
                "id": case_id,
                "output_hash": output_hash,
                "status": "CONVERGED",
            }
            for case_id, output_hash in zip(CASE_IDS, hashes)
        ],
    }
    output_directory.mkdir(parents=True, exist_ok=True)
    write_json(output_directory / "summary.json", summary)
    (output_directory / "stdout.txt").write_text(first, encoding="utf-8")
    print(first, end="")
    print(f"sealed={output_directory / 'summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
