#!/usr/bin/env python3
"""Compare old/new synthetic checkpoint roundtrip costs without G7 replay."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import time


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def invoke(binary: Path, mode: str, directory: Path, snapshot: Path | None = None):
    command = [str(binary), mode, "0", str(directory)]
    if snapshot is not None:
        command.append(str(snapshot))
    started = time.perf_counter_ns()
    result = subprocess.run(
        command,
        capture_output=True,
        timeout=30,
        env=dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1"),
    )
    return result, time.perf_counter_ns() - started


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--repeats", type=int, default=12)
    args = parser.parse_args()
    baseline = args.baseline.resolve(strict=True)
    candidate = args.candidate.resolve(strict=True)
    output = args.output.absolute()
    if args.repeats < 5:
        raise ValueError("at least five paired rounds are required")
    output.mkdir(parents=True, exist_ok=False)
    fixture = output / "synthetic-fixture"
    fixture.mkdir()
    attempt = fixture / "baseline-attempt"
    started, _ = invoke(baseline, "start", attempt)
    if started.returncode != 0 or started.stdout[:1] != b"P":
        raise RuntimeError(f"baseline failed to produce a synthetic checkpoint: {started.stderr!r}")
    snapshot = fixture / "checkpoint.bin"
    snapshot.write_bytes(started.stdout[1:])
    traces = sorted(attempt.glob("*.slqptrace"))
    if not traces:
        raise RuntimeError("synthetic checkpoint did not create a sealed trace")
    trace_hashes_before = {str(p): sha256(p) for p in traces}

    rows = []
    for repetition in range(args.repeats):
        labels = ["baseline", "candidate"] if repetition % 2 == 0 else ["candidate", "baseline"]
        for label in labels:
            binary = baseline if label == "baseline" else candidate
            output_dir = fixture / f"unused-{label}-{repetition:02d}"
            result, elapsed = invoke(binary, "roundtrip", output_dir, snapshot)
            if result.returncode != 0 or result.stdout != b"roundtrip verified\n":
                raise RuntimeError(f"{label} roundtrip failed: {result.stderr!r} {result.stdout!r}")
            if output_dir.exists():
                raise RuntimeError(f"roundtrip created an output directory before validation: {output_dir}")
            (output / f"{label}-{repetition:02d}.stdout").write_bytes(result.stdout)
            (output / f"{label}-{repetition:02d}.stderr").write_bytes(result.stderr)
            rows.append({"repetition": repetition, "binary": label, "elapsed_ns": elapsed})

    trace_hashes_after = {str(p): sha256(p) for p in traces}
    if trace_hashes_before != trace_hashes_after:
        raise RuntimeError("read-only roundtrip modified the source trace")
    by_name = {
        name: [row["elapsed_ns"] for row in rows if row["binary"] == name]
        for name in ("baseline", "candidate")
    }
    medians = {name: statistics.median(values) for name, values in by_name.items()}
    report = {
        "scope": "SYNTHETIC_CHECKPOINT_DECODE_AND_REENCODE_ONLY_NOT_G7",
        "source_run_or_optimizer_resume": "NOT_RUN",
        "fixture_created_by": "baseline executable; one synthetic stage only",
        "snapshot_sha256": sha256(snapshot),
        "trace_hashes_before_and_after_equal": True,
        "trace_sha256": trace_hashes_before,
        "binary_sha256": {"baseline": sha256(baseline), "candidate": sha256(candidate)},
        "roundtrip_output_identical": all((output / f"{name}-{i:02d}.stdout").read_bytes() == b"roundtrip verified\n"
            for name in ("baseline", "candidate") for i in range(args.repeats)),
        "timing_protocol": "fresh-process end-to-end wall clock under OMP/BLAS thread-limiting environment; includes process startup, decode, and checkpoint re-encode; alternating order",
        "threading_note": "thread-limit environment variables are requested controls, not an independent proof of every dependency's thread count",
        "paired_rounds": args.repeats,
        "elapsed_ns": by_name,
        "median_elapsed_ns": medians,
        "candidate_median_change_percent": (medians["candidate"] / medians["baseline"] - 1.0) * 100.0,
    }
    (output / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
