#!/usr/bin/env python3
"""Read-only audit of the frozen G7 v2 prefix; never invokes an optimizer."""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT / "python"), str(ROOT / "experiments" / "teapot_spin2014")]

from g7_transition_lineage import analyze_transition, fractional_counts_by_level
from make_it_stand.verification.source_checkpoint_store import _canonical, _fingerprints
from make_it_stand.verification.source_checkpoint_store_v2 import StageCheckpointStore

FROZEN = ROOT / "build" / "spin2014-v2-g7-20260919"
TRACE_READER = ROOT / "build" / "spin2014-g7-staged-20260916" / "native_trace_summary"

def read(path): return json.loads(Path(path).read_text())
def sha(path):
    with Path(path).open("rb") as stream: return hashlib.file_digest(stream, "sha256").hexdigest()

def readonly_store(directory, files):
    """The existing constructor may create directories; this reader never does."""
    store = object.__new__(StageCheckpointStore); store.directory = Path(directory)
    store._identity_files = {key: str(Path(value).absolute()) for key, value in files.items()}
    store._identity = _fingerprints(files)
    store._limits = dict(maximum_wall_ns=3_540_000_000_000, maximum_trace_bytes=256 * 1024 ** 3,
                         maximum_payload_bytes=64 * 1024 ** 2, maximum_identity_bytes=16 * 1024 ** 2)
    manifest = _canonical(dict(schema="source-stage-identity-v2", files=store._identity)); digest = hashlib.sha256(manifest).hexdigest()
    store._manifest_bytes = manifest
    store._manifest_reference = dict(name=f"identity-{digest}.json", sha256=digest, size=len(manifest), file_count=len(store._identity))
    return store

def verify_outer_receipts(run):
    job = run / "repeat-0"
    for path in (job / "config.json", job / "result.json", job / "complete.json"):
        if not path.is_file(): raise ValueError("MISSING_OUTER_RECEIPT:" + str(path))
    if read(job / "complete.json") != {"result_sha256": sha(job / "result.json")}: raise ValueError("OUTER_RECEIPT_HASH_MISMATCH")
    result = read(job / "result.json")
    if result.get("status") != "RESOURCE_LIMIT_WALL" or result.get("returncode") != -1: raise ValueError("UNEXPECTED_OUTER_JOB_STATUS")
    return job, result

def verify_source_and_checkpoints(run, job):
    spec_path = run / "spec.json"; spec = read(spec_path)
    if _fingerprints(spec["files"]) != spec["identity"]: raise ValueError("STALE_SOURCE_IDENTITY")
    files = dict(spec["files"], spec=str(spec_path)); store = readonly_store(job / "payload" / "checkpoints", files)
    # The immutable source carrier is mounted read-only; its existing lock is not
    # opened with O_CREAT. _latest_locked performs the full validation without a
    # write-capable lock acquisition, exactly as the frozen v2 verifier does.
    latest = store._latest_locked()  # validates every boundary, receipt, parent, artifact and identity
    if latest is None or latest.stage != 18: raise ValueError("EXPECTED_NINETEEN_CHECKPOINTS")
    return spec

def trace_paths(job):
    paths = []
    for index in range(19):
        path = job / "payload" / "attempts" / f"attempt-{index:08d}" / "payload" / "trace" / f"pass-{index:06d}.slqptrace"
        if not path.is_file(): raise ValueError("MISSING_COMPLETE_TRACE:" + str(path))
        paths.append(path)
    partial = job / "payload" / "attempts" / "attempt-00000019"
    if not partial.is_dir() or (partial / "result.json").exists(): raise ValueError("STAGE_TWENTY_PARTIAL_EVIDENCE_EXPECTED")
    return paths

def verify_complete_traces(paths):
    if not TRACE_READER.is_file(): raise ValueError("MISSING_FROZEN_TRACE_READER")
    rows = []
    for index, path in enumerate(paths, start=1):
        value = json.loads(subprocess.check_output([str(TRACE_READER), str(path)], text=True, timeout=180))
        if value.get("integrity") != "ALL_NATIVE_RECORD_CHECKSUMS_VERIFIED" or value.get("final_state") is not True: raise ValueError("INCOMPLETE_OR_UNVERIFIED_TRACE")
        rows.append(dict(stage=index, trace_sha256=sha(path), records=value["records"], final_iteration=value["iteration"], terminal_qp=value["qp_status"]))
    return rows

def _header(job, index):
    value = read(job / "payload" / "attempts" / f"attempt-{index:08d}" / "payload" / "header.json")
    if value.get("kind") != "header" or not isinstance(value.get("payload"), dict): raise ValueError("INVALID_STAGE_HEADER")
    return value["payload"]

def lineage_and_costs(job):
    headers = [_header(job, index) for index in range(19)]; epsilon = headers[0]["settings"]["solver_settings"]["fractional_epsilon"]
    rows = []
    for index, header in enumerate(headers):
        solve = header["solve"]; cells, beta = solve["final_partition"], solve["relaxed_fill"]
        if len(cells) != len(beta): raise ValueError("OUTPUT_PARTITION_FILL_MISMATCH")
        native_seconds = read(job / "payload" / "attempts" / f"attempt-{index:08d}" / "payload" / "stage.json")["native_seconds"]
        cost = dict(stage_worker_seconds=native_seconds, native_call_seconds=native_seconds, whole_job_seconds="AVAILABLE_AT_OUTER_RECEIPT_ONLY",
                    objective_derivative_seconds="UNAVAILABLE", qp_seconds="UNAVAILABLE", line_search_seconds="UNAVAILABLE",
                    split_merge_seconds="UNAVAILABLE", restore_seconds="UNAVAILABLE", trace_io_seconds="UNAVAILABLE",
                    replay_only_qp_seconds="NOT_RUN_NO_AUTHENTICATED_PER_ITERATION_QP_SELECTOR")
        if index == 0:
            row = dict(stage=1, input_partition=dict(recorded_hash=solve["input_hash"], verification="AUTHENTICATED_SOURCE_IDENTITY_INITIAL_PARTITION_SNAPSHOT_UNAVAILABLE"),
                       output_partition=dict(recorded_hash=solve["final_partition_hash"], verification="AUTHENTICATED_RECORD_ONLY_NATIVE_HASH"),
                       transition_evidence="UNAVAILABLE_INITIAL_BOUNDARY_NOT_SNAPSHOTTED", split=dict(parent_to_children={}, operation_order="UNAVAILABLE"),
                       merge=dict(parent_to_children={}, operation_order="UNAVAILABLE"),
                       fractional_counts_by_level=dict(input="UNAVAILABLE", output=fractional_counts_by_level(cells, beta, epsilon)),
                       moment_conservation_delta="UNAVAILABLE_INITIAL_PARTITION_SNAPSHOT",
                       moment_conservation_quantity="PARTITION_CELL_GEOMETRY_ONLY_NOT_MATERIAL_OR_CONSTRAINT",
                       historical_subtimers="UNAVAILABLE")
        else:
            prior = headers[index - 1]["solve"]
            record = dict(before=prior["final_partition"], after=cells, before_hash=prior["final_partition_hash"], after_hash=solve["final_partition_hash"],
                          hash_scheme="RECORDED_NATIVE_HASH", complete=True, stage=index + 1)
            row = analyze_transition(record, epsilon=epsilon, before_beta=prior["relaxed_fill"], after_beta=beta); row["stage"] = index + 1
        row["cost_accounting"] = cost; row["repeated_lineage"] = "UNAVAILABLE_NO_AUTHENTICATED_OPERATION_LOG"; rows.append(row)
    slowest = sorted(rows, key=lambda row: (-row["cost_accounting"]["native_call_seconds"], row["stage"]))[:2]
    return rows, [row["stage"] for row in slowest]

def diagnose(run, *, verify_traces=True):
    run = Path(run).resolve(); job, outer = verify_outer_receipts(run); verify_source_and_checkpoints(run, job); paths = trace_paths(job)
    traces = verify_complete_traces(paths) if verify_traces else "NOT_VERIFIED_BY_REQUEST"; stages, slowest = lineage_and_costs(job)
    return dict(schema="spin2014-g7-postrun-diagnostic-v1", scope="READ_ONLY_FROZEN_PREFIX_NO_NONLINEAR_SOLVE", source_run=str(run),
                outer_job_receipts="VERIFIED", source_identity="VERIFIED", checkpoint_boundaries="VERIFIED_19", complete_slqp_traces=traces,
                excluded_partial_stage=20, whole_job_seconds=outer["elapsed_ns"] / 1e9, mechanically_selected_slowest_stages=slowest,
                fixed_qp_replay=dict(status="REPORTED_SEPARATELY", detail="The controller verifies immutable source evidence only. The eight capped fixed-QP calls, hash-bound temporary build provenance, science equality, and KKT results are recorded in REPORT.md; no replay output is written into the frozen source run."), stages=stages)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("run", nargs="?", type=Path, default=FROZEN); parser.add_argument("--skip-trace-verification", action="store_true")
    args = parser.parse_args(); print(json.dumps(diagnose(args.run, verify_traces=not args.skip_trace_verification), indent=2, allow_nan=False))
