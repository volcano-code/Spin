"""Build source-bound V-HACD contact proxies for the sealed G7 teapot shell."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import trimesh

ROOT = Path(__file__).resolve().parents[2]
BATCH = ROOT / "build/spin2014-g7-cached-fingerprint-repeat-20260920"
INPUT = ROOT / "build/teapot-input-20260913-run2"
SEAL_DIR = ROOT / "build/spin2014-g7-offline-trace-audit-20260920/geometry-seal-job-retry-01/payload"
SEAL_BUNDLE = SEAL_DIR / "sealed-artifact"
JOB = BATCH / "collision-proxy-job-retry-01"


def read(path):
    return json.loads(Path(path).read_text())


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def inputs():
    from make_it_stand.verification.source_checkpoint_store import _fingerprints

    seal_report_path = SEAL_DIR / "geometry-postprocess-report.json"
    determinism_path = BATCH / "determinism-report-v2.json"
    audit_path = BATCH / "trace-audit-v2-job-retry-03/payload/report.json"
    seal_report, determinism, audit = map(read, (seal_report_path, determinism_path, audit_path))
    if (determinism.get("status") != "DETERMINISTIC_SCIENTIFIC_OUTPUT_VERIFIED" or
            determinism.get("terminal_geometry_files_match_sealed_artifact") is not True or
            audit.get("status") != "NATIVE_BINARY_TRACE_AUDIT_VERIFIED" or
            seal_report.get("independent_seal", {}).get("status") != "MESH_MATERIAL_VALIDATED"):
        raise ValueError("G7_REPEAT_AND_SEAL_GATES_REQUIRED")
    problem_payload = read(SEAL_BUNDLE / "problem.json")
    input_manifest = read(INPUT / "manifest.json")
    if input_manifest.get("problem_digest") != problem_payload.get("input_digest"):
        raise ValueError("COLLISION_SOURCE_PROBLEM_DIGEST_MISMATCH")
    bundle = read(SEAL_BUNDLE / "manifest.json")
    outer_rec = bundle["final"]["outer"]
    sealed_vertices = np.load(SEAL_BUNDLE / outer_rec["vertices"]["name"], allow_pickle=False)
    sealed_faces = np.load(SEAL_BUNDLE / outer_rec["faces"]["name"], allow_pickle=False)
    source_outer = trimesh.load(INPUT / "outer.obj", process=False, force="mesh")
    if (not np.array_equal(sealed_vertices, np.asarray(source_outer.vertices)) or
            not np.array_equal(sealed_faces, np.asarray(source_outer.faces))):
        raise ValueError("PURE_HOLLOW_OUTER_SURFACE_DIFFERS_FROM_COLLISION_SOURCE")
    if bundle.get("scientific_digest") != seal_report["independent_seal"]["scientific_digest"]:
        raise ValueError("SEALED_ARTIFACT_DIGEST_MISMATCH")

    import pybullet
    import scipy
    from make_it_stand import pybullet_runner
    paths = {
        "collision_runner": ROOT / "experiments/teapot_spin2014/collisions.py",
        "collision_job_script": Path(__file__).resolve(),
        "teapot_helper": ROOT / "experiments/teapot_spin2014/teapot.py",
        "stem_geometry": ROOT / "experiments/heart_stem_contact/stem_geometry.py",
        "pybullet_runner": Path(pybullet_runner.__file__).resolve(),
        "pybullet_binary": Path(pybullet.__file__).resolve(),
        "trimesh_entry": Path(trimesh.__file__).resolve(),
        "scipy_entry": Path(scipy.__file__).resolve(),
        "input_manifest": INPUT / "manifest.json",
        "input_outer": INPUT / "outer.obj",
        "input_teapot": INPUT / "teapot.obj",
        "input_rod": INPUT / "rod.obj",
        "sealed_bundle_manifest": SEAL_BUNDLE / "manifest.json",
        "sealed_problem_payload": SEAL_BUNDLE / "problem.json",
        "sealed_geometry_report": seal_report_path,
        "native_trace_audit_report": audit_path,
        "determinism_report": determinism_path,
    }
    for path in sorted(SEAL_BUNDLE.glob("mesh_*.npy")):
        paths[f"sealed_{path.name}"] = path
    return {key: str(Path(value).resolve()) for key, value in paths.items()}


def launch():
    from make_it_stand.verification.source_job import run_source_job

    command = ["/usr/bin/env", "LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6",
               sys.executable, str(Path(__file__).resolve()), "worker", "{job}/payload/collision"]
    environment = {
        "PYTHONPATH": f"{ROOT / 'python'}:{ROOT / 'experiments/teapot_spin2014'}:{ROOT / 'experiments/heart_stem_contact'}",
    }
    return run_source_job(
        JOB, command, inputs(), maximum_wall_ns=3_600_000_000_000,
        maximum_memory_bytes=24 * 1024**3, maximum_disk_bytes=80 * 1024**3,
        minimum_start_free_bytes=100 * 1024**3, minimum_free_bytes=20 * 1024**3,
        environment=environment)


def worker(output):
    output = Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, str(ROOT / "experiments/teapot_spin2014/collisions.py"),
               "--input", str(INPUT), "--output", str(output)]
    completed = subprocess.run(command, env=dict(os.environ), check=False)
    if completed.returncode:
        raise RuntimeError(f"V_HACD_COLLISION_GENERATION_FAILED:{completed.returncode}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        worker(sys.argv[2])
    else:
        job = launch()
        print(json.dumps({"status": job.status, "returncode": job.returncode,
                          "elapsed_seconds": job.elapsed_ns / 1e9,
                          "memory_peak": job.service.get("MemoryPeak")}, sort_keys=True), flush=True)
