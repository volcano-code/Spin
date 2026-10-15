"""Bounded free-table contact trials for the sealed reconstructed G7 teapot."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import trimesh

ROOT = Path(__file__).resolve().parents[2]
BATCH = ROOT / "build/spin2014-g7-cached-fingerprint-repeat-20260920"
INPUT = ROOT / "build/teapot-input-20260913-run2"
SEAL_DIR = ROOT / "build/spin2014-g7-offline-trace-audit-20260920/geometry-seal-job-retry-01/payload"
SEAL_BUNDLE = SEAL_DIR / "sealed-artifact"
COLLISION_JOB = BATCH / "collision-proxy-job-retry-01"
COLLISION_DIR = COLLISION_JOB / "payload/collision"
IDEAL_JOB = BATCH / "ideal-dynamics-job-retry-02"
JOB = BATCH / "free-contact-job-v1"


def read(path):
    return json.loads(Path(path).read_text())


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_new(path, value):
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()


def inputs():
    g7_spec = read(BATCH / "spec.json")
    collision_manifest_path = COLLISION_DIR / "manifest.json"
    collision = read(collision_manifest_path)
    bundle = read(SEAL_BUNDLE / "manifest.json")
    problem = read(SEAL_BUNDLE / "problem.json")
    dynamics = read(IDEAL_JOB / "payload/report.json")
    if (collision.get("input_digest") != problem.get("input_digest") or
            bundle.get("scientific_digest") != dynamics.get("sealed_artifact", {}).get("scientific_digest") or
            dynamics.get("sealed_artifact", {}).get("status") != "MESH_MATERIAL_VALIDATED" or
            dynamics.get("contact_simulation") != "NOT_RUN_COLLISION_PROXY_VALIDATION_REMAINS_NEXT_GATE"):
        raise ValueError("SEALED_IDEAL_AND_COLLISION_INPUTS_MUST_MATCH")
    from teapot import check_files
    check_files(collision["outer_hash"])
    check_files(collision["source_hashes"])
    for profile in collision["profiles"].values():
        check_files(profile["hashes"])
        if not profile["near_upright_passed"]:
            raise ValueError("V_HACD_PROXY_FAILED_STATIC_SUPPORT_ORACLE")
        near = [abs(item["error_m"]) for item in profile["support_errors"] if item["tilt"] <= 5]
        if not near or max(near) > 0.0001:
            raise ValueError("V_HACD_PROXY_STATIC_SUPPORT_ERROR_EXCEEDS_0_1MM")
    import pybullet
    paths = {
        "contact_script": Path(__file__).resolve(),
        "contact_integrator": ROOT / "experiments/heart_stem_contact/bullet_contact.py",
        "mass_integrator": ROOT / "experiments/reconstructed_heart_top/heart_top.py",
        "collision_builder": ROOT / "experiments/teapot_spin2014/collisions.py",
        "collision_manifest": collision_manifest_path,
        "collision_runner": ROOT / "python/make_it_stand/pybullet_runner.py",
        "pybullet_binary": Path(pybullet.__file__).resolve(),
        "source_helper": ROOT / "experiments/teapot_spin2014/teapot.py",
        "input_manifest": INPUT / "manifest.json",
        "input_outer": INPUT / "outer.obj",
        "input_teapot": INPUT / "teapot.obj",
        "input_rod": INPUT / "rod.obj",
        "sealed_geometry_report": SEAL_DIR / "geometry-postprocess-report.json",
        "sealed_bundle_manifest": SEAL_BUNDLE / "manifest.json",
        "sealed_problem": SEAL_BUNDLE / "problem.json",
        "ideal_dynamics_report": IDEAL_JOB / "payload/report.json",
        "collision_job_result": COLLISION_JOB / "result.json",
        "collision_job_config": COLLISION_JOB / "config.json",
        "collision_job_complete": COLLISION_JOB / "complete.json",
    }
    for item in sorted(SEAL_BUNDLE.glob("mesh_*.npy")):
        paths[f"sealed_{item.name}"] = item
    for item in sorted(COLLISION_DIR.rglob("*")):
        if item.is_file():
            paths[f"collision_{item.relative_to(COLLISION_DIR).as_posix()}"] = item
    return {key: str(Path(value).resolve()) for key, value in paths.items()}, collision


def worker(output: Path):
    from heart_top import integrate
    from bullet_contact import material_properties, simulate_contact

    paths, collision = inputs()
    bundle = read(SEAL_BUNDLE / "manifest.json")
    mat_record = bundle["final"]["material"]
    vertices = np.load(SEAL_BUNDLE / mat_record["vertices"]["name"], allow_pickle=False)
    faces = np.load(SEAL_BUNDLE / mat_record["faces"]["name"], allow_pickle=False)
    physical = material_properties(vertices * .025, faces, 1200.)
    outer_mesh = trimesh.load(collision["outer"], force="mesh", process=False)
    outer_vertices = np.asarray(outer_mesh.vertices, dtype=float)
    if sha256(collision["outer"]) != collision["outer_hash"][str(Path(collision["outer"]).resolve())]:
        raise ValueError("OUTER_VISUAL_SOURCE_HASH_MISMATCH")

    runs = []
    for profile in ("default", "robust"):
        for direction in range(8):
            print(f"CONTACT profile={profile} direction={direction}", flush=True)
            record = simulate_contact(
                Path(collision["outer"]), collision["profiles"][profile]["files"], physical,
                outer_vertices, hz=35., tilt=2., azimuth=direction * np.pi / 4,
                duration=8., dt=1 / 2000., clearance=.0002, gui=False)
            record.update({"profile": profile, "direction": direction,
                           "azimuth_deg": direction * 45., "artifact_digest": bundle["scientific_digest"],
                           "outer_collision_manifest_sha256": sha256(COLLISION_DIR / "manifest.json")})
            runs.append(record)

    nominal = {run["profile"]: run for run in runs if run["direction"] == 0}
    same_initial_oracle = all(run["collision_initial_support_passed"] for run in nominal.values())
    gap_delta = abs(nominal["default"]["initial_gap_m"] - nominal["robust"]["initial_gap_m"])
    max_tilt_delta = abs(nominal["default"]["maximum_tilt_degrees"] -
                         nominal["robust"]["maximum_tilt_degrees"])
    final_spin_delta = abs(nominal["default"]["final_spin_hz"] - nominal["robust"]["final_spin_hz"])
    report = {
        "schema": "spin2014-reconstructed-g7-free-contact-v1",
        "scope": "RECONSTRUCTED_G7_SAME_SEALED_MATERIAL_CONTACT_OBSERVATION",
        "evidence_level": "RECONSTRUCTED",
        "original_paper_case": False,
        "artifact_digest": bundle["scientific_digest"],
        "input_geometry": {
            "problem_digest": read(SEAL_BUNDLE / "problem.json")["input_digest"],
            "outer_surface_sha256": sha256(collision["outer"]),
            "collision_model_scope": collision["scope"],
            "default_static_support_oracle_passed": collision["profiles"]["default"]["near_upright_passed"],
            "robust_static_support_oracle_passed": collision["profiles"]["robust"]["near_upright_passed"],
        },
        "protocol": {
            "type": "PYBULLET_FREE_TABLE_CONTACT_NO_FIXED_CONSTRAINT",
            "duration_s": 8., "dt_s": 1 / 2000., "gravity_m_s2": 9.81,
            "spin_hz": 35., "initial_tilt_deg": 2., "directions": 8,
            "friction": .8, "rolling_friction": 0., "spinning_friction": 0.,
            "linear_and_angular_damping": 0., "maintained_torque": False,
            "classification": "OBSERVATIONAL_NOT_VALIDATED",
        },
        "nominal_default_vs_robust": {
            "initial_support_oracle_agreement": same_initial_oracle,
            "absolute_initial_gap_difference_m": gap_delta,
            "absolute_maximum_tilt_difference_deg": max_tilt_delta,
            "absolute_final_spin_difference_hz": final_spin_delta,
            "interpretation": "INITIAL_SUPPORT_COMPARED_TO_GEOMETRIC_ORACLE; DYNAMIC_DIFFERENCES_REPORTED_NOT_TUNED",
        },
        "ideal_fixed_pivot_reference": {
            "source_report": str(IDEAL_JOB / "payload/report.json"),
            "comparison_limit": "NO_CONTACT_THEORY_REFERENCE_NOT_A_CONTACT_GROUND_TRUTH",
        },
        "runs": runs,
        "completion_gates": {
            "contact_proxy_initial_support": "PASSED" if same_initial_oracle else "COLLISION_MODEL_BLOCKER",
            "physical_validation": "NOT_RUN",
            "paper_reproduction": "NOT_GRANTED",
            "viewer_or_video": "NOT_RUN",
        },
    }
    if not same_initial_oracle:
        raise ValueError("COLLISION_MODEL_BLOCKER_INITIAL_SUPPORT_ORACLE_DISAGREEMENT")
    output.mkdir(parents=True, exist_ok=False)
    write_new(output / "report.json", report)
    print(json.dumps({"status": "CONTACT_OBSERVATION_RECORDED",
                      "runs": len(runs), "max_nominal_tilt_delta_deg": max_tilt_delta,
                      "max_nominal_final_spin_delta_hz": final_spin_delta}, sort_keys=True), flush=True)


def launch():
    from make_it_stand.verification.source_job import run_source_job

    files, _ = inputs()
    command = ["/usr/bin/env", "LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6",
               sys.executable, str(Path(__file__).resolve()), "worker", "{job}/payload"]
    env = {"PYTHONPATH": f"{ROOT / 'python'}:{ROOT / 'experiments/teapot_spin2014'}:{ROOT / 'experiments/reconstructed_heart_top'}:{ROOT / 'experiments/heart_stem_contact'}"}
    return run_source_job(
        JOB, command, files, maximum_wall_ns=3_600_000_000_000,
        maximum_memory_bytes=24 * 1024**3, maximum_disk_bytes=80 * 1024**3,
        minimum_start_free_bytes=100 * 1024**3, minimum_free_bytes=20 * 1024**3,
        environment=env)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        worker(Path(sys.argv[2]))
    else:
        job = launch()
        print(json.dumps({"status": job.status, "returncode": job.returncode,
                          "elapsed_seconds": job.elapsed_ns / 1e9,
                          "memory_peak": job.service.get("MemoryPeak")}, sort_keys=True), flush=True)
