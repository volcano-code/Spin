"""Render a sealed reconstructed G7 artifact along a recorded Bullet trajectory.

The replay is not a live simulation. It consumes a completed contact report,
an independently sealed mesh bundle, and the independently integrated COM.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from scipy.spatial.transform import Rotation, Slerp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from make_it_stand.verification.source_job import run_source_job  # noqa: E402

BATCH = ROOT / "build/spin2014-g7-cached-fingerprint-repeat-20260920"
SEAL_JOB = ROOT / "build/spin2014-g7-offline-trace-audit-20260920/geometry-seal-job-retry-01"
SEAL_REPORT = SEAL_JOB / "payload/geometry-postprocess-report.json"
SEAL_BUNDLE = SEAL_JOB / "payload/sealed-artifact"
CONTACT_JOB = BATCH / "free-contact-job-v1"
CONTACT_REPORT = CONTACT_JOB / "payload/report.json"
IDEAL_JOB = BATCH / "ideal-dynamics-job-retry-02"
IDEAL_REPORT = IDEAL_JOB / "payload/report.json"
JOB = BATCH / "reconstructed-g7-contact-video-job-v2"


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_new(path: Path, value: dict) -> None:
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


def select_contact_run(report: dict, artifact_digest: str, *, profile: str, direction: int) -> dict:
    if report.get("schema") != "spin2014-reconstructed-g7-free-contact-v1":
        raise ValueError("UNEXPECTED_CONTACT_REPORT_SCHEMA")
    if (report.get("scope") != "RECONSTRUCTED_G7_SAME_SEALED_MATERIAL_CONTACT_OBSERVATION"
            or report.get("evidence_level") != "RECONSTRUCTED"
            or report.get("original_paper_case") is not False):
        raise ValueError("CONTACT_REPORT_EVIDENCE_SCOPE_INVALID")
    if report.get("artifact_digest") != artifact_digest:
        raise ValueError("ARTIFACT_DIGEST_MISMATCH")
    protocol = report.get("protocol", {})
    if protocol.get("type") != "PYBULLET_FREE_TABLE_CONTACT_NO_FIXED_CONSTRAINT":
        raise ValueError("UNEXPECTED_CONTACT_PROTOCOL")
    if direction < 0 or direction >= int(protocol.get("directions", 0)):
        raise ValueError("CONTACT_DIRECTION_OUT_OF_RANGE")
    runs = report.get("runs")
    if not isinstance(runs, list):
        raise ValueError("CONTACT_RUNS_MISSING")
    selected = [run for run in runs if run.get("profile") == profile and run.get("direction") == direction]
    if len(selected) != 1:
        raise ValueError("CONTACT_RUN_NOT_UNIQUE")
    record = selected[0]
    if (record.get("artifact_digest") != artifact_digest
            or record.get("protocol") != protocol["type"]):
        raise ValueError("CONTACT_RUN_ARTIFACT_OR_PROTOCOL_MISMATCH")
    samples = record.get("samples")
    if not isinstance(samples, list) or len(samples) < 2:
        raise ValueError("CONTACT_TRAJECTORY_TOO_SHORT")
    return record


def validate_seal_report_binding(report: dict, artifact_digest: str, problem_digest: str) -> bool:
    """Read the sealed digest from the report's actual independent-seal record."""
    independent = report.get("independent_seal", {})
    return (
        report.get("schema") == "spin2014-reconstructed-g7-sealed-geometry-v1"
        and report.get("source_binding") == "VERIFIED"
        and report.get("problem_digest") == problem_digest
        and independent.get("status") == "MESH_MATERIAL_VALIDATED"
        and independent.get("scientific_digest") == artifact_digest
    )


def load_sealed_mesh(bundle_dir: Path, manifest: dict, role: str) -> tuple[np.ndarray, np.ndarray]:
    try:
        record = manifest["final"][role]
        descriptors = (record["vertices"], record["faces"])
    except (KeyError, TypeError) as exc:
        raise ValueError("SEALED_MESH_ROLE_MISSING") from exc

    loaded = []
    for descriptor in descriptors:
        name = descriptor.get("name")
        if not isinstance(name, str) or Path(name).name != name:
            raise ValueError("SEALED_MESH_FILENAME_INVALID")
        path = Path(bundle_dir) / name
        if not path.is_file() or path.stat().st_size != descriptor.get("size"):
            raise ValueError("SEALED_MESH_FILE_SIZE_MISMATCH")
        if sha256(path) != descriptor.get("sha256"):
            raise ValueError("SEALED_MESH_FILE_HASH_MISMATCH")
        try:
            array = np.load(path, allow_pickle=False)
        except (OSError, ValueError) as exc:
            raise ValueError("SEALED_MESH_NPY_INVALID") from exc
        if list(array.shape) != descriptor.get("shape") or array.dtype.str != descriptor.get("dtype"):
            raise ValueError("SEALED_MESH_ARRAY_DESCRIPTOR_MISMATCH")
        loaded.append(array)

    vertices, faces = loaded
    if (vertices.ndim != 2 or vertices.shape[1] != 3 or not np.issubdtype(vertices.dtype, np.floating)
            or faces.ndim != 2 or faces.shape[1] != 3 or not np.issubdtype(faces.dtype, np.integer)
            or len(vertices) < 3 or len(faces) < 1 or not np.isfinite(vertices).all()
            or np.any(faces < 0) or np.any(faces >= len(vertices))
            or np.any(faces[:, 0] == faces[:, 1]) or np.any(faces[:, 1] == faces[:, 2])
            or np.any(faces[:, 2] == faces[:, 0])):
        raise ValueError("SEALED_MESH_TOPOLOGY_INVALID")
    return vertices, faces


def interpolate_recorded_poses(record: dict, requested_times: np.ndarray) -> dict[str, np.ndarray]:
    samples = record.get("samples")
    times = np.asarray([sample.get("time", np.nan) for sample in samples], dtype=float)
    positions = np.asarray([sample.get("position", [np.nan] * 3) for sample in samples], dtype=float)
    quaternions = np.asarray([sample.get("model_quaternion", [np.nan] * 4) for sample in samples], dtype=float)
    requested = np.asarray(requested_times, dtype=float)
    if (times.ndim != 1 or len(times) < 2 or positions.shape != (len(times), 3)
            or quaternions.shape != (len(times), 4) or not np.isfinite(times).all()
            or not np.isfinite(positions).all() or not np.isfinite(quaternions).all()
            or np.any(np.diff(times) <= 0)):
        raise ValueError("INVALID_RECORDED_TIME_SEQUENCE")
    norms = np.linalg.norm(quaternions, axis=1)
    if not np.allclose(norms, 1.0, rtol=0.0, atol=1e-7):
        raise ValueError("RECORDED_QUATERNION_NOT_UNIT_LENGTH")
    if (requested.ndim != 1 or not np.isfinite(requested).all()
            or np.any(np.diff(requested) < 0)
            or np.any(requested < times[0]) or np.any(requested > times[-1])):
        raise ValueError("POSE_TIME_OUTSIDE_RECORDED_INTERVAL")
    sampled_positions = np.column_stack([
        np.interp(requested, times, positions[:, axis]) for axis in range(3)
    ])
    sampled_rotations = Slerp(times, Rotation.from_quat(quaternions))(requested)
    return {"positions": sampled_positions, "quaternions": sampled_rotations.as_quat()}


def model_origins_from_com_poses(positions: np.ndarray, quaternions: np.ndarray,
                                 center_of_mass_m: np.ndarray) -> np.ndarray:
    positions = np.asarray(positions, dtype=float)
    quaternions = np.asarray(quaternions, dtype=float)
    center = np.asarray(center_of_mass_m, dtype=float)
    if (positions.ndim != 2 or positions.shape[1] != 3 or quaternions.shape != (len(positions), 4)
            or center.shape != (3,) or not np.isfinite(positions).all()
            or not np.isfinite(quaternions).all() or not np.isfinite(center).all()):
        raise ValueError("INVALID_COM_POSE_CONVERSION_INPUT")
    rotations = Rotation.from_quat(quaternions)
    return positions - rotations.apply(np.broadcast_to(center, positions.shape))


def _contact_job_inputs() -> tuple[dict, dict, dict, dict]:
    config = read_json(CONTACT_JOB / "config.json")
    result = read_json(CONTACT_JOB / "result.json")
    complete = read_json(CONTACT_JOB / "complete.json")
    report = read_json(CONTACT_REPORT)
    if (config.get("schema") != "source-job-v1"
            or result.get("schema") != "source-job-result-v1"
            or result.get("status") != "PROCESS_EXITED" or result.get("returncode") != 0
            or complete.get("result_sha256") != sha256(CONTACT_JOB / "result.json")):
        raise ValueError("CONTACT_SOURCE_JOB_RECEIPT_INVALID")
    for name, identity in config.get("identity", {}).items():
        path = Path(identity.get("path", ""))
        if not path.is_file() or path.stat().st_size != identity.get("size"):
            raise ValueError(f"CONTACT_SOURCE_IDENTITY_STALE:{name}")
        if sha256(path) != identity.get("sha256"):
            raise ValueError(f"CONTACT_SOURCE_IDENTITY_STALE:{name}")
    return config, result, complete, report


def inputs() -> tuple[dict[str, str], dict]:
    parent_config, _, _, report = _contact_job_inputs()
    bundle = read_json(SEAL_BUNDLE / "manifest.json")
    problem = read_json(SEAL_BUNDLE / "problem.json")
    seal_report = read_json(SEAL_REPORT)
    ideal_report = read_json(IDEAL_REPORT)
    if (bundle.get("schema") != "sealed-spin-bundle-v1"
            or bundle.get("validation_scope") != "HOLLOW_ONLY_MESH_MATERIAL_V1"
            or not validate_seal_report_binding(seal_report, bundle.get("scientific_digest"),
                                                problem.get("input_digest"))
            or read_json(SEAL_JOB / "complete.json").get("result_sha256")
            != sha256(SEAL_JOB / "result.json")
            or ideal_report.get("sealed_artifact", {}).get("scientific_digest") != bundle.get("scientific_digest")
            or ideal_report.get("sealed_artifact", {}).get("status") != "MESH_MATERIAL_VALIDATED"
            or problem.get("input_digest") != report.get("input_geometry", {}).get("problem_digest")
            or report.get("artifact_digest") != bundle.get("scientific_digest")):
        raise ValueError("SEALED_CONTACT_ARTIFACT_CHAIN_MISMATCH")
    select_contact_run(report, bundle["scientific_digest"], profile="default", direction=0)

    files: dict[str, str] = {
        "video_script": str(Path(__file__).resolve()),
        "contact_report": str(CONTACT_REPORT.resolve()),
        "contact_config": str((CONTACT_JOB / "config.json").resolve()),
        "contact_result": str((CONTACT_JOB / "result.json").resolve()),
        "contact_complete": str((CONTACT_JOB / "complete.json").resolve()),
        "seal_report": str(SEAL_REPORT.resolve()),
        "seal_manifest": str((SEAL_BUNDLE / "manifest.json").resolve()),
        "seal_problem": str((SEAL_BUNDLE / "problem.json").resolve()),
        "ideal_report": str(IDEAL_REPORT.resolve()),
        "ffmpeg": shutil.which("ffmpeg"),
    }
    if not files["ffmpeg"]:
        raise ValueError("FFMPEG_NOT_FOUND")
    for key, identity in parent_config["identity"].items():
        files[f"parent_contact_source_{key}"] = identity["path"]
    for name in sorted(SEAL_BUNDLE.glob("mesh_*.npy")):
        files[f"sealed_{name.name}"] = str(name.resolve())
    for module_name in ("pybullet", "scipy", "PIL", "numpy"):
        module = __import__(module_name)
        files[f"module_{module_name}"] = str(Path(module.__file__).resolve())
    return files, report


def worker(payload: Path) -> None:
    import pybullet as bullet
    from PIL import Image, ImageDraw, ImageFont

    parent_config, parent_result, parent_complete, contact_report = _contact_job_inputs()
    bundle = read_json(SEAL_BUNDLE / "manifest.json")
    ideal_report = read_json(IDEAL_REPORT)
    digest = bundle["scientific_digest"]
    run = select_contact_run(contact_report, digest, profile="default", direction=0)
    outer_v, outer_f = load_sealed_mesh(SEAL_BUNDLE, bundle, "outer")
    material_v, material_f = load_sealed_mesh(SEAL_BUNDLE, bundle, "material")
    removed_v, removed_f = load_sealed_mesh(SEAL_BUNDLE, bundle, "newly_carved")

    center_m = np.asarray(ideal_report["independent_mass_properties"]["center_of_mass_m"], dtype=float)
    scale = float(ideal_report["independent_mass_properties"]["scale_m_per_model_unit"])
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("SEALED_ARTIFACT_SCALE_INVALID")
    frame_times = np.linspace(0.0, float(run["samples"][-1]["time"]), 241)
    poses = interpolate_recorded_poses(run, frame_times)
    origins = model_origins_from_com_poses(poses["positions"], poses["quaternions"], center_m)

    payload.mkdir(parents=True, exist_ok=False)
    client = bullet.connect(bullet.DIRECT)
    if client < 0:
        raise RuntimeError("PYBULLET_DIRECT_CONNECT_FAILED")
    video_path = payload / "reconstructed-g7-contact-replay.mp4"
    try:
        floor_visual = bullet.createVisualShape(
            bullet.GEOM_BOX, halfExtents=[0.4, 0.4, 0.02], rgbaColor=[0.59, 0.55, 0.49, 1],
            physicsClientId=client)
        bullet.createMultiBody(baseMass=0, baseVisualShapeIndex=floor_visual,
                               basePosition=[0, 0, -0.02], physicsClientId=client)

        outer_shape = bullet.createVisualShape(
            bullet.GEOM_MESH, vertices=(outer_v * scale).tolist(), indices=outer_f.reshape(-1).tolist(),
            rgbaColor=[0.21, 0.48, 0.34, 1], physicsClientId=client)
        dynamic_body = bullet.createMultiBody(baseMass=0, baseVisualShapeIndex=outer_shape,
                                              basePosition=origins[0].tolist(),
                                              baseOrientation=poses["quaternions"][0].tolist(),
                                              physicsClientId=client)

        material_shape = bullet.createVisualShape(
            bullet.GEOM_MESH, vertices=(material_v * scale).tolist(), indices=material_f.reshape(-1).tolist(),
            rgbaColor=[0.91, 0.64, 0.15, 1], physicsClientId=client)
        material_body = bullet.createMultiBody(baseMass=0, baseVisualShapeIndex=material_shape,
                                               basePosition=[-0.045, 0, 0], physicsClientId=client)
        removed_shape = bullet.createVisualShape(
            bullet.GEOM_MESH, vertices=(removed_v * scale).tolist(), indices=removed_f.reshape(-1).tolist(),
            rgbaColor=[0.78, 0.24, 0.20, 1], physicsClientId=client)
        removed_body = bullet.createMultiBody(baseMass=0, baseVisualShapeIndex=removed_shape,
                                              basePosition=[0.055, 0, 0], physicsClientId=client)
        del material_body, removed_body

        main_view = bullet.computeViewMatrix([0.09, -0.14, 0.09], [0, 0, 0.02], [0, 0, 1])
        main_projection = bullet.computeProjectionMatrixFOV(36, 640 / 480, 0.001, 2.0)
        inset_view = bullet.computeViewMatrix([0.055, -0.25, 0.10], [0.055, 0, 0.025], [0, 0, 1])
        inset_projection = bullet.computeProjectionMatrixFOV(35, 320 / 480, 0.001, 2.0)

        def snap(view, projection, width, height):
            image = bullet.getCameraImage(width, height, viewMatrix=view, projectionMatrix=projection,
                lightDirection=[-2, -3, 5], shadow=1, renderer=bullet.ER_TINY_RENDERER,
                physicsClientId=client)
            return Image.fromarray(np.asarray(image[2], dtype=np.uint8).reshape(height, width, 4)[:, :, :3])

        try:
            font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 17)
            small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 13)
        except OSError:
            font, small = ImageFont.load_default(), ImageFont.load_default()

        log_path = payload / "ffmpeg.log"
        with log_path.open("x", encoding="utf-8") as log:
            encoder = subprocess.Popen([
                shutil.which("ffmpeg"), "-nostdin", "-y", "-f", "rawvideo", "-pixel_format", "rgb24",
                "-video_size", "960x560", "-framerate", "30", "-i", "pipe:0", "-an", "-c:v", "libx264",
                "-preset", "fast", "-crf", "19", "-pix_fmt", "yuv420p", str(video_path),
            ], stdin=subprocess.PIPE, stdout=log, stderr=log)
            try:
                for frame, time_s in enumerate(frame_times):
                    bullet.resetBasePositionAndOrientation(dynamic_body, origins[frame].tolist(),
                        poses["quaternions"][frame].tolist(), physicsClientId=client)
                    canvas = Image.new("RGB", (960, 560), "#f4f3ef")
                    draw = ImageDraw.Draw(canvas)
                    draw.text((14, 7), "Reconstructed G7 teapot | PyBullet free-contact replay", font=font, fill="#20252b")
                    draw.text((14, 30), f"default V-HACD | 35 rps | 2 deg | direction 0 | t={time_s:.2f}s",
                              font=small, fill="#43484d")
                    canvas.paste(snap(main_view, main_projection, 640, 480), (0, 55))
                    canvas.paste(snap(inset_view, inset_projection, 320, 480), (640, 55))
                    draw.rectangle((0, 55, 639, 80), fill="#e5e3dc")
                    draw.text((10, 59), "RECORDED CONTACT MOTION (1x)", font=small, fill="#20252b")
                    draw.rectangle((640, 55, 959, 80), fill="#e5e3dc")
                    draw.text((648, 59), "SEALED MATERIAL / REMOVED VOLUME", font=small, fill="#20252b")
                    draw.text((649, 518), "amber material | red void, offset for display",
                              font=small, fill="#30343a")
                    draw.text((14, 541),
                        "Recorded-trajectory replay, software renderer; NOT live dynamics, author-original, or physical validation.",
                        font=small, fill="#43484d")
                    if frame in (0, 120, 240):
                        canvas.save(payload / f"frame-{frame:03d}.png")
                    encoder.stdin.write(np.asarray(canvas, dtype=np.uint8).tobytes())
                    if frame % 60 == 0:
                        print(f"RENDER_FRAME {frame}/240", flush=True)
            finally:
                encoder.stdin.close()
            if encoder.wait(timeout=120) != 0 or not video_path.is_file() or video_path.stat().st_size == 0:
                raise RuntimeError("FFMPEG_VIDEO_ENCODING_FAILED")
    finally:
        bullet.disconnect(client)

    inputs_hashes = {name: {"path": path, "sha256": sha256(Path(path)), "size": Path(path).stat().st_size}
                     for name, path in inputs()[0].items()}
    video = {
        "schema": "reconstructed-g7-contact-video-v1",
        "status": "RECORDED_CONTACT_TRAJECTORY_RENDERED",
        "evidence_level": "RECONSTRUCTED",
        "scope": "RECONSTRUCTED_NOT_ORIGINAL_PAPER",
        "simulation": "REPLAY_NOT_LIVE_SIMULATION",
        "renderer": "PYBULLET_TINYRENDERER_SOFTWARE",
        "physical_validation": "NOT_RUN",
        "paper_reproduction": "NOT_GRANTED",
        "artifact_digest": digest,
        "contact_job": {
            "status": parent_result["status"],
            "result_sha256": sha256(CONTACT_JOB / "result.json"),
            "complete_sha256": sha256(CONTACT_JOB / "complete.json"),
            "contact_report_sha256": sha256(CONTACT_REPORT),
            "selection": {"profile": "default", "direction": 0, "spin_hz": 35.0},
        },
        "coordinate_mapping": {
            "recorded_position": "PyBullet center-of-mass world position in metres",
            "display_model_origin": "COM_position - quaternion.rotate(sealed_center_of_mass_m)",
            "mesh_scale_m_per_unit": scale,
        },
        "frames": len(frame_times), "fps": 30, "duration_s": 8.0, "playback_speed": 1.0,
        "trajectory_samples_used": len(run["samples"]),
        "source_input_hashes": inputs_hashes,
        "output_hashes": {p.name: {"sha256": sha256(p), "size": p.stat().st_size}
                          for p in sorted(payload.glob("frame-*.png")) + [video_path]},
    }
    write_new(payload / "video.json", video)
    print(json.dumps({"status": video["status"], "frames": len(frame_times),
                      "video_sha256": video["output_hashes"][video_path.name]["sha256"]}, sort_keys=True), flush=True)


def launch():
    files, _ = inputs()
    return run_source_job(
        JOB, [sys.executable, str(Path(__file__).resolve()), "worker", "{job}/payload"], files,
        maximum_wall_ns=900_000_000_000, maximum_memory_bytes=4 * 1024**3,
        maximum_disk_bytes=2 * 1024**3, minimum_start_free_bytes=5 * 1024**3,
        minimum_free_bytes=1 * 1024**3,
        environment={"PYTHONPATH": f"{ROOT / 'python'}:{ROOT / 'experiments/teapot_spin2014'}"})


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        worker(Path(sys.argv[2]))
    else:
        job = launch()
        print(json.dumps({"status": job.status, "returncode": job.returncode,
                          "elapsed_seconds": job.elapsed_ns / 1e9,
                          "memory_peak": job.service.get("MemoryPeak")}, sort_keys=True), flush=True)
