"""Source-bound headless admission check for the reconstructed G7 native viewer."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from make_it_stand.verification.source_job import run_source_job  # noqa: E402

BATCH = ROOT / "build/spin2014-g7-cached-fingerprint-repeat-20260920"
SEAL_JOB = ROOT / "build/spin2014-g7-offline-trace-audit-20260920/geometry-seal-job-retry-01"
BUNDLE = SEAL_JOB / "payload/sealed-artifact"
GEOMETRY_REPORT = SEAL_JOB / "payload/geometry-postprocess-report.json"
BUILD = ROOT / "build/reconstructed-g7-native-viewer-overlay-20260920c/cmake-build"
VIEWER = BUILD / "make-it-stand"
JOB = BATCH / "native-viewer-admission-job-v2"


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_viewer_admission_result(returncode: int, stdout: str, stderr: str, *,
                                     expected_artifact: str, expected_problem: str,
                                     expected_extension: str) -> dict:
    if returncode != 0:
        raise ValueError("NATIVE_VIEWER_CHECK_FAILED: " + stderr[-2000:])
    try:
        summary = json.loads(stdout)
    except json.JSONDecodeError as exc:
        raise ValueError("NATIVE_VIEWER_SUMMARY_INVALID_JSON") from exc
    if summary.get("artifact_digest") != expected_artifact:
        raise ValueError("NATIVE_VIEWER_ARTIFACT_DIGEST_MISMATCH")
    if summary.get("problem_digest") != expected_problem:
        raise ValueError("NATIVE_VIEWER_PROBLEM_DIGEST_MISMATCH")
    if summary.get("geometry_status") != "MESH_MATERIAL_VALIDATED":
        raise ValueError("NATIVE_VIEWER_GEOMETRY_STATUS_NOT_VALIDATED")
    if Path(summary.get("loaded_paths", {}).get("native", "")).resolve() != Path(expected_extension).resolve():
        raise ValueError("NATIVE_VIEWER_EXTENSION_PATH_MISMATCH")
    required_layers = {"original", "outer", "material", "cavity"}
    if not required_layers.issubset(summary.get("layers", {})):
        raise ValueError("NATIVE_VIEWER_REQUIRED_LAYERS_MISSING")
    if summary.get("disclosures", {}).get("project_status") != "PARTIAL_BLOCKED":
        raise ValueError("NATIVE_VIEWER_SCOPE_DISCLOSURE_MISSING")
    if summary.get("disclosures", {}).get("physical_validation") != "NOT_IN_SCOPE":
        raise ValueError("NATIVE_VIEWER_PHYSICAL_SCOPE_DISCLOSURE_MISSING")
    return summary


def inputs() -> tuple[dict[str, str], dict, str, str, str]:
    manifest = read_json(BUNDLE / "manifest.json")
    problem = read_json(BUNDLE / "problem.json")
    geometry = read_json(GEOMETRY_REPORT)
    cache = (BUILD / "CMakeCache.txt").read_text(encoding="utf-8")
    cache_key = "MAKE_IT_STAND_SEALED_SPIN_EXTENSION:FILEPATH="
    configured_extensions = [line[len(cache_key):] for line in cache.splitlines() if line.startswith(cache_key)]
    if len(configured_extensions) != 1:
        raise ValueError("VIEWER_EXACT_EXTENSION_CACHE_ENTRY_MISSING")
    extension = str(Path(configured_extensions[0]).resolve(strict=True))
    native_evidence = manifest.get("evidence", {}).get("files", {}).get("native", {})
    if (manifest.get("schema") != "sealed-spin-bundle-v1"
            or manifest.get("validation_scope") != "HOLLOW_ONLY_MESH_MATERIAL_V1"
            or manifest.get("scientific_digest") != geometry.get("independent_seal", {}).get("scientific_digest")
            or geometry.get("source_binding") != "VERIFIED"
            or geometry.get("independent_seal", {}).get("status") != "MESH_MATERIAL_VALIDATED"
            or problem.get("input_digest") != geometry.get("problem_digest")
            or native_evidence.get("sha256") != sha256(Path(extension))
            or native_evidence.get("size") != Path(extension).stat().st_size):
        raise ValueError("VIEWER_BUNDLE_OR_NATIVE_EXTENSION_BINDING_MISMATCH")
    digest = manifest["scientific_digest"]
    problem_digest = problem["input_digest"]

    paths: dict[str, str] = {
        "admission_script": str(Path(__file__).resolve()),
        "cmake_lists": str((ROOT / "CMakeLists.txt").resolve()),
        "viewer_main": str((ROOT / "main.cpp").resolve()),
        "sealed_view_source": str((ROOT / "src/viewer/SealedSpinView.cpp").resolve()),
        "viewer_binary": str(VIEWER.resolve(strict=True)),
        "viewer_build_cache": str((BUILD / "CMakeCache.txt").resolve()),
        "viewer_build_graph": str((BUILD / "build.ninja").resolve()),
        "viewer_core_archive": str((BUILD / "libmake_it_stand_core.a").resolve(strict=True)),
        "viewer_glfw_archive": str((BUILD / "glfw/src/libglfw3.a").resolve(strict=True)),
        "viewer_glad_archive": str((BUILD / "glad/libglad.a").resolve(strict=True)),
        "viewer_imgui_archive": str((BUILD / "imgui/libimgui.a").resolve(strict=True)),
        "sealed_native_extension": extension,
        "sealed_bundle_manifest": str((BUNDLE / "manifest.json").resolve()),
        "sealed_problem": str((BUNDLE / "problem.json").resolve()),
        "geometry_seal_report": str(GEOMETRY_REPORT.resolve()),
        "g7_determinism_report": str((BATCH / "determinism-report-v2.json").resolve()),
    }
    observed_names = {
        "spin_common": ROOT / "python/make_it_stand/verification/spin_common.py",
        "spin_problem": ROOT / "python/make_it_stand/verification/spin_problem.py",
        "spin_artifact": ROOT / "python/make_it_stand/verification/spin_artifact.py",
        "spin_bundle": ROOT / "python/make_it_stand/verification/spin_bundle.py",
        "mass_moments": ROOT / "python/make_it_stand/verification/mass_moments.py",
        "run_evidence": ROOT / "python/make_it_stand/verification/run_evidence.py",
        "spin_modes": ROOT / "python/make_it_stand/spin.py",
    }
    for key, source in observed_names.items():
        record = manifest["evidence"]["files"][key]
        if sha256(source) != record["sha256"] or source.stat().st_size != record["size"]:
            raise ValueError(f"STALE_SEALED_PYTHON_IMPLEMENTATION:{key}")
        paths[f"sealed_source_{key}"] = str(source.resolve())
    for item in sorted(BUNDLE.glob("mesh_*.npy")):
        paths[f"sealed_{item.name}"] = str(item.resolve())
    return paths, manifest, extension, digest, problem_digest


def worker(payload: Path, extension: str, digest: str, problem_digest: str) -> None:
    command = [
        "/usr/bin/env", "LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6",
        str(VIEWER), "--check-sealed-spin-bundle", str(BUNDLE), problem_digest,
    ]
    before = {path.name: sha256(path) for path in sorted(BUNDLE.iterdir()) if path.is_file()}
    result = subprocess.run(command, capture_output=True, text=True, timeout=840, check=False)
    summary = validate_viewer_admission_result(result.returncode, result.stdout, result.stderr,
        expected_artifact=digest, expected_problem=problem_digest, expected_extension=extension)
    after = {path.name: sha256(path) for path in sorted(BUNDLE.iterdir()) if path.is_file()}
    if before != after:
        raise ValueError("SEALED_BUNDLE_CHANGED_DURING_NATIVE_ADMISSION")
    payload.mkdir(parents=True, exist_ok=False)
    report = {
        "schema": "reconstructed-g7-native-viewer-admission-v1",
        "status": "NATIVE_VIEWER_SEALED_BUNDLE_ADMISSION_VERIFIED",
        "scope": "RECONSTRUCTED_G7_SEALED_GEOMETRY_ONLY",
        "evidence_level": "RECONSTRUCTED",
        "artifact_digest": digest,
        "problem_digest": problem_digest,
        "native_extension_path": extension,
        "native_extension_sha256": sha256(Path(extension)),
        "summary": summary,
        "bundle_before_sha256": before,
        "bundle_after_sha256": after,
        "bundle_unchanged": True,
        "paper_reproduction": "NOT_GRANTED",
        "physical_validation": "NOT_RUN",
        "viewer_interaction": "HEADLESS_ADMISSION_ONLY",
    }
    with (payload / "report.json").open("x", encoding="utf-8") as stream:
        json.dump(report, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    print(json.dumps({"status": report["status"], "artifact_digest": digest,
                      "geometry_status": summary["geometry_status"]}, sort_keys=True), flush=True)


def launch():
    files, _, extension, digest, problem_digest = inputs()
    return run_source_job(
        JOB, [sys.executable, str(Path(__file__).resolve()), "worker", "{job}/payload", extension,
              digest, problem_digest], files,
        maximum_wall_ns=900_000_000_000, maximum_memory_bytes=24 * 1024**3,
        maximum_disk_bytes=2 * 1024**3, minimum_start_free_bytes=5 * 1024**3,
        minimum_free_bytes=1 * 1024**3,
        environment={"PYTHONPATH": f"{ROOT / 'python'}:{ROOT / 'experiments/teapot_spin2014'}"})


if __name__ == "__main__":
    if len(sys.argv) == 6 and sys.argv[1] == "worker":
        worker(Path(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5])
    else:
        job = launch()
        print(json.dumps({"status": job.status, "returncode": job.returncode,
                          "elapsed_seconds": job.elapsed_ns / 1e9,
                          "memory_peak": job.service.get("MemoryPeak")}, sort_keys=True), flush=True)
