"""Physics-only validation of the already sealed reconstructed G7 teapot."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[2]
BATCH = ROOT / "build/spin2014-g7-cached-fingerprint-repeat-20260920"
FIRST = ROOT / "build/spin2014-g7-deferred-archive-20260920"
SEAL_DIR = ROOT / "build/spin2014-g7-offline-trace-audit-20260920/geometry-seal-job-retry-01/payload"
SEAL_BUNDLE = SEAL_DIR / "sealed-artifact"
SEAL_REPORT = SEAL_DIR / "geometry-postprocess-report.json"
AUDIT_REPORT = BATCH / "trace-audit-v2-job-retry-03/payload/report.json"
DETERMINISM_REPORT = BATCH / "determinism-report-v2.json"
JOB = BATCH / "ideal-dynamics-job-retry-02"
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "experiments/teapot_spin2014"))
sys.path.insert(0, str(ROOT / "experiments/reconstructed_heart_top"))


def read(path):
    return json.loads(Path(path).read_text())


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_new(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()


def physical_moments(model_moments, scale, density):
    powers = np.array([3, 4, 4, 4, 5, 5, 5, 5, 5, 5], dtype=int)
    return np.asarray(model_moments, dtype=float) * density * scale**powers


def normalize_moments(moments, frame):
    m = np.asarray(moments, dtype=float)
    origin = np.asarray(frame["origin"], dtype=float)
    length, density_scale = float(frame["L"]), float(frame["R"])
    first = m[1:4]
    q = np.array([[m[4], m[9], m[8]], [m[9], m[5], m[7]], [m[8], m[7], m[6]]])
    normalized_mass = m[0] / (density_scale * length**3)
    normalized_first = (first - origin * m[0]) / (density_scale * length**4)
    normalized_q = (q - np.outer(origin, first) - np.outer(first, origin) +
                    np.outer(origin, origin) * m[0]) / (density_scale * length**5)
    return np.array([normalized_mass, *normalized_first,
                     normalized_q[0, 0], normalized_q[1, 1], normalized_q[2, 2],
                     normalized_q[1, 2], normalized_q[0, 2], normalized_q[0, 1]])


def load_qualified_dynamics_bundle(bundle_directory, *, expected_problem_digest,
                                   expected_scientific_digest):
    """Recompute the sealed object; metadata flags alone cannot admit physics.

    The current seal checker still returns UNKNOWN for full feature thickness.
    Such objects remain blocked here, including analytic boxes with a verified
    outer-to-cavity gap. Kernel tests use explicit analytic inputs separately.
    """
    from make_it_stand.verification.spin_bundle import load_spin_artifact
    artifact = load_spin_artifact(bundle_directory,
                                  expected_problem_digest=expected_problem_digest)
    if artifact.scientific_digest != expected_scientific_digest:
        raise ValueError('SEALED_ARTIFACT_DIGEST_MISMATCH')
    facts = artifact.facts
    expected = {
        'overall': 'VERIFIED', 'same_component_feature': 'VERIFIED',
        'outer_to_newly_carved': ('NOT_APPLICABLE' if artifact.newly_carved.explicit_empty
                                  else 'VERIFIED'),
        'inter_cavity_web': ('NOT_APPLICABLE' if facts['bounded_void_components'] < 2
                             else 'VERIFIED'),
    }
    if (facts.get('geometry_status') != 'MESH_MATERIAL_VALIDATED'
            or facts.get('wall_thickness_audit') != expected):
        raise ValueError('COMPLETE_MATERIAL_GEOMETRY_REQUIRED')
    return artifact


def simulate_rotational_protocol(mass, center, central_inertia, omega, *, protocol,
                                 duration, tilt_degrees, azimuth, max_step):
    """Explicit inertia origins for two different rotational problems.

    This numerical kernel accepts analytic inputs for regression tests. It is
    not a material-admission API. The worker must load a qualified bundle first.
    Both inputs are body-frame quantities: center is measured from the chosen
    fixed pivot, while central_inertia is always measured about the mass center.
    """
    from heart_top import simulate
    c = np.asarray(center, float)
    Ic = np.asarray(central_inertia, float)
    if (c.shape != (3,) or Ic.shape != (3, 3) or not np.isfinite(c).all()
            or not np.isfinite(Ic).all() or not np.allclose(Ic, Ic.T, rtol=0., atol=1e-15)
            or np.linalg.eigvalsh(Ic).min() <= 0):
        raise ValueError('VALID_CENTRAL_INERTIA_REQUIRED')
    if protocol == 'FREE_CENTER_OF_MASS':
        origin, used_center, I, gravity = 'CENTER_OF_MASS', np.zeros(3), Ic, 0.
    elif protocol == 'IDEAL_FIXED_PIVOT_GRAVITY':
        origin, used_center, gravity = 'FIXED_PIVOT', c, 9.81
        I = Ic + mass*((c@c)*np.eye(3)-np.outer(c, c))
    else:
        raise ValueError('EXPLICIT_ROTATIONAL_PROTOCOL_REQUIRED')
    result = simulate(mass, used_center, I, omega, duration=duration,
        gravity=gravity, tilt_degrees=tilt_degrees, azimuth=azimuth, max_step=max_step)
    result['reference'] = dict(protocol=protocol, inertia_origin=origin,
        center_from_rotation_origin=used_center.tolist(), inertia=I.tolist(),
        gravity=gravity, translation='NOT_INTEGRATED', material_admission='NOT_GRANTED_BY_KERNEL')
    return result


def sealed_dynamics_runtime(bundle_directory, *, native_extension):
    """Bind an explicitly selected native binary to the seal's observation.

    The observation's `path` is the label 'native', not a filesystem path.
    Binary selection comes from the frozen run spec, never from that label.
    This prepares an environment only and grants no geometry/physics status.
    """
    from importlib.machinery import EXTENSION_SUFFIXES
    extension = Path(native_extension).resolve(strict=True)
    record = read(Path(bundle_directory)/'manifest.json')['evidence']['files']['native']
    if (extension.name not in {'_make_it_stand'+suffix for suffix in EXTENSION_SUFFIXES}
            or extension.stat().st_size != record['size'] or sha256(extension) != record['sha256']):
        raise ValueError('SEALED_NATIVE_RUNTIME_MISMATCH')
    runtime = Path('/usr/lib/x86_64-linux-gnu/libstdc++.so.6').resolve(strict=True)
    paths = [extension.parent, ROOT/'python', ROOT/'experiments/teapot_spin2014',
             ROOT/'experiments/reconstructed_heart_top']
    return {'command_prefix': ['/usr/bin/env', 'LD_PRELOAD='+str(runtime)],
            'environment': {'PYTHONPATH': ':'.join(map(str, paths))},
            'identity_files': {'sealed_native_extension': str(extension),
                               'native_cxx_runtime': str(runtime)}}


def source_files():
    from make_it_stand.verification.source_checkpoint_store import _fingerprints

    g7_spec_path = BATCH / "spec.json"
    g7_spec = read(g7_spec_path)
    if _fingerprints(g7_spec["files"]) != g7_spec["identity"]:
        raise ValueError("G7_REPEAT_SOURCE_IDENTITY_STALE")
    seal = read(SEAL_REPORT)
    deterministic = read(DETERMINISM_REPORT)
    audit = read(AUDIT_REPORT)
    if (deterministic.get("status") != "DETERMINISTIC_SCIENTIFIC_OUTPUT_VERIFIED" or
            deterministic.get("terminal_geometry_files_match_sealed_artifact") is not True or
            audit.get("status") != "NATIVE_BINARY_TRACE_AUDIT_VERIFIED" or
            audit.get("source_binding") != "VERIFIED" or audit.get("checkpoint_chain") != "VERIFIED"):
        raise ValueError("G7_REPEAT_OR_SEAL_GATE_NOT_PASSED")
    files = dict(g7_spec["files"])
    files.update({
        "g7_run_spec": str(g7_spec_path),
        "g7_repeat_result": str(BATCH / "repeat-0/payload/result.json"),
        "g7_repeat_outer_result": str(BATCH / "repeat-0/result.json"),
        "g7_first_result": str(FIRST / "repeat-0/payload/result.json"),
        "determinism_report": str(DETERMINISM_REPORT),
        "native_trace_audit_report": str(AUDIT_REPORT),
        "sealed_geometry_report": str(SEAL_REPORT),
        "simulation_script": str(Path(__file__).resolve()),
        "simulator_source": str(ROOT / "experiments/reconstructed_heart_top/heart_top.py"),
        "spin_artifact_source": str(ROOT / "python/make_it_stand/verification/spin_artifact.py"),
        "spin_bundle_source": str(ROOT / "python/make_it_stand/verification/spin_bundle.py"),
        "mass_moments_source": str(ROOT / "python/make_it_stand/verification/mass_moments.py"),
    })
    for base, label in ((SEAL_BUNDLE, "sealed"),):
        for path in sorted(base.rglob("*")):
            if path.is_file():
                files[f"{label}_{path.relative_to(base).as_posix()}"] = str(path)
    for run, label in ((FIRST, "first"), (BATCH, "repeat")):
        result = read(run / "repeat-0/payload/result.json")
        terminal = Path(result["terminal"]).resolve(strict=True)
        for name in ("header.json", "outer.npz", "material.npz", "cavity.npz", "trace-index.json"):
            files[f"{label}_terminal_{name}"] = str(terminal / name)
    input_dir = Path(read(g7_spec["files"]["config"]) ["input"])
    for name in ("outer.obj", "fixed.obj", "design.obj", "manifest.json"):
        files[f"original_input_{name.replace('.', '_')}"] = str(input_dir / name)
    runtime = sealed_dynamics_runtime(SEAL_BUNDLE, native_extension=g7_spec['files']['native'])
    files.update(runtime['identity_files'])
    for name in ('spin_common', 'spin_problem', 'run_evidence'):
        files[name+'_source'] = str(ROOT/'python/make_it_stand/verification'/f'{name}.py')
    return files, g7_spec, seal, audit, deterministic


def worker(output: Path):
    from heart_top import integrate, properties

    files, g7_spec, seal_report, audit_report, deterministic = source_files()
    from make_it_stand.verification.source_checkpoint_store import _fingerprints
    if _fingerprints(files) == {}:
        raise ValueError("EMPTY_PHYSICS_SOURCE_IDENTITY")
    bundle = read(SEAL_BUNDLE / "manifest.json")
    problem_path = SEAL_BUNDLE / "problem.json"
    problem_payload = read(problem_path)
    problem_record = bundle.get("problem", {})
    if (bundle.get("scientific_digest") != seal_report["independent_seal"]["scientific_digest"] or
            problem_record.get("name") != "problem.json" or
            problem_record.get("sha256") != sha256(problem_path) or
            problem_record.get("size") != problem_path.stat().st_size or
            problem_payload.get("input_digest") != seal_report["problem_digest"]):
        raise ValueError("SEALED_ARTIFACT_DIGEST_MISMATCH")
    facts = bundle.get("facts", {})
    if facts.get("geometry_status") != "MESH_MATERIAL_VALIDATED":
        raise ValueError("SEALED_GEOMETRY_NOT_VALIDATED")
    # Cheap rejection before expensive exact reintegration. This read-only
    # preflight can only reject; only the recomputing loader below may admit.
    if facts.get('wall_thickness_audit', {}).get('overall') != 'VERIFIED':
        raise ValueError('COMPLETE_MATERIAL_GEOMETRY_REQUIRED')
    qualified = load_qualified_dynamics_bundle(SEAL_BUNDLE,
        expected_problem_digest=seal_report['problem_digest'],
        expected_scientific_digest=seal_report['independent_seal']['scientific_digest'])
    facts = qualified.facts

    def load_sealed_mesh(record):
        arrays = []
        for key in ("vertices", "faces"):
            descriptor = record[key]
            path = SEAL_BUNDLE / descriptor["name"]
            if path.stat().st_size != descriptor["size"] or sha256(path) != descriptor["sha256"]:
                raise ValueError("SEALED_MESH_ARRAY_DIGEST_MISMATCH")
            array = np.load(path, allow_pickle=False)
            if list(array.shape) != descriptor["shape"] or array.dtype.str != descriptor["dtype"]:
                raise ValueError("SEALED_MESH_ARRAY_LAYOUT_MISMATCH")
            arrays.append(array)
        return arrays[0], arrays[1]

    sealed_outer = load_sealed_mesh(bundle["final"]["outer"])
    sealed_material = load_sealed_mesh(bundle["final"]["material"])
    sealed_cavity = load_sealed_mesh(bundle["final"]["newly_carved"])

    first_result = read(FIRST / "repeat-0/payload/result.json")
    repeat_result = read(BATCH / "repeat-0/payload/result.json")
    first_terminal = Path(first_result["terminal"]).resolve(strict=True)
    repeat_terminal = Path(repeat_result["terminal"]).resolve(strict=True)
    for name, record in seal_report["terminal_files"].items():
        if (sha256(first_terminal / name) != record["sha256"] or
                sha256(repeat_terminal / name) != record["sha256"]):
            raise ValueError("REPEAT_MATERIAL_DIFFERS_FROM_SEALED_GEOMETRY")

    def compare_terminal_mesh(filename, sealed):
        with np.load(repeat_terminal / filename, allow_pickle=False) as stored:
            if (not np.array_equal(stored["vertices"], sealed[0]) or
                    not np.array_equal(stored["faces"], sealed[1])):
                raise ValueError("TERMINAL_MESH_ARRAYS_DIFFER_FROM_SEALED_BUNDLE")
    compare_terminal_mesh("outer.npz", sealed_outer)
    compare_terminal_mesh("material.npz", sealed_material)
    compare_terminal_mesh("cavity.npz", sealed_cavity)

    vertices, faces = sealed_material
    moments_model = integrate(vertices, faces)
    normalized = normalize_moments(moments_model, facts["normalization"])
    sealed_normalized = np.asarray(facts["normalized"]["moments"], dtype=float)
    normalized_error = float(np.max(np.abs(normalized - sealed_normalized)))
    if normalized_error > 1e-6:
        raise ValueError("INDEPENDENT_MASS_MOMENTS_DISAGREE_WITH_SEAL")

    scale, density = .025, 1200.
    if abs(density * scale**3 - .01875) > 1e-15:
        raise ValueError("LOCKED_MODEL_UNIT_CONVERSION_MISMATCH")
    moments_si = physical_moments(moments_model, scale, density)
    mass, center, inertia_pivot = properties(moments_si)
    inertia_center = inertia_pivot - mass * ((center @ center) * np.eye(3) - np.outer(center, center))
    if np.linalg.eigvalsh(inertia_pivot).min() <= 0 or np.linalg.eigvalsh(inertia_center).min() <= 0:
        raise ValueError("NON_POSITIVE_INERTIA_FROM_SEALED_MATERIAL")

    duration, hz = 8., 35.
    steps = (1 / 240, 1 / 480, 1 / 960)
    omega = np.array([0., 0., 2 * np.pi * hz])
    free = {}
    for step in steps:
        free[str(step)] = simulate_rotational_protocol(mass, center, inertia_center, omega,
            protocol='FREE_CENTER_OF_MASS', duration=duration,
            tilt_degrees=0., azimuth=0., max_step=step)
    finest_q = Rotation.from_quat(np.asarray(free[str(steps[-1])]["quaternions"]))
    free_rows = []
    for step in steps:
        row = free[str(step)]
        q = Rotation.from_quat(np.asarray(row["quaternions"]))
        angle_error = float(np.max((finest_q.inv() * q).magnitude()) * 180 / np.pi)
        free_rows.append({"maximum_step_s": step,
                          "relative_energy_drift": row["relative_energy_drift"],
                          "relative_world_angular_momentum_drift": row["relative_momentum_drift"],
                          "maximum_orientation_error_to_1_over_960_deg": angle_error,
                          "function_evaluations": row["function_evaluations"]})

    gravity_rows = []
    for direction in range(8):
        azimuth = direction * np.pi / 4
        for step in steps:
            row = simulate_rotational_protocol(mass, center, inertia_center, omega,
                protocol='IDEAL_FIXED_PIVOT_GRAVITY', duration=duration,
                tilt_degrees=2., azimuth=azimuth, max_step=step)
            gravity_rows.append({
                "direction": direction, "azimuth_deg": direction * 45.,
                "maximum_step_s": step, "initial_tilt_deg": 2., "spin_hz": hz,
                "maximum_tilt_deg": row["max_tilt_degrees"],
                "final_tilt_deg": row["final_tilt_degrees"],
                "relative_energy_drift": row["relative_energy_drift"],
                "relative_vertical_angular_momentum_drift": row["relative_vertical_momentum_drift"],
                "function_evaluations": row["function_evaluations"],
            })

    central_eigenvalues = np.linalg.eigvalsh(inertia_center)
    symmetry_gap = float(abs(central_eigenvalues[1] - central_eigenvalues[0]) /
                         max(float(np.max(np.abs(central_eigenvalues))), 1e-30))
    transverse_com = float(np.linalg.norm(center[:2]))
    critical_speed = {
        "status": "NOT_APPLICABLE_NON_AXISYMMETRIC_OR_OFF_AXIS_RECONSTRUCTED_BODY",
        "relative_transverse_com_to_target_z": transverse_com,
        "relative_gap_between_two_smallest_central_moments": symmetry_gap,
        "formula_evaluated": False,
    }
    report = {
        "schema": "spin2014-reconstructed-g7-ideal-dynamics-v1",
        "scope": "RECONSTRUCTED_G7_SAME_SEALED_MATERIAL_PHYSICS_ONLY",
        "evidence_level": "RECONSTRUCTED",
        "original_paper_case": False,
        "source_binding": "VERIFIED_BY_SOURCE_JOB",
        "sealed_artifact": {
            "status": facts["geometry_status"],
            "scientific_digest": bundle["scientific_digest"],
            "problem_digest": problem_payload["input_digest"],
            "repeat_geometry_byte_identical": True,
            "native_trace_audit": audit_report["status"],
            "determinism": deterministic["status"],
        },
        "independent_mass_properties": {
            "integrator": "SIGNED_TETRAHEDRAL_INTEGRATION_ON_SHA256_VERIFIED_SEALED_FINAL_MATERIAL",
            "normalized_moment_max_abs_error": normalized_error,
            "model_moments_order": ["M", "LX", "LY", "LZ", "QXX", "QYY", "QZZ", "QYZ", "QXZ", "QXY"],
            "scale_m_per_model_unit": scale, "density_kg_per_m3": density,
            "mass_kg": float(mass), "center_of_mass_m": center.tolist(),
            "inertia_about_fixed_pivot_kg_m2": inertia_pivot.tolist(),
            "central_inertia_kg_m2": inertia_center.tolist(),
        },
        "torque_free_rotational_dynamics": {
            "protocol": "FREE_CENTER_OF_MASS_ROTATION_NO_EXTERNAL_TORQUE",
            "inertia_origin": "CENTER_OF_MASS",
            "translation": "NOT_INTEGRATED",
            "initial_spin_hz": hz, "duration_s": duration,
            "time_step_refinement": free_rows,
        },
        "ideal_fixed_pivot_gravity": {
            "protocol": "NO_CONTACT_OR_FRICTION_NO_MAINTAINED_TORQUE",
            "spin_hz": hz, "duration_s": duration,
            "directions": 8, "initial_tilt_deg": 2.,
            "maximum_step_s_values": list(steps), "runs": gravity_rows,
        },
        "symmetric_top_critical_speed": critical_speed,
        "contact_simulation": "NOT_RUN_COLLISION_PROXY_VALIDATION_REMAINS_NEXT_GATE",
        "video_or_desktop_viewer": "NOT_RUN",
    }
    output.mkdir(parents=True, exist_ok=False)
    write_new(output / "report.json", report)
    print(json.dumps({"status": "IDEAL_DYNAMICS_RECORDED",
                      "mass_kg": float(mass), "directions": len(gravity_rows) // len(steps),
                      "max_normalized_moment_error": normalized_error,
                      "artifact": bundle["scientific_digest"]}, sort_keys=True), flush=True)


def launch():
    from make_it_stand.verification.source_job import run_source_job
    files, _, _, _, _ = source_files()
    runtime = sealed_dynamics_runtime(SEAL_BUNDLE, native_extension=files['sealed_native_extension'])
    return run_source_job(
        JOB, [*runtime['command_prefix'], sys.executable, '-B',
              str(Path(__file__).resolve()), "worker", "{job}/payload"], files,
        maximum_wall_ns=3_600_000_000_000, maximum_memory_bytes=24 * 1024**3,
        maximum_disk_bytes=80 * 1024**3, minimum_start_free_bytes=100 * 1024**3,
        minimum_free_bytes=20 * 1024**3, environment=runtime['environment'])


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        worker(Path(sys.argv[2]))
    else:
        job = launch()
        print(json.dumps({"status": job.status, "returncode": job.returncode,
                          "elapsed_seconds": job.elapsed_ns / 1e9,
                          "service": job.service}, sort_keys=True), flush=True)
