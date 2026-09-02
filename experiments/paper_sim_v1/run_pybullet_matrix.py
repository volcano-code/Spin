#!/usr/bin/env python3
"""Run nominal plus eight 2-degree reseeded-pose drops for both V-HACD presets."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from make_it_stand.paper_evidence import capture_execution_fingerprint


AZIMUTHS = (0, 45, 90, 135, 180, 225, 270, 315)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--repeats", type=int, default=2)
    arguments = parser.parse_args()
    if arguments.repeats <= 0:
        parser.error("--repeats must be positive")
    manifest_path = arguments.manifest.expanduser().resolve()
    input_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if input_manifest.get("schema") != "make-it-stand-paper-sim-pybullet-input-v2":
        parser.error("--manifest is not a sealed PyBullet input v2 manifest")
    case_names = [str(item["name"]) for item in input_manifest.get("cases", [])]
    if len(case_names) != len(set(case_names)) or not case_names:
        parser.error("--manifest must contain unique, nonempty cases")
    scenarios = [("nominal", 0.0, 0.0)] + [
        (f"tilt-2deg-{azimuth:03d}", 2.0, float(azimuth))
        for azimuth in AZIMUTHS
    ]
    matrix = []
    all_stable = True
    for preset in ("default", "robust"):
        for name, tilt, azimuth in scenarios:
            for repetition in range(1, arguments.repeats + 1):
                output = (
                    arguments.output / preset / name / f"repeat-{repetition}"
                )
                command = [
                    sys.executable,
                    "-m",
                    "make_it_stand.pybullet_runner",
                    "batch",
                    "--manifest",
                    str(manifest_path),
                    "--output",
                    str(output.resolve()),
                    "--scene",
                    "table",
                    "--duration",
                    "8",
                    "--physics-hz",
                    "240",
                    "--friction",
                    "0.8",
                    "--preset",
                    preset,
                    "--initial-tilt-deg",
                    str(tilt),
                    "--initial-tilt-azimuth-deg",
                    str(azimuth),
                    "--no-render",
                ]
                completed = subprocess.run(command, check=False)
                summary_path = output / "summary.json"
                stable = False
                if completed.returncode == 0 and summary_path.is_file():
                    summary = json.loads(summary_path.read_text(encoding="utf-8"))
                    stable = all(case["stable"] for case in summary["cases"])
                all_stable = all_stable and stable
                matrix.append({
                    "preset": preset,
                    "scenario": name,
                    "repeat": repetition,
                    "perturbation_protocol": "reseeded_initial_pose_drop",
                    "initial_tilt_deg": tilt,
                    "initial_tilt_azimuth_deg": azimuth,
                    "process_returncode": completed.returncode,
                    "all_cases_stable": stable,
                    "summary_sha256": (
                        hashlib.sha256(summary_path.read_bytes()).hexdigest()
                        if summary_path.is_file()
                        else None
                    ),
                })
    repository_root = Path(__file__).resolve().parents[2]
    fingerprint = capture_execution_fingerprint(
        repository_root,
        source_paths=[
            repository_root / "python" / "make_it_stand" / "pybullet_runner.py",
            Path(__file__).resolve(),
        ],
        executable_paths=[Path(sys.executable).resolve()],
        configuration={
            "input_manifest_sha256": hashlib.sha256(
                manifest_path.read_bytes()
            ).hexdigest(),
            "case_names": case_names,
            "repeats": arguments.repeats,
        },
    )
    payload = {
        "schema": "make-it-stand-paper-sim-pybullet-matrix-v2",
        "schema_version": 2,
        "protocol": {
            "name": "reseeded_initial_pose_drop",
            "duration_seconds": 8,
            "frequency_hz": 240,
            "friction": 0.8,
            "final_tilt_limit_degrees": 5,
            "maximum_tilt_limit_degrees": 15,
            "unique_scenarios_per_case": 18,
            "repeats_per_scenario": arguments.repeats,
        },
        "freshness": {"status": "CURRENT", "mismatches": []},
        "execution_fingerprint": fingerprint,
        "input_manifest": str(manifest_path),
        "input_manifest_sha256": hashlib.sha256(
            manifest_path.read_bytes()
        ).hexdigest(),
        "case_names": case_names,
        "unique_scene_count": 18 * len(case_names),
        "simulation_count": 18 * len(case_names) * arguments.repeats,
        "all_stable": all_stable,
        "runs": matrix,
    }
    arguments.output.mkdir(parents=True, exist_ok=True)
    (arguments.output / "matrix.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0 if all_stable else 2


if __name__ == "__main__":
    raise SystemExit(main())
