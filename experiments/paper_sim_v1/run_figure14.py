#!/usr/bin/env python3
"""Run the full-resolution Spheres Figure 14 ablation under one sealed protocol."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from make_it_stand.paper_evidence import capture_execution_fingerprint


CONDITIONS = {
    "full": (),
    "no-carving": ("--no-enable-carving",),
    "no-laplacian": ("--no-enable-laplacian",),
    "neither": ("--no-enable-carving", "--no-enable-laplacian"),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--beta-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--shape-scale-lambda", type=float, default=20.0)
    arguments = parser.parse_args()
    root = arguments.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    source = arguments.beta_root.resolve() / "data/meshes/spheres.obj"
    records = []
    for name, switches in CONDITIONS.items():
        output = root / name
        command = [
            sys.executable,
            str(Path(__file__).resolve().with_name("run_official_case.py")),
            "--beta-root", str(arguments.beta_root.resolve()),
            "--case", "spheres",
            "--profile", "paper",
            "--shape-scale-lambda", str(arguments.shape_scale_lambda),
            "--maximum-carving-blocks", "400",
            "--maximum-deformation-iterations", "400",
            "--time-limit-seconds", "3600",
            "--memory-limit-gib", "24",
            "--repeat", "2",
            "--output", str(output),
            *switches,
        ]
        completed = subprocess.run(command, check=False)
        manifest_path = output / "manifest.json"
        manifest = (
            json.loads(manifest_path.read_text(encoding="utf-8"))
            if manifest_path.is_file()
            else {}
        )
        diagnostics_path = output / "repeat-1/diagnostics.json"
        diagnostics = (
            json.loads(diagnostics_path.read_text(encoding="utf-8"))
            if diagnostics_path.is_file()
            else {}
        )
        metrics = diagnostics.get("r10_metrics", {})
        records.append({
            "condition": name,
            "switches": list(switches),
            "process_returncode": completed.returncode,
            "classification": manifest.get("classification"),
            "deterministic_output_hashes": manifest.get(
                "deterministic_output_hashes"
            ),
            "repeat_count_requested": manifest.get("repeat_count_requested"),
            "completion": (
                manifest.get("runs", [{}])[0].get("completion")
                if manifest.get("runs")
                else "missing"
            ),
            "normalized_rms_vertex_displacement": metrics.get(
                "normalized_rms_vertex_displacement"
            ),
            "final_shape_energy": metrics.get("final_shape_energy"),
            "final_projected_target_energy": metrics.get(
                "final_projected_target_energy"
            ),
            "final_cavity_volume_fraction": metrics.get(
                "final_cavity_volume_fraction"
            ),
            "result_directory": str(output),
        })
    by_name = {record["condition"]: record for record in records}
    comparable = all(
        isinstance(record["normalized_rms_vertex_displacement"], (int, float))
        and isinstance(record["final_shape_energy"], (int, float))
        for record in records
    )
    protocol_complete = bool(
        comparable
        and all(
            record["repeat_count_requested"] == 2
            and record["deterministic_output_hashes"] is True
            for record in records
        )
    )
    qualitative_order = False
    if protocol_complete:
        full = by_name["full"]
        no_carving = by_name["no-carving"]
        no_laplacian = by_name["no-laplacian"]
        neither = by_name["neither"]
        qualitative_order = bool(
            no_carving["normalized_rms_vertex_displacement"]
            >= full["normalized_rms_vertex_displacement"]
            and no_laplacian["final_shape_energy"] >= full["final_shape_energy"]
            and neither["normalized_rms_vertex_displacement"]
            >= max(
                full["normalized_rms_vertex_displacement"],
                no_laplacian["normalized_rms_vertex_displacement"],
            )
        )
    repository_root = Path(__file__).resolve().parents[2]
    fingerprint = capture_execution_fingerprint(
        repository_root,
        source_paths=[
            repository_root / "CMakeLists.txt",
            repository_root / "bindings",
            repository_root / "include",
            repository_root / "src",
            repository_root / "python",
            Path(__file__).resolve(),
            Path(__file__).resolve().with_name("run_official_case.py"),
            *(
                path
                for name in CONDITIONS
                if (path := root / name / "manifest.json").is_file()
            ),
        ],
        executable_paths=[Path(sys.executable).resolve()],
        configuration={
            "shape_scale_lambda": arguments.shape_scale_lambda,
            "conditions": {name: list(value) for name, value in CONDITIONS.items()},
            "input_sha256": _sha256(source),
        },
    )
    payload = {
        "schema": "make-it-stand-figure14-reconstruction-v1",
        "schema_version": 1,
        "classification": "RECONSTRUCTED_PARAMETER",
        "freshness": {"status": "CURRENT", "mismatches": []},
        "execution_fingerprint": fingerprint,
        "claim_boundary": (
            "Official N80 Spheres cache; lambda and rendering camera are not "
            "archived, so this is a reconstructed protocol, not paper-exact imagery."
        ),
        "original_input": {
            "external_path": str(source),
            "sha256": _sha256(source),
            "copied_into_repository": False,
        },
        "shape_scale_lambda": arguments.shape_scale_lambda,
        "protocol_complete": protocol_complete,
        "qualitative_order_recovered": qualitative_order,
        "conditions": records,
    }
    (root / "figure14.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0 if qualitative_order else 2


if __name__ == "__main__":
    raise SystemExit(main())
