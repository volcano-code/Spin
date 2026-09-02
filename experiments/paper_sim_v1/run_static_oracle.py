#!/usr/bin/env python3
"""Evaluate fixed-support quasi-static tilt before any PyBullet drop test."""

from __future__ import annotations

from dataclasses import asdict
import argparse
import hashlib
import json
from pathlib import Path
import sys

from make_it_stand.paper_evidence import (
    audit_manifest_freshness,
    capture_execution_fingerprint,
)
from make_it_stand.sealed_artifact import load_sealed_pybullet_entry
from make_it_stand.stability_oracle import evaluate_tilt_azimuths


AZIMUTHS = tuple(range(0, 360, 45))


def _serializable_result(result):
    payload = asdict(result)
    payload["inside_support"] = bool(payload["inside_support"])
    return payload


def _accepted_required_angles(results: dict, tolerance: float) -> bool:
    return bool(all(
        float(item["signed_margin"]) >= -float(tolerance)
        for angle in ("0.0", "2.0", "5.0")
        for item in results[angle]
    ))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-run", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    repository_root = Path(__file__).resolve().parents[2]
    cases = []
    case_sources: list[tuple[str, Path, Path]] = []
    all_pass = True
    for raw_root in arguments.case_run:
        root = raw_root.expanduser().resolve()
        run_manifest = json.loads(
            (root / "manifest.json").read_text(encoding="utf-8")
        )
        audited = audit_manifest_freshness(run_manifest, repository_root)
        if audited["freshness"]["status"] != "CURRENT":
            raise ValueError(f"{run_manifest['case']}: STALE_RESULT")
        sealed_path = root / "repeat-1/sealed/sealed-artifact.json"
        sealed = json.loads(sealed_path.read_text(encoding="utf-8"))
        entry = load_sealed_pybullet_entry(sealed_path)
        support = sealed["support_polygon"]
        center = entry["center_of_mass_native"]
        gravity = entry["gravity_direction_native"]
        results = {
            str(angle): [
                _serializable_result(item)
                for item in evaluate_tilt_azimuths(
                    support,
                    center,
                    gravity,
                    tilt_degrees=angle,
                    azimuths_degrees=AZIMUTHS,
                )
            ]
            for angle in (0.0, 2.0, 5.0, 5.1)
        }
        safe_tolerance = float(sealed["safe_region"]["residual_limit"])
        accepted = _accepted_required_angles(results, safe_tolerance)
        all_pass = bool(all_pass and accepted)
        cases.append({
            "case": run_manifest["case"],
            "accepted": accepted,
            "protocol": "quasi_static_table_tilt",
            "contact_policy": "fixed_projected_support_polygon",
            "release_safe_region_tolerance": safe_tolerance,
            "support_polygon_sha256": sealed["support_polygon_sha256"],
            "mass_properties_sha256": sealed["mass_properties_sha256"],
            "angles": results,
        })
        case_sources.append((str(run_manifest["case"]), root / "manifest.json", sealed_path))
    payload = {
        "schema": "make-it-stand-quasi-static-tilt-matrix-v1",
        "schema_version": 1,
        "freshness": {"status": "CURRENT", "mismatches": []},
        "all_pass": all_pass,
        "cases": cases,
    }
    payload["execution_fingerprint"] = capture_execution_fingerprint(
        repository_root,
        source_paths=[
            repository_root / "python" / "make_it_stand" / "stability_oracle.py",
            repository_root / "python" / "make_it_stand" / "sealed_artifact.py",
            Path(__file__).resolve(),
            *(manifest for _, manifest, _ in case_sources),
            *(sealed for _, _, sealed in case_sources),
        ],
        executable_paths=[Path(sys.executable).resolve()],
        configuration={
            "angles_degrees": [0.0, 2.0, 5.0, 5.1],
            "azimuths_degrees": list(AZIMUTHS),
            "case_sha256": {
                name: hashlib.sha256(sealed.read_bytes()).hexdigest()
                for name, _, sealed in case_sources
            },
        },
    )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0 if all_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
