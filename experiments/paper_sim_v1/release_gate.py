#!/usr/bin/env python3
"""Fail-closed VERIFIED_PAPER_SIM release gate."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import subprocess

from make_it_stand.paper_evidence import audit_manifest_freshness


EXPECTED_CASES = {"toyrocking", "spheres", "armadillo1", "dinosaur", "horse"}
EXPECTED_CLASSIFICATIONS = {
    "spheres": "RECONSTRUCTED_PARAMETER",
    "armadillo1": "RECONSTRUCTED_PARAMETER",
    "dinosaur": "RECONSTRUCTED_PARAMETER",
    "horse": "RECONSTRUCTED_ASSET",
    "toyrocking": "RECONSTRUCTED_ASSET",
}


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _require_current(
    payload: dict,
    repository_root: Path,
    label: str,
    blockers: list[str],
) -> None:
    audited = audit_manifest_freshness(payload, repository_root)
    if audited["freshness"]["status"] != "CURRENT":
        blockers.append(f"{label}: STALE_RESULT")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table-case", action="append", default=[], type=Path)
    parser.add_argument("--figure14", type=Path)
    parser.add_argument("--static-oracle", type=Path)
    parser.add_argument("--pybullet-matrix", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    repository_root = Path(__file__).resolve().parents[2]
    blockers: list[str] = []
    case_evidence = {}
    for root in arguments.table_case:
        path = root.resolve() / "manifest.json"
        if not path.is_file():
            blockers.append(f"missing Table 1 manifest: {path}")
            continue
        manifest = _read(path)
        name = str(manifest.get("case"))
        case_evidence[name] = str(path)
        _require_current(manifest, repository_root, name, blockers)
        if manifest.get("classification") != EXPECTED_CLASSIFICATIONS.get(name):
            blockers.append(f"{name}: unexpected evidence classification")
        if manifest.get("repeat_count_requested") != 2:
            blockers.append(f"{name}: optimizer repeat count is not two")
        if manifest.get("deterministic_output_hashes") is not True:
            blockers.append(f"{name}: optimizer output hashes are not deterministic")
        runs = manifest.get("runs", [])
        if len(runs) != 2:
            blockers.append(f"{name}: missing two completed optimizer runs")
        for index, run in enumerate(runs, start=1):
            if run.get("completion") != "stable" or run.get("succeeded") is not True:
                blockers.append(f"{name} repeat {index}: optimizer not stable")
            first_safe = run.get("first_safe_event")
            if (
                isinstance(first_safe, bool)
                or not isinstance(first_safe, int)
                or first_safe < 0
            ):
                blockers.append(f"{name} repeat {index}: no safe committed event")
            residual = run.get("sealed_safe_region_residual")
            diagonal = run.get("sealed_bounding_box_diagonal")
            if (
                run.get("sealed_safe_region_state") != "SAFE"
                or isinstance(residual, bool)
                or not isinstance(residual, (int, float))
                or not math.isfinite(float(residual))
                or float(residual) < 0.0
                or isinstance(diagonal, bool)
                or not isinstance(diagonal, (int, float))
                or not math.isfinite(float(diagonal))
                or float(diagonal) <= 0.0
                or float(residual) > float(diagonal) * 1.0e-6
            ):
                blockers.append(
                    f"{name} repeat {index}: safe-region residual exceeds "
                    "normalized bounding-box diagonal x 1e-6"
                )
            if run.get("sealed_artifact_status") != "SEALED":
                blockers.append(f"{name} repeat {index}: final geometry not sealed")
            if run.get("sealed_outer_flattening_repair_status") not in {
                "NOT_REQUIRED",
                "APPLIED",
            }:
                blockers.append(
                    f"{name} repeat {index}: flattened support repair is "
                    "missing or failed"
                )
            if run.get("sealed_wall_thickness_status") != "PASS":
                blockers.append(
                    f"{name} repeat {index}: wall-thickness audit failed"
                )
            if run.get("sealed_containment_status") != "PASS":
                blockers.append(
                    f"{name} repeat {index}: cavity containment audit failed"
                )
    missing_cases = EXPECTED_CASES - set(case_evidence)
    blockers.extend(f"missing Table 1 case: {name}" for name in sorted(missing_cases))

    if arguments.figure14 is None or not arguments.figure14.is_file():
        blockers.append("missing Figure 14 report")
    else:
        figure14 = _read(arguments.figure14)
        _require_current(figure14, repository_root, "Figure 14", blockers)
        if figure14.get("protocol_complete") is not True:
            blockers.append("Figure 14 protocol is incomplete or nondeterministic")
        if figure14.get("qualitative_order_recovered") is not True:
            blockers.append("Figure 14 qualitative order not recovered")

    if arguments.static_oracle is None or not arguments.static_oracle.is_file():
        blockers.append("missing quasi-static tilt report")
    else:
        oracle = _read(arguments.static_oracle)
        _require_current(oracle, repository_root, "quasi-static oracle", blockers)
        if oracle.get("all_pass") is not True:
            blockers.append("quasi-static tilt gate failed")
        oracle_cases = oracle.get("cases", [])
        if {item.get("case") for item in oracle_cases} != EXPECTED_CASES:
            blockers.append("quasi-static tilt report does not contain all five cases")
        if any(item.get("accepted") is not True for item in oracle_cases):
            blockers.append("not every quasi-static case passed 0/2/5 degree tilts")

    if arguments.pybullet_matrix is None or not arguments.pybullet_matrix.is_file():
        blockers.append("missing PyBullet matrix")
    else:
        matrix = _read(arguments.pybullet_matrix)
        _require_current(matrix, repository_root, "PyBullet matrix", blockers)
        runs = matrix.get("runs", [])
        if matrix.get("all_stable") is not True:
            blockers.append("PyBullet matrix contains an unstable scenario")
        if len(runs) != 36:
            blockers.append("PyBullet matrix must contain 18 scenarios x 2 repeats")
        if set(matrix.get("case_names", [])) != EXPECTED_CASES:
            blockers.append("PyBullet matrix does not contain all five cases")
        if matrix.get("unique_scene_count") != 90:
            blockers.append("PyBullet matrix must contain 90 unique case-scenes")
        if matrix.get("simulation_count") != 180:
            blockers.append("PyBullet matrix must contain 180 repeated simulations")
        keys = {
            (item.get("preset"), item.get("scenario"), item.get("repeat"))
            for item in runs
        }
        if len(keys) != 36:
            blockers.append("PyBullet matrix contains duplicate or missing run keys")
        if any(item.get("all_cases_stable") is not True for item in runs):
            blockers.append("not every PyBullet batch reports all cases stable")

    test_results = {}
    for preset in ("headless", "paper-python", "python", "physics-python"):
        completed = subprocess.run(
            ("/usr/bin/ctest", "--preset", preset, "--output-on-failure"),
            cwd=repository_root,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        test_results[preset] = {
            "returncode": completed.returncode,
            "passed": completed.returncode == 0,
        }
        if completed.returncode != 0:
            blockers.append(f"CTest preset failed: {preset}")
        if any(
            marker in completed.stdout
            for marker in ("Not Run", "Skipped", "No tests were found")
        ):
            blockers.append(f"CTest preset has non-executed tests: {preset}")

    status = subprocess.run(
        ("git", "status", "--porcelain"),
        cwd=repository_root,
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout
    if status:
        blockers.append("worktree is not clean")
    if not subprocess.run(
        ("git", "config", "--get", "user.name"),
        cwd=repository_root,
        check=False,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout.strip() or not subprocess.run(
        ("git", "config", "--get", "user.email"),
        cwd=repository_root,
        check=False,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout.strip():
        blockers.append("Git identity is not configured; commit/tag forbidden")

    payload = {
        "schema": "make-it-stand-paper-sim-release-gate-v1",
        "schema_version": 1,
        "requested_release": "make-it-stand-paper-sim-v1",
        "claim": "VERIFIED_PAPER_SIM",
        "passed": not blockers,
        "blockers": blockers,
        "table_case_manifests": case_evidence,
        "ctest": test_results,
        "claim_boundary": (
            "VERIFIED_PAPER_SIM is not a claim of author-identical parameters, "
            "printing, or bitwise reproduction of unavailable assets."
        ),
        "spin_it_faster_unblocked": not blockers,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0 if not blockers else 2


if __name__ == "__main__":
    raise SystemExit(main())
