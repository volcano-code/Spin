#!/usr/bin/env python3
"""Source-bound local checks for selected Spin-It (2014) Givens equations.

Only Eqs. (1)--(3), (5)--(8), plus the displayed static objective are in
scope.  Equation (4), original models, and paper experiments are explicitly
outside this runner's grant.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import sys
import time
from typing import Any

import numpy as np

from make_it_stand.verification import spin2014_givens_reference as reference


ROOT = Path(__file__).resolve().parents[3]
EVIDENCE_SCHEMA = "spin2014-givens-equations-evidence-v1"
RESULT_SCHEMA = "spin2014-givens-equations-result-v1"
SCOPED_REQUIREMENTS = (
    "spin2014.eq.1", "spin2014.eq.2", "spin2014.eq.3", "spin2014.eq.5", "spin2014.eq.6",
    "spin2014.eq.7", "spin2014.eq.8", "spin2014.objective.static",
    "spin2014.objective.variants",
)


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path) -> dict[str, object]:
    path = path.resolve(strict=True)
    first, size = _sha256(path), path.stat().st_size
    second = _sha256(path)
    if first != second:
        raise ValueError("file changed while hashing: " + str(path))
    return {"path": str(path), "sha256": second, "size": size}


def _record_matches(record: object, label: str) -> list[str]:
    if not isinstance(record, dict):
        return [label + " record is malformed"]
    try:
        current = file_record(Path(record["path"]))
    except (KeyError, OSError, TypeError, ValueError):
        return [label + " path is unavailable"]
    return [] if current == record else [label + " record differs from current bytes"]


def _git_context() -> dict[str, object]:
    try:
        head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--verify", "HEAD"],
                              check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True).stdout.strip()
        status = subprocess.run(["git", "-C", str(ROOT), "status", "--short", "--untracked-files=all"],
                                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                text=True).stdout.splitlines()
        digests = {}
        for label, arguments in {
            "head_to_worktree": ("diff", "--binary", "HEAD"),
            "head_to_index": ("diff", "--binary", "--cached", "HEAD"),
            "index_to_worktree": ("diff", "--binary"),
        }.items():
            payload = subprocess.run(["git", "-C", str(ROOT), *arguments], check=True,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
            digests[label + "_sha256"] = hashlib.sha256(payload).hexdigest()
    except (OSError, subprocess.CalledProcessError) as error:
        raise ValueError("cannot capture Git context") from error
    return {"head": head, "changed_paths": [line[3:] for line in status if len(line) >= 4],
            "diff_sha256": digests}


def _native_module(directory: Path):
    directory = directory.resolve(strict=True)
    extension = next(iter(sorted(directory.glob("_make_it_stand*.so"))), None)
    if extension is None:
        raise ValueError("native module is missing")
    sys.path.insert(0, str(directory))
    try:
        native = importlib.import_module("_make_it_stand")
    finally:
        sys.path.pop(0)
    if Path(native.__file__).resolve() != extension.resolve():
        raise ValueError("loaded native module differs from declared module path")
    return native, extension


def _fixtures() -> tuple[dict[str, object], ...]:
    return (
        {"id": "asymmetric_nonzero_cross_terms",
         "raw": (4.0, 0.4, -0.8, 2.0, 10.0, 14.0, 6.0, 0.3, -0.4, 0.2),
         "inertia_weight": 2.5, "center_weight": 0.7, "phi": 0.37},
        {"id": "translated_moment_frame_nonzero_cross_terms",
         "raw": (5.0, -0.5, 0.7, -1.4, 8.0, 11.0, 9.0, -0.6, 0.8, -0.4),
         "inertia_weight": 1.3, "center_weight": 0.9, "phi": -0.29},
    )


def _stability_fixtures() -> tuple[dict[str, float | str], ...]:
    """Two non-degenerate Eq. (1) inputs, including non-unit gravity."""

    return (
        {"id": "positive_center_height", "mass": 4.0, "gravity": 9.81,
         "center_height": 0.5, "axial_inertia": 10.0, "middle_inertia": 7.0},
        {"id": "offset_heavier_top", "mass": 5.5, "gravity": 9.81,
         "center_height": 1.25, "axial_inertia": 20.0, "middle_inertia": 12.5},
    )


def _evaluate_native(native, objective: str, raw: np.ndarray, gamma_i: float,
                     gamma_c: float, phi: float) -> dict[str, Any]:
    enum = getattr(native.Spin2014Objective, objective)
    value = native.evaluate_spin_it_2014_givens(enum, raw, gamma_i, gamma_c, phi)
    if value.get("objective_token") != objective or value.get("formulation") != "PAPER_GIVENS":
        raise ValueError("native Givens identity differs")
    return value


def _differential_measurement(native, objective: str, raw: np.ndarray,
                              gamma_i: float, gamma_c: float, phi: float) -> dict[str, float]:
    state = np.concatenate((raw, (phi,)))
    exact = _evaluate_native(native, objective, raw, gamma_i, gamma_c, phi)
    gradient, hessian = np.asarray(exact["objective_gradient"], dtype=float), np.asarray(exact["objective_hessian"], dtype=float)
    active_dimension = 10 if objective == "STATIC_LOWERING" else 11
    maximum_gradient, maximum_hessian = 0.0, 0.0
    for step in (1e-4, 3e-5, 1e-5):
        finite_gradient, finite_hessian = np.zeros(11), np.zeros((11, 11))
        for column in range(active_dimension):
            plus, minus = state.copy(), state.copy()
            plus[column] += step; minus[column] -= step
            positive = _evaluate_native(native, objective, plus[:10], gamma_i, gamma_c, plus[10])
            negative = _evaluate_native(native, objective, minus[:10], gamma_i, gamma_c, minus[10])
            finite_gradient[column] = (positive["objective"] - negative["objective"]) / (2 * step)
            finite_hessian[:, column] = (np.asarray(positive["objective_gradient"]) -
                                         np.asarray(negative["objective_gradient"])) / (2 * step)
        maximum_gradient = max(maximum_gradient, float(np.max(np.abs(
            gradient[:active_dimension] - finite_gradient[:active_dimension]))))
        maximum_hessian = max(maximum_hessian, float(np.max(np.abs(
            hessian[:active_dimension, :active_dimension] -
            finite_hessian[:active_dimension, :active_dimension]))))
    if maximum_gradient > 2e-7 or maximum_hessian > 3e-6:
        raise ValueError("native Givens derivatives disagree with finite differences")
    return {"active_dimension": active_dimension,
            "gradient_max_absolute_error": maximum_gradient,
            "hessian_max_absolute_error": maximum_hessian}


def _scientific_data(native_directory: Path) -> tuple[dict[str, object], Path]:
    native, native_path = _native_module(native_directory)
    stability_checks = []
    for fixture in _stability_fixtures():
        expected = reference.stability_threshold_squared(
            fixture["mass"], fixture["gravity"], fixture["center_height"],
            fixture["axial_inertia"], fixture["middle_inertia"],
        )
        actual = native.spin_it_2014_stability_threshold_squared(
            fixture["mass"], fixture["gravity"], fixture["center_height"],
            fixture["axial_inertia"], fixture["middle_inertia"],
        )
        error = abs(float(actual) - expected)
        if error > 1e-13:
            raise ValueError("native Eq. (1) threshold differs from independent scalar formula")
        stability_checks.append({**fixture, "reference_omega_squared": expected,
                                 "native_omega_squared": float(actual),
                                 "max_absolute_error": error})
    checks = []
    for fixture in _fixtures():
        raw = np.asarray(fixture["raw"], dtype=float)
        gamma_i, gamma_c, phi = (float(fixture[key]) for key in ("inertia_weight", "center_weight", "phi"))
        objective_scopes = (
            ("YOYO", ("spin2014.eq.2", "spin2014.eq.5", "spin2014.eq.6",
                      "spin2014.objective.variants")),
            ("TOP_INERTIA_ONLY", ("spin2014.eq.3", "spin2014.eq.7", "spin2014.eq.8",
                                   "spin2014.objective.variants")),
            ("TOP_LOWERING_ONLY", ("spin2014.eq.3", "spin2014.eq.7", "spin2014.eq.8",
                                    "spin2014.objective.variants")),
            ("TOP_MASS_AND_LOWERING", ("spin2014.eq.3", "spin2014.eq.7", "spin2014.eq.8",
                                        "spin2014.objective.variants")),
            ("STATIC_LOWERING", ("spin2014.objective.static",)),
        )
        for objective, scope in objective_scopes:
            actual_phi = 0.0 if objective == "STATIC_LOWERING" else phi
            expected = reference.evaluate_paper_givens(objective, raw, gamma_i, gamma_c, actual_phi)
            actual = _evaluate_native(native, objective, raw, gamma_i, gamma_c, actual_phi)
            objective_error = abs(float(actual["objective"]) - float(expected["objective"]))
            constraint_error = float(np.max(np.abs(np.asarray(actual["constraints"]) - expected["constraints"]), initial=0.0))
            if objective_error > 1e-13 or constraint_error > 1e-13:
                raise ValueError("native Givens value or constraint differs from independent formula")
            checks.append({"fixture_id": fixture["id"], "objective": objective,
                           "requirements": list(scope), "raw_moments": raw.tolist(),
                           "inertia_weight": gamma_i, "center_weight": gamma_c, "phi": actual_phi,
                           "reference_objective": float(expected["objective"]),
                           "native_objective": float(actual["objective"]),
                           "reference_constraints": expected["constraints"].tolist(),
                           "native_constraints": np.asarray(actual["constraints"]).tolist(),
                           "objective_max_absolute_error": objective_error,
                           "constraint_max_absolute_error": constraint_error})
    differential = {
        objective: _differential_measurement(
            native, objective, np.asarray(_fixtures()[0]["raw"], dtype=float),
            2.5, 0.7, 0.0 if objective == "STATIC_LOWERING" else 0.37,
        )
        for objective, _scope in objective_scopes
    }
    return {"paper_locator": "Spin-It (2014), Eqs. (1)--(3), (5)--(8), static objective",
            "paper_source_form": "LOCAL_TEXT_EXTRACTION_PDF_NOT_PRESENT_IN_STRICT_INPUTS",
            "paper_text_sha256": reference.PAPER_TEXT_SHA256,
            "requirements": list(SCOPED_REQUIREMENTS), "stability_checks": stability_checks,
            "checks": checks,
            "native_differential_measurement": differential,
            "outcome": "REFERENCE_AND_NATIVE_GIVENS_VALUES_AGREE"}, native_path


def run(paper_text: Path, native_module_directory: Path, output_directory: Path) -> dict[str, object]:
    if output_directory.exists():
        raise ValueError("output directory already exists: " + str(output_directory))
    text_record = file_record(paper_text)
    if text_record["sha256"] != reference.PAPER_TEXT_SHA256:
        raise ValueError("Spin-It 2014 local text extraction SHA-256 mismatch")
    started, start = datetime.now(timezone.utc).isoformat(), time.perf_counter()
    scientific, native_path = _scientific_data(native_module_directory)
    scientific_digest = hashlib.sha256(canonical_json(scientific).encode("utf-8")).hexdigest()
    output_directory.mkdir(parents=True, exist_ok=False)
    result_path = output_directory / "result.json"
    result_path.write_text(canonical_json({"schema": RESULT_SCHEMA, "scientific_data": scientific,
                                           "scientific_data_sha256": scientific_digest}) + "\n", encoding="utf-8")
    build_root = native_path.resolve().parents[1]
    inputs = {"runner": file_record(Path(__file__)), "checker": file_record(Path(reference.__file__)),
              "paper_text": text_record, "native_module": file_record(native_path),
              "build_cmake_cache": file_record(build_root / "CMakeCache.txt"),
              "source_cmake_lists": file_record(ROOT / "CMakeLists.txt")}
    evidence = {"schema": EVIDENCE_SCHEMA, "inputs": inputs, "outputs": {"result": file_record(result_path)},
                "git": _git_context(), "runtime": {"python_executable": str(Path(sys.executable).resolve()),
                "python_version": sys.version, "started_utc": started, "elapsed_seconds": time.perf_counter() - start},
                "scientific_data_sha256": scientific_digest, "project_status": "PARTIAL_BLOCKED",
                "author_bitwise_exact": "NOT_VERIFIED", "physical_validation": "NOT_IN_SCOPE",
                "paper_source_admission": "LOCAL_TEXT_EXTRACTION_ONLY_NOT_ORIGINAL_PDF",
                "build_source_attestation": "NATIVE_BINARY_AND_BUILD_CONFIGURATION_BOUND"}
    evidence_path = output_directory / "evidence.json"
    evidence_path.write_text(canonical_json(evidence) + "\n", encoding="utf-8")
    return {"output_directory": str(output_directory.resolve()), "scientific_data_sha256": scientific_digest,
            "status": "CURRENT"}


def check_evidence_freshness(evidence_path: Path) -> dict[str, object]:
    try:
        evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("cannot load Givens equation evidence") from error
    mismatches = []
    if evidence.get("schema") != EVIDENCE_SCHEMA:
        mismatches.append("evidence schema differs")
    inputs = evidence.get("inputs")
    if not isinstance(inputs, dict):
        mismatches.append("input records are malformed")
    else:
        for label in ("runner", "checker", "paper_text", "native_module", "build_cmake_cache", "source_cmake_lists"):
            mismatches.extend(_record_matches(inputs.get(label), label))
        if inputs.get("paper_text", {}).get("sha256") != reference.PAPER_TEXT_SHA256:
            mismatches.append("paper text digest differs")
    outputs = evidence.get("outputs")
    if not isinstance(outputs, dict):
        mismatches.append("output records are malformed")
    else:
        mismatches.extend(_record_matches(outputs.get("result"), "result"))
    # ``git`` is retained as an observation for provenance, not as a scientific
    # dependency.  The actual calculation dependencies above are byte-bound;
    # an unrelated report or untracked file must not stale a numeric result.
    if not mismatches:
        result = json.loads(Path(outputs["result"]["path"]).read_text(encoding="utf-8"))
        digest = hashlib.sha256(canonical_json(result.get("scientific_data")).encode("utf-8")).hexdigest()
        if result.get("schema") != RESULT_SCHEMA or digest != result.get("scientific_data_sha256") or digest != evidence.get("scientific_data_sha256"):
            mismatches.append("scientific data digest differs")
    return {"status": "CURRENT" if not mismatches else "STALE_RESULT", "mismatches": mismatches}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--paper-text", type=Path, required=True)
    parser.add_argument("--native-module-directory", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    arguments = parser.parse_args()
    try:
        print(canonical_json(run(arguments.paper_text, arguments.native_module_directory, arguments.output_directory)))
    except (OSError, TypeError, ValueError, ImportError) as error:
        print("error: " + str(error), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
