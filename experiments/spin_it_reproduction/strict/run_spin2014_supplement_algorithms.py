#!/usr/bin/env python3
"""Create source-bound local evidence for Spin-It 2014 Supplement Algorithms 1--2.

This runner checks a direct Python transcription of the published surface
polynomials against the project's native implementation.  It deliberately
does not grant any original-model, geometry-legality, or release conclusion.
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

from make_it_stand.verification import spin2014_supplement_algorithms as supplement


ROOT = Path(__file__).resolve().parents[3]
EXPECTED_PDF_SHA256 = "f0acc779410bacac74ce9fbce305fe6dd9e1ab570885141644abac805a1396e5"
EXPECTED_MAIN_TEXT_SHA256 = "5cf5ab047189575df48b01796de73c72ceadc602a3380b74f33d2dda81459161"
EVIDENCE_SCHEMA = "spin2014-supplement-algorithms-evidence-v1"
RESULT_SCHEMA = "spin2014-supplement-algorithms-result-v1"


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
    before = _sha256(path)
    size = path.stat().st_size
    after = _sha256(path)
    if before != after:
        raise ValueError("file changed while hashing: " + str(path))
    return {"path": str(path), "sha256": after, "size": size}


def _record_matches(record: object, label: str) -> list[str]:
    if not isinstance(record, dict):
        return [label + " record is malformed"]
    try:
        path = Path(record["path"]).resolve(strict=True)
    except (KeyError, OSError, TypeError, ValueError):
        return [label + " path is unavailable"]
    current = file_record(path)
    return ([] if current == record else [label + " record differs from current bytes"])


def _fixture_data() -> list[dict[str, object]]:
    tetra_vertices = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                      (0.0, 2.0, 0.0), (0.0, 0.0, 3.0))
    tetra_faces = ((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3))
    translated = tuple((x + 3.5, y - 2.25, z + 1.75) for x, y, z in tetra_vertices)
    return [
        {"id": "asymmetric_tetrahedron", "vertices": tetra_vertices,
         "faces": tetra_faces, "density": 1.0},
        {"id": "density_scaled_asymmetric_tetrahedron", "vertices": tetra_vertices,
         "faces": tetra_faces, "density": 2.7},
        {"id": "translated_asymmetric_tetrahedron", "vertices": translated,
         "faces": tetra_faces, "density": 2.7},
    ]


def _equation4_checks(fixtures: list[dict[str, object]]) -> list[dict[str, object]]:
    """Check Eq. (4)'s exact monomial vector on a hand-derived tetrahedron.

    The hand-derived tetrahedron has volume one.  The explicit paper-order
    vector guards against treating the project's differently ordered raw
    storage as a source-free Equation (4) check.
    """

    base_paper = np.asarray((1.0, 0.25, 0.5, 0.75, 0.1, 0.3, 0.15, 0.1, 0.4, 0.9))
    paper_indices = np.asarray((0, 1, 2, 3, 9, 7, 8, 4, 5, 6), dtype=int)
    checks = []
    for fixture in fixtures:
        if fixture["fixture_id"] not in {
            "asymmetric_tetrahedron", "density_scaled_asymmetric_tetrahedron",
        }:
            continue
        density = float(fixture["density"])
        expected = density * base_paper
        reference = np.asarray(fixture["reference_moments"], dtype=float)[paper_indices]
        native = np.asarray(fixture["native_moments"], dtype=float)[paper_indices]
        reference_error = float(np.max(np.abs(reference - expected)))
        native_error = float(np.max(np.abs(native - expected)))
        cross_error = float(np.max(np.abs(native - reference)))
        if reference_error > 2e-12 or native_error > 2e-12 or cross_error > 2e-12:
            raise ValueError("Eq. (4) density-weighted monomial vector differs from hand-derived tetrahedron")
        checks.append({
            "fixture_id": fixture["fixture_id"], "density": density,
            "paper_term_order": list(supplement.PAPER_ORDER),
            "expected_paper_moments": _as_json(expected),
            "reference_paper_moments": _as_json(reference),
            "native_paper_moments": _as_json(native),
            "reference_exact_max_absolute_error": reference_error,
            "native_exact_max_absolute_error": native_error,
            "reference_native_max_absolute_error": cross_error,
        })
    return checks


def _native_module(native_directory: Path):
    directory = native_directory.resolve(strict=True)
    if not directory.is_dir():
        raise ValueError("native module directory is not a directory")
    extension = next(iter(sorted(directory.glob("_make_it_stand*.so"))), None)
    if extension is None:
        raise ValueError("native module is missing from " + str(directory))
    sys.path.insert(0, str(directory))
    try:
        module = importlib.import_module("_make_it_stand")
    finally:
        sys.path.pop(0)
    if Path(module.__file__).resolve() != extension.resolve():
        raise ValueError("loaded native module differs from declared module path")
    return module, extension


def _as_json(array: np.ndarray) -> list[object]:
    return np.asarray(array, dtype=float).tolist()


def _scientific_data(native_directory: Path) -> tuple[dict[str, object], Path]:
    native, native_path = _native_module(native_directory)
    fixtures = []
    for fixture in _fixture_data():
        vertices = np.asarray(fixture["vertices"], dtype=float)
        faces = fixture["faces"]
        density = float(fixture["density"])
        reference_moments = supplement.integrate_surface_moments(vertices, faces, density)
        reference_jacobian = supplement.differentiate_surface_moments(vertices, faces, density)
        native_result = native.differentiate_surface_raw_moments(vertices, faces, density)
        if tuple(native_result.get("order", ())) != supplement.ORDER or \
                native_result.get("algorithm") != supplement.ALGORITHM:
            raise ValueError("native order or algorithm identifier differs from published contract")
        native_moments = np.asarray(native_result["moments"], dtype=float)
        native_jacobian = np.asarray(native_result["jacobian"], dtype=float)
        if native_moments.shape != (10,) or native_jacobian.shape != reference_jacobian.shape:
            raise ValueError("native output shape differs from reference contract")
        moment_error = float(np.max(np.abs(reference_moments - native_moments)))
        derivative_error = float(np.max(np.abs(reference_jacobian - native_jacobian)))
        if moment_error > 2e-12 or derivative_error > 4e-12:
            raise ValueError("reference and native Supplement Algorithm outputs disagree")
        first = faces[0]
        trace = supplement.triangle_trace(*(vertices[index] for index in first))
        fixtures.append({
            "fixture_id": fixture["id"], "vertices": _as_json(vertices),
            "faces": [list(face) for face in faces], "density": density,
            "reference_moments": _as_json(reference_moments),
            "native_moments": _as_json(native_moments),
            "reference_jacobian": _as_json(reference_jacobian),
            "native_jacobian": _as_json(native_jacobian),
            "moment_max_absolute_error": moment_error,
            "derivative_max_absolute_error": derivative_error,
            "first_triangle_trace": {
                "normal": _as_json(trace["normal"]),
                "h": [_as_json(item) for item in trace["h"]],
                "paper_unscaled_moments": _as_json(trace["paper_unscaled_moments"]),
                "paper_unscaled_jacobian": _as_json(trace["paper_unscaled_jacobian"]),
            },
        })
    return {
        "paper_locator": "Spin-It (2014), Equation 4; Supplement, Algorithms 1 and 2",
        "algorithm": supplement.ALGORITHM,
        "order": list(supplement.ORDER), "paper_order": list(supplement.PAPER_ORDER),
        "tolerances": {"moment_max_absolute": 2e-12, "derivative_max_absolute": 4e-12},
        "fixtures": fixtures, "equation4_checks": _equation4_checks(fixtures),
        "outcome": "REFERENCE_AND_NATIVE_AGREE",
    }, native_path


def _git_context() -> dict[str, object]:
    try:
        head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--verify", "HEAD"],
                              check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True).stdout.strip()
        status = subprocess.run(["git", "-C", str(ROOT), "status", "--short", "--untracked-files=all"],
                                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                text=True).stdout.splitlines()
        diffs = {}
        for name, arguments in {
            "head_to_worktree": ("diff", "--binary", "HEAD"),
            "head_to_index": ("diff", "--binary", "--cached", "HEAD"),
            "index_to_worktree": ("diff", "--binary"),
        }.items():
            payload = subprocess.run(["git", "-C", str(ROOT), *arguments], check=True,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
            diffs[name + "_sha256"] = hashlib.sha256(payload).hexdigest()
    except (OSError, subprocess.CalledProcessError) as error:
        raise ValueError("cannot capture Git context") from error
    return {"head": head, "changed_paths": [line[3:] for line in status if len(line) >= 4],
            "diff_sha256": diffs}


def run(paper: Path, paper_text: Path, native_module_directory: Path,
        output_directory: Path) -> dict[str, object]:
    if output_directory.exists():
        raise ValueError("output directory already exists: " + str(output_directory))
    paper_record = file_record(paper)
    if paper_record["sha256"] != EXPECTED_PDF_SHA256:
        raise ValueError("published supplement PDF SHA-256 mismatch")
    paper_text_record = file_record(paper_text)
    if paper_text_record["sha256"] != EXPECTED_MAIN_TEXT_SHA256:
        raise ValueError("Spin-It 2014 local text extraction SHA-256 mismatch")
    started = datetime.now(timezone.utc).isoformat()
    start = time.perf_counter()
    scientific, native_path = _scientific_data(native_module_directory)
    scientific_digest = hashlib.sha256(canonical_json(scientific).encode("utf-8")).hexdigest()
    output_directory.mkdir(parents=True, exist_ok=False)
    result_path = output_directory / "result.json"
    result = {"schema": RESULT_SCHEMA, "scientific_data": scientific,
              "scientific_data_sha256": scientific_digest}
    result_path.write_text(canonical_json(result) + "\n", encoding="utf-8")
    build_root = native_path.resolve().parents[1]
    build_cache = build_root / "CMakeCache.txt"
    if not build_cache.is_file():
        raise ValueError("native module build directory lacks CMakeCache.txt")
    inputs = {
        "runner": file_record(Path(__file__)),
        "checker": file_record(Path(supplement.__file__)),
        "paper_pdf": paper_record,
        "paper_text": paper_text_record,
        "native_module": file_record(native_path),
        "build_cmake_cache": file_record(build_cache),
        "source_cmake_lists": file_record(ROOT / "CMakeLists.txt"),
    }
    evidence = {
        "schema": EVIDENCE_SCHEMA,
        "inputs": inputs,
        "outputs": {"result": file_record(result_path)},
        "git": _git_context(),
        "runtime": {"python_executable": str(Path(sys.executable).resolve()),
                    "python_version": sys.version, "started_utc": started,
                    "elapsed_seconds": time.perf_counter() - start},
        "scientific_data_sha256": scientific_digest,
        "project_status": "PARTIAL_BLOCKED",
        "author_bitwise_exact": "NOT_VERIFIED",
        "physical_validation": "NOT_IN_SCOPE",
        "build_source_attestation": "NOT_ESTABLISHED_NATIVE_BINARY_BOUND_ONLY",
    }
    evidence_path = output_directory / "evidence.json"
    evidence_path.write_text(canonical_json(evidence) + "\n", encoding="utf-8")
    return {"output_directory": str(output_directory.resolve()),
            "scientific_data_sha256": scientific_digest, "status": "CURRENT"}


def check_evidence_freshness(evidence_path: Path) -> dict[str, object]:
    try:
        evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("cannot load supplement evidence") from error
    mismatches: list[str] = []
    if evidence.get("schema") != EVIDENCE_SCHEMA:
        mismatches.append("evidence schema differs")
    inputs = evidence.get("inputs")
    if not isinstance(inputs, dict):
        mismatches.append("input records are malformed")
    else:
        for label in ("runner", "checker", "paper_pdf", "native_module",
                      "paper_text", "build_cmake_cache", "source_cmake_lists"):
            mismatches.extend(_record_matches(inputs.get(label), label))
        if inputs.get("paper_pdf", {}).get("sha256") != EXPECTED_PDF_SHA256:
            mismatches.append("published supplement PDF digest differs")
        if inputs.get("paper_text", {}).get("sha256") != EXPECTED_MAIN_TEXT_SHA256:
            mismatches.append("Spin-It 2014 local text extraction digest differs")
    outputs = evidence.get("outputs")
    if not isinstance(outputs, dict):
        mismatches.append("output records are malformed")
    else:
        mismatches.extend(_record_matches(outputs.get("result"), "result"))
    if not mismatches:
        # ``git`` is an observed provenance snapshot only.  Freshness is
        # governed by explicit byte-bound scientific inputs, so unrelated
        # documentation changes cannot stale numerical evidence.
        result_path = Path(outputs["result"]["path"])
        result = json.loads(result_path.read_text(encoding="utf-8"))
        digest = hashlib.sha256(canonical_json(result.get("scientific_data")).encode("utf-8")).hexdigest()
        if result.get("schema") != RESULT_SCHEMA or digest != result.get("scientific_data_sha256") or \
                digest != evidence.get("scientific_data_sha256"):
            mismatches.append("scientific data digest differs")
    return {"status": "CURRENT" if not mismatches else "STALE_RESULT", "mismatches": mismatches}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--paper", type=Path, required=True)
    parser.add_argument("--paper-text", type=Path, required=True)
    parser.add_argument("--native-module-directory", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    arguments = parser.parse_args()
    try:
        print(canonical_json(run(arguments.paper, arguments.paper_text,
                                 arguments.native_module_directory,
                                 arguments.output_directory)))
    except (OSError, TypeError, ValueError, ImportError) as error:
        print("error: " + str(error), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
