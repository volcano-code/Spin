#!/usr/bin/env python3
"""Capture a frozen full-cage adaptive run; failure is a scientific result.

Run under the project's external 24 GiB, swap=0, 60 minute process guard.
The fixed-partition hollowing branch is only a diagnostic comparator.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import importlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tarfile
import time

from make_it_stand.verification.mass_moments import TriangleMesh, compute_mass_properties
from make_it_stand.verification.run_evidence import (
    EXPLICIT_EFFECTIVE_SOURCE_POLICY, begin_run, check_freshness, finish_run,
)
from make_it_stand.verification.spin_common import canonical

ROOT = Path(__file__).resolve().parents[3]
BINARY = "spin-it-2014-source-cage-asymmetric-joint-regression"


def write(path, value):
    with Path(path).open("x") as stream:
        stream.write(canonical(value) + "\n")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def scientific_file_digests(native_directory: Path) -> tuple[dict[str, str], dict[str, str]]:
    """Separate solver science from replay-only geometry witness attachments.

    The witness preserves the exact rejected state and all BBW coefficients for
    diagnosis, but those coefficients are reconstructed from the cage and may
    differ by sub-ULP linear-algebra noise without changing the solver trace or
    material geometry.  Raw evidence still archives both groups; only the
    former is eligible for process-determinism comparison.
    """
    scientific, diagnostic = {}, {}
    for path in sorted(Path(native_directory).rglob("*")):
        if not path.is_file():
            continue
        name = path.name
        target = diagnostic if (
            name == "primitive_subset_witnesses.jsonl" or
            name.startswith("primitive_subset_callback_")
        ) else scientific
        target[str(path.relative_to(native_directory))] = digest(path)
    return scientific, diagnostic


def load(path):
    return json.loads(Path(path).read_text(), parse_constant=lambda value: (
        (_ for _ in ()).throw(ValueError("nonfinite JSON: " + value))))


def pack_raw_outputs(output):
    """Bind exact bytes even when a killed worker left empty/truncated files."""
    output = Path(output)
    archive = output / "raw_outputs.tar"
    paths = sorted(p for p in output.rglob("*") if p.is_file() and
                   p != archive and p.name != "begin.json")
    with tarfile.open(archive, mode="x") as stream:
        for path in paths:
            stream.add(path, arcname=str(path.relative_to(output)), recursive=False)
    return archive


def read_obj(path):
    vertices, faces = [], []
    with Path(path).open() as stream:
        for line in stream:
            fields = line.split()
            if fields and fields[0] == "v":
                vertices.append(tuple(map(float, fields[1:])))
            elif fields and fields[0] == "f":
                faces.append(tuple(int(value) - 1 for value in fields[1:]))
    return TriangleMesh(tuple(vertices), tuple(faces))


def independent_top_constraints(properties, phi):
    """Frozen input frame origin (.5,.5,.5), identity rotation, root side 1."""
    m, first, raw = properties.mass, properties.first_moment, properties.second_moment
    shifted_first = [x - .5 * m for x in first]
    q = [[raw[i][j] - .5 * first[i] - .5 * first[j] + .25 * m
          for j in range(3)] for i in range(3)]
    c, s = math.cos(phi), math.sin(phi)
    return [v / 1.3 for v in (shifted_first[0], shifted_first[1], q[0][2], q[1][2],
            c * s * (q[0][0] - q[1][1]) + (c * c - s * s) * q[0][1])]


def check_material_output(native):
    native = Path(native)
    result = load(native / "result.json")
    if result["binary_status"] != "BINARY_MATERIAL_VALIDATED_PRE_SEAL":
        return {"status": "NOT_VALIDATED", "native_status": result["binary_status"]}
    if not result.get("binary_candidate_available") or not result.get("binary_moments_available"):
        raise ValueError("validated status lacks admitted candidate and integrated moments")
    if result.get("adaptive_outcome") != 0 or result.get("solver_status") != "CONVERGED":
        raise ValueError("validated status lacks completed converged adaptive solve")
    if result.get("sealed") is not False or result.get("original_paper_case") is not False:
        raise ValueError("synthetic batch must not claim sealed material or an original case")
    for name, length in (("binary_moments", 10), ("binary_moment_residuals", 10),
                         ("binary_normalized_constraints", 5)):
        values = result.get(name)
        if not isinstance(values, list) or len(values) != length or any(
                type(v) not in (int, float) or not math.isfinite(v) for v in values):
            raise ValueError("missing or nonfinite " + name)
    if max(map(abs, result["binary_normalized_constraints"])) > 1e-6:
        raise ValueError("binary constraint gate failed")
    properties = compute_mass_properties(read_obj(native / "binary_material.obj"), density=1.3)
    q = properties.second_moment
    moments = [properties.mass, *properties.first_moment,
               q[0][0], q[1][1], q[2][2], q[1][2], q[0][2], q[0][1]]
    # Frozen input has root side 1 and density 1.3. Never fit scale to output.
    error = max(abs(a - b) / 1.3 for a, b in zip(moments, result["binary_moments"]))
    if error > 1e-6:
        raise ValueError("exported material disagrees with independent moments")
    constraints = independent_top_constraints(properties,
        result["final_state"][result["final_beta_count"]])
    if any(not math.isfinite(v) for v in constraints) or max(map(abs, constraints)) > 1e-6:
        raise ValueError("exported geometry fails independently recomputed axis constraints")
    return {"status": "INDEPENDENT_EXPORTED_MOMENTS_VALIDATED_PRE_SEAL",
            "normalized_moment_error": error, "mass_properties": asdict(properties),
            "independent_normalized_constraints": constraints}


def run(output, build, dependency_roots):
    output, build = Path(output).absolute(), Path(build).resolve()
    output.mkdir(exist_ok=False)
    binary = build / BINARY
    sys.path.insert(0, str(build / "python"))
    native = importlib.import_module("_make_it_stand")
    if Path(native.__file__).resolve().parent != build / "python":
        raise ValueError("wrong native extension loaded")
    sources = [ROOT / "tests/spin_it_2014_source_cage_asymmetric_joint_regression.cpp",
               ROOT / "tests/source_surface_test_utils.h", ROOT / "CMakeLists.txt"]
    # Explicit effective-source policy: include all native scientific sources,
    # including untracked headers; documentation/video changes do not stale it.
    for directory in ("src", "include"):
        sources.extend(sorted(p for p in (ROOT / directory).rglob("*")
                              if p.suffix in (".h", ".hpp", ".cpp", ".cc")))
    sources += [ROOT / "experiments/spin_it_reproduction/strict/source_cage_joint_batch.h",
                ROOT / "python/make_it_stand/verification/mass_moments.py",
                Path(sys.modules[begin_run.__module__].__file__).resolve(),
                Path(sys.modules[canonical.__module__].__file__).resolve()]
    begin = begin_run(repository_root=ROOT, runner_files=[Path(__file__)],
        checker_files=sources, asset_files=[binary],
        build_configuration_files=[ROOT / "CMakeLists.txt", build / "CMakeCache.txt"],
        parameters={"protocol": "source-cage-joint-batch-v1", "native_binary_sha256": digest(binary),
                    "classification": "SYNTHETIC_REGRESSION", "timeout_seconds": 3500},
        seed_protocol={"seed": 0, "protocol": "compiled frozen initial state; no restarts"},
        dependency_roots=dependency_roots, native_modules=[native], pure_python=False,
        source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY)
    write(output / "begin.json", begin.as_dict())
    started = time.monotonic()
    events = []
    try:
        with (output / "native.stdout").open("x") as transcript:
            completed = subprocess.run([str(binary), "--batch-run", str(output / "native")],
                stdout=transcript, stderr=subprocess.STDOUT, cwd=ROOT, timeout=3500)
        if completed.returncode:
            raise RuntimeError(f"native runner exited {completed.returncode}")
        result = load(output / "native/result.json")
        trace_rows = 0
        for path in sorted((output / "native").glob("*.jsonl")):
            with path.open() as stream:
                for line in stream:
                    json.loads(line)
                    trace_rows += 1
        checked = check_material_output(output / "native")
        write(output / "independent_check.json", checked)
        content, diagnostic_content = scientific_file_digests(output / "native")
        scientific_hash = hashlib.sha256(canonical(content).encode()).hexdigest()
        summary = {"status": result["binary_status"], "classification": "SYNTHETIC_REGRESSION",
            "project_status": "PARTIAL_BLOCKED", "solver_status": result.get("solver_status"),
            "scientific_sha256": scientific_hash, "scientific_files": content,
            "diagnostic_files": diagnostic_content,
            "parsed_event_rows": trace_rows, "independent_check": checked}
        write(output / "scientific.json", summary)
        events.append({"iteration": 0, "status": "ok" if result["binary_status"] ==
                       "BINARY_MATERIAL_VALIDATED_PRE_SEAL" else "failed", "data": summary})
    except Exception as error:
        summary = {"status": "EXECUTION_FAILED", "error": str(error)}
        write(output / "failure.json", summary)
        events.append({"iteration": 0, "status": "failed", "data": summary})
    archive = pack_raw_outputs(output)
    evidence = finish_run(begin, events=events, output_files=[archive],
        measurements={"elapsed_seconds": time.monotonic() - started})
    write(output / "evidence.json", evidence.as_dict())
    freshness = check_freshness(evidence)
    summary = {**summary, "freshness": freshness.status, "mismatches": list(freshness.mismatches)}
    write(output / "summary.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--build-directory", type=Path, required=True)
    parser.add_argument("--dependency-root", action="append", required=True)
    args = parser.parse_args()
    roots = {}
    for item in args.dependency_root:
        name, sep, value = item.partition("=")
        if not name or not sep or name in roots:
            parser.error("expected unique NAME=PATH")
        roots[name] = Path(value)
    result = run(args.output_directory, args.build_directory, roots)
    print(canonical({k: v for k, v in result.items() if k != "scientific_files"}))
    return int(result["status"] == "EXECUTION_FAILED" or result["freshness"] != "CURRENT")


if __name__ == "__main__":
    raise SystemExit(main())
