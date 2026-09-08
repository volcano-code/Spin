#!/usr/bin/env python3
"""Run the bounded exact checks for Faster Section 6.6."""

from __future__ import annotations

import argparse
import base64
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import time

from make_it_stand.verification import faster_section_6_6 as exact


ROOT = Path(__file__).resolve().parents[3]


def _git(*arguments: str) -> bytes:
    try:
        return subprocess.run(
            ["git", "-C", str(ROOT), *arguments],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
    except (OSError, subprocess.CalledProcessError) as error:
        raise ValueError("cannot capture Git context") from error


def _git_context() -> dict[str, object]:
    comparisons = {
        "head_to_worktree": ("diff", "--binary", "HEAD"),
        "head_to_index": ("diff", "--binary", "--cached", "HEAD"),
        "index_to_worktree": ("diff", "--binary"),
    }
    result: dict[str, object] = {
        "head": _git("rev-parse", "--verify", "HEAD").decode("ascii").strip(),
        "changed_paths": [
            line[3:] for line in _git("status", "--short", "--untracked-files=all")
            .decode("utf-8").splitlines() if len(line) >= 4
        ],
    }
    for label, command in comparisons.items():
        payload = _git(*command)
        result[label + "_diff_base64"] = base64.b64encode(payload).decode("ascii")
        result[label + "_diff_sha256"] = exact.scientific_sha256(
            {"base64": result[label + "_diff_base64"]}
        )
    return result


def _human_report(scientific: dict[str, object], spec_record: dict[str, object],
                  pdf_record: dict[str, object], digest: str) -> str:
    nonunique, full = scientific["cases"]
    conditions = []
    for case in (nonunique, full):
        conditions.extend(
            f"- `{case['source_case_id']}` `{name}`: `{str(value).lower()}`"
            for name, value in case["conditions"].items()
        )
    return "\n".join([
        "# Faster Section 6.6 exact formula evidence",
        "",
        "This is a pure-Python exact-rational check of only the two published one-dimensional examples.",
        "It is not a 3D proxy, Newton convergence run, CAD/physical validation, complete paper reproduction,",
        "or a classification of all problems with p=0.",
        "",
        "## Source binding and local comparison",
        "",
        f"- Local specification SHA256: `{spec_record['sha256']}` (required `{exact.EXPECTED_SPEC_SHA256}`; equal)",
        f"- Local published PDF SHA256: `{pdf_record['sha256']}` (required `{exact.EXPECTED_PDF_SHA256}`; equal)",
        f"- Locator: `{scientific['paper_locator']}`",
        "- The immutable historical specification remains NOT_VERIFIED / NOT_RUN; this is separate evidence.",
        "",
        "## Exact equations",
        "",
        "For Omega=(0,1), chi_tau(x)=1/2+tau(2x-1), tau in [-1/2,1/2].",
        f"Exact integration gives v(chi_tau)={nonunique['volume_polynomial_in_tau'][0]}+"
        f"{nonunique['volume_polynomial_in_tau'][1]} tau. The four bilinear rectangle corners are "
        f"{', '.join(nonunique['rectangle_corner_values'])}, proving 0<=chi<=1 by the corner-extremum argument.",
        f"Also integral_0^1 (2x-1)^2 dx={nonunique['parameter_derivative_l2_squared']}>0, so distinct real tau values",
        "give distinct functions modulo a.e. equality; the real interval therefore supplies uncountably many optima.",
        f"The two fixed binary indicators integrate to {', '.join(nonunique['binary_half_interval_integrals'])}, and chi=1/2",
        f"is nonbinary with integral {nonunique['constant_nonbinary_integral']}.",
        "For shifted binary intervals (t,t+1/2), 0<=t<=1/2. For 0<=s<t<=1/2, overlap length is",
        f"{nonunique['shifted_binary_overlap_length']} and symmetric-difference measure is "
        f"{nonunique['shifted_binary_symmetric_difference']}>0. This real interval is an uncountable binary family modulo a.e.",
        "",
        "For Omega=(0,1/2), chi=1 integrates to 1/2. At any optimum the nonnegative deficit 1-chi",
        "has integral |Omega|-1/2=0; the nonnegative-integrable-zero-integral lemma therefore gives chi=1 a.e.,",
        "which is uniqueness almost everywhere, not pointwise uniqueness.",
        "",
        "In both cases the squared objective has global lower bound 0, the witnesses attain 0,",
        "2(v-1/2)=0, and p=-2(v-1/2)=0. No normalized nonzero quadric is demanded at p=0.",
        "",
        "## Independently evaluated proof conditions",
        "",
        *conditions,
        "",
        "## Status and limits",
        "",
        "- Both scoped scientific outcomes: `VERIFIED_ANALYTIC_PROPERTIES`",
        "- Global project: `PARTIAL_BLOCKED`",
        "- Author bitwise exact: `NOT_VERIFIED`",
        "- Physical validation: `NOT_IN_SCOPE`",
        "- Build/source attestation: `NOT_APPLICABLE_PURE_PYTHON_FORMULA_PROOF`",
        "- Aggregate verification status unchanged; no manual positive release gate opened.",
        f"- Deterministic scientific-data SHA256: `{digest}`",
        "",
    ])


def run(spec_path: Path, paper_path: Path, output_directory: Path) -> dict[str, object]:
    if output_directory.exists():
        raise ValueError("output directory already exists: " + str(output_directory))
    spec_record = exact.file_record(spec_path)
    pdf_record = exact.file_record(paper_path)
    if spec_record["sha256"] != exact.EXPECTED_SPEC_SHA256:
        raise ValueError("specification SHA256 mismatch")
    if pdf_record["sha256"] != exact.EXPECTED_PDF_SHA256:
        raise ValueError("PDF SHA256 mismatch")

    started = datetime.now(timezone.utc).isoformat()
    start = time.perf_counter()
    spec = exact.load_spec(spec_path)
    scientific = exact.verify_published_section_6_6(
        spec, exact.published_witness_definitions()
    )
    scientific_digest = exact.scientific_sha256(scientific)

    output_directory.mkdir(parents=True, exist_ok=False)
    result_path = output_directory / "result.json"
    report_path = output_directory / "report.md"
    result = {
        "schema": "faster-section-6-6-exact-result-v1",
        "scientific_data": scientific,
        "scientific_data_sha256": scientific_digest,
        "source_binding": {
            "spec_sha256": spec_record["sha256"],
            "paper_pdf_sha256": pdf_record["sha256"],
        },
    }
    result_path.write_text(exact.canonical_json(result) + "\n", encoding="utf-8")
    report_path.write_text(
        _human_report(scientific, spec_record, pdf_record, scientific_digest),
        encoding="utf-8",
    )

    runner_record = exact.file_record(Path(__file__))
    checker_record = exact.file_record(Path(exact.__file__))
    cmake_record = exact.file_record(ROOT / "CMakeLists.txt")
    inputs = {
        "runner": runner_record,
        "checker": checker_record,
        "spec": spec_record,
        "paper_pdf": pdf_record,
        "cmake_configuration": cmake_record,
    }
    evidence = {
        "schema": "faster-section-6-6-exact-evidence-v1",
        "inputs": inputs,
        "source_content_sha256": exact.scientific_sha256(inputs),
        "outputs": {
            "result": exact.file_record(result_path),
            "human_report": exact.file_record(report_path),
        },
        "git": _git_context(),
        "runtime": {
            "python_executable": str(Path(sys.executable).resolve()),
            "python_version": sys.version,
            "native_extensions_required": [],
            "started_utc": started,
            "elapsed_seconds": time.perf_counter() - start,
        },
        "scientific_data_sha256": scientific_digest,
        "project_status": "PARTIAL_BLOCKED",
        "author_bitwise_exact": "NOT_VERIFIED",
        "physical_validation": "NOT_IN_SCOPE",
        "build_source_attestation": "NOT_APPLICABLE_PURE_PYTHON_FORMULA_PROOF",
    }
    evidence_path = output_directory / "evidence.json"
    evidence_path.write_text(exact.canonical_json(evidence) + "\n", encoding="utf-8")
    return {
        "output_directory": str(output_directory.resolve()),
        "scientific_data_sha256": scientific_digest,
        "freshness": exact.check_evidence_freshness(evidence_path)["status"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--paper", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    arguments = parser.parse_args()
    try:
        summary = run(arguments.spec, arguments.paper, arguments.output_directory)
    except (OSError, TypeError, ValueError) as error:
        print("error: " + str(error), file=sys.stderr)
        return 2
    print(exact.canonical_json(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
