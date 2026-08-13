"""Strict contract tests for the R10 presentation report generator."""

from __future__ import annotations

import csv
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

from make_it_stand.r10 import (
    CANONICALIZATION,
    MANIFEST_SCHEMA,
    RunRecorder,
    aggregate_run_records,
    bind_artifact,
    read_canonical_json,
    write_canonical_json,
)


REPOSITORY = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "r10_report", REPOSITORY / "python" / "r10_report.py"
)
assert SPEC is not None and SPEC.loader is not None
REPORT = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = REPORT
SPEC.loader.exec_module(REPORT)


def _hardware() -> dict[str, object]:
    return {
        "cpu_model": "Test CPU",
        "logical_cpu_count": 4,
        "memory_bytes": 8_589_934_592,
        "os": "TestOS 1",
        "machine": "x86_64",
    }


def _build() -> dict[str, object]:
    return {
        "build_type": "Release",
        "compiler": "GNU",
        "compiler_version": "13.3.0",
        "commit_sha": "a" * 40,
        "dirty_worktree": False,
        "flags": ["-O3", "-DNDEBUG"],
    }


def _make_case(evidence_root: Path, *, paper_metric: bool = True) -> Path:
    case_root = evidence_root / "synthetic-report-case"
    (case_root / "inputs").mkdir(parents=True)
    (case_root / "provenance").mkdir()
    (case_root / "inputs" / "fixture.obj").write_text(
        "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n", encoding="utf-8"
    )
    (case_root / "provenance" / "protocol.txt").write_text(
        "repository-defined synthetic protocol\n", encoding="utf-8"
    )
    references = []
    if paper_metric:
        references.append(
            {
                "name": "reported_iteration_count",
                "value": 987654.0,
                "unit": "iterations",
                "source": "paper_reference",
                "semantics": "Iteration count transcribed from the cited paper.",
                "paper_location": "Table 1, selected row",
                "comparability": "reference_only",
            }
        )
    definitions = [
        {
            "name": "projected_target_energy",
            "unit": "normalized_objective",
            "direction": "lower_is_better",
            "definition": "Final projected target energy.",
            "source": "local_measurement",
        },
        {
            "name": "removed_volume_fraction",
            "unit": "fraction",
            "direction": "higher_is_better",
            "definition": "Final removed volume fraction.",
            "source": "local_measurement",
        },
        {
            "name": "normalized_rms_displacement",
            "unit": "normalized_length",
            "direction": "lower_is_better",
            "definition": "Normalized RMS displacement.",
            "source": "local_measurement",
        },
        *(
            {
                "name": name,
                "unit": "indicator",
                "direction": "higher_is_better",
                "definition": f"Measured {name.replace('_', ' ')} indicator.",
                "source": "local_measurement",
            }
            for name in ("algorithm_success", "targets_reached", "manufacturing_ready")
        ),
    ]
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "schema_version": 1,
        "canonicalization": CANONICALIZATION,
        "experiment_id": "synthetic-report-case",
        "title": "Synthetic report contract case",
        "claim": {
            "level": "local_robustness",
            "statement": "This is a repository-defined synthetic test case.",
            "exact_reproduction_claim": False,
            "limitations": ["This fixture is not a published paper artifact."],
        },
        "paper_reference": {
            "citation_key": "prevost2013makeitstand",
            "scope": "Reference-only context; never pooled with local data.",
            "reported_metrics": references,
        },
        "protocol": {
            "variant": "full",
            "repeats": 2,
            "seed_policy": {"kind": "fixed_sequence", "first_seed": 1},
            "preprocessing": {"fixture": "synthetic"},
            "algorithm": {"carving": True, "deformation": True, "laplacian": True},
            "metric_definitions": definitions,
        },
        "provenance": {
            "repository_commit_sha": "a" * 40,
            "dirty_worktree": False,
            "implementation_profile": "release-headless-r10-test",
            "resources": [
                bind_artifact(
                    case_root,
                    "provenance/protocol.txt",
                    artifact_id="protocol",
                    role="protocol_reference",
                    media_type="text/plain",
                )
            ],
            "deviations": ["The fixture and grid are local rather than paper-identical."],
        },
        "inputs": [
            bind_artifact(
                case_root,
                "inputs/fixture.obj",
                artifact_id="fixture",
                role="model_input",
                media_type="model/obj",
            )
        ],
    }
    write_canonical_json(case_root / "manifest.json", manifest)
    records = []
    for index, energy in enumerate((1.0e-2, 5.0e-3)):
        run_root = case_root / "runs" / f"r{index:02d}"
        run_root.mkdir(parents=True)
        diagnostics = {
            "r10_events": [
                {"ordinal": 0, "projected_target_energy": 1.0e-1},
                {"ordinal": 1, "projected_target_energy": energy},
            ]
        }
        (run_root / "diagnostics.json").write_text(
            json.dumps(diagnostics, allow_nan=False, separators=(",", ":"), sort_keys=True),
            encoding="utf-8",
        )
        recorder = RunRecorder(
            manifest,
            run_id=f"synthetic-report-case-r{index:02d}",
            repeat_index=index,
            seed=index + 1,
            build=_build(),
            hardware=_hardware(),
            timing_scope="synthetic-test",
        )
        with recorder.phase("optimization"):
            sum(range(10))
        measurements = {
            "projected_target_energy": energy,
            "removed_volume_fraction": 0.2 + 0.1 * index,
            "normalized_rms_displacement": 0.01 + 0.01 * index,
            "algorithm_success": 1.0,
            "targets_reached": 1.0,
            "manufacturing_ready": float(index == 0),
        }
        semantics = {item["name"]: item["definition"] for item in definitions}
        for name, value in measurements.items():
            recorder.add_local_metric(name, value, semantics=semantics[name])
        recorder.add_artifact(
            bind_artifact(
                case_root,
                f"runs/r{index:02d}/diagnostics.json",
                artifact_id=f"diagnostics-r{index:02d}",
                role="diagnostics",
                media_type="application/json",
            )
        )
        record = recorder.finish(status="succeeded")
        write_canonical_json(case_root / "records" / f"r{index:02d}.json", record)
        records.append(record)
    write_canonical_json(case_root / "aggregate.json", aggregate_run_records(manifest, records))
    return case_root


class R10ReportTests(unittest.TestCase):
    def test_report_is_deterministic_portable_and_keeps_sources_separate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "evidence"
            _make_case(root)
            output = Path(temporary) / "report"
            self.assertEqual(
                REPORT.generate_report(root, output),
                ("summary.csv", "summary.md", "summary.svg"),
            )
            first = {path.name: path.read_bytes() for path in output.iterdir()}
            REPORT.generate_report(root, output)
            self.assertEqual(first, {path.name: path.read_bytes() for path in output.iterdir()})

            for payload in first.values():
                self.assertNotIn(str(root).encode(), payload)
            markdown = first["summary.md"].decode()
            self.assertIn("not an exact paper reproduction", markdown)
            self.assertIn("## Local measurements", markdown)
            self.assertIn("## Paper-reference values", markdown)
            self.assertIn("reference_only", markdown)

            rows = list(csv.DictReader(io.StringIO(first["summary.csv"].decode())))
            self.assertEqual({row["source"] for row in rows}, {"local_measurement", "paper_reference"})
            for row in rows:
                if row["source"] == "local_measurement":
                    self.assertNotEqual(row["median"], "")
                    self.assertEqual(row["reference_value"], "")
                    self.assertEqual(row["comparability"], "")
                else:
                    self.assertEqual(row["median"], "")
                    self.assertEqual(row["reference_value"], "987654")
                    self.assertEqual(row["comparability"], "reference_only")

            ET.fromstring(first["summary.svg"])
            svg = first["summary.svg"].decode()
            self.assertIn("local synthetic evidence", svg)
            self.assertNotIn("reported_iteration_count", svg)
            self.assertNotIn("987654", svg)

    def test_stale_but_canonical_aggregate_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "evidence"
            case_root = _make_case(root, paper_metric=False)
            aggregate = read_canonical_json(case_root / "aggregate.json")
            aggregate["local_metric_aggregates"][0]["median"] = 0.75
            aggregate["local_metric_aggregates"][0]["q1"] = 0.75
            aggregate["local_metric_aggregates"][0]["q3"] = 0.75
            aggregate["local_metric_aggregates"][0]["iqr"] = 0.0
            aggregate["local_metric_aggregates"][0]["minimum"] = 0.75
            aggregate["local_metric_aggregates"][0]["maximum"] = 0.75
            write_canonical_json(case_root / "aggregate.json", aggregate)
            with self.assertRaisesRegex(REPORT.R10ReportError, "stale or inconsistent"):
                REPORT.generate_report(root, Path(temporary) / "report")

    def test_tampered_bound_diagnostics_is_rejected_before_reporting(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "evidence"
            case_root = _make_case(root, paper_metric=False)
            (case_root / "runs" / "r00" / "diagnostics.json").write_text(
                '{"r10_events":[]}', encoding="utf-8"
            )
            with self.assertRaisesRegex(REPORT.R10ReportError, "SHA-256 mismatch"):
                REPORT.generate_report(root, Path(temporary) / "report")

    def test_empty_evidence_root_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(REPORT.R10ReportError, "no R10 case"):
                REPORT.generate_report(temporary, Path(temporary) / "report")


if __name__ == "__main__":
    unittest.main(verbosity=2)
