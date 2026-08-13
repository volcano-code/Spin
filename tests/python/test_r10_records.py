"""Deterministic tests for the R10 experiment evidence contract."""

from __future__ import annotations

import copy
import hashlib
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


def _metric(name: str, value: float) -> dict[str, object]:
    return {
        "name": name,
        "value": value,
        "unit": "iterations",
        "source": "paper_reference",
        "semantics": "iteration count printed in the cited paper table",
        "paper_location": "Table 1, selected row",
        "comparability": "reference_only",
    }


def _binding(identifier: str, path: str, digest: str = "1" * 64) -> dict[str, object]:
    return {
        "artifact_id": identifier,
        "relative_path": path,
        "sha256": digest,
        "byte_count": 12,
        "role": "model_input",
        "media_type": "model/obj",
    }


def _manifest(*, repeats: int = 2) -> dict[str, object]:
    return {
        "schema": "make-it-stand-r10-experiment-manifest-v1",
        "schema_version": 1,
        "canonicalization": "utf8-json-sorted-keys-compact-v1",
        "experiment_id": "synthetic-full-r10",
        "title": "Synthetic paper-governed full pipeline check",
        "claim": {
            "level": "paper_governed",
            "statement": "The numerical policy follows stated paper rules on a synthetic model.",
            "exact_reproduction_claim": False,
            "limitations": [
                "The synthetic model and current implementation do not reproduce the published artifact."
            ],
        },
        "paper_reference": {
            "citation_key": "prevost2013makeitstand",
            "scope": "Numerical-policy context only; values are not pooled with local data.",
            "reported_metrics": [_metric("reported_iteration_count", 128)],
        },
        "protocol": {
            "variant": "full",
            "repeats": repeats,
            "seed_policy": {"kind": "fixed_sequence", "first_seed": 31000},
            "preprocessing": {
                "normalization": "center bounding box and scale longest extent to 0.95"
            },
            "algorithm": {
                "carving": True,
                "deformation": True,
                "laplacian": True,
                "initial_step": 1.0,
                "warmup_runs_per_measured_repeat": 1,
                "native_extension_artifact_id": "native-extension",
            },
            "metric_definitions": [
                {
                    "name": "projected_com_error",
                    "unit": "normalized_length",
                    "direction": "lower_is_better",
                    "definition": "Maximum final normalized signed support violation.",
                    "source": "local_measurement",
                },
                {
                    "name": "final_cavity_volume_fraction",
                    "unit": "fraction",
                    "direction": "descriptive",
                    "definition": "Final cavity volume divided by final outer-envelope volume.",
                    "source": "local_measurement",
                },
            ],
        },
        "provenance": {
            "repository_commit_sha": "a" * 40,
            "dirty_worktree": False,
            "implementation_profile": "release-headless-r10",
            "resources": [
                {
                    "artifact_id": "paper-protocol-notes",
                    "relative_path": "provenance/protocol.txt",
                    "sha256": "2" * 64,
                    "role": "protocol_reference",
                    "media_type": "text/plain",
                },
                {
                    "artifact_id": "native-extension",
                    "relative_path": "provenance/native_extension.bin",
                    "sha256": "9" * 64,
                    "byte_count": 4096,
                    "role": "native_extension",
                    "media_type": "application/x-sharedlib",
                }
            ],
            "deviations": [
                "The adaptive local grid is not the publication's fixed cubic grid."
            ],
        },
        "inputs": [_binding("synthetic-model", "inputs/synthetic.obj")],
    }


def _hardware() -> dict[str, object]:
    return {
        "cpu_model": "Test CPU",
        "logical_cpu_count": 4,
        "memory_bytes": 8589934592,
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


def _run_record(api: object, manifest: dict[str, object], index: int, value: float) -> dict[str, object]:
    return {
        "schema": "make-it-stand-r10-run-record-v1",
        "schema_version": 1,
        "canonicalization": "utf8-json-sorted-keys-compact-v1",
        "experiment_id": manifest["experiment_id"],
        "run_id": f"synthetic-r{index:02d}",
        "repeat_index": index,
        "seed": 31000 + index,
        "manifest_sha256": api.document_sha256(manifest),
        "claim_level": "paper_governed",
        "status": "succeeded",
        "execution": {
            "hardware": _hardware(),
            "build": _build(),
            "timing": {
                "clock": "monotonic_performance_counter",
                "scope": "model_setup_start-through_run_record_finalization",
                "warmup_runs": 1,
                "phases": [{"name": "optimization", "seconds": 1.0 + index}],
                "total_seconds": 1.5 + index,
            },
        },
        "paper_reference_metrics": manifest["paper_reference"]["reported_metrics"],
        "local_metrics": [
            {
                "name": "projected_com_error",
                "value": value,
                "unit": "normalized_length",
                "source": "local_measurement",
                "semantics": "Maximum final normalized signed support violation.",
            },
            {
                "name": "final_cavity_volume_fraction",
                "value": 0.4 + value,
                "unit": "fraction",
                "source": "local_measurement",
                "semantics": "Final cavity volume divided by final outer-envelope volume.",
            },
        ],
        "artifacts": [
            copy.deepcopy(manifest["provenance"]["resources"][1]),
            {
                "artifact_id": f"diagnostics-r{index:02d}",
                "relative_path": f"runs/r{index:02d}/diagnostics.json",
                "sha256": f"{index + 3:x}" * 64,
                "byte_count": 123 + index,
                "role": "diagnostics",
                "media_type": "application/json",
            }
        ],
        "notes": [],
    }


class R10EvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.api = importlib.import_module("make_it_stand.r10")

    def test_schema_tokens_and_manifest_are_exact(self) -> None:
        value = _manifest()
        self.assertIs(self.api.validate_manifest(value), value)
        self.assertEqual(self.api.MANIFEST_SCHEMA, value["schema"])
        self.assertEqual(self.api.RUN_RECORD_SCHEMA, "make-it-stand-r10-run-record-v1")
        self.assertEqual(self.api.AGGREGATE_SCHEMA, "make-it-stand-r10-aggregate-v1")

    def test_every_claim_level_is_accepted_but_exact_reproduction_is_rejected(self) -> None:
        for level in sorted(self.api.CLAIM_LEVELS):
            with self.subTest(level=level):
                value = _manifest()
                value["claim"]["level"] = level
                self.api.validate_manifest(value)
        value = _manifest()
        value["claim"]["exact_reproduction_claim"] = True
        with self.assertRaisesRegex(self.api.R10RecordError, "exact paper reproduction"):
            self.api.validate_manifest(value)

    def test_manifest_requires_known_deviation_and_bound_input(self) -> None:
        for key in ("deviations",):
            value = _manifest()
            value["provenance"][key] = []
            with self.assertRaises(self.api.R10RecordError):
                self.api.validate_manifest(value)
        value = _manifest()
        value["inputs"] = []
        with self.assertRaisesRegex(self.api.R10RecordError, "SHA-256-bound input"):
            self.api.validate_manifest(value)

    def test_manifest_requires_one_byte_bound_declared_native_extension(self) -> None:
        value = _manifest()
        value["provenance"]["resources"].pop()
        with self.assertRaisesRegex(self.api.R10RecordError, "native extension"):
            self.api.validate_manifest(value)

        value = _manifest()
        value["provenance"]["resources"][1]["role"] = "protocol_reference"
        with self.assertRaisesRegex(self.api.R10RecordError, "native_extension"):
            self.api.validate_manifest(value)

        value = _manifest()
        del value["provenance"]["resources"][1]["byte_count"]
        with self.assertRaisesRegex(self.api.R10RecordError, "byte_count"):
            self.api.validate_manifest(value)

        for invalid in (-1, True):
            value = _manifest()
            value["protocol"]["algorithm"][
                "warmup_runs_per_measured_repeat"
            ] = invalid
            with self.assertRaisesRegex(self.api.R10RecordError, "warmup_runs"):
                self.api.validate_manifest(value)

    def test_paper_and_local_metric_semantics_cannot_be_swapped(self) -> None:
        value = _manifest()
        value["paper_reference"]["reported_metrics"][0]["source"] = "local_measurement"
        with self.assertRaisesRegex(self.api.R10RecordError, "paper_reference"):
            self.api.validate_manifest(value)

        manifest = _manifest()
        record = _run_record(self.api, manifest, 0, 0.1)
        record["local_metrics"][0]["source"] = "paper_reference"
        with self.assertRaisesRegex(self.api.R10RecordError, "local_measurement"):
            self.api.validate_run_record(record, manifest=manifest)

    def test_absolute_host_paths_and_identity_fields_are_rejected_recursively(self) -> None:
        bad_values = [
            ("inputs", "/home/example/model.obj"),
            ("preprocessing", "/workspace/private/model.obj"),
            ("windows", "C:\\Users\\example\\model.obj"),
        ]
        for kind, path in bad_values:
            with self.subTest(kind=kind):
                value = _manifest()
                if kind == "inputs":
                    value["inputs"][0]["relative_path"] = path
                else:
                    value["protocol"]["preprocessing"]["source"] = path
                with self.assertRaisesRegex(self.api.R10RecordError, "absolute|relative"):
                    self.api.validate_manifest(value)
        value = _manifest()
        value["protocol"]["algorithm"]["hostname"] = "private-machine"
        with self.assertRaisesRegex(self.api.R10RecordError, "host-specific"):
            self.api.validate_manifest(value)

    def test_relative_paths_reject_traversal_backslashes_and_noncanonical_forms(self) -> None:
        for path in ("../model.obj", "inputs/../model.obj", "inputs\\model.obj", "./model.obj"):
            with self.subTest(path=path):
                value = _manifest()
                value["inputs"][0]["relative_path"] = path
                with self.assertRaises(self.api.R10RecordError):
                    self.api.validate_manifest(value)

    def test_sha256_and_commit_formats_fail_closed(self) -> None:
        for field, bad in (("sha256", "A" * 64), ("sha256", "1" * 63)):
            value = _manifest()
            value["inputs"][0][field] = bad
            with self.assertRaisesRegex(self.api.R10RecordError, "SHA-256"):
                self.api.validate_manifest(value)
        value = _manifest()
        value["provenance"]["repository_commit_sha"] = "a" * 39
        with self.assertRaisesRegex(self.api.R10RecordError, "Git object"):
            self.api.validate_manifest(value)

    def test_canonical_bytes_are_stable_and_round_trip_only_canonical_json(self) -> None:
        first = _manifest()
        second = {key: first[key] for key in reversed(list(first))}
        self.assertEqual(
            self.api.canonical_json_bytes(first), self.api.canonical_json_bytes(second)
        )
        self.assertEqual(self.api.document_sha256(first), self.api.document_sha256(second))
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "manifest.json"
            self.api.write_canonical_json(path, first)
            self.assertFalse(path.read_bytes().endswith(b"\n"))
            self.assertEqual(self.api.read_canonical_json(path), first)
            path.write_text(json.dumps(first, indent=2), encoding="utf-8")
            with self.assertRaisesRegex(Exception, "not canonical"):
                self.api.read_canonical_json(path)

    def test_nonfinite_duplicate_and_non_json_values_are_rejected(self) -> None:
        value = _manifest()
        value["protocol"]["algorithm"]["bad"] = float("nan")
        with self.assertRaisesRegex(self.api.R10RecordError, "non-finite"):
            self.api.validate_manifest(value)
        with self.assertRaises(Exception):
            self.api.loads_json_strict('{"schema":1,"schema":1}')
        with self.assertRaisesRegex(self.api.R10RecordError, "unsupported JSON"):
            self.api.canonical_json_bytes({"set": {1, 2}})

    def test_bind_artifact_hashes_bytes_and_never_serializes_root(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifact = root / "inputs" / "model.obj"
            artifact.parent.mkdir()
            artifact.write_bytes(b"v 0 0 0\n")
            binding = self.api.bind_artifact(
                root,
                "inputs/model.obj",
                artifact_id="model",
                role="model_input",
            )
            self.assertEqual(binding["sha256"], hashlib.sha256(b"v 0 0 0\n").hexdigest())
            self.assertEqual(binding["byte_count"], 8)
            self.assertNotIn(str(root), self.api.canonical_json_bytes(binding).decode())

    @unittest.skipUnless(hasattr(os, "symlink"), "symbolic links unavailable")
    def test_bind_artifact_rejects_symlink_component(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            real = root / "real"
            real.mkdir()
            (real / "model.obj").write_text("mesh", encoding="utf-8")
            os.symlink(real, root / "inputs", target_is_directory=True)
            with self.assertRaisesRegex(self.api.R10RecordError, "symbolic-link"):
                self.api.bind_artifact(
                    root,
                    "inputs/model.obj",
                    artifact_id="model",
                    role="model_input",
                )

    def test_run_record_binds_manifest_claim_paper_values_and_local_definitions(self) -> None:
        manifest = _manifest()
        record = _run_record(self.api, manifest, 0, 0.1)
        self.api.validate_run_record(record, manifest=manifest)
        cases = [
            ("manifest_sha256", "0" * 64, "bind"),
            ("experiment_id", "different", "match"),
            ("claim_level", "smoke", "match"),
        ]
        for field, replacement, message in cases:
            with self.subTest(field=field):
                bad = copy.deepcopy(record)
                bad[field] = replacement
                with self.assertRaisesRegex(self.api.R10RecordError, message):
                    self.api.validate_run_record(bad, manifest=manifest)
        bad = copy.deepcopy(record)
        bad["paper_reference_metrics"][0]["value"] = 127
        with self.assertRaisesRegex(self.api.R10RecordError, "copy"):
            self.api.validate_run_record(bad, manifest=manifest)
        bad = copy.deepcopy(record)
        bad["local_metrics"][0]["unit"] = "millimeters"
        with self.assertRaisesRegex(self.api.R10RecordError, "unit"):
            self.api.validate_run_record(bad, manifest=manifest)
        bad = copy.deepcopy(record)
        bad["local_metrics"][0]["semantics"] = "A different measurement definition."
        with self.assertRaisesRegex(self.api.R10RecordError, "semantics"):
            self.api.validate_run_record(bad, manifest=manifest)
        bad = copy.deepcopy(record)
        bad["execution"]["build"]["commit_sha"] = "b" * 40
        with self.assertRaisesRegex(self.api.R10RecordError, "commit_sha"):
            self.api.validate_run_record(bad, manifest=manifest)
        bad = copy.deepcopy(record)
        bad["execution"]["build"]["dirty_worktree"] = True
        with self.assertRaisesRegex(self.api.R10RecordError, "dirty_worktree"):
            self.api.validate_run_record(bad, manifest=manifest)
        bad = copy.deepcopy(record)
        bad["execution"]["timing"]["warmup_runs"] = 2
        with self.assertRaisesRegex(self.api.R10RecordError, "warmup_runs"):
            self.api.validate_run_record(bad, manifest=manifest)
        bad = copy.deepcopy(record)
        bad["artifacts"][0]["sha256"] = "8" * 64
        with self.assertRaisesRegex(self.api.R10RecordError, "native extension"):
            self.api.validate_run_record(bad, manifest=manifest)
        bad = copy.deepcopy(record)
        bad["artifacts"].pop(0)
        with self.assertRaisesRegex(self.api.R10RecordError, "native extension"):
            self.api.validate_run_record(bad, manifest=manifest)

    def test_run_recorder_collects_explicit_build_hardware_timing_and_metrics(self) -> None:
        manifest = _manifest(repeats=1)
        recorder = self.api.RunRecorder(
            manifest,
            run_id="synthetic-r00",
            repeat_index=0,
            seed=31000,
            build=_build(),
            hardware=_hardware(),
            timing_scope="model_setup_start-through_run_record_finalization",
            warmup_runs=1,
        )
        with recorder.phase("optimization"):
            sum(range(20))
        recorder.add_local_metric(
            "projected_com_error",
            0.0,
            semantics="Maximum final normalized signed support violation.",
        )
        record = recorder.finish(status="succeeded")
        self.assertEqual(
            record["artifacts"][0],
            manifest["provenance"]["resources"][1],
        )
        self.api.validate_run_record(record, manifest=manifest)
        self.assertEqual(record["execution"]["hardware"], _hardware())
        self.assertEqual(record["execution"]["build"], _build())
        self.assertGreaterEqual(record["execution"]["timing"]["phases"][0]["seconds"], 0)
        with self.assertRaisesRegex(self.api.R10RecordError, "finalized"):
            with recorder.phase("optimization"):
                pass
        with self.assertRaisesRegex(self.api.R10RecordError, "finalized"):
            recorder.finish(status="succeeded")

    def test_timing_phases_cannot_exceed_the_recorded_total(self) -> None:
        manifest = _manifest()
        record = _run_record(self.api, manifest, 0, 0.1)
        record["execution"]["timing"]["phases"] = [
            {"name": "first", "seconds": 1.0},
            {"name": "second", "seconds": 1.0},
        ]
        record["execution"]["timing"]["total_seconds"] = 1.5
        with self.assertRaisesRegex(self.api.R10RecordError, "exceed"):
            self.api.validate_run_record(record, manifest=manifest)
        record = _run_record(self.api, manifest, 0, 0.1)
        record["execution"]["timing"]["scope"] = "optimization-only"
        with self.assertRaisesRegex(self.api.R10RecordError, "scope"):
            self.api.validate_run_record(record, manifest=manifest)

    def test_aggregate_uses_type7_median_iqr_min_max_and_does_not_pool_paper(self) -> None:
        manifest = _manifest(repeats=4)
        records = [
            _run_record(self.api, manifest, index, value)
            for index, value in enumerate((0.0, 1.0, 2.0, 9.0))
        ]
        aggregate = self.api.aggregate_run_records(manifest, records)
        self.api.validate_aggregate(aggregate, manifest=manifest)
        summary = next(
            item
            for item in aggregate["local_metric_aggregates"]
            if item["name"] == "projected_com_error"
        )
        self.assertEqual(summary["sample_count"], 4)
        self.assertEqual(summary["median"], 1.5)
        self.assertEqual(summary["q1"], 0.75)
        self.assertEqual(summary["q3"], 3.75)
        self.assertEqual(summary["iqr"], 3.0)
        self.assertEqual(summary["minimum"], 0.0)
        self.assertEqual(summary["maximum"], 9.0)
        self.assertFalse(aggregate["exact_reproduction_claim"])
        self.assertEqual(
            aggregate["paper_reference_metrics"],
            manifest["paper_reference"]["reported_metrics"],
        )
        self.assertFalse(
            any(item["name"] == "reported_iteration_count" for item in aggregate["local_metric_aggregates"])
        )

    def test_aggregate_requires_complete_unique_manifest_repeat_set(self) -> None:
        manifest = _manifest(repeats=2)
        first = _run_record(self.api, manifest, 0, 0.0)
        second = _run_record(self.api, manifest, 1, 1.0)
        with self.assertRaisesRegex(self.api.R10RecordError, "exactly 2"):
            self.api.aggregate_run_records(manifest, [first])
        duplicate = copy.deepcopy(second)
        duplicate["run_id"] = first["run_id"]
        with self.assertRaisesRegex(self.api.R10RecordError, "run_id"):
            self.api.aggregate_run_records(manifest, [first, duplicate])

    def test_schema_files_parse_and_cli_writes_canonical_bound_aggregate(self) -> None:
        repository = Path(__file__).resolve().parents[2]
        for name in (
            "r10-experiment-manifest-v1.schema.json",
            "r10-run-record-v1.schema.json",
            "r10-aggregate-v1.schema.json",
        ):
            with self.subTest(name=name):
                schema = json.loads((repository / "schemas" / name).read_text(encoding="utf-8"))
                self.assertEqual(schema["$schema"], "https://json-schema.org/draft/2020-12/schema")
        manifest = _manifest(repeats=2)
        records = [_run_record(self.api, manifest, index, float(index)) for index in range(2)]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest_path = root / "manifest.json"
            output_path = root / "aggregate.json"
            self.api.write_canonical_json(manifest_path, manifest)
            record_paths = []
            for index, record in enumerate(records):
                path = root / f"r{index}.json"
                self.api.write_canonical_json(path, record)
                record_paths.append(path)
            command = [
                sys.executable,
                str(repository / "python" / "r10_aggregate.py"),
                "--manifest",
                str(manifest_path),
            ]
            for path in record_paths:
                command.extend(["--run-record", str(path)])
            command.extend(["--output", str(output_path)])
            environment = dict(os.environ)
            environment["PYTHONPATH"] = str(repository / "python")
            completed = subprocess.run(
                command,
                cwd=root,
                env=environment,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            aggregate = self.api.read_canonical_json(output_path)
            self.api.validate_aggregate(aggregate, manifest=manifest)


if __name__ == "__main__":
    unittest.main(verbosity=2)
