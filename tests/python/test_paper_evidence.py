"""Evidence classification and executable/source freshness contracts."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from make_it_stand.paper_evidence import (
    EvidenceClass,
    FreshnessStatus,
    audit_manifest_freshness,
    assess_result_freshness,
    capture_execution_fingerprint,
    classify_evidence,
)


class PaperEvidenceTests(unittest.TestCase):
    def test_claim_classification_is_fail_closed(self) -> None:
        self.assertEqual(
            classify_evidence(complete_cache=True),
            EvidenceClass.EXACT_CACHE,
        )
        self.assertEqual(
            classify_evidence(complete_cache=True, parameter_override=True),
            EvidenceClass.RECONSTRUCTED_PARAMETER,
        )
        self.assertEqual(
            classify_evidence(complete_cache=True, reconstructed_surface=True),
            EvidenceClass.RECONSTRUCTED_ASSET,
        )
        self.assertEqual(
            classify_evidence(complete_cache=False),
            EvidenceClass.BLOCKED_MISSING_ASSET,
        )
        with self.assertRaisesRegex(ValueError, "missing assets"):
            classify_evidence(complete_cache=False, parameter_override=True)

    def test_source_or_executable_change_marks_result_stale(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "solver.py"
            executable = root / "solver"
            source.write_text("policy = 1\n", encoding="utf-8")
            executable.write_bytes(b"executable-v1")
            first = capture_execution_fingerprint(
                root,
                source_paths=[source],
                executable_paths=[executable],
                configuration={"profile": "paper", "lambda": 20.0},
            )
            self.assertEqual(
                assess_result_freshness(first, first).status,
                FreshnessStatus.CURRENT,
            )

            source.write_text("policy = 2\n", encoding="utf-8")
            changed_source = capture_execution_fingerprint(
                root,
                source_paths=[source],
                executable_paths=[executable],
                configuration={"profile": "paper", "lambda": 20.0},
            )
            source_assessment = assess_result_freshness(first, changed_source)
            self.assertEqual(source_assessment.status, FreshnessStatus.STALE_RESULT)
            self.assertIn("source_tree_sha256", source_assessment.mismatches)

            source.write_text("policy = 1\n", encoding="utf-8")
            executable.write_bytes(b"executable-v2")
            changed_executable = capture_execution_fingerprint(
                root,
                source_paths=[source],
                executable_paths=[executable],
                configuration={"profile": "paper", "lambda": 20.0},
            )
            executable_assessment = assess_result_freshness(
                first,
                changed_executable,
            )
            self.assertEqual(
                executable_assessment.status,
                FreshnessStatus.STALE_RESULT,
            )
            self.assertIn("executables", executable_assessment.mismatches)

    def test_fingerprint_is_json_serializable_and_configuration_bound(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "runner.py"
            source.write_text("pass\n", encoding="utf-8")
            fingerprint = capture_execution_fingerprint(
                root,
                source_paths=[source],
                configuration={"profile": "released-beta", "repeat": 2},
            )
            serialized = json.loads(json.dumps(fingerprint, sort_keys=True))
            self.assertEqual(serialized["schema_version"], 1)
            self.assertEqual(
                serialized["configuration"],
                {"profile": "released-beta", "repeat": 2},
            )
            changed = capture_execution_fingerprint(
                root,
                source_paths=[source],
                configuration={"profile": "released-beta", "repeat": 1},
            )
            assessment = assess_result_freshness(fingerprint, changed)
            self.assertEqual(assessment.status, FreshnessStatus.STALE_RESULT)
            self.assertIn("configuration_sha256", assessment.mismatches)

    def test_runtime_caches_do_not_change_source_fingerprint(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            package = root / "package"
            package.mkdir()
            (package / "module.py").write_text("VALUE = 1\n", encoding="utf-8")
            first = capture_execution_fingerprint(root, source_paths=[package])
            cache = package / "__pycache__"
            cache.mkdir()
            (cache / "module.cpython-313.pyc").write_bytes(b"runtime-cache")
            second = capture_execution_fingerprint(root, source_paths=[package])
            self.assertEqual(
                first["source_tree_sha256"],
                second["source_tree_sha256"],
            )

    def test_manifest_audit_marks_legacy_or_changed_result_stale(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "solver.py"
            source.write_text("VERSION = 1\n", encoding="utf-8")
            fingerprint = capture_execution_fingerprint(
                root,
                source_paths=[source],
                configuration={"case": "armadillo1"},
            )
            manifest = {"execution_fingerprint": fingerprint}
            current = audit_manifest_freshness(manifest, root)
            self.assertEqual(current["freshness"]["status"], "CURRENT")

            source.write_text("VERSION = 2\n", encoding="utf-8")
            stale = audit_manifest_freshness(manifest, root)
            self.assertEqual(stale["freshness"]["status"], "STALE_RESULT")
            self.assertIn(
                "source_tree_sha256",
                stale["freshness"]["mismatches"],
            )

            legacy = audit_manifest_freshness({"case": "legacy"}, root)
            self.assertEqual(legacy["freshness"]["status"], "STALE_RESULT")
            self.assertIn(
                "missing_execution_fingerprint",
                legacy["freshness"]["mismatches"],
            )


if __name__ == "__main__":
    unittest.main()
