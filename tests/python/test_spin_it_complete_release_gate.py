from __future__ import annotations

import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
BUILDER = ROOT / "experiments/spin_it_reproduction/build_release_gate.py"


def load_builder():
    spec = importlib.util.spec_from_file_location("spin_it_release_gate", BUILDER)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load Spin-It release gate builder")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SpinItCompleteReleaseGateTests(unittest.TestCase):
    def test_synthetic_all_true_payload_cannot_publish(self) -> None:
        builder = load_builder()
        faster = {"cases": [{"passed": True}], "author_exact": False, "physical_validation": False}
        spin2014 = {
            "cases": [{"passed": True}],
            "author_exact": False,
            "physical_validation": False,
            "paper_resolution_reproduced": True,
            "paper_numerics_reproduced": True,
            "cage_deformation_executed": True,
        }
        figures = {
            "figure11": {"result": {"front_audit": {"finite": True, "nondominated": True}}},
            "figure12": {"result": {"front_audit": {"finite": True, "nondominated": True}}},
            "figure13": {"result": {"rows": [{"status": "CONVERGED"}]}},
        }
        gate = builder.aggregate(faster, spin2014, figures, tests_passed=True)
        self.assertFalse(gate["public_complete_reproduction"])
        self.assertFalse(gate["author_exact"])
        self.assertFalse(gate["physical_validation"])
        self.assertFalse(gate["formula_coverage"])
        self.assertFalse(gate["experiment_coverage"])
        self.assertEqual(gate["status"], "SUPERSEDED_BLOCKED")
        self.assertIn("strict requirements audit", " ".join(gate["limitations"]))
        report = builder.report(gate)
        self.assertIn("superseded legacy compatibility report", report)
        self.assertIn("cannot grant strict reproduction claims", report)

    def test_legacy_cli_rejects_removed_tests_passed_switch(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            result = subprocess.run(
                [
                    sys.executable,
                    str(BUILDER),
                    "--tests-passed",
                    "--faster", str(root / "faster.json"),
                    "--spin2014", str(root / "spin2014.json"),
                    "--figures", str(root / "figures.json"),
                    "--output", str(root / "gate.json"),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("unrecognized arguments: --tests-passed", result.stderr)

    def test_missing_claim_dimensions_fail_closed(self) -> None:
        builder = load_builder()
        with self.assertRaises(ValueError):
            builder.validate_gate({"public_complete_reproduction": True})


if __name__ == "__main__":
    unittest.main()
