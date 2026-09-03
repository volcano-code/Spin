from __future__ import annotations

import importlib.util
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
    def test_public_gate_requires_both_computational_matrices(self) -> None:
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
        self.assertTrue(gate["public_complete_reproduction"])
        self.assertFalse(gate["author_exact"])
        self.assertFalse(gate["physical_validation"])
        self.assertTrue(gate["formula_coverage"])
        self.assertTrue(gate["experiment_coverage"])

        faster["cases"][0]["passed"] = False
        self.assertFalse(
            builder.aggregate(faster, spin2014, figures, tests_passed=True)[
                "public_complete_reproduction"
            ]
        )

        faster["cases"][0]["passed"] = True
        spin2014["paper_resolution_reproduced"] = False
        blocked = builder.aggregate(faster, spin2014, figures, tests_passed=True)
        self.assertFalse(blocked["public_complete_reproduction"])
        self.assertFalse(blocked["requirements"]["spin_it_2014_paper_resolution"])

        spin2014["paper_resolution_reproduced"] = True
        spin2014["paper_numerics_reproduced"] = True
        spin2014["cage_deformation_executed"] = False
        blocked = builder.aggregate(faster, spin2014, figures, tests_passed=True)
        self.assertFalse(blocked["public_complete_reproduction"])
        self.assertFalse(blocked["requirements"]["spin_it_2014_cage_deformation"])

    def test_missing_claim_dimensions_fail_closed(self) -> None:
        builder = load_builder()
        with self.assertRaises(ValueError):
            builder.validate_gate({"public_complete_reproduction": True})


if __name__ == "__main__":
    unittest.main()
