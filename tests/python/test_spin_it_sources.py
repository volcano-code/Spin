from __future__ import annotations

import json
from pathlib import Path
import subprocess
import unittest


class SpinItSourceContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repository = Path(__file__).resolve().parents[2]
        cls.manifest_path = (
            cls.repository
            / "experiments/spin_it_reproduction/sources/manifest.json"
        )
        cls.equations_path = cls.repository / "docs/spin_it_faster_equations.md"

    def test_primary_source_manifest_is_exact_and_license_safe(self) -> None:
        manifest = json.loads(self.manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(manifest["schema"], "spin-it-faster-sources-v1")
        self.assertEqual(manifest["doi"], "10.1145/3658194")
        self.assertEqual(manifest["author_code_status"], "NOT_PUBLICLY_RELEASED")
        self.assertEqual(manifest["maximum_claim"], "PAPER_FORMULA_REPRODUCTION")
        sources = manifest["sources"]
        self.assertEqual([item["id"] for item in sources], ["main", "supplement"])
        self.assertEqual(
            [item["sha256"] for item in sources],
            [
                "e7fc6395362ef7e63a703851176afd81d8e953bae4e8080279d7e244b3c51fcc",
                "3e18e0e784bb917897f01d1ce182e55dbc3bee5e055aca0f71582673f4431200",
            ],
        )
        for source in sources:
            self.assertTrue(source["url"].startswith("https://research-explorer.ista.ac.at/"))
            self.assertRegex(source["sha256"], r"^[0-9a-f]{64}$")
        completed = subprocess.run(
            ["git", "-C", str(self.repository), "ls-files", "*.pdf"],
            check=True, capture_output=True, text=True,
        )
        committed_pdf_names = {Path(path).name for path in completed.stdout.splitlines()}
        self.assertNotIn("sif-final.pdf", committed_pdf_names)
        self.assertNotIn("sif-supp-final.pdf", committed_pdf_names)
        ignored = subprocess.run(
            ["git", "-C", str(self.repository), "check-ignore", "-q",
             "build/strict-paper-sources/sif-final.pdf"],
            check=False,
        )
        self.assertEqual(ignored.returncode, 0)

    def test_equation_contract_records_paper_order_and_solver_constants(self) -> None:
        text = self.equations_path.read_text(encoding="utf-8")
        for required in (
            "(v, l1, l2, l3, Q11, Q22, Q33, Q23, Q13, Q12)",
            "(a, b1, b2, b3, A11, A22, A33, A23, A13, A12)",
            "beta(x) > 0",
            "alpha = 1/2",
            "t_max = 1/2",
            "1e-12",
            "60 function evaluations",
        ):
            self.assertIn(required, text)


if __name__ == "__main__":
    unittest.main()
