from __future__ import annotations

import hashlib
import json
import pathlib
import unittest
from urllib.parse import urlparse


ROOT = pathlib.Path(__file__).resolve().parents[2]
INVENTORY = ROOT / "experiments" / "spin_it_reproduction" / "paper_inventory.json"


class SpinItPaperInventoryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.payload = json.loads(INVENTORY.read_text(encoding="utf-8"))

    def test_schema_and_claim_dimensions_are_explicit(self) -> None:
        self.assertEqual(self.payload["schema"], "spin-it-public-paper-inventory-v1")
        self.assertEqual(self.payload["schema_version"], 1)
        self.assertEqual(
            self.payload["claim_dimensions"],
            [
                "formula_coverage",
                "experiment_coverage",
                "asset_class",
                "author_exact",
                "physical_validation",
            ],
        )
        self.assertFalse(self.payload["aggregate_claim"]["author_exact"])
        self.assertFalse(self.payload["aggregate_claim"]["physical_validation"])

    def test_both_papers_and_public_attachments_are_identified(self) -> None:
        papers = {paper["id"]: paper for paper in self.payload["papers"]}
        self.assertEqual(set(papers), {"spin-it-2014", "spin-it-faster-2024"})
        self.assertEqual(papers["spin-it-2014"]["doi"], "10.1145/2601097.2601157")
        self.assertEqual(papers["spin-it-faster-2024"]["doi"], "10.1145/3658194")
        faster_kinds = {source["kind"] for source in papers["spin-it-faster-2024"]["sources"]}
        self.assertEqual(faster_kinds, {"PAPER", "SUPPLEMENT", "VIDEO"})
        for paper in papers.values():
            for source in paper["sources"]:
                parsed = urlparse(source["url"])
                self.assertEqual(parsed.scheme, "https")
                self.assertTrue(parsed.netloc)
                if "sha256" in source:
                    self.assertEqual(len(source["sha256"]), 64)
                    int(source["sha256"], 16)

    def test_published_experiment_rows_are_complete(self) -> None:
        rows = self.payload["experiments"]
        ids = {row["id"] for row in rows}
        required = {
            "sif-static-tilted-vase",
            "sif-floating-tea-infuser",
            "sif-floating-island-cat",
            "sif-floating-island-dog-human",
            "sif-heart-top",
            "sif-noodles-top",
            "sif-jellyfish-yoyo",
            "sif-flower-fidget",
            "sif-quadratic-convergence",
            "sif-static-pareto",
            "sif-top-pareto",
            "sif-near-infeasible",
            "sif-degenerate-problems",
            "si14-ellipsoid",
            "si14-heart",
            "si14-armadillo-hand",
            "si14-armadillo-shell",
            "si14-teapot-variants",
            "si14-cuboid-yoyo",
            "si14-woven-ring-yoyo",
            "si14-dancing-couple",
            "si14-sumone-dual-material",
            "si14-trex-static",
            "si14-elephant",
            "si14-cuts-and-voids",
        }
        self.assertTrue(required.issubset(ids), sorted(required - ids))
        for row in rows:
            self.assertIn(row["paper"], {"spin-it-2014", "spin-it-faster-2024"})
            self.assertIn(row["asset_class"], {"AUTHOR_ASSET", "RECONSTRUCTED_ASSET", "BLOCKED_MISSING_ASSET"})
            self.assertIsInstance(row["formula_coverage"], bool)
            self.assertIsInstance(row["experiment_coverage"], bool)
            self.assertIsInstance(row["author_exact"], bool)
            self.assertIsInstance(row["physical_validation"], bool)
            self.assertFalse(row["author_exact"])

    def test_inventory_hash_is_content_addressed(self) -> None:
        declared = self.payload["inventory_sha256"]
        canonical = dict(self.payload)
        canonical["inventory_sha256"] = ""
        encoded = json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode()
        self.assertEqual(declared, hashlib.sha256(encoded).hexdigest())


if __name__ == "__main__":
    unittest.main()
