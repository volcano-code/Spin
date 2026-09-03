from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
GENERATOR = (
    ROOT
    / "experiments"
    / "spin_it_reproduction"
    / "faster"
    / "generate_paper_assets.py"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_faces(path: Path) -> list[tuple[int, int, int]]:
    faces: list[tuple[int, int, int]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if fields and fields[0] == "f":
            faces.append(tuple(int(value) - 1 for value in fields[1:4]))
    return faces


def assert_closed_oriented(test: unittest.TestCase, path: Path) -> None:
    directed: dict[tuple[int, int], int] = {}
    for a, b, c in read_faces(path):
        for edge in ((a, b), (b, c), (c, a)):
            directed[edge] = directed.get(edge, 0) + 1
    test.assertTrue(directed, path)
    for (a, b), count in directed.items():
        test.assertEqual(count, 1, (path, a, b, count))
        test.assertEqual(directed.get((b, a)), 1, (path, b, a))


class SpinItFasterPaperAssetTests(unittest.TestCase):
    def test_generator_is_deterministic_and_inventory_is_complete(self) -> None:
        with tempfile.TemporaryDirectory() as first_dir, tempfile.TemporaryDirectory() as second_dir:
            first = Path(first_dir)
            second = Path(second_dir)
            subprocess.run(
                ["python3", str(GENERATOR), "--output", str(first)],
                check=True,
            )
            subprocess.run(
                ["python3", str(GENERATOR), "--output", str(second)],
                check=True,
            )
            manifest_path = first / "paper_cases.json"
            self.assertEqual(manifest_path.read_bytes(), (second / "paper_cases.json").read_bytes())
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["schema"], "spin-it-faster-paper-cases-v1")
            self.assertEqual(payload["asset_class"], "RECONSTRUCTED_ASSET")
            cases = payload["cases"]
            self.assertEqual(
                {case["id"] for case in cases},
                {
                    "static-tilted-vase",
                    "floating-tea-infuser",
                    "floating-island-cat",
                    "floating-island-dog-human",
                    "heart-spinning-top",
                    "noodles-spinning-top",
                    "jellyfish-yoyo",
                    "flower-fidget-spinner",
                },
            )
            self.assertEqual(len(cases), 8)
            self.assertEqual(sum(case["mode"] == "BUOYANT_STABILITY" for case in cases), 3)
            for case in cases:
                self.assertFalse(case["author_exact"])
                self.assertFalse(case["physical_validation"])
                self.assertEqual(case["units"], "normalized_unit")
                self.assertEqual(len(case["target_axes"]), 3)
                for key in ("outer", "design_domain"):
                    path = first / case[key]["path"]
                    counterpart = second / case[key]["path"]
                    self.assertEqual(sha256(path), case[key]["sha256"])
                    self.assertEqual(path.read_bytes(), counterpart.read_bytes())
                    assert_closed_oriented(self, path)
                if case["mode"] == "BUOYANT_STABILITY":
                    self.assertGreater(case["buoyancy"]["liquid_density"], 0.0)
                    self.assertGreater(case["buoyancy"]["displaced_liquid_volume"], 0.0)


if __name__ == "__main__":
    unittest.main()
