from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
GENERATOR = ROOT / "experiments/spin_it_reproduction/spin2014/generate_paper_assets.py"
EXPECTED = {
    "ellipsoid-top", "heart-top", "armadillo-hand", "armadillo-shell",
    "teapot-inertia", "teapot-lowering", "teapot-mass-reduction",
    "cuboid-yoyo", "woven-ring-yoyo", "dancing-couple-hollowing",
    "dancing-couple-deformation", "sumone-dual-density", "trex-static",
    "elephant-accessibility",
}


def load_generator():
    spec = importlib.util.spec_from_file_location("spin2014_assets", GENERATOR)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load Spin-It 2014 asset generator")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SpinIt2014PaperAssetTests(unittest.TestCase):
    def test_generator_is_deterministic_and_inventory_is_complete(self) -> None:
        generator = load_generator()
        with tempfile.TemporaryDirectory() as first_dir, tempfile.TemporaryDirectory() as second_dir:
            subprocess.run(["python3", str(GENERATOR), "--output", first_dir], check=True)
            subprocess.run(["python3", str(GENERATOR), "--output", second_dir], check=True)
            first_root, second_root = Path(first_dir), Path(second_dir)
            first = json.loads((first_root / "paper_cases.json").read_text(encoding="utf-8"))
            second = json.loads((second_root / "paper_cases.json").read_text(encoding="utf-8"))
            self.assertEqual(first, second)
            self.assertEqual(first["asset_class"], "RECONSTRUCTED_ASSET")
            self.assertFalse(first["author_exact"])
            self.assertFalse(first["physical_validation"])
            self.assertEqual({case["id"] for case in first["cases"]}, EXPECTED)
            for case in first["cases"]:
                self.assertIn(case["mode"], {"STATIC_BALANCE", "YO_YO", "SPINNING_TOP"})
                self.assertEqual(len(case["target_axis"]), 3)
                self.assertGreaterEqual(case["expected_void_components"], 1)
                self.assertGreaterEqual(case["manual_planar_cuts"], 1)
                self.assertEqual(case["maximum_octree_level"], 4)
                self.assertEqual(case["paper_maximum_octree_level"], 9)
                self.assertGreater(case["binary_constraint_tolerance"], case["paper_binary_constraint_tolerance"])
                self.assertIn("accessibility", case)
                asset = first_root / case["occupancy"]["path"]
                mirror = second_root / case["occupancy"]["path"]
                self.assertEqual(asset.read_bytes(), mirror.read_bytes())
                self.assertEqual(generator.sha256(asset), case["occupancy"]["sha256"])
                payload = json.loads(asset.read_text(encoding="utf-8"))
                if case["mode"] == "YO_YO":
                    self.assertEqual(payload["minimum"][2], -1.0)
                    self.assertEqual(payload["maximum"][2], 1.0)
                else:
                    self.assertEqual(payload["minimum"][2], 0.0)
                    self.assertEqual(payload["maximum"][2], 2.0)
                material = generator.decode_rle(payload["material_rle"])
                carveable = generator.decode_rle(payload["carveable_rle"])
                self.assertEqual(len(material), 16 ** 3)
                self.assertTrue(any(material))
                self.assertTrue(any(carveable))
                self.assertTrue(all(not value or material[i] for i, value in enumerate(carveable)))


if __name__ == "__main__":
    unittest.main()
