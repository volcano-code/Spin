from __future__ import annotations

import importlib.util
import inspect
import json
import math
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "experiments/spin_it_reproduction/spin2014/run_paper_cases.py"


def load_runner():
    spec = importlib.util.spec_from_file_location("spin_it_2014_paper_runner", RUNNER)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load Spin-It 2014 paper runner")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SpinIt2014PaperRunnerTests(unittest.TestCase):
    def test_rle_decode_and_component_count_are_exact(self) -> None:
        runner = load_runner()
        self.assertEqual(runner.decode_rle([[0, 2], [1, 3], [0, 1]], 6), (0, 0, 1, 1, 1, 0))
        with self.assertRaises(ValueError):
            runner.decode_rle([[1, 2]], 3)
        mask = [0] * 27
        mask[0] = mask[1] = 1
        mask[-1] = 1
        self.assertEqual(runner.count_components(mask, (3, 3, 3)), 2)

    def test_all_objective_metrics_follow_independent_equations(self) -> None:
        runner = load_runner()
        self.assertIn("objective", inspect.signature(runner.paper_metrics).parameters)
        moments = (4.0, 0.0, 0.0, 2.0, 10.0, 10.0, 1.0, 0.0, 0.0, 0.0)
        # Central inertia is diag(10, 10, 20), so f_yo at gamma_I=2 is 1.
        inertia = runner.paper_metrics(
            "SPINNING_TOP", moments, 3.0,
            objective="TOP_INERTIA_ONLY", inertia_weight=2.0,
        )
        lowering = runner.paper_metrics(
            "SPINNING_TOP", moments, 3.0,
            objective="TOP_LOWERING_ONLY", inertia_weight=2.0,
        )
        full = runner.paper_metrics(
            "SPINNING_TOP", moments, 3.0,
            objective="TOP_MASS_AND_LOWERING", inertia_weight=2.0,
        )
        yoyo = runner.paper_metrics(
            "YO_YO", moments, 99.0,
            objective="YOYO", inertia_weight=2.0,
        )
        static = runner.paper_metrics(
            "STATIC_BALANCE", moments, 3.0,
            objective="STATIC_LOWERING", inertia_weight=99.0,
        )
        self.assertAlmostEqual(full["center_height"], 0.5)
        self.assertAlmostEqual(inertia["objective"], 1.0)
        self.assertAlmostEqual(lowering["objective"], 1.75)
        self.assertAlmostEqual(full["objective"], 13.0)
        self.assertAlmostEqual(yoyo["objective"], 1.0)
        self.assertAlmostEqual(static["objective"], 0.75)
        self.assertAlmostEqual(full["axis_alignment_error"], 0.0)
        self.assertTrue(full["axial_is_maximum_principal_moment"])
        self.assertTrue(math.isfinite(full["principal_moment_ratio"]))

    def test_legacy_teapot_cases_resolve_to_three_distinct_objectives(self) -> None:
        runner = load_runner()
        self.assertTrue(hasattr(runner, "objective_for_case"))
        manifest = json.loads(runner.DEFAULT_CASES.read_text(encoding="utf-8"))
        teapots = {
            case["id"]: runner.objective_for_case(case)
            for case in manifest["cases"]
            if case["id"].startswith("teapot-")
        }
        self.assertEqual(teapots, {
            "teapot-inertia": "TOP_INERTIA_ONLY",
            "teapot-lowering": "TOP_LOWERING_ONLY",
            "teapot-mass-reduction": "TOP_MASS_AND_LOWERING",
        })

    def test_scientific_hash_ignores_only_timings(self) -> None:
        runner = load_runner()
        first = runner.scientific_hash({"output_hash": "a", "elapsed_seconds": 1.0})
        second = runner.scientific_hash({"output_hash": "a", "elapsed_seconds": 9.0})
        third = runner.scientific_hash({"output_hash": "b", "elapsed_seconds": 1.0})
        self.assertEqual(first, second)
        self.assertNotEqual(first, third)


class SpinIt2014PaperRunnerNativeTests(unittest.TestCase):
    def test_ellipsoid_two_repeat_miniature_is_sealed_and_deterministic(self) -> None:
        runner = load_runner()
        manifest = json.loads(runner.DEFAULT_CASES.read_text(encoding="utf-8"))
        case = next(row for row in manifest["cases"] if row["id"] == "ellipsoid-top")
        result = runner.run_case(case, repeats=2, maximum_level=4)
        self.assertEqual(len(result["repeats"]), 2)
        self.assertTrue(result["deterministic"])
        self.assertTrue(result["passed"])
        self.assertTrue(all(row["geometry_sealed"] for row in result["repeats"]))
        self.assertTrue(all(row["metrics"]["mass"] > 0.0 for row in result["repeats"]))


if __name__ == "__main__":
    unittest.main()
