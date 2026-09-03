from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "experiments/spin_it_reproduction/faster/run_paper_figures.py"


def load_runner():
    spec = importlib.util.spec_from_file_location("spin_it_faster_figures", RUNNER)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load Faster paper-figure runner")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SpinItFasterPaperFigureTests(unittest.TestCase):
    def test_figure10_pairs_successful_monotone_residuals_by_application(self) -> None:
        runner = load_runner()
        table = {
            "cases": [
                {
                    "id": "static",
                    "mode": "STATIC_BALANCE",
                    "attempts": [{
                        "status": "CONVERGED",
                        "trace": [
                            {"normalized_equation_residual": 1e-1},
                            {"normalized_equation_residual": 1e-2},
                            {"normalized_equation_residual": 1e-4},
                        ],
                    }],
                },
                {
                    "id": "failed",
                    "mode": "SPINNING_TOP",
                    "attempts": [{"status": "RESOURCE_LIMIT_REACHED", "trace": []}],
                },
            ]
        }
        grouped = runner.quadratic_convergence_pairs(table)
        self.assertEqual(grouped["static_stability"], [(1e-1, 1e-2), (1e-2, 1e-4)])
        self.assertEqual(grouped["spinning_tops"], [])

    def test_pareto_weight_grid_and_front_invariants(self) -> None:
        runner = load_runner()
        weights = runner.pareto_weights(11)
        self.assertEqual(weights[0], 0.0)
        self.assertEqual(weights[-1], 1.0)
        self.assertTrue(all(a < b for a, b in zip(weights, weights[1:])))
        front = runner.validate_pareto_front([
            {"weight": 0.0, "first_objective": 4.0, "center_height": 0.2},
            {"weight": 0.5, "first_objective": 3.0, "center_height": 0.3},
            {"weight": 1.0, "first_objective": 2.0, "center_height": 0.4},
        ])
        self.assertTrue(front["finite"])
        self.assertTrue(front["nondominated"])

    def test_figure13_binary_search_and_distances_match_paper_definition(self) -> None:
        runner = load_runner()
        critical = runner.binary_search_critical_distance(
            lambda distance: distance <= 1.0,
            lower=0.0,
            upper=2.0,
            tolerance=1e-9,
        )
        self.assertLessEqual(abs(critical - 1.0), 1e-9)
        distances = runner.feasibility_distances(critical)
        self.assertEqual(len(distances), 9)
        for exponent, distance in enumerate(distances, 1):
            self.assertAlmostEqual(critical - distance, 10.0 ** (-exponent))

    def test_degenerate_one_dimensional_examples_are_exact(self) -> None:
        runner = load_runner()
        examples = runner.degenerate_one_dimensional_examples()
        self.assertEqual(examples[0]["optimizer_cardinality"], "UNCOUNTABLY_INFINITE")
        self.assertEqual(examples[1]["optimizer_cardinality"], "ONE")
        self.assertEqual(examples[1]["unique_fill"], 1.0)
        self.assertTrue(all(math.isclose(row["gradient_at_optimum"], 0.0) for row in examples))

    def test_final_payload_has_timing_independent_scientific_hash(self) -> None:
        runner = load_runner()
        first = runner.finalize_payload({"schema": "x", "elapsed_seconds": 1.0})
        second = runner.finalize_payload({"schema": "x", "elapsed_seconds": 9.0})
        self.assertEqual(first["scientific_hash"], second["scientific_hash"])


class SpinItFasterPaperFigureNativeTests(unittest.TestCase):
    def test_three_point_static_and_top_pareto_sweeps(self) -> None:
        runner = load_runner()
        manifest = json.loads(
            (ROOT / "experiments/spin_it_reproduction/faster/paper_cases.json").read_text(
                encoding="utf-8"
            )
        )
        cases = {case["id"]: case for case in manifest["cases"]}
        static = runner.run_pareto_sweep(
            cases["static-tilted-vase"],
            [0.0, 0.5, 1.0],
        )
        top = runner.run_pareto_sweep(
            cases["heart-spinning-top"],
            [0.0, 0.5, 1.0],
        )
        for sweep in (static, top):
            self.assertEqual(sweep["reference_status"], "CONVERGED")
            self.assertTrue(all(row["status"] == "CONVERGED" for row in sweep["rows"]))
            self.assertTrue(sweep["front_audit"]["finite"])
            self.assertTrue(sweep["front_audit"]["nondominated"])
        self.assertEqual(
            {row["interface_subtype"] for row in static["rows"]},
            {"PLANE"},
        )
        self.assertGreaterEqual(
            len({row["interface_subtype"] for row in top["rows"]}),
            2,
        )

    def test_tilted_cube_feasibility_smoke_uses_unperturbed_warm_start(self) -> None:
        runner = load_runner()
        result = runner.run_tilted_cube_feasibility(tolerance=0.05, upper=0.8)
        self.assertEqual(result["asset_class"], "RECONSTRUCTED_ASSET")
        self.assertEqual(result["unperturbed_status"], "CONVERGED")
        self.assertGreater(result["critical_distance"], 0.1)
        self.assertLess(result["critical_distance"], 0.4)
        self.assertEqual([row["m"] for row in result["rows"]], list(range(1, 10)))
        self.assertTrue(all(row["status"] == "CONVERGED" for row in result["rows"]))
        self.assertTrue(all(row["trace"] for row in result["rows"]))


if __name__ == "__main__":
    unittest.main()
