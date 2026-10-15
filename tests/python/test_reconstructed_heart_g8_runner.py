import importlib.util
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
RUNNER_PATH = ROOT / "experiments/spin_it_reproduction/strict/run_reconstructed_heart_g8.py"
SPEC = importlib.util.spec_from_file_location("run_reconstructed_heart_g8", RUNNER_PATH)
RUNNER = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(RUNNER)


class ReconstructedHeartG8RunnerTests(unittest.TestCase):
    def test_observed_levels_reads_native_count_arrays_without_confusing_budget_with_observation(self):
        result = RUNNER._observed_levels({
            "solve": {"passes": [
                {"level_counts": [0, 0, 3, 0],
                 "output_level_counts": [0, 1, 2, 0]},
                {"level_counts": [0, 0, 0, 4],
                 "output_fractional_level_counts": [0, 0, 0, 2]},
            ]}
        })

        self.assertEqual(result["maximum_observed_level"], 3)
        self.assertEqual(result["level_counts_union"], {"1": 1, "2": 3, "3": 4})

    def test_evidence_parameters_are_normalized_to_json_arrays(self):
        payload = RUNNER.evidence_parameters(
            {"bounds": ((-1.0, 1.0),)}, {"center": (0.0, 0.0, 0.95)})

        self.assertEqual(payload["problem_identity"]["bounds"], [[-1.0, 1.0]])
        self.assertEqual(payload["config"]["center"], [0.0, 0.0, 0.95])

    def test_runner_declares_its_pybind11_build_dependency_root(self):
        roots = RUNNER.declared_dependency_roots()

        self.assertEqual(set(roots), {"pybind11"})
        self.assertTrue(roots["pybind11"].is_dir())
        self.assertTrue((roots["pybind11"] / "pybind11" / "pybind11.h").is_file())

    def test_g8_protocol_changes_only_the_declared_depth_caps(self):
        # The runner is launched by file path from an arbitrary working
        # directory, so it must resolve the frozen Heart generator itself.
        previous_path = sys.path[:]
        try:
            sys.path[:] = [entry for entry in sys.path
                           if Path(entry or ".").resolve() != ROOT]
            problem, source_settings, solver_settings, _ = RUNNER.build_protocol()
        finally:
            sys.path[:] = previous_path

        self.assertEqual(source_settings.maximum_depth, 8)
        self.assertEqual(solver_settings.maximum_optimization_level, 8)
        self.assertEqual(solver_settings.initial_level, 3)
        self.assertEqual(solver_settings.maximum_leaf_count, 4000)
        self.assertEqual(solver_settings.objective, "TOP_MASS_AND_LOWERING")
        self.assertEqual(solver_settings.formulation, "PAPER_GIVENS")
        self.assertEqual(solver_settings.binary_constraint_tolerance, 1e-6)

        prior = __import__("json").loads(
            (ROOT / "build/reconstructed-heart-current-20260923-run1/input.json").read_text()
        )
        actual = problem.identity
        for region, label in (("original", "original"),):
            self.assertEqual(
                actual[region]["vertices"]["sha256"],
                prior[label]["vertices"]["sha256"],
            )
            self.assertEqual(
                actual[region]["faces"]["sha256"],
                prior[label]["faces"]["sha256"],
            )

        self.assertEqual(actual["design_regions"][0]["mesh"]["vertices"]["sha256"],
                         prior["design_regions"][0]["mesh"]["vertices"]["sha256"])
        self.assertEqual(actual["design_regions"][0]["mesh"]["faces"]["sha256"],
                         prior["design_regions"][0]["mesh"]["faces"]["sha256"])
        self.assertEqual(actual["fixed_regions"][0]["mesh"]["vertices"]["sha256"],
                         prior["fixed_regions"][0]["mesh"]["vertices"]["sha256"])
        self.assertEqual(actual["fixed_regions"][0]["mesh"]["faces"]["sha256"],
                         prior["fixed_regions"][0]["mesh"]["faces"]["sha256"])


if __name__ == "__main__":
    unittest.main()
