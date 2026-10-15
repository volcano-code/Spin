import importlib.util
from pathlib import Path
import unittest


class RoundingBreakdownTests(unittest.TestCase):
    def test_rounding_half_removes_material_and_preserves_sign(self):
        path = Path(__file__).resolve().parents[2] / "experiments/spin_it_reproduction/strict/rounding_moment_breakdown.py"
        self.assertTrue(path.is_file(), "a per-cell rounding decomposition is required")
        spec = importlib.util.spec_from_file_location("rounding_breakdown", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        moment = [2., 1., 0., 0., 1., 1., 1., 0., 0., 0.]
        leaves = [{"level": 3, "coordinates": [0,0,0], "moments": moment}]
        output = module.decompose(leaves, [.5], [1.], 2.)
        self.assertEqual(output["moment_delta"], [-2.,-1.,0.,0.,-1.,-1.,-1.,0.,0.,0.])
        self.assertEqual(output["by_level"]["3"]["moment_delta"], output["moment_delta"])
        with self.assertRaises(ValueError):
            module.decompose(leaves, [.5], [0.], 2.)
        with self.assertRaises(ValueError):
            module.decompose(leaves, [float("nan")], [1.], 2.)


if __name__ == "__main__":
    unittest.main()
