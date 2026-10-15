import importlib
from pathlib import Path
import unittest
import numpy as np

HERE = Path(__file__).resolve().parent
BUNDLE = HERE/'paper-flower-method-v2/table-scene-v1'


class MuJoCoTests(unittest.TestCase):
    def adapter(self):
        self.assertIsNotNone(importlib.util.find_spec('flower_mujoco'), 'MuJoCo adapter missing')
        return importlib.import_module('flower_mujoco')

    def test_actual_loaded_full_inertia_and_unpowered_balance(self):
        a = self.adapter()
        packet = a.load_input(BUNDLE)
        for variant in ('original', 'optimized'):
            result = a.simulate(packet, variant, 'ideal', .00025, .05, 0)
            np.testing.assert_allclose(result['loaded_inertia'], packet['properties'][variant]['inertia_kg_m2'], rtol=1e-12, atol=1e-18)
            self.assertEqual(result['actuator_count'], 0)
            self.assertLess(result['energy_drift'], 1e-10)
            self.assertAlmostEqual(result['force_mean'], result['analytic_force'], delta=1e-9)

    def test_flexible_support_shared_between_variants(self):
        a = self.adapter(); packet = a.load_input(BUNDLE)
        results = [a.simulate(packet, v, 'flexible', .00025, .1, 0) for v in ('original','optimized')]
        self.assertEqual(results[0]['support'], results[1]['support'])
        self.assertGreater(results[0]['peak_displacement'], 1e-5)
        self.assertLess(results[1]['peak_displacement'], 1e-8)
        for r in results:
            self.assertLess(r['energy_balance_error'], 1e-4)

    def test_nonpositive_inertia_rejected(self):
        a = self.adapter(); packet = a.load_input(BUNDLE)
        packet['properties']['original']['inertia_kg_m2'][0][0] = -1
        with self.assertRaises(ValueError):
            a.simulate(packet, 'original', 'ideal', .00025, .01, 0)


if __name__ == '__main__': unittest.main()
