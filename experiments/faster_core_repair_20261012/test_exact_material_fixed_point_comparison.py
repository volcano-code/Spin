"""Existing saved actual evaluations, not a new geometry calculation."""
import json
from pathlib import Path
import unittest


class ExactMaterialFixedPointComparison(unittest.TestCase):
    def test_same_geometry_moments_derivatives_changed_policy_identity(self):
        here=Path(__file__).resolve().parent
        old=json.loads((here/'heart-exact-material-gate-diagnosis-job-01/output/evaluation.json').read_text())
        new=json.loads((here/'heart-exact-material-gate-diagnosis-job-02/output/evaluation.json').read_text())
        self.assertEqual(set(old),set(new))
        self.assertNotEqual(old['scientific_hash'],new['scientific_hash'])
        for key in old:
            if key != 'scientific_hash':
                self.assertEqual(old[key],new[key],key)
        self.assertTrue(new['derivative_defined'])
