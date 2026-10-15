"""Admission composes existing evidence without rerunning its numerical jobs."""
import copy
import unittest
from pathlib import Path
from noodles_core_run import load_evidence, check_evidence

HERE = Path(__file__).resolve().parent


class CoreAdmissionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.evidence = load_evidence(HERE)

    def test_full_frozen_protocol_and_seed_only(self):
        result = check_evidence(self.evidence, 0)
        self.assertEqual(result['request']['settings']['maximum_restarts'], 8)
        self.assertEqual(result['request']['settings']['maximum_function_evaluations_per_restart'], 60)
        self.assertTrue(result['evaluator']['algebraic_finite_map'])
        self.assertNotIn('initial_quadrics', result)
        self.assertEqual(check_evidence(self.evidence, 1)['request']['settings']['deterministic_seed'], 1)

    def test_partial_certificate_rejected(self):
        evidence = copy.deepcopy(self.evidence)
        evidence['branch']['passed_indices'].pop()
        with self.assertRaisesRegex(ValueError, 'FINITE81'):
            check_evidence(evidence, 0)

    def test_mismatched_build_rejected(self):
        evidence = copy.deepcopy(self.evidence)
        evidence['native_sha256'] = 'wrong-build'
        with self.assertRaisesRegex(ValueError, 'BUILD'):
            check_evidence(evidence, 0)

    def test_old_request_and_failed_derivative_rejected(self):
        for key in ('request_sha256', 'manifest_sha256'):
            evidence = copy.deepcopy(self.evidence)
            evidence[key] = 'stale'
            with self.assertRaisesRegex(ValueError, 'IDENTITY'):
                check_evidence(evidence, 0)
        evidence = copy.deepcopy(self.evidence)
        evidence['numerical']['numerical_checks_passed'] = False
        with self.assertRaisesRegex(ValueError, 'NUMERICAL'):
            check_evidence(evidence, 0)


if __name__ == '__main__':
    unittest.main()
