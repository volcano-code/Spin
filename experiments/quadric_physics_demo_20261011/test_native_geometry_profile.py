"""Native timing diagnosis must observe real execution and preserve failures."""
import importlib.util
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest


class NativeProfileTests(unittest.TestCase):
    def module(self):
        self.assertIsNotNone(importlib.util.find_spec('native_geometry_profile'),
                             'native geometry stack profiler has not been implemented')
        import native_geometry_profile
        return native_geometry_profile

    def test_caller_identifies_exact_scan_even_inside_gmp(self):
        module = self.module()
        self.assertEqual(module.classify_stack(
            '^done,stack=[frame={func="__gmpq_mul"},frame={func="verification::PairScan::scan"}]'),
            'BOUNDARY_INTERSECTION_SCAN')
        self.assertEqual(module.classify_stack(
            '^done,stack=[frame={func="CGAL::SNC_constructor::build"}]'), 'NEF_CONSTRUCTION')
        self.assertEqual(module.classify_stack('^done,stack=[]'), 'UNRESOLVED')

    def test_real_child_completes_and_diagnostic_grants_no_admission(self):
        module = self.module()
        with TemporaryDirectory() as folder:
            report = module.profile_command(
                [sys.executable, '-c', 'import time; time.sleep(.25); print("observed-child")'],
                Path(folder), interval=.04, maximum_samples=2, timeout=15)
            self.assertEqual(report['inferior_exit_code'], 0)
            self.assertGreater(len(report['samples']), 0)
            self.assertFalse(report['geometry_admission'])
            self.assertFalse(report['optimization_admission'])
            self.assertIn('observed-child', (Path(folder) / 'gdb-mi.log').read_text())
            self.assertEqual(json.loads((Path(folder) / 'profile.json').read_text())['inferior_exit_code'], 0)

    def test_real_child_failure_is_not_promoted_to_success(self):
        module = self.module()
        with TemporaryDirectory() as folder:
            report = module.profile_command([sys.executable, '-c', 'raise SystemExit(7)'],
                Path(folder), interval=.04, maximum_samples=1, timeout=15)
            self.assertEqual(report['inferior_exit_code'], 7)
            self.assertFalse(report['optimization_admission'])


if __name__ == '__main__':
    unittest.main()
