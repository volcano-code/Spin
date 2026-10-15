"""Prevent a frozen-input diagnostic from silently accepting another problem."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
import audit_lateral_capacity as audit


class CapacityScopeTests(unittest.TestCase):
    def setUp(self):
        self.identity = dict(
            target=dict(mode='SPINNING_TOP', origin=[0., 0., 0.], axis=[0., 0., 1.]),
            fixed_regions=[dict(density=1.)], design_regions=[dict(density=1.)])

    def test_rejects_changed_origin_axis_or_density(self):
        check = getattr(audit, 'validate_diagnostic_assumptions', None)
        self.assertTrue(callable(check), 'capacity assumptions must be enforced')
        check(self.identity)
        for field, value in [('origin', [0.1, 0., 0.]), ('axis', [1., 0., 0.])]:
            edited = copy.deepcopy(self.identity)
            edited['target'][field] = value
            with self.assertRaises(ValueError):
                check(edited)
        for regions in ('fixed_regions', 'design_regions'):
            for density in (2., float('nan'), 0.):
                edited = copy.deepcopy(self.identity)
                edited[regions][0]['density'] = density
                with self.assertRaises(ValueError):
                    check(edited)

    def test_digest_mismatch_is_rejected_without_python_assertions(self):
        check = getattr(audit, 'require_partition_identity', None)
        self.assertTrue(callable(check), 'partition mismatch needs an explicit exception')
        check('same', 'same')
        with self.assertRaises(ValueError):
            check('changed', 'frozen')

    def test_entrypoint_checks_assumptions_before_reading_geometry(self):
        edited = copy.deepcopy(self.identity)
        edited['target']['origin'] = [0.1, 0., 0.]
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / 'input-admission.json').write_text(
                json.dumps(dict(problem_identity=edited)))
            previous = audit.CASE
            try:
                audit.CASE = root
                with self.assertRaisesRegex(ValueError, 'zero origin'):
                    audit.main()
            finally:
                audit.CASE = previous


if __name__ == '__main__':
    unittest.main()
