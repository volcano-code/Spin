"""Actual exact-polyhedral checks of the fixed shell, not a full wall metric."""
import importlib.util
import unittest

import numpy as np

from make_it_stand.verification.spin_problem import SpinArtifactProblem
from test_exact_material_region import box, combine


class FixedBoundaryShellTests(unittest.TestCase):
    def setUp(self):
        name = 'make_it_stand.verification.fixed_boundary_shell'
        self.assertIsNotNone(importlib.util.find_spec(name), 'fixed boundary-shell checker missing')
        self.check = importlib.import_module(name).check_fixed_boundary_shell

    def fixture(self, gap=.25, wall=.25):
        outer = box()
        interior = box((gap,)*3, (1-gap,)*3)
        shell = combine(outer, (interior[0], interior[1][:, ::-1]))
        problem = SpinArtifactProblem(original=outer,
            fixed_regions=[('shell', shell, 1.)], design_regions=[('interior', interior, 1.)],
            length_unit='m', mass_unit='kg', minimum_outer_wall=wall,
            target={'mode': 'SPINNING_TOP', 'axis': [0.,0.,1.],
                    'origin': [0.,0.,0.], 'support_points': []})
        return problem, outer, interior, shell

    def test_exact_threshold_and_preserved_shell_do_not_grant_full_wall(self):
        p, outer, interior, shell = self.fixture()
        result = self.check(p, outer=outer, design_results=[('interior', None, interior)],
                            final_material=shell, newly_carved=interior)
        self.assertEqual(result['status'], 'FIXED_BOUNDARY_SHELL_VALIDATED')
        self.assertEqual(result['domain_clearances'][0]['clearance']['squared_distance_exact'], '1/16')
        self.assertEqual(result['full_feature_thickness'], 'UNKNOWN')
        self.assertEqual(result['material_admission'], 'NOT_GRANTED')
        self.assertEqual(result['original_paper_reproduction'], 'NOT_GRANTED')
        again = self.check(p, outer=outer, design_results=[('interior', None, interior)],
                           final_material=shell, newly_carved=interior)
        self.assertEqual(result, again)

    def test_empty_carving_cannot_hide_too_thin_input_shell(self):
        p, outer, interior, _ = self.fixture(gap=.125)
        with self.assertRaisesRegex(ValueError, 'design-domain clearance'):
            self.check(p, outer=outer, design_results=[('interior', interior, None)],
                       final_material=outer, newly_carved=None)

    def test_next_float_above_true_clearance_is_rejected(self):
        p, outer, interior, shell = self.fixture(wall=np.nextafter(.25, 1.))
        with self.assertRaisesRegex(ValueError, 'design-domain clearance'):
            self.check(p, outer=outer, design_results=[('interior', None, interior)],
                       final_material=shell, newly_carved=interior)

    def test_missing_fixed_material_cannot_pass(self):
        p, outer, interior, _ = self.fixture()
        with self.assertRaisesRegex(ValueError, 'partition'):
            self.check(p, outer=outer, design_results=[('interior', interior, None)],
                       final_material=interior, newly_carved=None)

    def test_changed_outer_or_incomplete_or_duplicate_results_rejected(self):
        p, outer, interior, _ = self.fixture()
        for changed_outer, rows in (
            ((outer[0] + [0.,0.,.125], outer[1]), [('interior', interior, None)]),
            (outer, []),
            (outer, [('interior', interior, None)]*2),
            (outer, [('unknown', interior, None)]),
        ):
            with self.subTest(rows=rows[0][0] if rows else 'missing'):
                with self.assertRaises(ValueError):
                    self.check(p, outer=changed_outer, design_results=rows,
                               final_material=outer, newly_carved=None)

    def test_removal_outside_design_domain_rejected(self):
        p, outer, interior, _ = self.fixture()
        with self.assertRaisesRegex(ValueError, 'partition'):
            self.check(p, outer=outer, design_results=[('interior', None, outer)],
                       final_material=interior, newly_carved=outer)

    def test_nontrivial_split_and_wrong_aggregate_removal(self):
        p, outer, _, shell = self.fixture()
        kept = box((.25,)*3, (.5,.75,.75))
        removed = box((.5,.25,.25), (.75,)*3)
        # Use the outer-minus-removed boundary, not touching shell/kept meshes.
        material = combine(outer, (removed[0], removed[1][:, ::-1]))
        rows = [('interior', kept, removed)]
        result = self.check(p, outer=outer, design_results=rows,
                            final_material=material, newly_carved=removed)
        self.assertEqual(result['status'], 'FIXED_BOUNDARY_SHELL_VALIDATED')
        for wrong in (None, kept):
            with self.subTest(wrong_is_empty=wrong is None):
                with self.assertRaisesRegex(ValueError, 'partition'):
                    self.check(p, outer=outer, design_results=rows,
                               final_material=material, newly_carved=wrong)

    def test_overlapping_retained_and_carved_rejected(self):
        p, outer, interior, shell = self.fixture()
        with self.assertRaisesRegex(ValueError, 'partition'):
            self.check(p, outer=outer, design_results=[('interior', interior, interior)],
                       final_material=shell, newly_carved=interior)


if __name__ == '__main__':
    unittest.main()
