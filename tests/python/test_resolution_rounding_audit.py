"""Diagnostic algebra only: these stencils are not paper model substitutes."""
import importlib.util
import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'experiments/teapot_spin2014'))


def fixture():
    outer=np.array([3.,.251,0,1,1,1,1,0,0,0])
    cells=np.array([[1,.2,0,.3,.1,.1,.1,0,.1,0],
                    [1,.3,0,.4,.2,.2,.2,0,.2,0]])
    return dict(settings=dict(source_settings=dict(root_side=2.),density=1.,
                solver_settings=dict(objective='TOP_INERTIA_ONLY')),
        source=dict(outer_surface_moments=outer.tolist()),solve=dict(
            status='BINARY_VALIDATION_FAILED',final_phi=0.,
            relaxed_fill=[.49,.51],binary_fill=[0.,1.],
            binary_moments=(outer-cells[1]).tolist(),
            final_partition=[dict(level=3,moments=m.tolist()) for m in cells]))


class RoundingAuditTests(unittest.TestCase):
    def audit(self,record):
        self.assertIsNotNone(importlib.util.find_spec('resolution_rounding_audit'),
                             'missing independent rounding decomposition')
        from resolution_rounding_audit import rounding_diagnostics
        return rounding_diagnostics(record)

    def test_equal_removed_mass_can_still_move_com(self):
        r=self.audit(fixture())
        self.assertTrue(r['available'])
        self.assertAlmostEqual(r['raw_moment_delta'][0],0.)
        self.assertAlmostEqual(r['raw_moment_delta'][1],-.049)
        self.assertAlmostEqual(r['constraint_delta_normalized'][0],-.049/16)
        self.assertAlmostEqual(r['constraint_delta_normalized'][2],-.049/32)
        self.assertAlmostEqual(r['relaxed_constraints_normalized'][0],0.)
        self.assertEqual(r['rounding_changed_cells'],2)
        self.assertEqual(r['by_level']['3']['cells'],2)

    def test_reject_wrong_moments_and_wrong_rounding(self):
        r=fixture();r['solve']['binary_moments'][1]+=.01
        with self.assertRaisesRegex(ValueError,'MOMENT_MISMATCH'):self.audit(r)
        r=fixture();r['solve']['binary_fill'][0]=1.
        with self.assertRaisesRegex(ValueError,'ROUNDING_MISMATCH'):self.audit(r)

    def test_unknown_is_not_diagnosed_as_rounding_failure(self):
        r=fixture();r['solve']['status']='RELAXATION_UNKNOWN';r['solve']['binary_fill']=[]
        self.assertFalse(self.audit(r)['available'])

    def test_givens_cross_constraint_uses_actual_final_angle(self):
        r=fixture();r['solve']['final_phi']=np.pi/4
        r['source']['outer_surface_moments'][4]+=.4
        r['solve']['binary_moments'][4]+=.4
        self.assertAlmostEqual(self.audit(r)['binary_constraints_normalized'][4],.2/32)

    def test_reject_nonfinite_and_out_of_range_fill(self):
        for bad in (float('nan'),float('inf'),-.01,1.01):
            r=fixture();r['solve']['relaxed_fill'][0]=bad
            with self.assertRaises(ValueError):self.audit(r)

    def test_reject_overflow_from_finite_trace_numbers(self):
        r=fixture();r['settings']['density']=4.
        r['source']['outer_surface_moments'][0]=1e308
        r['solve']['final_partition'][0]['moments'][0]=1e308
        r['solve']['final_partition'][1]['moments'][0]=1e308
        # Scale the other claimed moments to the new physical density.
        r['solve']['binary_moments']=(4*np.array(r['solve']['binary_moments'])).tolist()
        with np.errstate(over='ignore',invalid='ignore'):
            with self.assertRaises(ValueError):self.audit(r)


if __name__=='__main__':unittest.main()
