"""Approved FD arithmetic/layout, NOT production map or optimizer evidence."""
from fractions import Fraction as F
import unittest

import numpy as np

from make_it_stand.verification.surface_derivatives import (
    five_point_actual_weights, quadric_five_point_layout,
    compare_quadric_five_point_differential,require_complete_quadric_trial_ledger)


class QuadricFivePointProtocol(unittest.TestCase):
    def test_long_authority_jacobian_and_report_keep_exact_fraction(self):
        # Catch conversion through Fraction(str) or decimal formatting that
        # trips Python's digit guard on legitimate rational geometry values.
        from make_it_stand.verification.rational_text import fraction_to_rational_text, fraction_from_rational_text
        slope=F(1,3)+F(1,10**5000)
        nodes=(-2,-1,0,1,2)
        coefficients=[[F(x)]+[F(0)]*9 for x in nodes]
        moments=[[slope*x]+[F(0)]*9 for x in nodes]
        jacobian=[['0/1']*10 for _ in range(10)]
        jacobian[0][0]=fraction_to_rational_text(slope)
        result=compare_quadric_five_point_differential(nodes,0,coefficients,moments,jacobian)
        self.assertEqual(fraction_from_rational_text(result['prediction_authority'][0]),slope)
        self.assertTrue(result['within_tolerance'])

    def test_new_authority_protocol_preserves_frozen_source_recipe(self):
        # Catch a CLI default silently replacing the frozen compact-source
        # setting; this tests configuration behavior, not source text.
        import importlib.util
        from pathlib import Path
        from make_it_stand.spin import SpinEvaluatorConfig
        path=Path(__file__).resolve().parents[2]/'.superpowers/sdd/2026-10-07-spin2014-objectives-quadric/validate_production_five_point_map.py'
        spec=importlib.util.spec_from_file_location('five_point_config_runner',path)
        runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
        frozen=SpinEvaluatorConfig(derivative_method='DISCRETE_FULL_RANK_QUADRIC',rational_source_coordinates=True,
            compact_rational_source_coordinates=True,icosphere_subdivisions=2,
            maximum_icosphere_subdivisions=2)
        actual=runner.production_evaluator_config(frozen,compact_source=False,
                                                  coefficient_authority=True)
        self.assertEqual(actual.identity_dict(),frozen.identity_dict())

    def test_authority_jacobian_is_not_rounded_before_difference_prediction(self):
        nodes=(-2,-1,0,1,2)
        slope=F(1,3)+F(1,10**25)
        coefficients=[[F(x)]+[F(0)]*9 for x in nodes]
        moments=[[slope*x]+[F(0)]*9 for x in nodes]
        jacobian=[['0/1']*10 for _ in range(10)]
        jacobian[0][0]=f'{slope.numerator}/{slope.denominator}'
        result=compare_quadric_five_point_differential(nodes,0,coefficients,moments,jacobian)
        self.assertEqual(F(result['prediction_authority'][0]),slope)
        self.assertEqual(F(result['finite_difference_authority'][0]),slope)
        self.assertEqual(F(result['per_moment_error_authority'][0]),0)
        self.assertTrue(result['within_tolerance'])

    def test_error_just_above_gate_is_not_rounded_into_a_pass(self):
        nodes=(-2,-1,0,1,2)
        error=F(1,10**6)+F(1,10**40)
        coefficients=[[F(x)]+[F(0)]*9 for x in nodes]
        moments=[[error*x]+[F(0)]*9 for x in nodes]
        result=compare_quadric_five_point_differential(nodes,0,coefficients,moments,
                                                     [['0/1']*10 for _ in range(10)])
        self.assertEqual(F(result['per_moment_error_authority'][0]),error)
        self.assertFalse(result['within_tolerance'])

    def test_missing_duplicate_or_stale_trial_cannot_pass_complete_ledger(self):
        trials,_=quadric_five_point_layout((.3,)*10,(2e-5,)*10)
        # Failure-ledger transport fixtures, NOT fabricated material results.
        records=[dict(index=i,raw_coefficients=row,exception='fixture-only failure')
                 for i,row in enumerate(trials)]
        require_complete_quadric_trial_ledger(trials,records)
        with self.assertRaises(ValueError):
            require_complete_quadric_trial_ledger(trials,records[:-1])
        wrong=[dict(row) for row in records];wrong[-1]['index']=0
        with self.assertRaises(ValueError):
            require_complete_quadric_trial_ledger(trials,wrong)
        wrong=[dict(row) for row in records];wrong[1]['raw_coefficients']=trials[2]
        with self.assertRaises(ValueError):
            require_complete_quadric_trial_ledger(trials,wrong)

    def test_actual_abscissae_weights_reproduce_zero_to_four(self):
        self.assertEqual(five_point_actual_weights((-2.,-1.,0.,1.,2.),0.),
                         (F(1,12),F(-2,3),F(0),F(2,3),F(-1,12)))
        base=.3
        nodes=tuple(base+x*2e-5 for x in (-1.,-.5,0.,.5,1.))
        weights=five_point_actual_weights(nodes,base)
        offsets=tuple(F(x)-F(base) for x in nodes)
        for degree in range(5):
            self.assertEqual(sum(w*x**degree for w,x in zip(weights,offsets)),
                             1 if degree==1 else 0)
        for bad in ((-2.,-1.,0.,0.,2.),(-2.,-1.,float('nan'),1.,2.),(-1.,0.,1.)):
            with self.assertRaises(ValueError):
                five_point_actual_weights(bad,0.)

    def test_layout_has_81_unique_trials_and_preserves_original_three_steps(self):
        base=(.3,)*10;steps=(2e-5,)*10
        trials,rows=quadric_five_point_layout(base,steps)
        self.assertEqual(len(trials),81)
        self.assertEqual(len(set(trials)),81)
        self.assertEqual(len(rows),30)
        for column,factor,indices in rows:
            self.assertEqual(trials[indices[2]],base)
            for index,offset in zip(indices,(-1.,-.5,0.,.5,1.)):
                expected=list(base);expected[column]+=steps[column]*factor*offset
                self.assertEqual(trials[index],tuple(expected))
        with self.assertRaises(ValueError):
            quadric_five_point_layout((1e20,)*10,(1.,)*10)

    def test_same_weights_compare_real_values_and_reject_wrong_jacobian(self):
        nodes=(-2.,-1.,0.,1.,2.)
        # Exact polynomial map purely tests FD arithmetic. It is deliberately
        # NOT relabeled as normalized coefficients or a production Boolean.
        coefficients=[tuple(F(i+1)*F(x)+F(x)**4 for i in range(10)) for x in nodes]
        independent=[tuple(2*y for y in row) for row in coefficients]
        jacobian=2*np.eye(10)
        checked=compare_quadric_five_point_differential(nodes,0.,coefficients,independent,jacobian)
        self.assertEqual(checked['maximum_component_error'],0.)
        self.assertTrue(checked['within_tolerance'])
        jacobian[0,0]+=1e-3
        wrong=compare_quadric_five_point_differential(nodes,0.,coefficients,independent,jacobian)
        self.assertFalse(wrong['within_tolerance'])
        self.assertGreater(wrong['maximum_component_error'],1e-6)
        with self.assertRaises(ValueError):
            compare_quadric_five_point_differential(nodes,0.,coefficients[:-1],independent,jacobian)
        jacobian[0,0]=np.nan
        with self.assertRaises(ValueError):
            compare_quadric_five_point_differential(nodes,0.,coefficients,independent,jacobian)


if __name__=='__main__':
    unittest.main()
