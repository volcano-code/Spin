"""Exact ideal-model criteria checked against independent Euler linearization."""
from fractions import Fraction as Q
import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'experiments/reconstructed_heart_top'))
from upright_linearization import upright_matrix


class FixedPivotCriticalOracleTests(unittest.TestCase):
    def setUp(self):
        name = 'make_it_stand.verification.fixed_pivot_critical'
        self.assertIsNotNone(importlib.util.find_spec(name), 'fixed-pivot critical oracle missing')
        self.check = importlib.import_module(name).principal_axis_critical_report

    def test_asymmetric_threshold_is_about_pivot(self):
        for omega, expected in ((8, 'EXPONENTIALLY_UNSTABLE'), (12, 'SPECTRALLY_STABLE_SIMPLE')):
            row = self.check(1, '0.1', ['0.04','0.06','0.08'], omega)
            self.assertEqual(row['pivot_principal_moments'], ['1/20','7/100','2/25'])
            self.assertEqual(row['asymmetric_sufficient_omega_squared'], '981/10')
            self.assertEqual(row['linear_state'], expected)
            self.assertEqual(row['material_admission'], 'NOT_GRANTED')

    def test_central_inertia_false_positive_is_rejected_by_pivot_analysis(self):
        # Wrong central-origin substitution says 225 > 196.2. Actual fixed
        # pivot tensor is diag(.06,.07,.04) and has growing modes.
        row = self.check(1, '0.2', ['0.02','0.03','0.04'], 15)
        self.assertIsNone(row['asymmetric_sufficient_omega_squared'])
        self.assertEqual(row['asymmetric_condition'], 'NOT_APPLICABLE_AXIS_NOT_MAXIMAL')
        self.assertEqual(row['linear_state'], 'EXPONENTIALLY_UNSTABLE')
        self.assertLess(Q(row['polynomial']['discriminant']), 0)
        values = np.linalg.eigvals(upright_matrix(1., .2, np.diag([.06,.07,.04]), 15.))
        self.assertGreater(values.real.max(), 2.9)

    def test_polynomial_matches_independent_euler_matrix(self):
        for m,h,ic,w in ((1, Q(1,10), [Q(4,100),Q(6,100),Q(8,100)],12),
                         (1, Q(1,5), [Q(2,100),Q(3,100),Q(4,100)],15)):
            row = self.check(m,h,ic,w)
            p = [float(Q(x)) for x in row['pivot_principal_moments']]
            actual = np.poly(upright_matrix(float(m),float(h),np.diag(p),float(w)))
            expected = [1.,0.,float(Q(row['polynomial']['a'])),0.,float(Q(row['polynomial']['b']))]
            np.testing.assert_allclose(actual, expected, rtol=1e-11, atol=1e-9)

    def test_axisymmetric_formula_is_separate_and_equality_is_not_strict_pass(self):
        # Pivot A=B=3, C=4, mg*h=12: classical threshold squared=9.
        for omega, expected in ((2,'EXPONENTIALLY_UNSTABLE'),
                                (3,'MARGINAL_REQUIRES_FURTHER_ANALYSIS'),
                                (4,'SPECTRALLY_STABLE_SIMPLE')):
            row = self.check(1, 1, [2,2,4], omega, gravity=12)
            self.assertEqual(row['axisymmetric_linear_omega_squared'], '9')
            self.assertIsNone(row['asymmetric_sufficient_omega_squared'])
            self.assertEqual(row['linear_state'], expected)

    def test_boundary_of_asymmetric_sufficient_condition_is_not_pass(self):
        # J=(5,6,7), mgh=1, so omega^2=1 is exact sufficient-boundary.
        row = self.check(1,1,[4,5,7],1,gravity=1)
        self.assertEqual(row['asymmetric_condition'], 'NOT_SATISFIED')
        self.assertEqual(row['linear_state'], 'MARGINAL_REQUIRES_FURTHER_ANALYSIS')

    def test_density_and_length_scaling(self):
        a = self.check(1,'0.1',['0.04','0.06','0.08'],12)
        density = self.check(3,'0.1',['0.12','0.18','0.24'],12)
        doubled = self.check(8,'0.2',['1.28','1.92','2.56'],12)
        key = 'asymmetric_sufficient_omega_squared'
        self.assertEqual(a[key], density[key])
        self.assertEqual(Q(doubled[key]), Q(a[key])/2)

    def test_below_sufficient_bound_is_not_an_instability_proof(self):
        row = self.check(1,1,[4,5,7],'0.7',gravity=1)
        self.assertEqual(row['asymmetric_sufficient_omega_squared'], '1')
        self.assertEqual(row['asymmetric_condition'], 'NOT_SATISFIED')
        self.assertEqual(row['linear_state'], 'SPECTRALLY_STABLE_SIMPLE')

    def test_spin_reversal_preserves_stability_criteria(self):
        positive = self.check(1,'0.1',['0.04','0.06','0.08'],12)
        negative = self.check(1,'0.1',['0.04','0.06','0.08'],-12)
        for key in ('polynomial','linear_state','asymmetric_condition',
                    'asymmetric_sufficient_omega_squared','axisymmetric_linear_omega_squared'):
            self.assertEqual(positive[key], negative[key])

    def test_nonphysical_nonfinite_or_ambiguous_inputs_rejected(self):
        for values in ([1,1,3], [0,1,1], [1,1], [1,float('nan'),1], [True,1,1]):
            with self.subTest(values=values):
                with self.assertRaises(ValueError): self.check(1,1,values,10)
        for mass,height,omega in ((0,1,1),(1,0,1),(True,1,1),(1,1,float('inf'))):
            with self.assertRaises(ValueError): self.check(mass,height,[2,3,4],omega)


if __name__ == '__main__': unittest.main()
