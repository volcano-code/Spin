import unittest
import numpy as np
from aggregate_rounded_production import checked_rows
from make_it_stand.verification.surface_derivatives import quadric_five_point_layout


class RegisteredDifferentialTests(unittest.TestCase):
    def fixture(self):
        trials,rows=quadric_five_point_layout([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.],[.01]*10)
        registration=dict(base=list(trials[0]),points=list(trials),rows=[dict(column=c,factor=f,indices=ids) for c,f,ids in rows])
        records=[dict(index=i,derivative_defined=True,moment_errors=[0]*10,material_facts={'valid':True},
            normalized_coefficients=list(q),independent_total_moments=list(q),branch_identity='same') for i,q in enumerate(trials)]
        return registration,records,np.eye(10).tolist()

    def test_correct_matrix_passes_all_thirty_rows(self):
        registration,records,jacobian=self.fixture()
        self.assertTrue(all(r['passed'] for r in checked_rows(registration,records,jacobian)))

    def test_wrong_jacobian_and_unproved_branch_are_rejected(self):
        registration,records,jacobian=self.fixture();jacobian[0][0]=2
        self.assertFalse(all(r['passed'] for r in checked_rows(registration,records,jacobian)))
        registration,records,jacobian=self.fixture();records[1]['branch_identity']='changed'
        self.assertFalse(all(r['passed'] for r in checked_rows(registration,records,jacobian)))

    def test_missing_point_cannot_grant_batch(self):
        registration,records,jacobian=self.fixture()
        with self.assertRaisesRegex(ValueError,'ALL_81'):
            checked_rows(registration,records[:-1],jacobian)
