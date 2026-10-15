import unittest
import numpy as np
from independent_eq14 import formulas, dual


class Eq14Tests(unittest.TestCase):
    def test_isotropic_centered_reference(self):
        f, gradient, hessian, constraints, jacobian, ch = formulas()
        r = [1,0,0,0,1,1,1,0,0,0]
        self.assertAlmostEqual(float(f(*r)),2.)
        np.testing.assert_array_equal(np.asarray(constraints(*r)).ravel(),np.zeros(4))
        self.assertAlmostEqual(float(np.asarray(gradient(*r)).ravel()[3]),1.)

    def test_cross_terms_are_doubled_in_dual_coordinates(self):
        q=np.ones(10)/np.sqrt(10)
        p=dual(q)
        self.assertAlmostEqual(p[7]/p[0],2.)


if __name__ == '__main__': unittest.main()
