import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'experiments/reconstructed_heart_top'))
from upright_linearization import upright_matrix


class LinearizationTests(unittest.TestCase):
    def test_axisymmetric_critical_speed(self):
        m,h,i,k=1.,.1,.02,.01
        critical=np.sqrt(4*i*m*9.81*h)/k
        low=np.linalg.eigvals(upright_matrix(m,h,np.diag([i,i,k]),.95*critical))
        high=np.linalg.eigvals(upright_matrix(m,h,np.diag([i,i,k]),1.05*critical))
        self.assertGreater(low.real.max(),.1)
        self.assertLess(np.abs(high.real).max(),1e-10)

    def test_matches_differentiated_full_euler_and_vertical(self):
        m,h,I,w=1.3,.2,np.diag([.04,.055,.02]),40.
        def rhs(x):
            n=np.r_[x[:2],np.sqrt(1-x[:2]@x[:2])]; omega=np.r_[x[2:],w]
            ndot=np.cross(n,omega)
            torque=np.cross([0,0,h],-m*9.81*n)
            acceleration=np.linalg.solve(I,torque-np.cross(omega,I@omega))
            return np.r_[ndot[:2],acceleration[:2]]
        eps=1e-6
        actual=np.column_stack([(rhs(eps*e)-rhs(-eps*e))/(2*eps) for e in np.eye(4)])
        np.testing.assert_allclose(upright_matrix(m,h,I,w),actual,rtol=1e-10,atol=1e-10)


if __name__=='__main__': unittest.main()
