"""Independent simplex derivatives for exact boundaries with collapsed displays."""
from fractions import Fraction
import unittest
import numpy as np

from make_it_stand.verification import cherchi_motion_differential as transport

class RationalBoundaryTransportTests(unittest.TestCase):
    def test_streaming_exact_sum_preserves_cancellation_and_odd_tail(self):
        accumulator=getattr(transport,'_ExactFractionSum',None)
        self.assertTrue(callable(accumulator),'streaming exact accumulator missing')
        total=accumulator()
        for value in (Fraction(2**200),Fraction(1,3),-Fraction(2**200),Fraction(1,7),Fraction(1,11)):
            total.add(value)
        self.assertEqual(total.total(),Fraction(131,231))
        self.assertEqual(total.total(),Fraction(131,231))
        total.add(-Fraction(131,231))
        self.assertEqual(total.total(),0)

    def test_streaming_exact_sum_matches_fraction_reference_with_long_denominators(self):
        accumulator=getattr(transport,'_ExactFractionSum',None)
        self.assertTrue(callable(accumulator),'streaming exact accumulator missing')
        total=accumulator();reference=Fraction(0)
        self.assertEqual(total.total(),0)
        for index in range(257):
            value=Fraction((-1)**index,2**(index+10)+1)
            total.add(value);reference+=value
        self.assertEqual(total.total(),reference)

    def test_primary_surface_tetrahedron_has_exact_density_and_mixed_moments(self):
        vertices=np.asarray(((0.,0.,0.),(1.,0.,0.),(0.,1.,0.),(0.,0.,1.)))
        faces=np.asarray(((0,2,1),(0,1,3),(0,3,2),(1,2,3)))
        exact=tuple(tuple(Fraction(float(x)) for x in row) for row in vertices)
        actual=transport.integrate_rational_polyhedral_raw_moments(vertices,faces,density=2.,rational_vertices=exact)
        expected=(Fraction(1,3),*(Fraction(1,12) for _ in range(3)),
                  *(Fraction(1,30) for _ in range(3)),*(Fraction(1,60) for _ in range(3)))
        np.testing.assert_array_equal(actual,np.asarray(expected,dtype=float))

    def fixture(self):
        epsilon = Fraction(1, 2**60)
        exact = ((Fraction(0),Fraction(0),Fraction(1)),
                 (Fraction(1,2),Fraction(0),1+epsilon),
                 (Fraction(1),Fraction(0),Fraction(1)),
                 (Fraction(0),Fraction(1),Fraction(0)))
        display = np.asarray(exact, dtype=float)
        faces = np.asarray(((0,2,1),(0,1,3),(0,3,2),(1,2,3)))
        return epsilon, exact, display, faces

    def test_exact_orientation_preserves_positive_sliver_and_rejects_float_only(self):
        epsilon, exact, display, faces = self.fixture()
        with self.assertRaises(transport.CherchiDifferentialUnavailable):
            transport.normalize_closed_component_orientation(display, faces)
        oriented, audit = transport.normalize_closed_component_orientation(
            display, faces, rational_vertices=exact)
        np.testing.assert_array_equal(oriented, faces)
        self.assertEqual(audit['coordinate_authority'], 'RATIONAL_BOUNDARY')
        self.assertEqual(Fraction(audit['components'][0]['exact_signed_volume']), epsilon/6)
        self.assertEqual(audit['components'][0]['orientation_sign'], 1)
        self.assertFalse(audit['floating_export_approved'])
        self.assertEqual(audit['near_degenerate_predicate'], 'exact_rational_area_nonzero_v1')

    def test_reynolds_derivative_uses_exact_sliver_not_collapsed_display(self):
        epsilon, exact, display, faces = self.fixture()
        motion = np.zeros((4,3,1)); motion[1,2,0] = 1
        jac = transport.differentiate_polyhedral_raw_moments(
            display, faces, motion, rational_vertices=exact)
        expected = (Fraction(1,6),Fraction(1,16),Fraction(1,24),Fraction(1,8)+epsilon/12,
            Fraction(7,240),Fraction(1,60),Fraction(1,10)+2*epsilon/15+epsilon**2/20,
            Fraction(1,40)+epsilon/60,Fraction(1,20)+epsilon/30,Fraction(1,80))
        np.testing.assert_allclose(jac[:,0], np.asarray(expected,dtype=float), rtol=1e-14, atol=1e-16)
        # The sign/density contract must remain linear for inward cavity shells.
        inward = transport.differentiate_polyhedral_raw_moments(
            display, faces[:,(0,2,1)], motion, density=2.5, rational_vertices=exact)
        np.testing.assert_allclose(inward[:,0], -2.5*np.asarray(expected,dtype=float), rtol=1e-14, atol=1e-16)

    def test_primary_surface_moments_preserve_sliver_and_match_simplex_integrals(self):
        e,exact,display,faces=self.fixture()
        primary=getattr(transport,'integrate_rational_polyhedral_raw_moments',None)
        self.assertTrue(callable(primary),'primary surface-integral map is missing')
        actual=primary(display,faces,density=2.5,rational_vertices=exact)
        expected=(e/6,e/16,e/24,e*(3+e)/24,7*e/240,e/60,
            e/10+e**2/15+e**3/60,e/40+e**2/120,e/20+e**2/60,e/80)
        np.testing.assert_allclose(actual,2.5*np.asarray(expected,dtype=float),rtol=1e-14,atol=0.)
        with self.assertRaises(transport.CherchiDifferentialUnavailable):
            primary(display,faces[:3],density=2.5,rational_vertices=exact)
        with self.assertRaises(transport.CherchiDifferentialUnavailable):
            primary(display,faces[:,(0,2,1)],density=2.5,rational_vertices=exact)

    def test_rational_input_rejects_zero_area_and_mismatched_display(self):
        _, exact, display, faces = self.fixture()
        collinear = tuple(tuple(Fraction(float(value)) for value in row) for row in display)
        for function, args in ((transport.normalize_closed_component_orientation, (display,faces)),
            (transport.differentiate_polyhedral_raw_moments, (display,faces,np.zeros((4,3,1))))):
            with self.subTest(function=function.__name__), self.assertRaises(transport.CherchiDifferentialUnavailable):
                function(*args, rational_vertices=collinear)
            bad_display=display.copy(); bad_display[0,0]=0.1
            with self.subTest(mismatch=function.__name__), self.assertRaises(transport.CherchiDifferentialUnavailable):
                function(bad_display,*args[1:], rational_vertices=exact)

    def test_exact_component_extent_is_not_lost_when_all_display_points_coincide(self):
        e=Fraction(1,2**60)
        rows=((Fraction(1),)*3,(1+e,Fraction(1),Fraction(1)),
              (Fraction(1),1+e,Fraction(1)),(Fraction(1),Fraction(1),1+e))
        _,_,_,faces=self.fixture()
        display=np.asarray(rows,dtype=float)
        self.assertTrue(np.array_equal(display,np.ones((4,3))))
        oriented,audit=transport.normalize_closed_component_orientation(
            display,faces,rational_vertices=rows)
        np.testing.assert_array_equal(oriented,faces)
        self.assertEqual(Fraction(audit['components'][0]['exact_signed_volume']), e**3/6)

    def test_rational_reynolds_derivative_matches_independent_three_step_integrals(self):
        from make_it_stand.spin import integrate_raw_mass_moments
        # This is an independent signed-tetrahedron integrator, not the
        # Reynolds/barycentric coefficients used by the staged derivative.
        exact=((Fraction(1,8),Fraction(-1,4),Fraction(1,2)),
               (Fraction(9,8),Fraction(1,8),Fraction(5,8)),
               (Fraction(-1,8),Fraction(5,4),Fraction(3,8)),
               (Fraction(3,8),Fraction(-1,8),Fraction(3,2)))
        points=np.asarray(exact,dtype=float);_,_,_,faces=self.fixture()
        motion=np.asarray((((.1,-.2),(.2,.1),(.3,-.1)),
                           ((-.2,.1),(.4,-.3),(.1,.2)),
                           ((.3,.2),(-.1,.1),(.2,.4)),
                           ((.2,-.1),(.1,.3),(-.2,.1))))
        jac=transport.differentiate_polyhedral_raw_moments(
            points,faces,motion,density=2.5,rational_vertices=exact)
        primary=getattr(transport,'integrate_rational_polyhedral_raw_moments',None)
        self.assertTrue(callable(primary))
        np.testing.assert_allclose(primary(points,faces,density=2.5,rational_vertices=exact),
            integrate_raw_mass_moments(points.tolist(),faces.tolist(),2.5),rtol=1e-14,atol=1e-14)
        errors=[]
        for step in (1e-3,1e-4,1e-5):
            step_errors=[]
            for direction in range(2):
                plus=integrate_raw_mass_moments((points+step*motion[:,:,direction]).tolist(),faces.tolist(),2.5)
                minus=integrate_raw_mass_moments((points-step*motion[:,:,direction]).tolist(),faces.tolist(),2.5)
                estimate=(np.asarray(plus)-np.asarray(minus))/(2*step)
                np.testing.assert_allclose(jac[:,direction],estimate,rtol=1e-6,atol=1e-8)
                step_errors.append(float(np.max(np.abs(jac[:,direction]-estimate))))
                h=Fraction(str(step))
                shifted=[]
                for sign in (1,-1):
                    rows=tuple(tuple(exact[i][c]+sign*h*Fraction(float(motion[i,c,direction]))
                        for c in range(3)) for i in range(4))
                    shifted.append(primary(np.asarray(rows,dtype=float),faces,density=2.5,rational_vertices=rows))
                estimate_primary=(shifted[0]-shifted[1])/(2*step)
                np.testing.assert_allclose(jac[:,direction],estimate_primary,rtol=1e-6,atol=1e-8)
            errors.append(max(step_errors))
        self.assertLess(errors[1],errors[0]/20)
        self.assertLess(errors[2],errors[0]/20)

if __name__=='__main__':
    unittest.main()
