"""Physical checks, not screenshot checks, for the approved flexible support."""
import unittest
import numpy as np


class FlexibleBearingTests(unittest.TestCase):
    def test_centered_rotor_does_not_excite_support_or_use_a_motor(self):
        from flexible_bearing import support_protocol, simulate
        protocol = support_protocol(.005, .002)
        result = simulate(.005, [0, 0, 0], np.diag([1e-6, 1e-6, 2e-6]),
                          protocol=protocol, dt=.00025, duration=.03, hz=35, phase=0)
        self.assertLess(result['maximum_support_displacement_m'], 1e-12)
        self.assertAlmostEqual(result['records'][-1]['omega'], 2*np.pi*35, places=7)
        self.assertFalse(result['rotor_motor_enabled'])

    def test_eccentric_rotor_excites_support_and_preserves_full_inertia(self):
        from flexible_bearing import support_protocol, simulate
        inertia = np.array([[1.5e-6, .2e-6, .1e-6], [.2e-6, 1e-6, 0], [.1e-6, 0, 2e-6]])
        protocol = support_protocol(.005, .002)
        result = simulate(.005, [.002, 0, 0], inertia,
                          protocol=protocol, dt=.00025, duration=.08, hz=35, phase=0)
        self.assertGreater(result['maximum_support_displacement_m'], 1e-5)
        np.testing.assert_allclose(result['loaded_inertia_kg_m2'], inertia, rtol=1e-12, atol=1e-18)
        self.assertLess(result['maximum_relative_energy_balance_error'], .02)

    def test_changed_rotor_mass_does_not_retune_support(self):
        from flexible_bearing import support_protocol, simulate
        protocol = support_protocol(.005, .002)
        before = simulate(.005, [.001, 0, 0], np.diag([1e-6, 1e-6, 2e-6]),
                          protocol=protocol, dt=.0005, duration=.01, hz=35, phase=0)
        after = simulate(.003, [0, 0, 0], np.diag([.8e-6, .8e-6, 1.5e-6]),
                         protocol=protocol, dt=.0005, duration=.01, hz=35, phase=0)
        self.assertEqual(before['support_protocol'], after['support_protocol'])
        self.assertEqual(before['support_protocol']['frequency_hz'], 10.)
        self.assertEqual(before['support_protocol']['damping_ratio'], .1)

    def test_nonfinite_and_indefinite_inputs_are_rejected(self):
        from flexible_bearing import support_protocol, simulate
        for inertia in (np.diag([-1., 1., 1.]), np.eye(3)*float('nan')):
            with self.assertRaises(ValueError):
                simulate(.005, [0, 0, 0], inertia, protocol=support_protocol(.005, .002),
                         dt=.0005, duration=.01, hz=35, phase=0)


if __name__ == '__main__':
    unittest.main()
