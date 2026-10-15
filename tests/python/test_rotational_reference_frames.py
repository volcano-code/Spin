"""Regression for central-inertia free rotation versus fixed-pivot rotation."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    'g7_frame_subject', ROOT/'experiments/teapot_spin2014/reconstructed_g7_dynamics.py')
SUBJECT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SUBJECT)
from heart_top import simulate


class RotationalFrameTests(unittest.TestCase):
    def run_protocol(self, center, inertia=None, protocol='FREE_CENTER_OF_MASS'):
        self.assertTrue(hasattr(SUBJECT, 'simulate_rotational_protocol'),
                        'separate origin-aware rotational protocols are required')
        return SUBJECT.simulate_rotational_protocol(2., np.asarray(center, float),
            np.diag([2., 3., 4.]) if inertia is None else inertia,
            np.array([1.1, .7, 5.]), protocol=protocol,
            duration=.125, tilt_degrees=2., azimuth=.7, max_step=1/240)

    def test_free_orientation_is_invariant_to_model_origin(self):
        first = self.run_protocol([0., 0., 0.])
        shifted = self.run_protocol([1., 2., 3.])
        np.testing.assert_array_equal(first['quaternions'], shifted['quaternions'])
        np.testing.assert_array_equal(first['angular_velocity_body'], shifted['angular_velocity_body'])
        self.assertEqual(shifted['reference']['inertia_origin'], 'CENTER_OF_MASS')
        self.assertEqual(shifted['reference']['center_from_rotation_origin'], [0., 0., 0.])

    def test_fixed_pivot_uses_full_parallel_axis_tensor(self):
        actual = self.run_protocol([1., 2., 3.], protocol='IDEAL_FIXED_PIVOT_GRAVITY')
        # Hand-computed Ip = diag(2,3,4) + 2*(14*E - c*c^T).
        Ip = np.array([[28., -4., -6.], [-4., 23., -12.], [-6., -12., 14.]])
        expected = simulate(2., [1., 2., 3.], Ip, [1.1, .7, 5.], gravity=9.81,
                            duration=.125, tilt_degrees=2., azimuth=.7, max_step=1/240)
        np.testing.assert_array_equal(actual['quaternions'], expected['quaternions'])
        np.testing.assert_array_equal(actual['reference']['inertia'], Ip)
        self.assertEqual(actual['reference']['inertia_origin'], 'FIXED_PIVOT')

    def test_old_gravity_zero_pivot_motion_is_not_free_center_rotation(self):
        correct = self.run_protocol([1., 2., 3.])
        old_joint = simulate(2., [1., 2., 3.],
            [[28., -4., -6.], [-4., 23., -12.], [-6., -12., 14.]], [1.1, .7, 5.],
            gravity=0., duration=.125, tilt_degrees=2., azimuth=.7, max_step=1/240)
        error = (Rotation.from_quat(correct['quaternions']).inv() *
                 Rotation.from_quat(old_joint['quaternions'])).magnitude()
        self.assertGreater(np.max(error), .001)

    def test_unrecognized_protocol_and_invalid_central_inertia_are_rejected(self):
        for protocol in ('TORQUE_FREE_UNSPECIFIED_ORIGIN', '', None):
            with self.subTest(protocol=protocol), self.assertRaises(ValueError):
                self.run_protocol([1., 2., 3.], protocol=protocol)
        # A parallel-axis addition must not hide a negative central inertia.
        with self.assertRaises(ValueError):
            self.run_protocol([1., 2., 3.], np.diag([-.1, 3., 4.]),
                              protocol='IDEAL_FIXED_PIVOT_GRAVITY')


if __name__ == '__main__':
    unittest.main()
