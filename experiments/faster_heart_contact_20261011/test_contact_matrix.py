import copy
import json
from pathlib import Path
import unittest
import numpy as np
from scipy.spatial.transform import Rotation

from contact_matrix import scenarios, physical_properties


class ContactMatrixTests(unittest.TestCase):
    def setUp(self):
        self.protocol = json.loads((Path(__file__).parent/'rounded-asymmetric-v2/contact-protocol.json').read_text())

    def test_all_predefined_scenarios_and_video_direction(self):
        cases = scenarios(self.protocol)
        self.assertEqual(len({x['case_id'] for x in cases}), 64)
        self.assertEqual(sum(x['azimuth_degrees'] == self.protocol['video_direction_degrees'] for x in cases), 8)
        for field, value in [('azimuth_degrees', [0]), ('rolling_friction', .001),
                             ('persistent_drive_torque', True), ('unique_scenarios', 63)]:
            changed = copy.deepcopy(self.protocol); changed[field] = value
            with self.assertRaises(ValueError):
                scenarios(changed)

    def test_full_inertia_units_and_density_applied_once(self):
        # Independent analytic covariance of a rotated rectangular solid.
        mass_model = 3.
        center_model = np.array([.2, -.1, 1.1])
        rotation = Rotation.from_euler('xyz', [12, 23, 34], degrees=True).as_matrix()
        covariance = rotation @ np.diag(np.array([1., 2., 3.])**2/12) @ rotation.T
        q = mass_model*(covariance+np.outer(center_model, center_model))
        raw = [mass_model, *(mass_model*center_model), *np.diag(q), q[1, 2], q[0, 2], q[0, 1]]
        m, c, inertia, eigenvalues, axes = physical_properties(raw, mass_unit_kg=.01875, length_unit_m=.025)
        expected = mass_model*.01875*.025**2*(np.trace(covariance)*np.eye(3)-covariance)
        self.assertAlmostEqual(m, .05625)
        np.testing.assert_allclose(c, center_model*.025, atol=1e-16)
        np.testing.assert_allclose(inertia, expected, rtol=1e-13, atol=1e-18)
        np.testing.assert_allclose(axes@np.diag(eigenvalues)@axes.T, expected, atol=1e-18)
        self.assertGreater(abs(inertia[0, 2]), 1e-9)
        doubled = physical_properties(np.array(raw)*2, mass_unit_kg=.01875, length_unit_m=.025)
        self.assertAlmostEqual(doubled[0], 2*m)
        np.testing.assert_allclose(doubled[1], c, atol=1e-16)
        np.testing.assert_allclose(doubled[2], 2*inertia, atol=1e-18)

    def test_invalid_and_nonphysical_states_are_rejected(self):
        raw = np.array([1., 0, 0, 0, 1, 1, 1, 0, 0, 0])
        for index, value in [(0, 0.), (0, -1.), (8, float('nan')), (4, -3.)]:
            changed = raw.copy(); changed[index] = value
            with self.assertRaises(ValueError):
                physical_properties(changed, mass_unit_kg=.01875, length_unit_m=.025)


if __name__ == '__main__':
    unittest.main()
