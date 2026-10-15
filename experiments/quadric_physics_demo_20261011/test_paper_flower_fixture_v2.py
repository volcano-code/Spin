"""The editable middle must not secretly protect every hole boundary."""
from pathlib import Path
from fractions import Fraction
import unittest
import paper_flower_fixture as legacy

HERE = Path(__file__).resolve().parent


class OpenMiddleFixtureTests(unittest.TestCase):
    def fixture(self):
        try:
            from paper_flower_fixture_v2 import build_fixture
        except ImportError:
            build_fixture = legacy.build_fixture
        return build_fixture(HERE / 'paper-outline-flower-preview-v3')

    def test_entire_middle_is_editable_with_no_fixed_hole_web(self):
        problem, record, geometry = self.fixture()
        self.assertEqual(geometry['original'].symmetric_difference(geometry['design']).area, 0)
        self.assertTrue(geometry['fixed_middle'].is_empty)
        self.assertNotIn('fixed-middle-boundary-rim', [x['id'] for x in record['layers']])
        self.assertEqual(len(geometry['design'].interiors), 48)

    def test_real_fixed_mass_keeps_layers_and_bearing_density(self):
        problem, record, geometry = self.fixture()
        self.assertEqual([x['id'] for x in record['layers']],
                         ['bottom', 'top', 'bearing-seat', 'outer-race'])
        mass = sum((legacy.exact_extrusion_moments(x)[0] for x in record['layers']), Fraction(0))
        self.assertEqual(str(mass), record['fixed_exact_rational_moments'][0])
        self.assertEqual(record['layers'][-1]['volume_density'], 6.5)
        self.assertEqual(geometry['seat'].intersection(geometry['design']).area, 0)
        self.assertEqual(set(record['excluded_from_rotor']),
                         {'upper_cap', 'lower_cap', 'shaft', 'inner_race'})

    def test_no_saved_circle_or_cut_location_drives_initialization(self):
        problem, record, geometry = self.fixture()
        self.assertIsNone(record['target_cut_position'])
        self.assertFalse(record['postsolve_rotation'])
        self.assertNotIn('initial_center', problem)
        self.assertNotIn('initial_radius', problem)
        self.assertEqual(record['protocol']['seeds'], [0, 1])


if __name__ == '__main__':
    unittest.main()
