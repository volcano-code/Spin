import unittest
from compare_pilot_repeats import scientific


class ScientificRecordTests(unittest.TestCase):
    def test_only_known_timing_fields_are_removed(self):
        self.assertEqual(scientific({'objective': 2, 'evaluation_payload': {'evaluator_wall_seconds': 9}}),
                         {'objective': 2, 'evaluation_payload': {}})
        self.assertNotEqual(scientific({'objective': 2}), scientific({'objective': 3}))
        self.assertNotEqual(scientific({'reason': 'failed'}), scientific({'reason': 'accepted'}))
        self.assertNotEqual(scientific({'material_geometry_hash': 'a'}),
                            scientific({'material_geometry_hash': 'b'}))


if __name__ == '__main__': unittest.main()
