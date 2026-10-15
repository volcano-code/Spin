"""Synthetic accounting fixtures; no Heart physics success is asserted."""
import copy
import importlib.util
import json
from pathlib import Path
import unittest
from contact_matrix import scenarios


class ContactSummaryTests(unittest.TestCase):
    def setUp(self):
        self.protocol = json.loads((Path(__file__).parent/'rounded-asymmetric-v2/contact-protocol.json').read_text())
        self.cases = scenarios(self.protocol)
        self.records = [dict(scene=s, status='FAILED', exception_type='FixtureFailure',
            exception_message='synthetic accounting fixture') for s in self.cases]

    def summarize(self, records):
        path = Path(__file__).with_name('contact_summary.py')
        self.assertTrue(path.exists(), 'full-matrix contact accounting is missing')
        spec = importlib.util.spec_from_file_location('contact_summary', path)
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        return module.summarize(self.protocol, records)

    def test_failed_scenes_remain_in_denominator_and_video_is_fixed(self):
        result = self.summarize(self.records)
        self.assertEqual(result['denominator'], 64)
        self.assertEqual(result['failed_scene_count'], 64)
        self.assertEqual(result['recorded_scene_count'], 0)
        self.assertEqual(len(result['video_case_ids']), 8)
        self.assertTrue(all('dir000' in x for x in result['video_case_ids']))
        self.assertFalse(result['physical_validation_passed'])

    def test_missing_duplicate_and_unregistered_records_are_rejected(self):
        for records in (self.records[:-1], self.records+[self.records[0]],
                        self.records[:-1]+[self.records[0]]):
            with self.assertRaises(ValueError):
                self.summarize(records)
        changed = copy.deepcopy(self.records); changed[0]['scene']['dt_seconds'] = .001
        with self.assertRaises(ValueError):
            self.summarize(changed)

    def test_nan_and_incomplete_trajectory_cannot_count_as_recorded(self):
        for end, tilt in ((7.9, 4.), (8., float('nan'))):
            changed = copy.deepcopy(self.records)
            changed[0].update(status='PHYSICAL_TRAJECTORY_RECORDED_NOT_CALIBRATED',
                duration_seconds=8., actual_steps=16000, maximum_tilt_degrees=tilt,
                maximum_penetration_m=0., final_tilt_degrees=4., final_spin_hz=20.,
                records=[dict(time_s=i*.005) for i in range(1601)])
            changed[0]['records'][-1]['time_s'] = end
            reason = 'INCOMPLETE_CONTACT_TRAJECTORY' if end != 8. else 'NONFINITE_CONTACT_METRIC'
            with self.assertRaisesRegex(ValueError, reason):
                self.summarize(changed)

    def test_completed_pair_differences_do_not_discard_failed_cases(self):
        changed = copy.deepcopy(self.records)
        for record in changed[:2]:
            record.update(status='PHYSICAL_TRAJECTORY_RECORDED_NOT_CALIBRATED',
                duration_seconds=8., actual_steps=round(8./record['scene']['dt_seconds']),
                maximum_tilt_degrees=4., maximum_penetration_m=0., final_tilt_degrees=3.,
                final_spin_hz=20., records=[dict(time_s=i*.005) for i in range(1601)])
        changed[1]['final_spin_hz'] = 19.
        result = self.summarize(changed)
        self.assertEqual(result['recorded_scene_count'], 2)
        self.assertEqual(result['failed_scene_count'], 62)
        measured = [r for r in result['comparisons'] if r['status']=='MEASURED_DIFFERENCE']
        self.assertEqual(len(measured), 1)
        self.assertEqual(measured[0]['kind'], 'TIME_STEP')
        self.assertEqual(measured[0]['absolute_differences']['final_spin_hz'], 1.)
        self.assertFalse(result['physical_validation_passed'])


if __name__ == '__main__':
    unittest.main()
