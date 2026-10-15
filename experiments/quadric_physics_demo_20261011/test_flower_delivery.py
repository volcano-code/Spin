"""Final assembly may not silently refill its cavity with diagnostic solids."""
from pathlib import Path
import json
import unittest

HERE=Path(__file__).resolve().parent


class DeliverySelectionTests(unittest.TestCase):
    def selection(self, record):
        try:
            from flower_delivery import final_parts
        except ImportError:
            return list(record['files'])
        return final_parts(record)

    def test_delivered_parts_exclude_all_preoptimization_and_cavity_solids(self):
        record=json.loads((HERE/'paper-flower-method-v1/prepare-job-03/output/output-preparation.json').read_text())
        names=self.selection(record)
        self.assertEqual(set(names),{'retained',*record['rotor_layers'],*record['clamp_layers']})
        self.assertTrue(set(names).isdisjoint({'original-flower','original-middle','design','cavity','bearing-bore-region'}))

    def test_diagnostic_object_in_physical_roles_is_rejected(self):
        record=dict(files={'retained':{},'cavity':{}},rotor_layers={'cavity':1},clamp_layers={})
        with self.assertRaisesRegex(ValueError,'DIAGNOSTIC_SOLID_IN_FINAL_ASSEMBLY'):
            self.selection(record)


if __name__=='__main__':unittest.main()
