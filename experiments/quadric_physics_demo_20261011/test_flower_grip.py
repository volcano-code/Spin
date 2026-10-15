"""Prevent clamping the rotating flower or assigning caps to rotating mass."""
import importlib.util
from pathlib import Path
import unittest

class GripTests(unittest.TestCase):
    def test_caps_clear_flower_and_remain_stationary(self):
        path=Path(__file__).with_name('flower_grip_preview.py')
        self.assertTrue(path.exists(),'grip assembly is not implemented')
        spec=importlib.util.spec_from_file_location('grip',path)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        parts=module.assembly_spec()
        for name in ('upper_cap','lower_cap','shaft','inner_race'):
            self.assertEqual(parts[name]['motion'],'FIXED_CLAMP')
        self.assertEqual(parts['outer_race']['motion'],'ROTATING_WITH_FLOWER')
        self.assertGreater(parts['upper_cap']['z'][0],.08)
        self.assertLess(parts['lower_cap']['z'][1],-.08)
        self.assertLess(parts['shaft']['radius'],parts['inner_race']['inner_radius'])
        self.assertLess(parts['inner_race']['radius'],parts['outer_race']['inner_radius'])

if __name__=='__main__':unittest.main()
