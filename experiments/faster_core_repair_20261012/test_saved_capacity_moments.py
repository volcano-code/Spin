"""Catch target-frame loss or a bounds check promoted to feasibility."""
import unittest
from pathlib import Path
import sys
import trimesh
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'faster_heart_contact_20261011'))

class SavedCapacityMoments(unittest.TestCase):
    def function(self):
        try:
            from saved_input_capacity import capacity_from_fixed_moments
        except ImportError:
            capacity_from_fixed_moments=None
        self.assertTrue(callable(capacity_from_fixed_moments),'saved exact moment capacity adapter missing')
        return capacity_from_fixed_moments

    def test_off_axis_fixed_mass_is_excluded_by_exact_bound(self):
        mesh=trimesh.creation.box(extents=(2,2,2))
        result=self.function()([1,3,0,0,0,0,0,0,0,0],mesh.vertices,mesh.faces,[0,0,-2])
        self.assertEqual(result['x']['total_min_exact'],'1')
        self.assertEqual(result['x']['total_max_exact'],'5')
        self.assertEqual(result['status'],'EXCLUDED_BY_FIRST_MOMENT_BOUND')
        self.assertFalse(result['feasible_witness'])

    def test_origin_translation_changes_both_fixed_and_design_moments(self):
        mesh=trimesh.creation.box(extents=(2,2,2))
        result=self.function()([1,3,0,0,0,0,0,0,0,0],mesh.vertices,mesh.faces,[3,0,0])
        self.assertEqual(result['x']['fixed_first_exact'],'0')
        self.assertEqual(result['x']['total_min_exact'],'-24')
        self.assertEqual(result['x']['total_max_exact'],'0')
        self.assertEqual(result['status'],'UNKNOWN_NOT_EXCLUDED_BY_FIRST_MOMENT_BOUNDS')
        self.assertFalse(result['feasible_witness'])

if __name__=='__main__':unittest.main()
