"""Geometry/density regressions for the frozen paper-image flower."""
from pathlib import Path
import json
import unittest
from fractions import Fraction
from paper_flower_fixture import build_fixture, exact_profile_moments, exact_extrusion_moments, boundary_loops

HERE=Path(__file__).resolve().parent


class FlowerFixtureTests(unittest.TestCase):
    def test_polygon_moments_against_asymmetric_rectangle(self):
        loop=[[1.,2.],[4.,2.],[4.,6.],[1.,6.]]
        actual=exact_profile_moments([loop])
        self.assertEqual(actual,[Fraction(12),Fraction(30),Fraction(48),
                                 Fraction(84),Fraction(208),Fraction(120)])

    def test_all_holes_and_true_metric_partition(self):
        problem,record,g=build_fixture(HERE/'paper-outline-flower-preview-v3')
        self.assertEqual(len(g['original'].interiors),48)
        self.assertEqual(record['admission']['design_component_count'],2)
        self.assertEqual(g['design'].intersection(g['fixed_middle']).area,0)
        self.assertEqual(g['original'].symmetric_difference(g['design'].union(g['fixed_middle'])).area,0)
        self.assertAlmostEqual(g['design'].boundary.distance(g['original'].boundary),.015,places=12)
        self.assertIsNone(record['target_cut_position'])
        self.assertNotIn('initial_center',problem)

    def test_mass_density_and_rotating_bearing_role(self):
        problem,record,g=build_fixture(HERE/'paper-outline-flower-preview-v3')
        layers=record['layers']
        self.assertEqual(layers[-1]['volume_density'],6.5)
        fixed_mass=sum(exact_extrusion_moments(layer)[0] for layer in layers)
        self.assertEqual(str(fixed_mass),record['fixed_exact_rational_moments'][0])
        wrong=dict(layers[-1],volume_density=1.)
        self.assertGreater(exact_extrusion_moments(layers[-1])[0],exact_extrusion_moments(wrong)[0])
        self.assertEqual(set(record['excluded_from_rotor']),{'upper_cap','lower_cap','shaft','inner_race'})
        self.assertEqual(g['seat'].intersection(g['outer_race']).area,0)
        self.assertEqual(g['seat'].intersection(g['original']).area,0)

    def test_frozen_objs_reload_and_match_independent_polygon_volume(self):
        import trimesh
        folder=HERE/'paper-flower-method-v1'
        record=json.loads((folder/'input.json').read_text())
        for layer in record['layers']:
            mesh=trimesh.load(folder/(layer['id']+'.obj'),force='mesh',process=False)
            exact=exact_profile_moments(layer['loops'])[0]*Fraction(layer['z_interval'][1]-layer['z_interval'][0])
            expected=float(exact)*record['meters_per_model_unit']**3
            self.assertTrue(mesh.is_watertight)
            self.assertTrue(mesh.is_winding_consistent)
            self.assertLess(abs(mesh.volume-expected)/expected,1e-10)


if __name__=='__main__':unittest.main()
