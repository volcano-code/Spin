"""Input contract: deleting a fixed layer or filling a hole must fail."""
import importlib.util
from pathlib import Path
import unittest

class FlowerInputTests(unittest.TestCase):
    def test_original_material_contains_all_three_layers(self):
        path = Path(__file__).with_name('perforated_flower_physics.py')
        self.assertTrue(path.exists(), 'three-layer physics adapter is missing')
        spec = importlib.util.spec_from_file_location('flower_physics', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        from perforated_flower import build_input
        problem, _ = build_input()
        parts = module.original_entities(problem)
        self.assertEqual(len(parts), 5)
        self.assertEqual(parts[0]['loops'], problem['design_region_loops'])
        self.assertEqual(parts[0]['z_interval'], [-.07, .07])
        self.assertEqual([x['z_interval'] for x in parts[1:3]], [[-.10,-.07],[.07,.10]])

    def test_multilayer_perforated_partition(self):
        path = Path(__file__).with_name('perforated_flower.py')
        self.assertTrue(path.exists(), 'paper-structure flower builder is not implemented')
        spec = importlib.util.spec_from_file_location('flower_input', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        from shapely.geometry import Polygon
        problem, geometry = module.build_input()
        def polygon(loops):
            return Polygon(loops[0], loops[1:])
        original = polygon(geometry['original_loops'])
        design = polygon(problem['design_region_loops'])
        self.assertTrue(original.is_valid and design.is_valid)
        self.assertEqual(len(original.interiors), 13)  # 12 petal holes + bearing bore
        self.assertLess(original.difference(original.buffer(0)).area, 1e-12)
        layers = problem['fixed_extrusion_layers']
        self.assertEqual(len(layers), 4)  # top, bottom, middle rim, bearing ring
        self.assertEqual(layers[0]['z_interval'], [-.10, -.07])
        self.assertEqual(layers[1]['z_interval'], [.07, .10])
        middle = polygon(layers[2]['loops']).union(polygon(layers[3]['loops']))
        self.assertLess(middle.intersection(design).area, 1e-12)
        self.assertLess(middle.union(design).symmetric_difference(original).area, 1e-12)
        self.assertGreater(abs(original.centroid.x)+abs(original.centroid.y), 1e-4)
        self.assertEqual(problem['design_areal_density'], .14)
        self.assertNotIn('initial_center', problem)

if __name__ == '__main__':
    unittest.main()
