import unittest,importlib.util
from pathlib import Path
import numpy as np

class SolidViewTests(unittest.TestCase):
    def test_extrusion_preserves_hole_volume_and_closed_faces(self):
        path=Path(__file__).with_name('solid_view.py')
        self.assertTrue(path.exists(),'Missing real 3D extrusion')
        spec=importlib.util.spec_from_file_location('sv',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
        mesh=m.extrude([[[0,0],[4,0],[4,4],[0,4]],[[1,1],[1,3],[3,3],[3,1]]],-.1,.1)
        self.assertTrue(mesh.is_watertight)
        self.assertTrue(mesh.is_winding_consistent)
        self.assertAlmostEqual(mesh.volume,2.4,places=12)
        self.assertEqual(len(mesh.split()),1)

if __name__=='__main__':unittest.main()
