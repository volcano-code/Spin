"""Reject contour loss and incorrect deprojection; no optimization targets."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
import cv2

class OutlineTests(unittest.TestCase):
    def test_rectification_preserves_holes_and_circular_bearing(self):
        path=Path(__file__).with_name('paper_outline_flower.py')
        self.assertTrue(path.exists(), 'reference contour reconstruction is missing')
        spec=importlib.util.spec_from_file_location('outline',path)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        image=np.full((240,400,3),220,np.uint8)
        cv2.ellipse(image,(200,120),(170,85),0,0,360,(40,40,40),-1)
        cv2.ellipse(image,(200,120),(32,16),0,0,360,(220,220,220),-1)
        cv2.ellipse(image,(290,120),(20,10),0,0,360,(220,220,220),-1)
        loops,record=module.trace_outline(image,(200,120))
        self.assertEqual(len(loops),3)
        bore=np.asarray(loops[record['bearing_loop_index']])
        radii=np.linalg.norm(bore,axis=1)
        self.assertLess(np.std(radii)/np.mean(radii),.04)
        # Simplification samples the boundary nonuniformly: vertex averages
        # are not area centroids. Check the geometric center instead.
        moments=cv2.moments(bore.astype(np.float32))
        centroid=np.array([moments['m10'],moments['m01']])/moments['m00']
        self.assertLess(np.linalg.norm(centroid),.025)
        self.assertGreater(module.signed_area(loops[0]),0)
        self.assertTrue(all(module.signed_area(h)<0 for h in loops[1:]))

if __name__=='__main__':unittest.main()
