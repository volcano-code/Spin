import unittest
import exact_planar_patches as exact


class ExactIncidenceTests(unittest.TestCase):
    def inside(self,p):
        fn=getattr(exact,'point_on_triangle',None)
        self.assertTrue(callable(fn),'exact triangle incidence missing')
        return fn(p,[(0,0,0),(1,0,0),(0,1,0)])

    def test_actual_triangle_not_infinite_plane(self):
        self.assertTrue(self.inside(('1/4','1/4','0/1')))
        self.assertFalse(self.inside((2,2,0)))

    def test_no_near_plane_epsilon(self):
        self.assertFalse(self.inside(('1/4','1/4','1/'+str(10**60))))

    def test_boundary_point_is_incident(self):
        self.assertTrue(self.inside(('1/2','1/2',0)))

    def test_exact_collinear_subdivision_is_not_a_patch_corner(self):
        fn=getattr(exact,'coalesce_collinear_cycle',None)
        self.assertTrue(callable(fn),'exact patch subdivision certificate missing')
        points=[(0,0,0),('1/2',0,0),(1,0,0),(1,1,0),(0,1,0)]
        self.assertEqual(fn(points,[0,1,2,3,4]),[0,2,3,4])

    def test_near_collinear_corner_is_preserved(self):
        fn=getattr(exact,'coalesce_collinear_cycle',None)
        self.assertTrue(callable(fn),'exact patch subdivision certificate missing')
        points=[(0,0,0),('1/2','1/'+str(10**60),0),(1,0,0),(1,1,0),(0,1,0)]
        self.assertEqual(fn(points,[0,1,2,3,4]),[0,1,2,3,4])
