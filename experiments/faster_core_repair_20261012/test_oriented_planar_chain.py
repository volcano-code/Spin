import unittest
import exact_planar_patches as exact


class PlanarChainTests(unittest.TestCase):
    def signature(self,points,faces):
        fn=getattr(exact,'oriented_planar_chain_signature',None)
        self.assertTrue(callable(fn),'exact oriented planar-chain certificate missing')
        return fn(points,faces)['signature']

    def test_exact_diagonal_flip_same_oriented_surface(self):
        points=[(0,0,0),(1,0,0),(1,1,0),(0,1,0)]
        self.assertEqual(self.signature(points,[[0,1,2],[0,2,3]]),
                         self.signature(points,[[0,1,3],[1,2,3]]))

    def test_null_triangle_only_in_certificate_does_not_change_surface(self):
        points=[(0,0,0),(1,0,0),(1,1,0),(0,1,0),('1/2',0,0)]
        faces=[[0,1,2],[0,2,3]]
        self.assertEqual(self.signature(points,faces),self.signature(points,faces+[[0,4,1]]))

    def test_reversed_surface_is_not_equal(self):
        points=[(0,0,0),(1,0,0),(0,1,0)]
        self.assertNotEqual(self.signature(points,[[0,1,2]]),self.signature(points,[[0,2,1]]))

    def test_tiny_nonplanarity_is_not_erased(self):
        points=[(0,0,0),(1,0,0),(1,1,'1/'+str(10**60)),(0,1,0)]
        self.assertNotEqual(self.signature(points,[[0,1,2],[0,2,3]]),
                            self.signature(points,[[0,1,3],[1,2,3]]))
