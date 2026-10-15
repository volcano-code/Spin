"""Waterplane integrals may not depend on Boolean triangle enumeration."""
from fractions import Fraction as F
import unittest
import numpy as np
from make_it_stand.spin import compute_hydrostatic_geometry
from test_spin_it_faster import box


def exact_horizontal_cap(geometry, origin):
    """Independent rational triangle oracle over the actual returned boundary."""
    area, lx, ly, xx, yy, xy = (F(0) for _ in range(6))
    for indices in geometry.submerged_faces:
        points = [geometry.submerged_vertices[i] for i in indices]
        if not all(abs(p[2]-origin[2]) < 1e-9 for p in points):
            continue
        (a,b),(c,d),(e,f) = [(F(p[0]-origin[0]),F(p[1]-origin[1])) for p in points]
        weight=((c-a)*(f-b)-(d-b)*(e-a))/2
        assert weight > 0
        area+=weight; lx+=weight*(a+c+e)/3; ly+=weight*(b+d+f)/3
        xx+=weight*(a*a+c*c+e*e+a*c+c*e+e*a)/6
        yy+=weight*(b*b+d*d+f*f+b*d+d*f+f*b)/6
        xy+=weight*(2*(a*b+c*d+e*f)+a*d+c*b+c*f+e*d+e*b+a*f)/12
    return float(area), (float(lx/area),float(ly/area)), np.array([
        [float(xx-lx*lx/area),float(xy-lx*ly/area)],
        [float(xy-lx*ly/area),float(yy-ly*ly/area)]])


class HydrostaticDeterminismTests(unittest.TestCase):
    def test_displaced_volume_and_buoyancy_center_ignore_triangle_enumeration(self):
        vertices,faces=box((-1.,-8.,1.),(7.,0.,9.))
        z=7.399999999347529
        expected_volume=float(F(64)*(F(z)-1))
        expected_center=(3.,-4.,float((F(z)+1)/2))
        for repeat in range(96):
            ordered=np.roll(faces,repeat % len(faces),axis=0).tolist()
            g=compute_hydrostatic_geometry(vertices,ordered,waterline_origin=(3.,-4.,z),
                                          waterline_normal=(0,0,1))
            self.assertEqual(g.displaced_volume,expected_volume)
            self.assertEqual(g.center_of_buoyancy,expected_center)

    def test_shifted_square_cap_has_exact_zero_cross_moment_for_all_face_orders(self):
        vertices, faces=box((-1.,-8.,1.),(7.,0.,9.))
        origin=(3.,-4.,7.399999999347529)
        signatures=set()
        for repeat in range(24):
            # Cyclic face permutations preserve the same oriented input entity.
            ordered=np.roll(np.asarray(faces),repeat % len(faces),axis=0).tolist()
            g=compute_hydrostatic_geometry(vertices,ordered,waterline_origin=origin,
                                           waterline_normal=(0.,0.,1.))
            area,center,matrix=exact_horizontal_cap(g,origin)
            self.assertEqual(matrix[0,1],0.)
            self.assertEqual(g.waterline_central_second_moment[0][1],0.)
            self.assertEqual(g.waterline_area,area)
            np.testing.assert_allclose(g.waterline_central_second_moment,matrix,rtol=1e-15,atol=0)
            signatures.add((g.waterline_area,g.waterline_centroid,g.waterline_central_second_moment))
        self.assertEqual(len(signatures),1)

    def test_nonzero_cross_moment_is_integrated_not_snapped_to_zero(self):
        vertices,faces=box((-1.,-1.,-1.),(1.,1.,1.))
        v=np.asarray(vertices);v[:,0]+=0.375*v[:,1];v+=np.array([3.,-4.,5.])
        origin=(3.,-4.,5.3)
        for repeat in range(4):
            g=compute_hydrostatic_geometry(v,np.roll(faces,repeat,axis=0).tolist(),
                                          waterline_origin=origin,waterline_normal=(0,0,1))
            area,center,matrix=exact_horizontal_cap(g,origin)
            self.assertAlmostEqual(matrix[0,1],0.5,places=14)
            np.testing.assert_allclose(g.waterline_central_second_moment,matrix,rtol=1e-14,atol=1e-15)
            self.assertEqual(g.waterline_area,area)
