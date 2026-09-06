"""Real native exact represented-boundary regressions; fixtures are not paper cases."""
import hashlib
import importlib
import json
import struct
import unittest
from fractions import Fraction

import numpy as np


V = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
F = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]], dtype=np.int64)


class ExactBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.native = importlib.import_module('_make_it_stand')

    def audit(self, v=V, f=F, orientation='OutwardPositive'):
        self.assertTrue(hasattr(self.native, 'audit_exact_triangle_boundary'),
                        'common exact native boundary API is missing')
        return self.native.audit_exact_triangle_boundary(v, f, orientation)

    def test_shared_corner_analytical_short_crossing(self):
        for delta in [1e-4, 1e-5, 1e-7, 1e-8]:
            v = np.array([[0., 0., 0.], [1, 0, 0], [0, 1, 0],
                          [1, 1, -1], [-1, -1, 1 + delta]])
            f = np.array([[0, 1, 2], [0, 3, 4]])
            # The exact represented segment ends at (d/(2+d),d/(2+d),0).
            d = Fraction(float(v[4, 2])) - 1
            self.assertGreater(d / (2 + d), 0)
            with self.subTest(delta=delta):
                out = self.audit(v, f)
                self.assertTrue(out['self_intersection_completed'])
                self.assertFalse(out['self_intersection_free'])
                self.assertEqual(out['intersecting_pairs'], [[0, 1]])
                self.assertFalse(out['boundary_valid'])
                self.assertFalse(out['closed'])

    def test_literal_volume_orientation_and_input_identity(self):
        v, f = V.copy(), F.copy()
        out = self.audit(v, f)
        self.assertTrue(out['boundary_valid'], out)
        self.assertEqual(out['component_volume_signs'], [1])
        self.assertEqual([Fraction(x) for x in out['component_signed_volumes_exact']], [1])
        self.assertEqual(out['face_components'], [0, 0, 0, 0])
        self.assertEqual(self.audit(v, f), out)
        np.testing.assert_array_equal(v, V)
        np.testing.assert_array_equal(f, F)
        # Independently specify canonical bytes: version, mode, dimensions,
        # row-major IEEE754 binary64 bit patterns and signed index values as u64 BE.
        canonical = b'exact-boundary-input-v1\0' + bytes([0])
        canonical += struct.pack('>QQ', 4, 4)
        canonical += b''.join(struct.pack('>d', x) for x in v.flat)
        canonical += b''.join(struct.pack('>Q', int(x)) for x in f.flat)
        self.assertEqual(out['input_digest'], hashlib.sha256(canonical).hexdigest())
        self.assertNotEqual(self.audit(orientation='AnyNonzero')['input_digest'], out['input_digest'])
        self.assertEqual(out['scope'], 'single-represented-boundary-only')
        self.assertEqual(out['material_domain_status'], 'NOT_VERIFIED')
        json.dumps(out, allow_nan=False)

    def test_winding_and_explicit_inward(self):
        rev = F[:, ::-1].copy()
        self.assertFalse(self.audit(V, rev)['boundary_valid'])
        inward = self.audit(V, rev, 'InwardNegative')
        self.assertTrue(inward['boundary_valid'])
        self.assertEqual(inward['component_volume_signs'], [-1])
        self.assertEqual([Fraction(x) for x in inward['component_signed_volumes_exact']], [-1])
        bad = F.copy()
        bad[0] = bad[0, ::-1]
        self.assertFalse(self.audit(V, bad)['consistently_oriented'])

    def test_indexed_simplex_contact_rules(self):
        cases = [
            # Legal shared point; legal coplanar opposite-side diagonal.
            ([[0,0,0],[1,0,0],[0,1,0],[-1,0,1],[0,-1,1]], [[0,1,2],[0,3,4]], []),
            ([[0,0,0],[1,0,0],[0,1,0],[1,1,0]], [[0,1,2],[1,3,2]], []),
            # Coplanar shared-edge positive area overlap.
            ([[0,0,0],[1,0,0],[0,1,0],[.2,.2,0]], [[0,1,2],[0,1,3]], [[0,1]]),
            # Distinct indices: coincident point, edge and area all illegal.
            ([[0,0,0],[1,0,0],[0,1,0],[0,0,0],[-1,0,0],[0,-1,0]], [[0,1,2],[3,4,5]], [[0,1]]),
            ([[0,0,0],[1,0,0],[0,1,0],[0,0,0],[1,0,0],[0,-1,0]], [[0,1,2],[3,4,5]], [[0,1]]),
            ([[0,0,0],[1,0,0],[0,1,0],[0,0,0],[1,0,0],[0,1,0]], [[0,1,2],[3,4,5]], [[0,1]]),
            # Six-vertex polygon overlap exercises CGAL polygon variant.
            ([[-2,-1,0],[2,-1,0],[0,2,0],[-2,1,0],[2,1,0],[0,-2,0]], [[0,1,2],[3,4,5]], [[0,1]]),
            # Shared indexed origin with four-vertex polygon overlap, so the
            # polygon must be tested against the indexed point, not excused.
            ([[0,0,0],[2,0,0],[0,2,0],[1,.25,0],[.25,3,0]], [[0,1,2],[0,3,4]], [[0,1]]),
            # All three indices shared is always duplicate geometry.
            ([[0,0,0],[1,0,0],[0,1,0]], [[0,1,2],[2,1,0]], [[0,1]]),
        ]
        for v, f, pairs in cases:
            with self.subTest(v=v):
                out = self.audit(np.array(v, dtype=float), np.array(f))
                self.assertTrue(out['self_intersection_completed'])
                self.assertEqual(out['intersecting_pairs'], pairs)

    def test_open_degenerate_unreferenced_and_vertex_link(self):
        self.assertFalse(self.audit(V, F[:-1])['closed'])
        collinear = self.audit(np.array([[0.,0,0],[1,1,1],[2,2,2]]), np.array([[0,1,2]]))
        self.assertEqual(collinear['degenerate_faces'], [0])
        self.assertFalse(collinear['self_intersection_completed'])
        self.assertFalse(collinear['boundary_valid'])
        self.assertEqual(self.audit(np.vstack([V, [9,9,9]]))['unreferenced_vertices'], [4])
        # Two tetrahedra meeting at one indexed vertex have 2 faces per edge,
        # but a disconnected vertex link (two cycles).
        v = np.vstack([V, -V[1:]])
        f2 = np.where(F == 0, 0, F + 3)[:, ::-1]
        out = self.audit(v, np.vstack([F, f2]))
        self.assertTrue(out['edge_manifold'])
        self.assertFalse(out['vertex_manifold'])
        self.assertEqual(out['non_manifold_vertices'], [0])
        self.assertTrue(out['self_intersection_free'])
        self.assertFalse(out['boundary_valid'])
        triple = self.audit(V, np.vstack([F, F[:1]]))
        self.assertFalse(triple['edge_manifold'])

    def test_multiple_nested_components_only_boundary_facts(self):
        inner = V / 8 + [1/8, 1/8, 1/8]
        v = np.vstack([V, inner])
        f = np.vstack([F, F[:, ::-1] + 4])
        out = self.audit(v, f, 'AnyNonzero')
        self.assertTrue(out['boundary_valid'])
        self.assertEqual(out['component_volume_signs'], [1, -1])
        self.assertEqual([Fraction(x) for x in out['component_signed_volumes_exact']], [1, -Fraction(1,512)])
        self.assertEqual(out['face_components'], [0]*4 + [1]*4)
        self.assertEqual(out['material_domain_status'], 'NOT_VERIFIED')
        self.assertFalse(self.audit(v, f)['component_orientations_valid'])

    def test_translated_scaled_large_offset_and_overflow_metadata(self):
        for scale, shift in [(2., 7.), (1., 2.**40), (2.**-200, 0.),
                             (2.**-1070, 0.), (2.**400, 0.)]:
            out = self.audit(V * scale + shift)
            self.assertTrue(out['boundary_valid'], out)
            self.assertEqual(Fraction(out['component_signed_volumes_exact'][0]), Fraction(scale)**3)
            json.dumps(out, allow_nan=False)

    def test_closed_crossing_components(self):
        out = self.audit(np.vstack([V, V + [.2,.2,.2]]), np.vstack([F,F+4]))
        self.assertTrue(out['closed'])
        self.assertFalse(out['self_intersection_free'])
        self.assertFalse(out['boundary_valid'])

    def test_single_closed_connected_self_crossing(self):
        # Inverted lower octahedron apex: the W->B edge crosses the upper
        # E->T edge at the literal midpoint (.5,0,.5), away from shared indices.
        v = np.array([[1.,0,0],[0,1,0],[-1,0,0],[0,-1,0],[0,0,1],[2,0,1]])
        f = np.array([[0,1,4],[1,2,4],[2,3,4],[3,0,4],
                      [1,0,5],[2,1,5],[3,2,5],[0,3,5]])
        out = self.audit(v, f, 'AnyNonzero')
        self.assertTrue(out['closed'])
        self.assertTrue(out['edge_manifold'])
        self.assertTrue(out['vertex_manifold'])
        self.assertEqual(set(out['face_components']), {0})
        self.assertIn([0, 5], out['intersecting_pairs'])
        self.assertFalse(out['boundary_valid'])

    def test_python_conversion_cannot_change_represented_coordinates(self):
        # Binary64 API must not silently round high precision or integer inputs.
        v = V.astype(np.int64)
        v[1,0] = 2**53 + 1
        with self.assertRaises(ValueError):
            self.audit(v)
        if np.dtype(np.longdouble).itemsize > 8:
            with self.assertRaises(TypeError):
                self.audit(V.astype(np.longdouble))
        self.assertTrue(self.audit(V.astype(np.float32))['boundary_valid'])
        self.assertTrue(self.audit(V.astype(np.int64))['boundary_valid'])
        for dtype, value in [(np.int64, 2**63-1), (np.uint64, 2**64-1), (np.uint64, 2**53+1)]:
            large=V.astype(dtype)
            large[1,0]=value
            with self.assertRaises(ValueError):
                self.audit(large)

    def test_digest_binds_row_order_coordinates_indices_and_signed_zero(self):
        original=self.audit()['input_digest']
        self.assertNotEqual(self.audit(V, F[::-1])['input_digest'], original)
        perm=np.array([1,0,2,3])
        self.assertNotEqual(self.audit(V[perm],perm[F])['input_digest'], original)
        changed=V.copy()
        changed[1,0]=np.nextafter(1.,2.)
        self.assertNotEqual(self.audit(changed)['input_digest'],original)
        changed=V.copy()
        changed[0,0]=-0.0
        self.assertNotEqual(self.audit(changed)['input_digest'],original)

    def test_malformed_inputs_fail_before_predicates(self):
        cases = [(V, F.astype(float)), (V, F.astype(bool)), (V, F + 99),
                 (V, F - 99), (V, np.array([[0,0,1]])), (V[:, :2], F),
                 (V, F[:, :2]), (np.empty((0,3)), np.empty((0,3),dtype=int)),
                 (np.full_like(V, float('nan')), F), (np.full_like(V, float('inf')), F)]
        for v, f in cases:
            with self.subTest(v=v, f=f):
                with self.assertRaises((ValueError, TypeError)):
                    self.audit(v, f)
        with self.assertRaises(ValueError):
            self.audit(orientation='caller-passed')


if __name__ == '__main__':
    unittest.main()
