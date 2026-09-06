"""Analytic represented-set fixtures only; no paper asset admission."""
import importlib
import json
import unittest
import numpy as np
from test_exact_boundary_audit import V, F


def box(lo=(0, 0, 0), hi=(1, 1, 1), inward=False):
    v = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],
                  [0,0,1],[1,0,1],[1,1,1],[0,1,1]], dtype=float)
    f = np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],
                  [0,1,5],[0,5,4],[1,2,6],[1,6,5],
                  [2,3,7],[2,7,6],[3,0,4],[3,4,7]])
    return v * (np.array(hi)-lo) + lo, f[:, ::-1] if inward else f


def combine(*meshes):
    vertices, faces, offset = [], [], 0
    for v, f in meshes:
        vertices.append(v)
        faces.append(f + offset)
        offset += len(v)
    return np.vstack(vertices), np.vstack(faces)


def l_prism():
    # CCW polygon, area 5: the missing upper-right 2x2 square is exterior.
    p = np.array([[0,0],[3,0],[3,1],[1,1],[1,3],[0,3]], dtype=float)
    v = np.vstack([np.column_stack([p, np.zeros(6)]),
                   np.column_stack([p, np.ones(6)])])
    triangles = [[0,1,3],[1,2,3],[0,3,5],[3,4,5]]
    f = [t[::-1] for t in triangles] + [[i+6 for i in t] for t in triangles]
    for i in range(6):
        j = (i+1) % 6
        f.extend([[i,j,j+6],[i,j+6,i+6]])
    return v, np.array(f)


class ExactMaterialTests(unittest.TestCase):
    def setUp(self):
        self.n = importlib.import_module('_make_it_stand')
        self.assertTrue(hasattr(self.n, 'make_exact_material_region'),
                        'common exact material native API is missing')

    def region(self, mesh):
        return self.n.make_exact_material_region(*mesh)

    def valid(self, mesh, parents, depths, signs, material, void):
        r = self.region(mesh)
        f = r.facts
        self.assertTrue(f['completed'], f)
        self.assertTrue(f['valid'], f)
        self.assertFalse(f['explicit_empty'])
        self.assertEqual(f['parents'], parents)
        self.assertEqual(f['depths'], depths)
        self.assertEqual(f['original_signs'], signs)
        self.assertEqual(f['checker_reversed_components'], [i for i,s in enumerate(signs) if s < 0])
        self.assertEqual((f['material_components'], f['bounded_void_components']), (material, void))
        json.dumps(f, allow_nan=False)
        return r

    def relation(self, a, b, expected):
        r = self.n.relate_exact_material_regions(a, b)
        self.assertTrue(r['completed'], r)
        self.assertEqual(r['error'], '')
        self.assertEqual(tuple(r[k] for k in ('intersection_nonempty', 'positive_volume_overlap',
                         'a_subset_b', 'b_subset_a', 'equal')), expected)
        self.assertEqual(r['a_digest'], a.facts['region_digest'])
        self.assertEqual(r['b_digest'], b.facts['region_digest'])
        json.dumps(r, allow_nan=False)
        return r

    def test_positive_tetra_and_reversed_root(self):
        self.valid((V,F), [-1], [0], [1], 1, 0)
        bad = self.region((V,F[:,::-1]))
        self.assertFalse(bad.facts['valid'])
        self.assertIn('orientation', bad.facts['error'])

    def test_topology_hollow_disjoint_and_nested_island_void(self):
        outer = box((0,0,0),(8,8,8))
        hole = box((1,1,1),(7,7,7), True)
        island = box((2,2,2),(6,6,6))
        inner_void = box((3,3,3),(5,5,5), True)
        self.valid(combine(outer,hole), [-1,0], [0,1], [1,-1], 1,1)
        self.valid(combine(outer,hole,island), [-1,0,1], [0,1,2], [1,-1,1], 2,1)
        self.valid(combine(outer,hole,island,inner_void), [-1,0,1,2], [0,1,2,3], [1,-1,1,-1], 2,2)
        self.valid(combine(box(),box((2,0,0),(3,1,1))), [-1,-1], [0,0], [1,1], 2,0)

    def test_boundary_valid_does_not_authenticate_region(self):
        cases = [combine(box((0,0,0),(4,4,4)),box((1,1,1),(2,2,2))),
                 combine(box((0,0,0),(4,4,4),True),box((1,1,1),(2,2,2),True)),
                 combine(box(),box((2,0,0),(3,1,1),True)),
                 combine(box((0,0,0),(8,8,8)),box((1,1,1),(7,7,7),True),
                         box((2,2,2),(6,6,6),True))]
        for mesh in cases:
            self.assertTrue(self.n.audit_exact_triangle_boundary(*mesh,'AnyNonzero')['boundary_valid'])
            bad = self.region(mesh)
            self.assertFalse(bad.facts['valid'])
            self.assertTrue(bad.facts['error'])
            with self.assertRaises(ValueError):
                self.n.relate_exact_material_regions(bad, self.region(box()))

    def test_nonconvex_exact_set_not_bbox_or_centroid(self):
        l = self.valid(l_prism(), [-1], [0], [1], 1,0)
        missing = box((1.5,1.5,.25),(2.5,2.5,.75))
        self.relation(l,self.region(missing),(False,False,False,False,False))
        # A negative shell in the notch is an invalid root despite bbox inclusion.
        v,f = missing
        self.assertFalse(self.region(combine(l_prism(),(v,f[:,::-1]))).facts['valid'])
        self.relation(l,self.region(box((.25,1.5,.25),(.75,2.5,.75))),
                      (True,True,False,True,False))

    def test_relation_contact_volume_containment_and_equality(self):
        a = self.region(box())
        cases = [((2,0,0),(3,1,1),(False,False,False,False,False)),
                 ((0,0,0),(1,1,1),(True,True,True,True,True)),
                 ((.25,.25,.25),(.75,.75,.75),(True,True,False,True,False)),
                 ((.5,0,0),(1.5,1,1),(True,True,False,False,False)),
                 ((1,0,0),(2,1,1),(True,False,False,False,False)),
                 ((1,1,0),(2,2,1),(True,False,False,False,False)),
                 ((1,1,1),(2,2,2),(True,False,False,False,False)),
                 ((0,0,0),(.5,1,1),(True,True,False,True,False)),
                 ((1-2**-40,0,0),(2,1,1),(True,True,False,False,False))]
        for lo, hi, expected in cases:
            self.relation(a,self.region(box(lo,hi)),expected)

    def test_disjoint_multicomponent_bbox_and_cavity_interface(self):
        a = self.region(combine(box(),box((4,0,0),(5,1,1))))
        b = self.region(combine(box((2,0,0),(3,1,1)),box((6,0,0),(7,1,1))))
        self.relation(a,b,(False,False,False,False,False))
        hollow = self.region(combine(box((0,0,0),(4,4,4)),box((1,1,1),(3,3,3),True)))
        self.relation(hollow,self.region(box((1,1,1),(3,3,3))),
                      (True,False,False,False,False))

    def test_empty_is_explicit_and_algebraic_only(self):
        e = self.n.empty_exact_material_region()
        self.assertTrue(e.facts['valid'])
        self.assertTrue(e.facts['explicit_empty'])
        self.assertEqual((e.facts['material_components'],e.facts['bounded_void_components']),(0,0))
        a = self.region(box())
        self.relation(e,a,(False,False,True,False,False))
        self.relation(a,e,(False,False,False,True,False))
        self.relation(e,self.n.empty_exact_material_region(),(False,False,True,True,True))
        with self.assertRaises(ValueError):
            self.region((np.empty((0,3)),np.empty((0,3),dtype=int)))
        with self.assertRaises(TypeError):
            self.n.ExactMaterialRegion()

    def test_identity_order_immutability_and_topology_reordering(self):
        outer, hole = box((0,0,0),(4,4,4)),box((1,1,1),(3,3,3),True)
        mesh = combine(outer,hole)
        saved = tuple(x.copy() for x in mesh)
        a = self.valid(mesh,[-1,0],[0,1],[1,-1],1,1)
        b = self.valid(combine(hole,outer),[1,-1],[1,0],[-1,1],1,1)
        self.assertNotEqual(a.facts['region_digest'],b.facts['region_digest'])
        self.assertEqual(a.facts,self.region(mesh).facts)
        self.assertEqual(a.facts['input_digest'],self.n.audit_exact_triangle_boundary(*mesh,'AnyNonzero')['input_digest'])
        for old, current in zip(saved,mesh):
            np.testing.assert_array_equal(old,current)
        facts = a.facts
        facts['valid'] = False
        self.assertTrue(a.facts['valid'])
        with self.assertRaises(AttributeError):
            a.facts = facts
        ab = self.relation(a,b,(True,True,True,True,True))
        self.assertEqual(ab,self.n.relate_exact_material_regions(a,b))
        ba = self.relation(b,a,(True,True,True,True,True))
        self.assertNotEqual(ab['relation_digest'],ba['relation_digest'])
        with self.assertRaises(TypeError):
            self.n.relate_exact_material_regions(a.facts,b)
        mesh[0][:] = 123
        self.assertEqual(ab,self.n.relate_exact_material_regions(a,b))

    def test_invalid_boundaries_and_lossless_conversion(self):
        bad_meshes = [(V,F[:-1]),combine(box(),box()),
                      combine(box(),box((.5,0,0),(1.5,1,1))),
                      combine(box(),box((1,0,0),(2,1,1)))]
        for mesh in bad_meshes:
            r = self.region(mesh)
            self.assertFalse(r.facts['valid'])
            self.assertTrue(r.facts['error'])
            with self.assertRaises(ValueError):
                self.n.relate_exact_material_regions(r,self.n.empty_exact_material_region())
        for mesh in [(V,F+99),(V,F.astype(float)),(V,np.array([[0,0,1]])),
                     (V,np.array([[0,1,2**63-1]])),(np.full_like(V,np.nan),F),
                     (np.full_like(V,np.inf),F),(V[:,:2],F)]:
            with self.assertRaises((ValueError,TypeError)):
                self.region(mesh)
        integer = V.astype(np.int64)
        integer[1,0] = 2**53+1
        with self.assertRaises(ValueError):
            self.region((integer,F))
        if np.dtype(np.longdouble).itemsize > 8:
            with self.assertRaises(TypeError):
                self.region((V.astype(np.longdouble),F))
        self.assertTrue(self.region((V.astype(np.float32),F)).facts['valid'])

    def test_incomplete_boundary_processing_keeps_counts_unknown(self):
        degenerate = self.region((np.array([[0.,0,0],[1,1,1],[2,2,2]]),np.array([[0,1,2]])))
        self.assertFalse(degenerate.facts['completed'])
        self.assertFalse(degenerate.facts['valid'])
        self.assertIsNone(degenerate.facts['material_components'])
        self.assertIsNone(degenerate.facts['bounded_void_components'])
        self.assertIn('degenerate_faces=1',degenerate.facts['error'])
        with self.assertRaises(ValueError):
            self.n.relate_exact_material_regions(degenerate,self.region(box()))

    def test_exact_translation_scale_and_thin_overlap(self):
        for scale, offset in [(1.,0.),(8.,32.),(.5,2**20)]:
            va,fa = box()
            vb,fb = box((1-2**-20,0,0),(2,1,1))
            a = self.region((va*scale+offset,fa))
            b = self.region((vb*scale+offset,fb))
            self.relation(a,b,(True,True,False,False,False))


if __name__ == '__main__':
    unittest.main()
