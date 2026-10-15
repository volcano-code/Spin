"""Public exact checks with analytic fixtures and an independent Fraction oracle.

These tests catch volume/bbox substitutes, omitted nearest features, rounded
pruning, unstable witness ties, incomplete evidence and unauthenticated inputs.
"""
import importlib
import json
import unittest
from fractions import Fraction as Q
import numpy as np
from test_exact_material_region import box, combine, l_prism
from test_exact_boundary_audit import F


def sub(a, b):
    return tuple(x-y for x,y in zip(a,b))


def add(a, b):
    return tuple(x+y for x,y in zip(a,b))


def mul(a, t):
    return tuple(x*t for x in a)


def dot(a, b):
    return sum(x*y for x,y in zip(a,b))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def point_segment(p, a, b):
    u = sub(b,a)
    t = max(Q(0), min(Q(1), dot(sub(p,a),u)/dot(u,u)))
    d = sub(p,add(a,mul(u,t)))
    return dot(d,d)


def in_triangle(p, t):
    a,b,c = t
    n = cross(sub(b,a),sub(c,a))
    return all(dot(cross(sub(y,x),sub(p,x)),n) >= 0
               for x,y in ((a,b),(b,c),(c,a)))


def point_triangle(p, t):
    a,b,c = t
    n = cross(sub(b,a),sub(c,a))
    h = dot(sub(p,a),n)
    projection = sub(p,mul(n,h/dot(n,n)))
    if in_triangle(projection,t):
        return h*h/dot(n,n)
    return min(point_segment(p,x,y) for x,y in ((a,b),(b,c),(c,a)))


def segment_segment(a,b,c,d):
    u,v,w = sub(b,a),sub(d,c),sub(a,c)
    aa,bb,cc,dd,ee = dot(u,u),dot(u,v),dot(v,v),dot(u,w),dot(v,w)
    candidates = [point_segment(a,c,d),point_segment(b,c,d),
                  point_segment(c,a,b),point_segment(d,a,b)]
    determinant = aa*cc-bb*bb
    if determinant:
        s,t = (bb*ee-cc*dd)/determinant,(aa*ee-bb*dd)/determinant
        if 0 <= s <= 1 and 0 <= t <= 1:
            delta = sub(add(w,mul(u,s)),mul(v,t))
            candidates.append(dot(delta,delta))
    return min(candidates)


def triangle_distance(a,b):
    # Independent rational plane crossing detects edge/face intersections.
    # Coplanar intersections are covered by vertex/face and edge/edge minima.
    for t,other in ((a,b),(b,a)):
        n = cross(sub(other[1],other[0]),sub(other[2],other[0]))
        for i in range(3):
            p,q = t[i],t[(i+1)%3]
            hp,hq = dot(sub(p,other[0]),n),dot(sub(q,other[0]),n)
            if hp != hq:
                s = hp/(hp-hq)
                if 0 <= s <= 1 and in_triangle(add(p,mul(sub(q,p),s)),other):
                    return Q(0)
    return min([point_triangle(p,b) for p in a] + [point_triangle(p,a) for p in b] +
               [segment_segment(a[i],a[(i+1)%3],b[j],b[(j+1)%3])
                for i in range(3) for j in range(3)])


def exhaustive(mesh_a, mesh_b):
    meshes = []
    for v,f in (mesh_a,mesh_b):
        points = [tuple(Q(float(x)) for x in p) for p in v]
        meshes.append([tuple(points[int(i)] for i in face) for face in f])
    return min((triangle_distance(a,b), i,j)
               for i,a in enumerate(meshes[0]) for j,b in enumerate(meshes[1]))


class ExactPartitionClearanceTests(unittest.TestCase):
    def setUp(self):
        self.n = importlib.import_module('_make_it_stand')
        missing = [api for api in ('check_exact_material_partition','check_exact_boundary_clearance')
                   if not hasattr(self.n,api)]
        self.assertEqual(missing,[], f'common native APIs missing: {missing}')

    def region(self, mesh):
        r = self.n.make_exact_material_region(*mesh)
        self.assertTrue(r.facts['valid'],r.facts)
        return r

    def partition(self, whole, parts, expected):
        r = self.n.check_exact_material_partition(whole,parts)
        self.assertTrue(r['completed'],r)
        self.assertEqual(r['error'],'')
        self.assertEqual(tuple(r[k] for k in ('union_equals_whole','no_positive_volume_overlap',
                         'missing_whole_content','excess_union_content','complete_partition')),expected)
        self.assertEqual(r['whole_digest'],whole.facts['region_digest'])
        self.assertEqual(r['part_digests'],[p.facts['region_digest'] for p in parts])
        self.assertEqual(r['whole_explicit_empty'],whole.facts['explicit_empty'])
        self.assertEqual(r['part_explicit_empty'],[p.facts['explicit_empty'] for p in parts])
        self.assertTrue(r['algorithm'] and r['backend'] and r['semantics'] and r['empty_protocol'])
        self.assertEqual(r['scope'],'exact-represented-polyhedral-partition-only')
        self.assertEqual(len(r['partition_digest']),64)
        json.dumps(r,allow_nan=False)
        return r

    def test_partition_cover_gap_excess_overlap_duplicate_and_subset_only(self):
        whole = self.region(box())
        left = self.region(box(hi=(.5,1,1)))
        right = self.region(box(lo=(.5,0,0)))
        r = self.partition(whole,[left,right],(True,True,False,False,True))
        self.assertEqual(r['part_subset_whole'],[True,True])
        self.assertEqual(r['positive_volume_overlap_pairs'],[])
        reordered = self.partition(whole,[right,left],(True,True,False,False,True))
        self.assertNotEqual(r['partition_digest'],reordered['partition_digest'])
        larger = self.partition(self.region(box(hi=(2,1,1))),[left,right],(False,True,True,False,False))
        self.assertNotEqual(r['partition_digest'],larger['partition_digest'])
        for parts,expected,subsets,pairs in [
            ([left],(False,True,True,False,False),[True],[]),
            ([left,self.region(box(lo=(.75,0,0)))],(False,True,True,False,False),[True,True],[]),
            ([left,self.region(box(lo=(.5,0,0),hi=(1.25,1,1)))],
             (False,True,False,True,False),[True,False],[]),
            ([left,self.region(box(lo=(.25,0,0)))],(True,False,False,False,False),[True,True],[[0,1]]),
            ([whole,whole],(True,False,False,False,False),[True,True],[[0,1]])]:
            result = self.partition(whole,parts,expected)
            self.assertEqual(result['part_subset_whole'],subsets)
            self.assertEqual(result['positive_volume_overlap_pairs'],pairs)

    def test_partition_contact_codimension_and_underflow_scale_overlap(self):
        whole = self.region(box(hi=(3,3,3)))
        first = self.region(box())
        for lo,hi in [((1,0,0),(2,1,1)),((1,1,0),(2,2,1)),((1,1,1),(2,2,2))]:
            r = self.partition(whole,[first,self.region(box(lo,hi))],(False,True,True,False,False))
            self.assertEqual(r['positive_volume_overlap_pairs'],[])
        # Exact missing sliver and overlapping sliver next to zero, with a gap
        # much smaller than the smallest normal binary64 squared length.
        tiny = 2.**-600
        outer = self.region(box((-1,0,0),(1,1,1)))
        left = self.region(box((-1,0,0),(0,1,1)))
        self.partition(outer,[left,self.region(box((tiny,0,0),(1,1,1)))],(False,True,True,False,False))
        self.partition(outer,[left,self.region(box((-tiny,0,0),(1,1,1)))],(True,False,False,False,False))
        r = self.partition(first,[first,first,first],(True,False,False,False,False))
        self.assertEqual(r['positive_volume_overlap_pairs'],[[0,1],[0,2],[1,2]])

    def test_partition_hollow_nonconvex_nested_island_multicomponent(self):
        outer = box(hi=(4,4,4))
        cavity = box((1,1,1),(3,3,3))
        v,f = cavity
        hollow = self.region(combine(outer,(v,f[:,::-1])))
        self.partition(self.region(outer),[hollow,self.region(cavity)],(True,True,False,False,True))
        self.partition(self.region(l_prism()),[self.region(box(hi=(3,1,1))),
                       self.region(box((0,1,0),(1,3,1)))],(True,True,False,False,True))
        island = box((1.5,1.5,1.5),(2.5,2.5,2.5))
        remote = box((5,0,0),(6,1,1))
        nested = self.region(combine(outer,(v,f[:,::-1]),island,remote))
        self.partition(nested,[hollow,self.region(island),self.region(remote)],(True,True,False,False,True))

    def test_partition_explicit_empty_missing_evidence_and_authentication(self):
        a,e = self.region(box()),self.n.empty_exact_material_region()
        self.partition(a,[e,a,e],(True,True,False,False,True))
        self.partition(e,[e,e],(True,True,False,False,True))
        r = self.partition(e,[a],(False,True,False,True,False))
        self.assertEqual(r['part_subset_whole'],[False])
        self.partition(a,[e],(False,True,True,False,False))
        for whole,parts in ((a,[]),(e,[])):
            with self.assertRaises(ValueError):
                self.n.check_exact_material_partition(whole,parts)
        invalid = self.n.make_exact_material_region(*box(inward=True))
        for whole,parts in ((invalid,[a]),(a,[invalid])):
            with self.assertRaises(ValueError):
                self.n.check_exact_material_partition(whole,parts)
        for whole,parts in ((a.facts,[a]),(a,[a.facts]),(None,[a])):
            with self.assertRaises(TypeError):
                self.n.check_exact_material_partition(whole,parts)

    def distance(self, ma, mb, threshold, squared, oracle=False):
        a,b = self.region(ma),self.region(mb)
        r = self.n.check_exact_boundary_clearance(a,b,threshold)
        self.assertTrue(r['completed'],r)
        self.assertTrue(r['applicable'])
        self.assertEqual(r['state'],'COMPLETE')
        self.assertEqual(r['error'],'')
        q = Q(squared)
        self.assertEqual(Q(r['squared_distance_exact']),q)
        self.assertEqual(Q(r['minimum_clearance_exact']),Q(float(threshold)))
        self.assertEqual(r['at_least_minimum'],q >= Q(float(threshold))**2)
        self.assertEqual(r['contact'],q == 0)
        self.assertEqual((r['a_digest'],r['b_digest']),(a.facts['region_digest'],b.facts['region_digest']))
        self.assertEqual((r['a_input_digest'],r['b_input_digest']),(a.facts['input_digest'],b.facts['input_digest']))
        self.assertEqual((r['a_explicit_empty'],r['b_explicit_empty']),(False,False))
        self.assertTrue(r['algorithm'] and r['backend'] and r['semantics'])
        self.assertEqual(r['scope'],'exact-represented-boundary-distance-only')
        self.assertEqual(len(r['clearance_digest']),64)
        i,j = r['witness_faces']
        self.assertTrue(0 <= i < len(ma[1]) and 0 <= j < len(mb[1]))
        if oracle:
            self.assertEqual((q,i,j),exhaustive(ma,mb))
        json.dumps(r,allow_nan=False)
        return r

    def test_clearance_boxes_nested_not_volume_distance_or_containment(self):
        self.distance(box(),box((3,0,0),(4,1,1)),2,4,True)
        self.distance(box(hi=(4,4,4)),box((1,1,1),(3,3,3)),1,1,True)
        hollow = combine(box(hi=(4,4,4)),box((1,1,1),(3,3,3),True))
        self.distance(hollow,box((1.5,1.5,1.5),(2.5,2.5,2.5)),.5,Q(1,4))
        self.distance(l_prism(),box((1.5,1.5,.25),(2.5,2.5,.75)),.5,Q(1,4),True)
        outside = self.distance(box(),box((2,2,2),(3,3,3)),1,3)
        for unsupported in ('contained','valid_shell','sealed','thickness_passed','feasible'):
            self.assertNotIn(unsupported,outside)

    def test_clearance_face_edge_vertex_contact_and_crossing(self):
        for lo,hi in [((1,0,0),(2,1,1)),((1,1,0),(2,2,1)),((1,1,1),(2,2,2)),
                      ((.5,.5,.5),(1.5,1.5,1.5)),((1-2**-40,0,0),(2,1,1))]:
            r = self.distance(box(),box(lo,hi),0,0,True)
            self.assertTrue(r['at_least_minimum'] and r['contact'])

    def test_clearance_dyadic_sub_epsilon_gap_and_threshold(self):
        # Adjacent to zero: the squared value underflows binary64, not rationals.
        gap = 2.**-600
        a = box((-1,0,0),(0,1,1))
        b = box((gap,0,0),(1,1,1))
        hashes = []
        for threshold in (gap/2,gap,gap*2):
            r = self.distance(a,b,threshold,Q(1,2**1200))
            hashes.append(r['clearance_digest'])
        self.assertEqual(len(set(hashes)),3)

    def test_clearance_asymmetric_tetra_nearest_features_and_covariance(self):
        # A top edge runs x at z=0; B bottom edge runs y at z=2.
        # Their interior projections cross; every vertex pair is farther away.
        a = (np.array([[-2,0,0],[3,0,0],[0,1,-2],[0,-2,-3.]]),F[:,::-1])
        b = (np.array([[0,-3,2],[0,2,2],[1,0,4],[-2,0,5.]]),F[:,::-1])
        self.distance(a,b,2,4,True)
        # Nearest upper vertex projects inside lower tetra's horizontal face.
        c = (np.array([[0,0,0],[4,0,0],[0,3,0],[0,0,-2.]]),F[:,::-1])
        d = (np.array([[1,1,2],[2,1,4],[1,2,5],[0,0,6.]]),F)
        self.distance(c,d,2,4,True)
        for scale,offset in ((8.,2.**40),(.5,32.)):
            av,af = a; bv,bf = b
            self.distance((av*scale+offset,af),(bv*scale+offset,bf),2*scale,Q(4)*Q(scale)**2,True)

    def test_clearance_acceleration_nonaxis_oracle_and_ties(self):
        ma,mb = box(),box((2,0,0),(3,1,1))
        # Integral shear keeps the fixture exact and its surfaces non-axis-aligned.
        shear = np.array([[1,1,0],[0,1,1],[0,0,1.]])
        a,b = (ma[0]@shear,ma[1]),(mb[0]@shear,mb[1])
        q,i,j = exhaustive(a,b)
        r = self.distance(a,b,0,q,True)
        self.assertEqual(r['witness_faces'],[i,j])
        # Face order changes indices: the lexicographically first minimizer must
        # still match exhaustive enumeration, independent of BVH traversal.
        self.distance((a[0],a[1][::-1]),(b[0],b[1][::-1]),0,q,True)

    def test_clearance_empty_and_invalid_thresholds_operands(self):
        a,e = self.region(box()),self.n.empty_exact_material_region()
        for x,y in ((a,e),(e,a),(e,e)):
            r = self.n.check_exact_boundary_clearance(x,y,1.)
            self.assertTrue(r['completed'])
            self.assertFalse(r['applicable'])
            self.assertEqual(r['state'],'NOT_APPLICABLE')
            for k in ('squared_distance_exact','witness_faces','at_least_minimum','contact'):
                self.assertIsNone(r[k])
            json.dumps(r,allow_nan=False)
        for bad in (True,np.bool_(False),-1.,float('nan'),float('inf'),2**53+1,
                    np.uint64(2**64-1),np.longdouble(1),[1.],None,'1',Q(1,3)):
            with self.subTest(threshold=bad), self.assertRaises((ValueError,TypeError,OverflowError)):
                self.n.check_exact_boundary_clearance(a,a,bad)
        for good in (0,1,np.float32(.5),np.int64(2),-0.):
            self.assertTrue(self.n.check_exact_boundary_clearance(a,a,good)['completed'])
        invalid = self.n.make_exact_material_region(*box(inward=True))
        for x,y in ((invalid,a),(a,invalid)):
            with self.assertRaises(ValueError):
                self.n.check_exact_boundary_clearance(x,y,0)
        for x,y in ((a.facts,a),(None,a)):
            with self.assertRaises(TypeError):
                self.n.check_exact_boundary_clearance(x,y,0)
        with self.assertRaises(TypeError):
            self.n.check_exact_boundary_clearance(a,a)

    def test_owned_operands_and_repeated_identities(self):
        ma,mb = box(),box((2,0,0),(3,1,1))
        a,b = self.region(ma),self.region(mb)
        before = tuple(x.copy() for x in ma)
        r = self.n.check_exact_boundary_clearance(a,b,.5)
        p = self.n.check_exact_material_partition(a,[a])
        for x,y in zip(before,ma):
            np.testing.assert_array_equal(x,y)
        ma[0][:] = 42; ma[1][:] = 0
        a.facts['valid'] = False
        r['contact'] = True; p['complete_partition'] = False
        fresh = self.n.check_exact_boundary_clearance(a,b,.5)
        self.assertFalse(fresh['contact'])
        self.assertEqual(fresh,self.n.check_exact_boundary_clearance(a,b,.5))
        self.assertTrue(self.n.check_exact_material_partition(a,[a])['complete_partition'])
        self.assertEqual(p['partition_digest'],self.n.check_exact_material_partition(a,[a])['partition_digest'])
        self.assertNotEqual(fresh['clearance_digest'],self.n.check_exact_boundary_clearance(b,a,.5)['clearance_digest'])


    def test_rational_coordinate_clearance_uses_exact_constructed_faces(self):
        unit = box()
        faces = unit[1]
        a_vertices = [[('1/3' if value else '0/1') for value in row]
                      for row in unit[0].astype(int)]
        b_vertices = [[('1/1' if axis == 0 and value else
                        '2/3' if axis == 0 else
                        '1/3' if value else '0/1')
                       for axis, value in enumerate(row)]
                      for row in unit[0].astype(int)]
        a = self.n.make_exact_material_region_from_rationals(a_vertices, faces)
        b = self.n.make_exact_material_region_from_rationals(b_vertices, faces)
        result = self.n.check_exact_boundary_clearance(a, b, 1/3)
        self.assertTrue(result['completed'], result)
        self.assertEqual(result['state'], 'COMPLETE')
        self.assertEqual(Q(result['squared_distance_exact']), Q(1, 9))
        self.assertTrue(result['at_least_minimum'])
        self.assertFalse(result['contact'])
        self.assertEqual(result['a_digest'], a.facts['region_digest'])
        self.assertEqual(result['b_digest'], b.facts['region_digest'])

        # Mixed input must convert the binary64 operand's exact represented
        # coordinates into the same rational kernel without rounding the 1/3.
        binary = self.region(box((1, 0, 0), (2, 1, 1)))
        mixed = self.n.check_exact_boundary_clearance(a, binary, 2/3)
        self.assertTrue(mixed['completed'], mixed)
        self.assertEqual(Q(mixed['squared_distance_exact']), Q(4, 9))
        self.assertTrue(mixed['at_least_minimum'])


if __name__ == '__main__':
    unittest.main()
