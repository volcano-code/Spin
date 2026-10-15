"""Reuse requires exact unchanged triangles and a native immutable valid proof."""
from copy import deepcopy
import numpy as np
import pytest
from make_it_stand import spin

V = [['0/1','0/1','0/1'], ['1/1','0/1','0/1'], ['0/1','1/1','0/1'], ['0/1','0/1','1/1']]
F = np.asarray([[1,2,3],[0,3,2],[0,1,3],[0,2,1]], dtype=np.int32)


def api():
    native = spin._native_module()
    factory = getattr(native, 'verify_fixed_rational_boundary', None)
    assert callable(factory), 'immutable verified fixed boundary API is missing'
    checker = getattr(native, 'audit_rational_boundary_with_fixed_proof', None)
    assert callable(checker), 'exact fixed-face reuse audit API is missing'
    return native, factory, checker


def same_geometric_facts(first, second):
    excluded = {'algorithm','backend','exact_pair_tests','reused_fixed_pair_tests','fixed_proof_identity'}
    assert {k:v for k,v in first.items() if k not in excluded} == {
        k:v for k,v in second.items() if k not in excluded}


def test_exact_unchanged_faces_reuse_proof_under_vertex_renumbering():
    native, factory, checker = api()
    proof = factory(V, F)
    permutation = [2,0,3,1]; inverse = np.argsort(permutation)
    vertices = [V[i] for i in permutation]; faces = inverse[F]
    before = native.audit_exact_rational_triangle_boundary(vertices, faces, 'OutwardPositive')
    after = checker(vertices, faces, 'OutwardPositive', proof)
    same_geometric_facts(before, after)
    assert after['boundary_valid'] and after['reused_fixed_pair_tests'] > 0
    assert after['exact_pair_tests'] < before['exact_pair_tests']


def test_new_faces_intersecting_original_faces_remain_rejected():
    native, factory, checker = api()
    proof = factory(V, F)
    shifted = [['1/4','1/4','1/4'], ['5/4','1/4','1/4'],
               ['1/4','5/4','1/4'], ['1/4','1/4','5/4']]
    vertices = V + shifted; faces = np.vstack([F, F+4])
    before = native.audit_exact_rational_triangle_boundary(vertices, faces, 'OutwardPositive')
    after = checker(vertices, faces, 'OutwardPositive', proof)
    same_geometric_facts(before, after)
    assert not after['boundary_valid'] and not after['self_intersection_free']


def test_coincident_vertices_with_different_indices_cannot_borrow_original_adjacency():
    native, factory, checker = api()
    proof = factory(V, F)
    vertices = V + deepcopy(V); faces = np.vstack([F, F+4])
    before = native.audit_exact_rational_triangle_boundary(vertices, faces, 'OutwardPositive')
    after = checker(vertices, faces, 'OutwardPositive', proof)
    same_geometric_facts(before, after)
    assert not after['boundary_valid'] and after['reused_fixed_pair_tests'] == 0


def test_duplicate_faces_are_not_certified_by_a_single_original_face():
    native, factory, checker = api()
    proof = factory(V, F)
    faces = np.vstack([F, F[0]])
    before = native.audit_exact_rational_triangle_boundary(V, faces, 'OutwardPositive')
    after = checker(V, faces, 'OutwardPositive', proof)
    same_geometric_facts(before, after)
    assert not after['boundary_valid']


def test_invalid_reference_is_rejected_and_proof_does_not_alias_input():
    native, factory, checker = api()
    with pytest.raises(ValueError): factory(V, F[:3])
    vertices = deepcopy(V); faces = F.copy()
    proof = factory(vertices, faces)
    vertices[0][0] = '100/1'; faces[0] = [0,1,2]
    before = native.audit_exact_rational_triangle_boundary(V, F, 'OutwardPositive')
    after = checker(V, F, 'OutwardPositive', proof)
    same_geometric_facts(before, after)
    assert after['boundary_valid'] and after['reused_fixed_pair_tests'] > 0
