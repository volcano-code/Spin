"""An explicitly filtered exact build must identify its predicate backend."""
import numpy as np
from fractions import Fraction
from make_it_stand import spin


def test_filtered_build_declares_exact_predicate_backend():
    vertices=[['0/1','0/1','0/1'],['1/1','0/1','0/1'],['0/1','1/1','0/1'],['0/1','0/1','1/1']]
    faces=np.asarray([[0,2,1],[0,1,3],[0,3,2],[1,2,3]],dtype=np.int64)
    facts=spin._native_module().audit_exact_rational_triangle_boundary(vertices,faces,'AnyNonzero')
    assert facts['boundary_valid']
    assert 'Filtered_kernel<Simple_cartesian<Gmpq>>' in facts['backend'], \
        'explicit filtered candidate is not identified in the actual audit'


def test_uncertain_interval_coordinates_keep_exact_nonzero_volume():
    # All displayed z coordinates coincide at 1; the exact tetra is nevertheless
    # closed and nondegenerate. A double-only fallback or epsilon rejects it.
    epsilon=Fraction(1,2**1200)
    vertices=[['0/1','0/1','1/1'],['1/1','0/1','1/1'],['0/1','1/1','1/1'],
              ['0/1','0/1',str(1+epsilon)]]
    faces=np.asarray([[0,2,1],[0,1,3],[0,3,2],[1,2,3]],dtype=np.int64)
    facts=spin._native_module().audit_exact_rational_triangle_boundary(vertices,faces,'OutwardPositive')
    assert facts['boundary_valid'] and facts['degenerate_faces']==[]
    assert Fraction(facts['component_signed_volumes_exact'][0])==epsilon/6


def test_sub_binary64_shared_corner_overlap_is_not_excused():
    # Triangle intersection extends beyond the common indexed origin by the
    # exact amount d/(2+d)>0, although converting 1+d to binary64 gives 1.
    delta=Fraction(1,2**200)
    vertices=[['0/1','0/1','0/1'],['1/1','0/1','0/1'],['0/1','1/1','0/1'],
              ['1/1','1/1','-1/1'],['-1/1','-1/1',str(1+delta)]]
    faces=np.asarray([[0,1,2],[0,3,4]],dtype=np.int64)
    facts=spin._native_module().audit_exact_rational_triangle_boundary(vertices,faces,'AnyNonzero')
    assert facts['self_intersection_completed']
    assert facts['intersecting_pairs']==[[0,1]]
    assert not facts['boundary_valid']
