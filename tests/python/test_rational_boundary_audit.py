"""Expose the existing exact rational checker without changing its predicates."""
from fractions import Fraction
import numpy as np
import pytest
from make_it_stand import spin


def api():
    native=spin._native_module()
    assert hasattr(native,'audit_exact_rational_triangle_boundary'), 'exact rational boundary audit binding missing'
    return native


def tetra():
    vertices=[[Fraction(0),Fraction(0),Fraction(0)], [Fraction(1),Fraction(0),Fraction(0)],
              [Fraction(0),Fraction(1),Fraction(0)], [Fraction(0),Fraction(0),Fraction(1)]]
    faces=np.asarray([[0,2,1],[0,1,3],[0,3,2],[1,2,3]],dtype=np.int64)
    return vertices,faces


def tokens(vertices):
    return [[str(x.numerator)+'/'+str(x.denominator) for x in row] for row in vertices]


def test_valid_single_shell_matches_full_region_acceptance():
    vertices,faces=tetra();native=api()
    audit=native.audit_exact_rational_triangle_boundary(tokens(vertices),faces,'AnyNonzero')
    region=native.make_exact_material_region_from_rationals(tokens(vertices),faces)
    assert audit['boundary_valid'] and region.facts['valid']
    assert audit['component_volume_signs']==[1]
    assert audit['self_intersection_completed'] and audit['self_intersection_free']


def test_rational_sliver_is_checked_before_binary64_collapse():
    vertices,faces=tetra()
    vertices=[[p[0],p[1],Fraction(1)+p[2]*Fraction(1,2**100)] for p in vertices]
    audit=api().audit_exact_rational_triangle_boundary(tokens(vertices),faces,'AnyNonzero')
    assert audit['boundary_valid'] and not audit['degenerate_faces']
    assert audit['component_signed_volumes_exact']==[str(Fraction(1,6*2**100))]


@pytest.mark.parametrize('fault',['open','wrong_orientation','self_intersection'])
def test_invalid_boundaries_are_rejected(fault):
    vertices,faces=tetra();mode='AnyNonzero'
    if fault=='open':faces=faces[:3]
    elif fault=='wrong_orientation':faces=faces[:,::-1];mode='OutwardPositive'
    else:
        vertices=vertices+[[x+Fraction(1,4) for x in p] for p in vertices]
        faces=np.vstack([faces,faces+4])
    audit=api().audit_exact_rational_triangle_boundary(tokens(vertices),faces,mode)
    assert not audit['boundary_valid']
    if fault=='self_intersection':assert not audit['self_intersection_free']
