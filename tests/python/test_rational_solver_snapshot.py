"""Single rational shells retain exact checking; multiple shells retain Nef."""
from fractions import Fraction
import numpy as np
import pytest
from make_it_stand.verification import spin_problem
from make_it_stand.verification.faster_material import _raw_moments


def api():
    assert hasattr(spin_problem,'rational_solver_candidate_snapshot'), 'rational solver snapshot missing'
    return spin_problem.rational_solver_candidate_snapshot


def tetra():
    vertices=[[Fraction(0),Fraction(0),Fraction(0)], [Fraction(1),Fraction(0),Fraction(0)],
              [Fraction(0),Fraction(1),Fraction(0)], [Fraction(0),Fraction(0),Fraction(1)]]
    faces=np.asarray([[0,2,1],[0,1,3],[0,3,2],[1,2,3]],dtype=np.int64)
    return vertices,faces


def tokens(vertices):
    return [[str(x.numerator)+'/'+str(x.denominator) for x in row] for row in vertices]


def test_single_shell_is_immutable_and_matches_independent_region_moments():
    vertices,faces=tetra();rows=tokens(vertices)
    snapshot,mode=api()(rows,faces)
    full=spin_problem.MeshSnapshot.from_rationals(rows,faces)
    assert mode=='EXACT_RATIONAL_SINGLE_OUTWARD_BOUNDARY'
    assert not isinstance(snapshot,spin_problem.MeshSnapshot)
    np.testing.assert_array_equal(_raw_moments(snapshot,1.2),_raw_moments(full,1.2))
    identity=snapshot.identity
    rows[0][0]='1/1';faces[0,0]=3
    assert snapshot.identity==identity
    view,_=snapshot.arrays();view[0,0]=9.
    assert snapshot.arrays()[0][0,0]==0.


def test_exact_sliver_stays_authoritative_when_display_collapses():
    vertices,faces=tetra()
    vertices=[[p[0],p[1],Fraction(1)+p[2]*Fraction(1,2**100)] for p in vertices]
    snapshot,mode=api()(tokens(vertices),faces)
    assert mode=='EXACT_RATIONAL_SINGLE_OUTWARD_BOUNDARY'
    assert len(set(snapshot.arrays()[0][:,2]))==1
    assert _raw_moments(snapshot,1.)[0]==float(Fraction(1,6*2**100))


def test_multiple_shell_cavity_still_uses_full_region_validation():
    vertices,faces=tetra()
    inner=[[Fraction(1,8)+x/4 for x in p] for p in vertices]
    combined=vertices+inner
    triangles=np.vstack([faces,(faces+4)[:,::-1]])
    snapshot,mode=api()(tokens(combined),triangles)
    assert mode=='EXACT_NEF_SIGNED_REGION'
    assert isinstance(snapshot,spin_problem.MeshSnapshot)
    full=spin_problem.MeshSnapshot.from_rationals(tokens(combined),triangles)
    assert snapshot.identity==full.identity


def test_open_boundary_is_rejected():
    vertices,faces=tetra()
    with pytest.raises(ValueError,match='invalid exact rational solver candidate boundary'):
        api()(tokens(vertices),faces[:3])
