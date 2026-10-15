"""Exact outside-box witnesses distinguish disjoint roots from nesting.

Full embedded-boundary checks run first. This optimization applies only to
trial snapshots and must not create a final material/Nef certificate.
"""
from fractions import Fraction

import numpy as np
import pytest

from test_exact_material_region import box, combine, l_prism
from make_it_stand.verification.spin_problem import MeshSnapshot, rational_solver_candidate_snapshot


def tokens(vertices):
    return [[f'{q.numerator}/{q.denominator}' for q in map(lambda x: Fraction(float(x)), row)]
            for row in vertices]


def snapshot(mesh, enabled=True):
    return rational_solver_candidate_snapshot(tokens(mesh[0]), mesh[1], allow_separated_roots=enabled)


def test_three_disjoint_positive_roots_get_bound_trial_proof():
    mesh = combine(box(), box((2, 0, 0), (3, 1, 1)), box((0, 2, 0), (1, 3, 1)))
    result, mode = snapshot(mesh)
    assert mode == 'EXACT_RATIONAL_SEPARATED_OUTWARD_BOUNDARIES'
    assert not isinstance(result, MeshSnapshot)
    with pytest.raises(ValueError, match='mesh must be'):
        MeshSnapshot(result)
    assert result.identity['scope'] == mode
    proof = result.separation_proof
    assert proof['parents'] == [-1, -1, -1]
    assert len(proof['outside_witnesses']) == 6
    proof['parents'][0] = 1
    assert result.separation_proof['parents'] == [-1, -1, -1]
    np.testing.assert_array_equal(result.arrays()[1], mesh[1])


def test_overlapping_boxes_can_still_have_proven_disjoint_regions():
    # The second solid is outside the L, although their bounding boxes overlap.
    mesh = combine(l_prism(), box((1.5, 1.5, .2), (3.5, 2.5, .8)))
    result, mode = snapshot(mesh)
    assert mode == 'EXACT_RATIONAL_SEPARATED_OUTWARD_BOUNDARIES'
    assert len(result.separation_proof['outside_witnesses']) == 2


def test_nested_positive_shell_is_rejected_not_blessed_by_boundary_validity():
    mesh = combine(box((0, 0, 0), (4, 4, 4)), box((1, 1, 1), (2, 2, 2)))
    with pytest.raises(ValueError, match='orientation'):
        snapshot(mesh)


def test_signed_cavity_keeps_full_material_region_path():
    cavity = box((.2, .2, .2), (.8, .8, .8))
    result, mode = snapshot(combine(box(), (cavity[0], cavity[1][:, ::-1])))
    assert mode == 'EXACT_NEF_SIGNED_REGION'
    assert isinstance(result, MeshSnapshot)
    assert result.exact_region().facts['bounded_void_components'] == 1


def test_touching_roots_rejected_by_complete_boundary_check():
    with pytest.raises(ValueError, match='invalid exact rational solver candidate'):
        snapshot(combine(box(), box((1, 0, 0), (2, 1, 1))))


def test_old_default_keeps_full_multiple_shell_check():
    mesh = combine(box(), box((2, 0, 0), (3, 1, 1)))
    result, mode = rational_solver_candidate_snapshot(tokens(mesh[0]), mesh[1])
    assert mode == 'EXACT_NEF_SIGNED_REGION'
    assert isinstance(result, MeshSnapshot)


def test_inconclusive_box_witness_retains_exact_nef_classification():
    # This disjoint root lies inside the L's bounding box, so the cheap
    # certificate cannot prove exteriority. The full region check must run.
    mesh = combine(l_prism(), box((1.5, 1.5, .2), (2.5, 2.5, .8)))
    result, mode = snapshot(mesh)
    assert mode == 'EXACT_NEF_SIGNED_REGION'
    assert result.exact_region().facts['parents'] == [-1, -1]


def test_nonboolean_opt_in_cannot_select_a_validation_path():
    with pytest.raises(ValueError, match='must be boolean'):
        snapshot(box(), enabled='true')
