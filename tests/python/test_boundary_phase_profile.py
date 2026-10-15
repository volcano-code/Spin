"""Diagnostics must measure the real checker without changing its decisions."""
import math

import numpy as np
import pytest

from make_it_stand import spin

V = [['0/1','0/1','0/1'], ['1/1','0/1','0/1'],
     ['0/1','1/1','0/1'], ['0/1','0/1','1/1']]
F = np.asarray([[1,2,3],[0,3,2],[0,1,3],[0,2,1]], dtype=np.int32)


def measured(vertices, faces):
    n = spin._native_module()
    proof = n.verify_fixed_rational_boundary(V, F)
    plain = n.audit_rational_boundary_with_fixed_proof(vertices, faces, 'AnyNonzero', proof)
    result = n.audit_rational_boundary_with_fixed_proof(
        vertices, faces, 'AnyNonzero', proof, profile=True)
    diagnostics = result.pop('performance')
    assert result == plain  # Missing checks or time leaking into facts breaks this.
    assert diagnostics['vertex_count'] == len(vertices)
    assert diagnostics['face_count'] == len(faces)
    assert diagnostics['exact_pair_tests'] == plain['exact_pair_tests']
    assert diagnostics['reused_pair_tests'] == plain['reused_fixed_pair_tests']
    phases = diagnostics['exclusive_seconds']
    for name in ('index_validation', 'rational_parse', 'canonical_digest',
                 'original_face_matching', 'topology_and_degeneracy',
                 'component_volumes', 'bvh_build', 'intersection_query', 'finalize'):
        assert name in phases
        assert math.isfinite(phases[name]) and phases[name] >= 0
    assert sum(phases.values()) <= diagnostics['total_seconds']
    assert diagnostics['total_seconds'] > 0
    return result, diagnostics


def test_profile_preserves_valid_boundary_and_reports_reused_work():
    facts, diagnostics = measured(V, F)
    assert facts['boundary_valid']
    assert diagnostics['matched_original_faces'] == 4
    assert diagnostics['reused_pair_tests'] == 6


def test_profile_preserves_actual_intersection_rejection():
    shifted = [['1/4','1/4','1/4'], ['5/4','1/4','1/4'],
               ['1/4','5/4','1/4'], ['1/4','1/4','5/4']]
    facts, diagnostics = measured(V + shifted, np.vstack([F, F+4]))
    assert not facts['boundary_valid']
    assert not facts['self_intersection_free']
    assert diagnostics['exact_pair_tests'] > 0


def test_profile_does_not_bypass_invalid_indices():
    n = spin._native_module()
    proof = n.verify_fixed_rational_boundary(V, F)
    bad = F.copy(); bad[0, 1] = bad[0, 0]
    with pytest.raises(ValueError, match='distinct'):
        n.audit_rational_boundary_with_fixed_proof(V, bad, 'AnyNonzero', proof, profile=True)
