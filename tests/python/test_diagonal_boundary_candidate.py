"""Removing the production diagonal filter must lose measurable pair pruning."""
import numpy as np
from make_it_stand import spin


def test_separated_closed_components_prune_pairs_without_losing_contacts():
    vertices = np.array([[0,0,0],[2,0,0],[0,2,0],[0,0,1],
                         [2,2,0],[2,.5,0],[.5,2,0],[2,2,1]], dtype=float)
    faces = np.array([[1,2,3],[0,3,2],[0,1,3],[0,2,1]], dtype=np.int32)
    faces = np.vstack([faces, faces+4])
    bounds = vertices[faces]
    low, high = bounds.min(axis=1), bounds.max(axis=1)
    coordinate_candidates = sum(np.all(high[i]>=low[j]) and np.all(high[j]>=low[i])
                                for i in range(8) for j in range(i+1,8))
    n = spin._native_module()
    coordinates = [[f'{float(x).as_integer_ratio()[0]}/{float(x).as_integer_ratio()[1]}'
                    for x in row] for row in vertices]
    proof = n.verify_fixed_rational_boundary(coordinates[:4], faces[:4])
    result = n.audit_rational_boundary_with_fixed_proof(coordinates, faces, 'AnyNonzero', proof)
    assert result['boundary_valid']
    # Six pairs of the first tetrahedron are covered by its immutable proof.
    assert result['exact_pair_tests'] + result['reused_fixed_pair_tests'] < coordinate_candidates
    assert result['reused_fixed_pair_tests'] == 6
