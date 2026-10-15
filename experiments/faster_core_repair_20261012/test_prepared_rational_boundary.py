"""Catch reparsed/mutable geometry and weakened validation in native reuse.

The new object is arithmetic-only: it must not assert self-intersection or
source certification, which remain the production caller's responsibility.
"""
import numpy as np
import pytest
import _cherchi_rational_moments as arithmetic


def tetra():
    return ([['0', '0', '0'], ['1', '0', '0'],
             ['0', '1', '0'], ['0', '0', '1']],
            np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]],
                     dtype=np.int64))


def prepare(vertices, faces):
    assert hasattr(arithmetic, 'PreparedBoundary'), 'immutable native boundary not implemented'
    return arithmetic.PreparedBoundary(vertices, faces)


def test_owned_geometry_and_analytic_integrals_survive_caller_mutation():
    vertices, faces = tetra()
    owned = prepare(vertices, faces)
    vertices[1][0] = '9'
    faces[:] = 0
    # Unit simplex integrals, independently known for volume, first and second moments.
    expected = np.array([1/6, 1/24, 1/24, 1/24,
                         1/60, 1/60, 1/60, 1/120, 1/120, 1/120])
    np.testing.assert_array_equal(owned.integrate_moments(1.), expected)
    np.testing.assert_array_equal(owned.integrate_moments(2.), expected * 2)
    returned = owned.integrate_moments(1.)
    returned[:] = -1
    np.testing.assert_array_equal(owned.integrate_moments(1.), expected)


def test_selected_transport_matches_legacy_and_rejects_invalid_selections():
    vertices, faces = tetra()
    owned = prepare(vertices, faces)
    expected = arithmetic.transport_coefficients(vertices, faces, 1.25)
    np.testing.assert_array_equal(
        owned.transport_coefficients_selected(1.25, [1, 3]), expected[[1, 3]])
    assert owned.transport_coefficients_selected(1., []).shape == (0, 3, 10, 3)
    for indices in ([1, 1], [3, 1], [-1], [4]):
        with pytest.raises(ValueError):
            owned.transport_coefficients_selected(1., indices)
    for density in (0., -1., float('nan'), float('inf')):
        with pytest.raises(ValueError):
            owned.integrate_moments(density)
        with pytest.raises(ValueError):
            owned.transport_coefficients_selected(density, [])


def test_all_faces_validated_before_any_selection():
    vertices, faces = tetra()
    faces[-1] = faces[0]
    with pytest.raises(ValueError):
        prepare(vertices, faces)
    vertices, faces = tetra()
    vertices[3] = ['0', '0', '0']
    with pytest.raises(ValueError):
        prepare(vertices, faces)
