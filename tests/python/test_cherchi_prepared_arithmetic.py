"""One owned exact boundary must preserve the existing arithmetic bitwise."""
import numpy as np
import pytest


def api():
    import _cherchi_rational_moments as arithmetic
    assert hasattr(arithmetic,'PreparedBoundary'), 'prepared exact arithmetic is unavailable'
    return arithmetic


def tetrahedron():
    return [['0','0','0'],['1/3','0','0'],['0','2/7','0'],['0','0','5/11']], np.array(
        [[0,2,1],[0,1,3],[1,2,3],[2,0,3]],dtype=np.int64)


def test_owned_prepared_boundary_preserves_integrals_and_transport_exactly():
    a = api(); vertices,faces = tetrahedron()
    prepared = a.PreparedBoundary(vertices,faces)
    for density in (.7,1.,1200.):
        np.testing.assert_array_equal(prepared.integrate_moments(density),a.integrate_moments(vertices,faces,density))
        np.testing.assert_array_equal(prepared.transport_coefficients_selected(density,[1,3]),
            a.transport_coefficients_selected(vertices,faces,density,[1,3]))


def test_prepared_boundary_owns_coordinates_and_connectivity():
    a = api(); vertices,faces = tetrahedron()
    expected = a.integrate_moments(vertices,faces,1.)
    prepared = a.PreparedBoundary(vertices,faces)
    vertices[0][0] = '200'; faces[0,0] = 999
    np.testing.assert_array_equal(prepared.integrate_moments(1.),expected)


def test_prepared_boundary_rejects_invalid_geometry_and_duplicate_selection():
    a = api(); vertices,faces = tetrahedron()
    with pytest.raises(ValueError):
        a.PreparedBoundary(vertices,faces[:-1])
    prepared = a.PreparedBoundary(vertices,faces)
    with pytest.raises(ValueError):
        prepared.transport_coefficients_selected(1.,[0,0])
    with pytest.raises(ValueError):
        prepared.transport_coefficients_selected(1.,[3,1])
