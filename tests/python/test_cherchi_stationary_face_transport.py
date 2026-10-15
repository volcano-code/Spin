"""Exact-zero motion faces contribute no Reynolds transport; tiny motion remains."""
from fractions import Fraction
import numpy as np
import pytest
from test_cherchi_layered_validation import actual_producer, actual_material


@pytest.mark.parametrize('speed', [1., 2.**-200])
def test_stationary_faces_skip_work_without_discarding_tiny_motion(actual_producer, monkeypatch, speed):
    *_, producer = actual_producer
    points = np.asarray([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    faces = np.asarray([[1,2,3],[0,3,2],[0,1,3],[0,2,1]])
    rational = tuple(tuple(Fraction(float(x)) for x in p) for p in points)
    moving = np.zeros((4,3,1)); moving[0,0,0] = speed
    original = getattr(producer.arithmetic, 'transport_coefficients_selected', None)
    assert callable(original), 'closed-boundary selected-face transport interface is missing'
    observed = []
    def observe(tokens, full_faces, density, selected):
        observed.append((len(full_faces), list(selected)))
        return original(tokens, full_faces, density, selected)
    monkeypatch.setattr(producer.arithmetic, 'transport_coefficients_selected', observe)
    derivative = producer._differentiate(points, faces, moving, rational)
    # V=(1-x)/6; first/second moments follow tetrahedron monomial integrals.
    expected = speed*np.asarray([-1/6,0.,-1/24,-1/24,0.,-1/60,-1/60,-1/120,0.,0.])
    np.testing.assert_allclose(derivative[:,0], expected, atol=0., rtol=1e-14)
    assert observed == [(4, [1,2,3])], 'full boundary must be checked and only moving faces computed'


def test_entirely_stationary_material_has_exact_zero_transport(actual_producer, monkeypatch):
    *_, producer = actual_producer
    points = np.asarray([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    faces = np.asarray([[1,2,3],[0,3,2],[0,1,3],[0,2,1]])
    rational = tuple(tuple(Fraction(float(x)) for x in p) for p in points)
    calls = []
    original = getattr(producer.arithmetic, 'transport_coefficients_selected', None)
    assert callable(original), 'closed-boundary selected-face transport interface is missing'
    def observe(tokens, full_faces, density, selected):
        calls.append((len(full_faces), list(selected)))
        return original(tokens, full_faces, density, selected)
    monkeypatch.setattr(producer.arithmetic, 'transport_coefficients_selected', observe)
    result = producer._differentiate(points, faces, np.zeros((4,3,10)), rational)
    np.testing.assert_array_equal(result, np.zeros((10,10)))
    assert calls == [(4, [])], 'full stationary boundary must be validated without computing coefficients'


def test_selected_transport_rejects_invalid_full_boundary_and_indices(actual_producer):
    *_, producer = actual_producer
    selected = getattr(producer.arithmetic, 'transport_coefficients_selected', None)
    assert callable(selected), 'closed-boundary selected-face transport interface is missing'
    vertices = [['0/1','0/1','0/1'],['1/1','0/1','0/1'],['0/1','1/1','0/1'],['0/1','0/1','1/1']]
    faces = np.asarray([[1,2,3],[0,3,2],[0,1,3],[0,2,1]], dtype=np.int64)
    with pytest.raises(ValueError): selected(vertices, faces[:3], 1., [])
    for indices in ([1,1], [3,1], [-1], [4]):
        with pytest.raises(ValueError): selected(vertices, faces, 1., indices)
