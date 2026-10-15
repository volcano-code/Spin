"""Native arithmetic must preserve the registered exact surface formulas."""
import importlib.util
from fractions import Fraction
from math import fsum

import numpy as np
import pytest
from make_it_stand.verification.cherchi_motion_differential import (
    integrate_rational_polyhedral_raw_moments, differentiate_polyhedral_raw_moments)


def engine():
    assert importlib.util.find_spec('_cherchi_rational_moments') is not None, 'native exact arithmetic module missing'
    import _cherchi_rational_moments
    return _cherchi_rational_moments


def tetra():
    v = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    f = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]], dtype=np.int64)
    return v, f


def rationals(v):
    return [[Fraction(float(x)) for x in p] for p in v]


def strings(q):
    return [[str(x) for x in p] for p in q]


def test_surface_moments_match_exact_production_formula():
    v, f = tetra()
    v = v * np.array([1.25, .75, 2.]) + [.125, -.25, .5]
    q = rationals(v)
    expected = integrate_rational_polyhedral_raw_moments(v, f, 1.2, rational_vertices=q)
    np.testing.assert_array_equal(engine().integrate_moments(strings(q), f, 1.2), expected)


def test_native_coefficients_preserve_actual_reynolds_differential():
    v, f = tetra()
    q = rationals(v)
    motion = np.arange(24, dtype=float).reshape(4, 3, 2) / 32.
    expected = differentiate_polyhedral_raw_moments(v, f, motion, 1.2, rational_vertices=q)
    coefficients = engine().transport_coefficients(strings(q), f, 1.2)
    terms = np.einsum('fjmc,fjcd->fmd', coefficients, motion[f])
    actual = np.array([[fsum(terms[:, m, d]) for d in range(2)] for m in range(10)])
    np.testing.assert_array_equal(actual, expected)


def test_tiny_nonzero_faces_remain_exact():
    v, f = tetra()
    v[:, 1:] *= 2. ** -200
    q = rationals(v)
    expected = integrate_rational_polyhedral_raw_moments(v, f, rational_vertices=q)
    np.testing.assert_array_equal(engine().integrate_moments(strings(q), f, 1.), expected)
    assert expected[0] > 0.


@pytest.mark.parametrize('value', [Fraction(1) + Fraction(1, 2**53),
    Fraction(1) + Fraction(3, 2**53), Fraction(3, 2**1075),
    Fraction(3 * 2**120 - 1, 2**1195), -Fraction(3 * 2**120 - 1, 2**1195)])
def test_nearest_even_rounding_including_subnormal_values(value):
    assert engine().round_to_binary64(str(value)) == float(value)


@pytest.mark.parametrize('fault', ['open', 'zero_area', 'bad_density', 'bad_denominator'])
def test_invalid_input_is_rejected(fault):
    v, f = tetra()
    q = strings(rationals(v))
    density = 1.
    if fault == 'open':
        f = np.vstack([f[:3], f[:1]])
    elif fault == 'zero_area':
        q[3] = ['0', '0', '0']
    elif fault == 'bad_density':
        density = float('nan')
    else:
        q[0][0] = '1/0'
    with pytest.raises(ValueError):
        engine().integrate_moments(q, f, density)
