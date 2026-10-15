"""Interval acceleration must preserve exact nearest-even material moments."""
from fractions import Fraction

import numpy as np
import pytest
import _cherchi_rational_moments as arithmetic


FACES = np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]], dtype=np.int64)


def certified(vertices, density):
    owned = arithmetic.PreparedBoundary(vertices, FACES)
    assert hasattr(owned, 'integrate_moments_certified'), 'certified rounding path is not implemented'
    result, evidence = owned.integrate_moments_certified(density)
    np.testing.assert_array_equal(result, owned.integrate_moments(density))
    assert evidence['enclosure_certified']
    return result, evidence


def test_unit_simplex_density_and_exact_zero_components():
    vertices = [['0','0','0'],['1','0','0'],['0','1','0'],['0','0','1']]
    result, _ = certified(vertices, 1.25)
    expected = [Fraction(5,24), *([Fraction(5,96)]*3),
                *([Fraction(1,48)]*3), *([Fraction(1,96)]*3)]
    np.testing.assert_array_equal(result, np.array([float(x) for x in expected]))
    shifted = [[str(Fraction(v[0])-Fraction(1,4)),v[1],v[2]] for v in vertices]
    centered, _ = certified(shifted, 1.25)
    assert centered[1] == 0.


def test_near_halfway_mass_requires_exact_fallback():
    mass = Fraction(1)+Fraction(1,2**53)+Fraction(1,2**600)
    vertices = [['0','0','0'],[str(6*mass),'0','0'],['0','1','0'],['0','0','1']]
    result, evidence = certified(vertices, 1.)
    assert result[0] == float(mass) == np.nextafter(1.,2.)
    assert evidence['fallback_components'] >= 1


def test_large_translation_cancellation_matches_exact_gmp():
    offset = Fraction(2**200)
    vertices = [[str(offset+int(i==j)) for j in range(3)] for i in (-1,0,1,2)]
    result, _ = certified(vertices, .75)
    assert result[0] == .125


def test_certified_transport_matches_exact_and_keeps_selection_validation():
    vertices = [['0','0','0'],['1','0','0'],['0','1','0'],['0','0','1']]
    owned = arithmetic.PreparedBoundary(vertices,FACES)
    assert hasattr(owned,'transport_coefficients_selected_certified'), 'certified transport is not implemented'
    for selection in ([],[0,2],list(range(4))):
        actual,evidence=owned.transport_coefficients_selected_certified(1.25,selection)
        np.testing.assert_array_equal(actual,owned.transport_coefficients_selected(1.25,selection))
        assert evidence['enclosure_certified']
    for selection in ([1,1],[2,1],[-1],[4]):
        with pytest.raises(ValueError):owned.transport_coefficients_selected_certified(1.,selection)


def test_production_adapter_consumes_certified_integral_evidence():
    from types import SimpleNamespace
    from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator
    vertices=[['0','0','0'],['1','0','0'],['0','1','0'],['0','0','1']]
    owned=arithmetic.PreparedBoundary(vertices,FACES)
    adapter=object.__new__(InProcessCherchiEvaluator)
    adapter.request=SimpleNamespace(material_density=1.)
    adapter.arithmetic=arithmetic
    adapter.integral_profile='CERTIFIED_MPFR256_GMP_FALLBACK_V1'
    result=adapter._integrate(None,FACES,None,prepared_boundary=owned)
    np.testing.assert_array_equal(result,owned.integrate_moments(1.))
    assert getattr(adapter,'last_integral_evidence',{}).get('enclosure_certified') is True
    adapter._integrate(None,FACES,[[Fraction(x) for x in v] for v in vertices])


def test_formal_entry_rejects_unknown_integral_profile_before_problem_loading(tmp_path):
    from make_it_stand.verification.cherchi_inprocess import solve_cherchi_inprocess_and_seal_faster
    with pytest.raises(ValueError,match='integral profile'):
        solve_cherchi_inprocess_and_seal_faster(None,evaluator_settings=None,
            solver_settings=None,output_directory=tmp_path/'unused',integral_profile='BOGUS')
    assert not (tmp_path/'unused').exists()
