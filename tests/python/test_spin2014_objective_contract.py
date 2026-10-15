"""Catch objective substitution and unit loss before formal solver integration."""
from dataclasses import replace
import pytest
from make_it_stand.spin import Spin2014QuadricObjectiveConfig
from make_it_stand.verification.spin2014_objective_contract import Spin2014ObjectiveContract


@pytest.mark.parametrize('selector,want,count', [
    ('TOP_INERTIA_ONLY', '25/8', 4),
    ('TOP_LOWERING_ONLY', '57/8', 4),
    ('TOP_MASS_AND_LOWERING', '153/8', 4),
    ('YOYO', '25/8', 5),
    ('STATIC_LOWERING', '4/1', 2),
])
def test_original_formula_and_constraint_count(selector, want, count):
    # M=2, COM=(0,0,2), centered Q=diag(2,2,3): I=(5,5,4).
    contract = Spin2014ObjectiveContract(Spin2014QuadricObjectiveConfig(objective=selector))
    result = contract.independent_reference(('2/1','0/1','0/1','4/1','2/1','2/1','11/1','0/1','0/1','0/1'))
    assert result['objective'] == want
    assert len(result['constraints']) == count
    assert contract.mode == ('STATIC_BALANCE' if count == 2 else 'YO_YO' if count == 5 else 'SPINNING_TOP')


def test_units_are_bound_and_applied_once():
    config = Spin2014QuadricObjectiveConfig(objective='STATIC_LOWERING')
    original = Spin2014ObjectiveContract(config)
    scaled = Spin2014ObjectiveContract(replace(config, solver_length_unit_in_declared_units=3))
    moments = ('2/1','0/1','0/1','4/1','2/1','2/1','11/1','0/1','0/1','0/1')
    assert scaled.independent_reference(moments)['objective'] == '36/1'
    assert scaled.identity != original.identity
    assert scaled.identity == Spin2014ObjectiveContract(scaled.settings).identity


def test_mode_mismatch_is_rejected():
    contract = Spin2014ObjectiveContract(Spin2014QuadricObjectiveConfig(objective='YOYO'))
    with pytest.raises(ValueError, match='mode'):
        contract.require_mode('SPINNING_TOP')


def test_public_external_entry_rejects_original_mode_before_native_loading():
    from test_spin2014_quadric_interface import Quadric2014Interface
    from make_it_stand.spin import solve_spin_it_faster_with_external_moment_evaluator
    request = Quadric2014Interface().request()
    with pytest.raises(ValueError, match='mode'):
        solve_spin_it_faster_with_external_moment_evaluator(request, lambda q: None,
            objective_settings=Spin2014QuadricObjectiveConfig(objective='YOYO'))
