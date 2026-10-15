"""Comparison colors encode physical roles, not optimization status."""
import importlib
import pytest


def test_comparison_uses_same_polymer_for_design_and_retained():
    policy = importlib.import_module('flower_render_material_policy')
    assert policy.comparison_material_role('design', ['top-cap']) == 'polymer'
    assert policy.comparison_material_role('retained', ['top-cap']) == 'polymer'
    assert policy.comparison_material_role('top', ['top-cap']) == 'fixed'
    assert policy.comparison_material_role('bottom', ['top-cap']) == 'fixed'
    assert policy.comparison_material_role('top-cap', ['top-cap']) == 'metal'
    assert policy.comparison_material_role('outer-race', ['top-cap']) == 'metal'


@pytest.mark.parametrize('profile', ['TYPO', '', None])
def test_unknown_material_profile_is_rejected(profile):
    policy = importlib.import_module('flower_render_material_policy')
    with pytest.raises(ValueError, match='UNKNOWN_RENDER_MATERIAL_PROFILE'):
        policy.require_material_profile(profile)
