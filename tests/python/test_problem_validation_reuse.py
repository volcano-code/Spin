"""Reuse only the exact constructor-checked live input, never caller flags."""
import gc
from unittest.mock import patch
import weakref

import pytest

from make_it_stand.verification import spin_problem as p
from make_it_stand.verification import spin_artifact as a
from test_exact_material_region import box


def problem():
    outer = box()
    return p.SpinArtifactProblem(original=outer, fixed_regions=[], design_regions=[('d', outer, 1.)],
        length_unit='m', mass_unit='kg', minimum_outer_wall=.1,
        target=dict(mode='SPINNING_TOP', axis=[0., 0., 1.], origin=[0., 0., 0.], support_points=[]))


def test_opt_in_uses_actual_constructor_partition_once():
    value = problem()
    with patch.object(p, 'partition', side_effect=AssertionError('duplicate partition')):
        assert a._recheck_problem(value, reuse_validation=True) is value


def test_default_still_reconstructs_and_rechecks_input():
    value = problem()
    with patch.object(p, 'partition', wraps=p.partition) as partition:
        fresh = a._recheck_problem(value)
    assert fresh is not value and fresh.input_digest == value.input_digest
    assert partition.call_count == 1


@pytest.mark.parametrize('fault', ['density', 'bool_density', 'mesh_bytes', 'target', 'normalization', 'digest'])
def test_reuse_rejects_live_dependency_mutation(fault):
    value = problem()
    if fault in ('density', 'bool_density'):
        object.__setattr__(value.design_regions[0], 'density', 2. if fault == 'density' else True)
    elif fault == 'mesh_bytes':
        mesh = value.design_regions[0].mesh
        object.__setattr__(mesh, '_vertices', b'\x01'+mesh._vertices[1:])
    else:
        key = dict(target='_target', normalization='_normalization', digest='input_digest')[fault]
        object.__setattr__(value, key, 'changed')
    with pytest.raises(ValueError):
        a._recheck_problem(value, reuse_validation=True)


def test_copied_or_self_reported_problem_cannot_borrow_constructor_proof():
    value = problem()
    forged = object.__new__(p.SpinArtifactProblem)
    for name, item in vars(value).items():
        object.__setattr__(forged, name, item)
    for candidate in (forged, dict(validated=True, problem=value)):
        with pytest.raises(ValueError):
            a._recheck_problem(candidate, reuse_validation=True)


def test_reuse_rejects_changed_implementation_observation():
    value = problem()
    with patch('make_it_stand.verification.spin_common.implementation_observation', return_value={'changed': True}):
        with pytest.raises(ValueError, match='implementation'):
            a._recheck_problem(value, reuse_validation=True)


def test_constructor_registry_does_not_keep_problem_alive():
    value = problem()
    reference, key = weakref.ref(value), id(value)
    registry = getattr(p, '_VALIDATED_PROBLEMS', None)
    assert registry is not None and key in registry
    del value
    gc.collect()
    assert reference() is None and key not in registry
