"""Run-scoped preparation must preserve real material/frame validation."""
from dataclasses import replace
from copy import deepcopy
import json
import numpy as np
import pytest
from make_it_stand.verification import faster_material as fm
from test_cherchi_rational_boundary_profile import actual_material
from make_it_stand import spin
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator
from test_cherchi_rational_boundary_profile import CORE


def test_fixed_preparation_reused_without_reintegration(actual_material):
    problem, request, q, material = actual_material
    prepare = getattr(fm, '_prepare_external_fixed_material', None)
    assert callable(prepare), 'run-scoped fixed material preparation is missing'
    stages = []
    observer = lambda name, op: (stages.append(name), op())[1]
    prepared = prepare(problem, request, stage_observer=observer)
    assert sum(s.startswith('map_fixed_material_integration_') for s in stages) == 1
    expected = fm._external_geometry_moment_map(problem, request, material, expected_coefficients=q)
    stages.clear()
    for _ in range(2):
        actual = fm._external_geometry_moment_map(problem, request, material,
            expected_coefficients=q, prepared_fixed_material=prepared, stage_observer=observer)
        np.testing.assert_array_equal(actual['raw_moments'], expected['raw_moments'])
        np.testing.assert_array_equal(actual['moment_jacobian'], expected['moment_jacobian'])
    assert not any(s.startswith('map_fixed_material_integration_') for s in stages)
    assert stages.count('map_variable_independent_integration') == 2


def test_fixed_preparation_rejects_other_target_request(actual_material):
    problem, request, q, material = actual_material
    prepare = getattr(fm, '_prepare_external_fixed_material', None)
    assert callable(prepare), 'run-scoped fixed material preparation is missing'
    prepared = prepare(problem, request)
    changed = replace(request, target_origin=(.01, 0., 0.))
    with pytest.raises(ValueError, match='fixed material preparation'):
        fm._external_geometry_moment_map(problem, changed, material,
            expected_coefficients=q, prepared_fixed_material=prepared)


@pytest.fixture
def actual_producer(actual_material):
    problem, request, q, _ = actual_material
    config = json.loads((CORE/'core-benchmark-inputs-job-01/output/small.json').read_text())
    producer = InProcessCherchiEvaluator(request, spin.SpinEvaluatorConfig(**config['evaluator']),
        normalization=problem.normalization,
        design_vertices_physical=config['request']['design_vertices'],
        rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1')
    producer(q)
    return problem, request, q, producer.material(q), producer


def test_layered_map_preserves_values_and_keeps_boundary_checks(actual_producer):
    problem, request, q, material, producer = actual_producer
    validator = getattr(producer, 'validate_same_call_material', None)
    assert callable(validator), 'same-call production material validation is missing'
    prepared = fm._prepare_external_fixed_material(problem, request)
    expected = fm._external_geometry_moment_map(problem, request, material, expected_coefficients=q)
    stages = []
    result = fm._external_geometry_moment_map(problem, request, material,
        expected_coefficients=q, prepared_fixed_material=prepared,
        validation_profile='LAYERED_SAME_CALL_V1', material_validator=validator,
        solver_boundary_profile='EXACT_RATIONAL_BOUNDARY_V1',
        stage_observer=lambda name, op: (stages.append(name), op())[1])
    for field in ('raw_moments', 'variable_raw_moments', 'moment_jacobian'):
        np.testing.assert_array_equal(result[field], expected[field])
    assert 'map_variable_mesh_snapshot' in stages
    assert 'map_variable_independent_integration' not in stages
    assert 'map_same_call_material_validation' in stages


@pytest.mark.parametrize('field', ['variable_raw_moments', 'moment_jacobian', 'variable_vertices'])
def test_layered_map_rejects_modified_same_call_attachment(actual_producer, field):
    problem, request, q, material, producer = actual_producer
    validator = getattr(producer, 'validate_same_call_material', None)
    assert callable(validator), 'same-call production material validation is missing'
    changed = deepcopy(material)
    if field == 'variable_raw_moments': changed[field][0] += .001
    else: changed[field][0][0] += .001
    with pytest.raises(ValueError, match='SAME_CALL_MATERIAL'):
        fm._external_geometry_moment_map(problem, request, changed,
            expected_coefficients=q, validation_profile='LAYERED_SAME_CALL_V1',
            material_validator=validator)


def test_layered_map_requires_same_call_producer(actual_material):
    problem, request, q, material = actual_material
    with pytest.raises(ValueError, match='same-call material validator'):
        fm._external_geometry_moment_map(problem, request, material,
            expected_coefficients=q, validation_profile='LAYERED_SAME_CALL_V1')


def test_layered_map_rejects_producer_from_other_density(actual_producer):
    problem, request, q, material, producer = actual_producer
    changed = replace(request, material_density=2 * request.material_density)
    with pytest.raises(ValueError, match='SAME_CALL_REQUEST'):
        fm._external_geometry_moment_map(problem, changed, material,
            expected_coefficients=q, validation_profile='LAYERED_SAME_CALL_V1',
            material_validator=producer.validate_same_call_material)


def test_layered_map_rejects_unbound_validator_that_only_claims_success(actual_material):
    problem, request, q, material = actual_material
    with pytest.raises(ValueError, match='bound production material validator'):
        fm._external_geometry_moment_map(problem, request, material,
            expected_coefficients=q, validation_profile='LAYERED_SAME_CALL_V1',
            material_validator=lambda *args: None)
