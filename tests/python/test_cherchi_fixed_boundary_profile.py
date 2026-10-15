"""The formal cached profile must consume the producer's actual fixed proof."""
import inspect
import numpy as np
import pytest
from make_it_stand.verification import faster_material as fm
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator
from test_cherchi_layered_validation import actual_material, actual_producer


def cached_producer(actual):
    problem, request, q, _, producer = actual
    assert 'fixed_boundary_reuse' in inspect.signature(InProcessCherchiEvaluator).parameters, \
        'formal fixed-boundary preparation option is missing'
    candidate = InProcessCherchiEvaluator(request, producer.settings,
        normalization=problem.normalization, design_vertices_physical=producer.design,
        rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1', fixed_boundary_reuse=True)
    candidate(q)
    return candidate


def test_fixed_boundary_profile_preserves_actual_moments_and_jacobian(actual_producer):
    problem, request, q, material, _ = actual_producer
    candidate = cached_producer(actual_producer)
    expected = fm._external_geometry_moment_map(problem, request, material, expected_coefficients=q)
    actual = fm._external_geometry_moment_map(problem, request, candidate.material(q),
        expected_coefficients=q, validation_profile='LAYERED_SAME_CALL_V1',
        material_validator=candidate.validate_same_call_material,
        solver_boundary_profile='EXACT_RATIONAL_FIXED_SUBSET_V1')
    for key in ('raw_moments','variable_raw_moments','moment_jacobian'):
        np.testing.assert_array_equal(actual[key], expected[key])
    assert candidate.call_counts['fixed_boundary_proof'] == 1
    candidate(q)
    assert candidate.call_counts['fixed_boundary_proof'] == 1


def test_cached_profile_requires_prepared_actual_native_proof(actual_producer):
    problem, request, q, material, producer = actual_producer
    with pytest.raises(ValueError, match='prepared native fixed boundary proof'):
        fm._external_geometry_moment_map(problem, request, material, expected_coefficients=q,
            validation_profile='LAYERED_SAME_CALL_V1',material_validator=producer.validate_same_call_material,
            solver_boundary_profile='EXACT_RATIONAL_FIXED_SUBSET_V1')
