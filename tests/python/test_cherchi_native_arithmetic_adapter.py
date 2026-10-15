"""Explicit arithmetic selection must preserve actual Boolean material."""
from dataclasses import replace
import inspect
import json
from pathlib import Path

import numpy as np
import pytest
from make_it_stand import spin
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator

CORE = Path(__file__).resolve().parents[2] / 'experiments/faster_core_repair_20261012'


def inputs():
    config = json.loads((CORE / 'projective-formal-pipeline-job-01/output/registered-inputs.json').read_text())
    payload = dict(config['request'])
    payload['settings'] = spin.SpinSolverConfig(**payload['settings'])
    return spin.SpinItFasterRequest(**payload), spin.SpinEvaluatorConfig(**config['settings']), config['states']


def evaluator(request, settings, backend):
    assert 'rational_arithmetic_backend' in inspect.signature(InProcessCherchiEvaluator).parameters, \
        'formal adapter has no explicit native arithmetic selection'
    return InProcessCherchiEvaluator(request, settings, rational_arithmetic_backend=backend)


@pytest.mark.parametrize('state', ['ellipsoid', 'one_sheet', 'two_sheet'])
def test_actual_material_and_jacobian_match_with_explicit_native_arithmetic(state):
    request, settings, states = inputs()
    q = np.asarray(states[state]); q /= np.linalg.norm(q)
    reference = evaluator(request, settings, 'PYTHON_FRACTION_V1')
    candidate = evaluator(request, settings, 'NATIVE_GMP_SURFACE_V1')
    before, after = reference(q), candidate(q)
    for field in ('raw_moments', 'variable_raw_moments', 'moment_jacobian'):
        np.testing.assert_array_equal(before[field], after[field])
    for field in ('geometry_hash', 'topology_hash'):
        assert before[field] == after[field]
    assert before['evaluation_map_identity'] != after['evaluation_map_identity']
    old, new = reference.material(q), candidate.material(q)
    for field in ('variable_vertices', 'variable_rational_vertices', 'variable_faces', 'material_proof'):
        assert old[field] == new[field]
    assert new['rational_arithmetic_identity']['backend'] == 'NATIVE_GMP_SURFACE_V1'
    assert len(new['rational_arithmetic_identity']['module_sha256']) == 64
    assert candidate.call_counts['source'] == candidate.call_counts['boolean'] == 1


def test_default_stays_equivalent_to_explicit_python_arithmetic():
    request, settings, _ = inputs()
    default = InProcessCherchiEvaluator(request, settings)
    explicit = evaluator(request, settings, 'PYTHON_FRACTION_V1')
    assert default.context_identity == explicit.context_identity


def test_unsupported_arithmetic_is_rejected():
    request, settings, _ = inputs()
    with pytest.raises(ValueError, match='UNSUPPORTED_RATIONAL_ARITHMETIC'):
        evaluator(request, settings, 'UNREGISTERED')


def test_native_arithmetic_cache_does_not_cross_input_density():
    request, settings, states = inputs()
    changed = replace(request, material_density=request.material_density * 2.)
    q = np.asarray(states['ellipsoid']); q /= np.linalg.norm(q)
    first = evaluator(request, settings, 'NATIVE_GMP_SURFACE_V1')
    second = evaluator(changed, settings, 'NATIVE_GMP_SURFACE_V1')
    a, b = first(q), second(q)
    assert first.context_identity != second.context_identity
    np.testing.assert_array_equal(b['variable_raw_moments'], 2. * np.asarray(a['variable_raw_moments']))
    np.testing.assert_array_equal(b['moment_jacobian'], 2. * np.asarray(a['moment_jacobian']))


def test_same_boolean_prepares_exact_arithmetic_once_for_moments_and_derivative(monkeypatch):
    import _cherchi_rational_moments as arithmetic
    assert hasattr(arithmetic,'PreparedBoundary')
    original = arithmetic.PreparedBoundary
    calls=[]
    class ObservedBoundary:
        def __init__(self,*args):
            calls.append('prepare')
            self.owned=original(*args)
        def integrate_moments(self,*args):
            calls.append('integrate')
            return self.owned.integrate_moments(*args)
        def transport_coefficients_selected(self,*args):
            calls.append('transport')
            return self.owned.transport_coefficients_selected(*args)
    monkeypatch.setattr(arithmetic,'PreparedBoundary',ObservedBoundary)
    request,settings,states=inputs()
    q=np.asarray(states['ellipsoid']);q/=np.linalg.norm(q)
    candidate=evaluator(request,settings,'NATIVE_GMP_SURFACE_V1')
    candidate(q)
    assert calls==['prepare','transport','integrate']
