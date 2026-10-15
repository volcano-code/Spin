"""Diagnostic-only callbacks must reach the ordinary canonical transport."""
import json
from pathlib import Path
import numpy as np
from make_it_stand import spin


def test_canonical_transport_forwards_observer_without_changing_jacobian():
    core=Path(__file__).resolve().parents[2]/'experiments/faster_core_repair_20261012'
    frozen=json.loads((core/'projective-formal-pipeline-job-01/output/registered-inputs.json').read_text())
    payload=dict(frozen['request'])
    payload['settings']=spin.SpinSolverConfig(**payload['settings'])
    request=spin.SpinItFasterRequest(**payload)
    settings=spin.SpinEvaluatorConfig(**frozen['settings'])
    q=np.asarray(frozen['states']['ellipsoid'])
    stages=[]
    observed=spin.evaluate_spin_it_faster_fullrank_quadric(request,q,
        evaluator_settings=settings,diagnostic_stage_observer=stages.append)
    assert observed.derivative_defined
    begin=stages.index('canonical-material-transport-begin')
    end=stages.index('canonical-material-transport-end')
    assert 'transport-source-boundary-audit-begin' in stages[begin+1:end], \
        'canonical overload lost its actual transport observer'
    plain=spin.evaluate_spin_it_faster_fullrank_quadric(request,q,evaluator_settings=settings)
    assert plain.derivative_defined
    np.testing.assert_array_equal(observed.variable_raw_moments,plain.variable_raw_moments)
    np.testing.assert_array_equal(observed.coefficient_jacobian,plain.coefficient_jacobian)


def test_canonical_inner_observer_failure_aborts_evaluation():
    core=Path(__file__).resolve().parents[2]/'experiments/faster_core_repair_20261012'
    frozen=json.loads((core/'projective-formal-pipeline-job-01/output/registered-inputs.json').read_text())
    payload=dict(frozen['request'])
    payload['settings']=spin.SpinSolverConfig(**payload['settings'])
    request=spin.SpinItFasterRequest(**payload)
    active=False
    def observer(stage):
        nonlocal active
        if stage=='canonical-material-transport-begin':active=True
        if active and stage=='transport-source-boundary-audit-begin':
            raise RuntimeError('INDEPENDENT_OBSERVER_ABORT')
    import pytest
    with pytest.raises(RuntimeError,match='INDEPENDENT_OBSERVER_ABORT'):
        spin.evaluate_spin_it_faster_fullrank_quadric(request,frozen['states']['ellipsoid'],
            evaluator_settings=spin.SpinEvaluatorConfig(**frozen['settings']),diagnostic_stage_observer=observer)
