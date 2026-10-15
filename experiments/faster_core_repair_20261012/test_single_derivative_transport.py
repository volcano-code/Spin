"""Canonical evaluator must assemble its final Jacobian only once."""
import json
from pathlib import Path
import numpy as np
from make_it_stand import spin
HERE=Path(__file__).resolve().parent


def test_canonical_fullrank_has_one_transport():
    frozen=json.loads((HERE/'same-input-job-03/output/frozen-inputs.json').read_text())
    request=dict(frozen['request']);request['settings']=spin.SpinSolverConfig(**request['settings'])
    request=spin.SpinItFasterRequest(**request)
    settings_values=dict(frozen['settings'])
    # The exact-material production profile exposes both transport stages.
    settings_values['exact_binary64_material_boundary']=True
    settings=spin.SpinEvaluatorConfig(**settings_values)
    q=np.asarray(frozen['states']['ellipsoid']);q/=np.linalg.norm(q)
    events=[]
    result=spin.evaluate_spin_it_faster_fullrank_quadric(request,q,evaluator_settings=settings,
        diagnostic_stage_observer=events.append)
    assert result.derivative_defined,result.derivative_diagnostic
    transport=[x for x in events if x in ('material-transport-begin','canonical-material-transport-begin')]
    assert len(transport)==1,transport
