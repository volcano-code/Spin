"""Registered two-component projective source must have a consistent selector."""
import json
from pathlib import Path
import numpy as np
import pytest
from make_it_stand import spin

HERE=Path(__file__).resolve().parent


@pytest.mark.parametrize('sign', (1., -1.))
def test_registered_projective_source_component_selector(sign):
    registration=json.loads((HERE/'noodles-first-core-failures-registration-01/registration.json').read_text())
    config=json.loads(Path(registration['mapping_config']).read_text())
    request=config['request'];settings=spin.SpinEvaluatorConfig(**config['evaluator'])
    native=spin._native_module();tri=native.SpinTriangulationSettings()
    for key in ('icosphere_subdivisions','maximum_icosphere_subdivisions','maximum_radius_doublings',
                'retain_fixed_normal_provenance','truncation_radius','paper_projective_hyperboloid','algebraic_finite_map'):
        setattr(tri,key,getattr(settings,key))
    vertices=(np.asarray(request['design_vertices'])-np.asarray(request['target_origin']))@np.asarray(request['target_axes'])
    q=sign*np.asarray(registration['failures']['source_side']['record']['coefficients'])
    result=native.evaluate_spin_it_faster_source_motion(q,vertices,tri)
    assert result['defined'],result['diagnostic']
    assert result['side_selection_defined'],result.get('side_selection_diagnostic')
    assert result['interior_is_positive'] is (sign > 0)
