"""The trial boundary profile keeps material and frame checks in force."""
from copy import deepcopy
import inspect
import json
from pathlib import Path

import numpy as np
import pytest
from make_it_stand import spin
from make_it_stand.verification.spin_problem import SpinArtifactProblem
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator
from make_it_stand.verification.faster_material import _prepare, _external_geometry_moment_map

CORE=Path(__file__).resolve().parents[2]/'experiments/faster_core_repair_20261012'


@pytest.fixture
def actual_material():
    config=json.loads((CORE/'core-benchmark-inputs-job-01/output/small.json').read_text())
    r=config['request'];fixed=config['fixed_material']
    problem=SpinArtifactProblem(original=(r['outer_vertices'],r['outer_faces']),
        fixed_regions=[('shell',(fixed['vertices'],fixed['faces']),fixed['density'])],
        design_regions=[('design',(r['design_vertices'],r['design_faces']),r['material_density'])],
        length_unit='m',mass_unit='kg',minimum_outer_wall=.05,
        target=dict(mode='SPINNING_TOP',axis=[0.,0.,1.],origin=r['target_origin'],support_points=[]))
    settings=spin.SpinEvaluatorConfig(**config['evaluator'])
    problem,request,_=_prepare(problem,settings,spin.SpinSolverConfig(**r['settings']))
    q=json.loads((CORE/'core-benchmark-small-seed-0-job-01/output/run/terminal-material-reference.json').read_text())['coefficients']
    evaluator=InProcessCherchiEvaluator(request,settings,normalization=problem.normalization,
        design_vertices_physical=r['design_vertices'],rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1')
    evaluator(q)
    return problem,request,q,evaluator.material(q)


def mapped(problem,request,q,material,profile):
    assert 'solver_boundary_profile' in inspect.signature(_external_geometry_moment_map).parameters, \
        'external map has no explicit rational trial-boundary profile'
    return _external_geometry_moment_map(problem,request,material,expected_coefficients=q,
        solver_boundary_profile=profile)


def test_rational_trial_boundary_preserves_actual_residual_and_jacobian(actual_material):
    problem,request,q,material=actual_material
    full=mapped(problem,request,q,material,'EXACT_NEF_REGION_V1')
    candidate=mapped(problem,request,q,material,'EXACT_RATIONAL_BOUNDARY_V1')
    for field in ('raw_moments','variable_raw_moments','moment_jacobian'):
        np.testing.assert_array_equal(full[field],candidate[field])
    for field in ('geometry_hash','topology_hash','evaluation_map_identity'):
        assert full[field]==candidate[field]
    assert candidate['variable_geometry_validation']=='EXACT_RATIONAL_SINGLE_OUTWARD_BOUNDARY'


@pytest.mark.parametrize('fault',['wrong_moments','wrong_frame','wrong_display','bad_jacobian','wrong_backend'])
def test_existing_material_checks_still_reject_bad_trial_data(actual_material,fault):
    problem,request,q,material=actual_material
    changed=deepcopy(material)
    if fault=='wrong_moments':changed['variable_raw_moments'][0]+=.001
    elif fault=='wrong_frame':changed['variable_vertices_model'][0][0]+=.01
    elif fault=='wrong_display':changed['variable_vertices'][0][0]+=.01
    elif fault=='bad_jacobian':changed['moment_jacobian']=changed['moment_jacobian'][:-1]
    else:changed['boolean_backend_identity']='cgal'
    with pytest.raises(ValueError):
        mapped(problem,request,q,changed,'EXACT_RATIONAL_BOUNDARY_V1')


def test_unknown_boundary_profile_is_rejected(actual_material):
    problem,request,q,material=actual_material
    with pytest.raises(ValueError,match='unknown solver boundary profile'):
        mapped(problem,request,q,material,'UNREGISTERED')


def test_separated_root_profile_preserves_same_call_material_and_jacobian(actual_material):
    problem, request, q, _ = actual_material
    config = json.loads((CORE/'core-benchmark-inputs-job-01/output/small.json').read_text())
    producer = InProcessCherchiEvaluator(request, spin.SpinEvaluatorConfig(**config['evaluator']),
        normalization=problem.normalization,
        design_vertices_physical=problem.design_regions[0].mesh.arrays()[0],
        rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1', fixed_boundary_reuse=True)
    producer(q)
    material = producer.material(q)
    def evaluate(profile, value):
        return _external_geometry_moment_map(problem, request, value, expected_coefficients=q,
            solver_boundary_profile=profile, validation_profile='LAYERED_SAME_CALL_V1',
            material_validator=producer.validate_same_call_material)
    baseline = evaluate('EXACT_RATIONAL_FIXED_SUBSET_V1', material)
    candidate = evaluate('EXACT_RATIONAL_FIXED_SEPARATED_ROOTS_V1', material)
    for field in ('raw_moments', 'variable_raw_moments', 'moment_jacobian'):
        np.testing.assert_array_equal(baseline[field], candidate[field])
    changed = deepcopy(material)
    changed['variable_faces'][0][0] = changed['variable_faces'][0][1]
    with pytest.raises(ValueError):
        evaluate('EXACT_RATIONAL_FIXED_SEPARATED_ROOTS_V1', changed)


def test_signed_classification_profile_preserves_full_region_map(actual_material):
    problem, request, q, material = actual_material
    config = json.loads((CORE/'core-benchmark-inputs-job-01/output/small.json').read_text())
    producer = InProcessCherchiEvaluator(request, spin.SpinEvaluatorConfig(**config['evaluator']),
        normalization=problem.normalization,
        design_vertices_physical=problem.design_regions[0].mesh.arrays()[0],
        rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1', fixed_boundary_reuse=True)
    producer(q)
    material = producer.material(q)
    baseline = mapped(problem, request, q, material, 'EXACT_NEF_REGION_V1')
    candidate = _external_geometry_moment_map(problem, request, material, expected_coefficients=q,
        solver_boundary_profile='EXACT_RATIONAL_FIXED_CLASSIFIED_V1',
        validation_profile='LAYERED_SAME_CALL_V1', material_validator=producer.validate_same_call_material)
    for field in ('raw_moments', 'variable_raw_moments', 'moment_jacobian'):
        np.testing.assert_array_equal(baseline[field], candidate[field])
    for field in ('geometry_hash', 'topology_hash', 'evaluation_map_identity'):
        assert baseline[field] == candidate[field]
    assert candidate['variable_geometry_validation'] == \
        'EXACT_RATIONAL_SIGNED_BOUNDARY_CLASSIFICATION'


def test_signed_classification_profile_keeps_same_call_material_guards(actual_material):
    problem, request, q, material = actual_material
    changed = deepcopy(material)
    changed['variable_raw_moments'][0] += .001
    with pytest.raises(ValueError):
        mapped(problem, request, q, changed, 'EXACT_RATIONAL_FIXED_CLASSIFIED_V1')
