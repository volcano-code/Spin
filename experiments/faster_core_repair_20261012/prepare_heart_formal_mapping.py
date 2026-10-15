"""Formal input/request and two saved-state comparisons, not optimization.

Keep frozen model coordinates: length label explicitly denotes 25 mm and mass
label 18.75 g. SI conversion belongs to the later verified export contract.
"""
from dataclasses import asdict, replace
from pathlib import Path
import hashlib
import json
import sys
import time
import numpy as np
from make_it_stand import spin
from make_it_stand.verification.spin_problem import SpinArtifactProblem
from make_it_stand.verification.faster_material import (
    _prepare, _prepare_external_fixed_material, _external_geometry_moment_map)
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator, _json_value
from noodles_formal_case import compare_requests
from case_coordinate_contract import normalized_faster_weight

HERE = Path(__file__).resolve().parent
HEART = HERE.parent / 'faster_heart_contact_20261011'


def main(output):
    load = lambda p: json.loads(p.read_text())
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    frozen_path = HEART / 'rounded-clearance-v3-faster-frozen-v2/request.json'
    frozen = load(frozen_path)
    for filename, expected in frozen['input_sources'].items():
        if sha(Path(filename)) != expected:
            raise ValueError('FROZEN_INPUT_CHANGED')
    clearance = load(HERE / 'heart-derived-clearance-transfer-01.json')
    if not clearance['at_least_required'] or clearance['request_sha256'] != sha(frozen_path):
        raise ValueError('FROZEN_CLEARANCE_NOT_BOUND')
    if frozen['units'] != dict(model_length_m=.025, model_mass_kg=.01875,
            density_kg_m3=1200., objective_coordinate_contract='INPUT_MODEL_UNITS_TARGET_ORIGIN_ZERO_V1'):
        raise ValueError('UNREGISTERED_MODEL_UNITS')
    paths = {Path(p).name: Path(p) for p in frozen['input_sources']}
    meshes = {r: load(paths[r + '.binary64.json']) for r in ('outer', 'fixed', 'design')}
    registration = load(HERE / 'heart-cherchi-coordinate-map-job-01/output/registration.json')
    if registration['request_sha256'] != sha(frozen_path):
        raise ValueError('COORDINATE_REFERENCE_NOT_BOUND')
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    mesh = lambda r: (meshes[r]['vertices'], meshes[r]['faces'])
    problem = SpinArtifactProblem(original=mesh('outer'),
        fixed_regions=[('fixed-shell-handle-pin', mesh('fixed'), 1.)],
        design_regions=[('design', mesh('design'), 1.)],
        length_unit='MODEL_UNIT_0.025_M', mass_unit='MODEL_UNIT_0.01875_KG',
        minimum_outer_wall=.015,
        target=dict(mode='SPINNING_TOP', axis=[0., 0., 1.], origin=[0., 0., 0.], support_points=[]))
    print('LOCKED_PROBLEM_READY', time.monotonic() - started, flush=True)
    evaluator = spin.SpinEvaluatorConfig(**registration['settings'])
    solver = spin.SpinSolverConfig(**frozen['requests'][0]['settings'])
    _, request, _ = _prepare(problem, evaluator, solver)
    request = replace(request,
        center_of_mass_weight=normalized_faster_weight(frozen['requests'][0]['center_of_mass_weight'], problem.normalization['L']),
        spinning_top_objective=frozen['requests'][0]['spinning_top_objective'])
    comparison = compare_requests(registration['normalized'], asdict(request), .015/problem.normalization['L'])
    (output / 'prepared-request.json').write_text(json.dumps(asdict(request)))
    mapping = InProcessCherchiEvaluator(request, evaluator, normalization=problem.normalization,
        design_vertices_physical=meshes['design']['vertices'],
        rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1', fixed_boundary_reuse=True)
    fixed = _prepare_external_fixed_material(problem, request)
    rows = []
    for state in ('heart-first-cap-failure', 'heart-first-small-face-rejection'):
        path = HERE / f'heart-cherchi-coordinate-map-job-01/output/{state}-normalized.json'
        saved = load(path)
        mapped = mapping(saved['quadric_coefficients'])
        material = mapping.material(saved['quadric_coefficients'])
        # Same callback used by formal Newton: include actual boundary admission,
        # fixed material and same-call validation, without another Boolean.
        production = _external_geometry_moment_map(problem, request, material,
            expected_coefficients=saved['quadric_coefficients'], prepared_fixed_material=fixed,
            solver_boundary_profile='EXACT_RATIONAL_FIXED_SUBSET_V1',
            validation_profile='LAYERED_SAME_CALL_V1', material_validator=mapping.validate_same_call_material)
        (output / (state + '-material.json')).write_text(json.dumps(material, default=_json_value))
        exact_equal = all(material[k] == saved[k] for k in (
            'variable_rational_vertices', 'variable_faces', 'quadric_source_vertices', 'quadric_source_faces'))
        moment_error = float(np.max(np.abs(np.asarray(mapped['raw_moments'])-saved['raw_moments'])))
        jacobian_error = float(np.max(np.abs(np.asarray(mapped['moment_jacobian'])-saved['moment_jacobian'])))
        row = dict(state=state, material_geometry_equal=exact_equal,
                   moment_error=moment_error, jacobian_error=jacobian_error, reference_sha256=sha(path),
                   production_boundary_validation=production['variable_geometry_validation'])
        rows.append(row)
        if not exact_equal or moment_error > 1e-10 or jacobian_error > 1e-10:
            raise ValueError('FORMAL_MAPPING_DIFFERS_FROM_REGISTERED_STATE')
    report = dict(passed=True, request_comparison=comparison, rows=rows,
        input_digest=problem.input_digest, request_sha256=sha(frozen_path),
        seconds=time.monotonic()-started, optimization_invocations=0, case_admission=False,
        scope='FORMAL_INPUT_AND_FIXED_STATE_MAPPING_ONLY_NOT_DERIVATIVE_OR_SOLVE_ADMISSION')
    (output / 'report.json').write_text(json.dumps(report))
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main(Path(sys.argv[1]))
