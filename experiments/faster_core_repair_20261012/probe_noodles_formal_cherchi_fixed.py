"""One saved input/state through formal Cherchi; never grants case admission.

Uses the registered algebraic finite construction, not a new initial guess.
Independent saved-boundary integration and CGAL comparisons are separate facts.
"""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
import trimesh
from make_it_stand import spin
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator, _json_value
from independent_tetra_moments import integrate

HERE = Path(__file__).resolve().parent


def main(output):
    frozen_path = HERE / 'noodles-balanced-request-freeze-job-01/output/request.json'
    frozen = json.loads(frozen_path.read_text())
    for name, digest in frozen['input_sources'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != digest:
            raise ValueError('FROZEN_INPUT_CHANGED')
    prior = HERE / 'noodles-algebraic-production-moments-job-01/output'
    manifest = json.loads((prior/'manifest.json').read_text())
    if manifest['request_sha256'] != hashlib.sha256(frozen_path.read_bytes()).hexdigest():
        raise ValueError('DIAGNOSTIC_REQUEST_MISMATCH')
    settings = dict(frozen['evaluator_settings'], algebraic_finite_map=True)
    if settings != manifest['settings']:
        raise ValueError('REGISTERED_CONSTRUCTION_MISMATCH')
    payload = dict(frozen['requests'][0])
    payload['settings'] = spin.SpinSolverConfig(**payload['settings'])
    request = spin.SpinItFasterRequest(**payload)
    design_path = next(Path(p) for p in frozen['input_sources'] if Path(p).name=='design-si-candidate.obj')
    physical = trimesh.load(design_path, force='mesh', process=False)
    if not np.array_equal(physical.faces, request.design_faces):
        raise ValueError('PHYSICAL_FACE_IDENTITY_MISMATCH')
    output.mkdir(parents=True, exist_ok=False)
    q = np.asarray(manifest['samples'][0])
    (output/'registration.json').write_text(json.dumps(dict(
        request_sha256=manifest['request_sha256'], settings=settings,
        normalization=frozen['normalization'], coefficients=q.tolist(),
        case_admission=False, optimization_invocations=0)))
    start = time.monotonic()
    mapping = InProcessCherchiEvaluator(request, spin.SpinEvaluatorConfig(**settings),
        normalization=frozen['normalization'], design_vertices_physical=physical.vertices,
        rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1')
    mapped = mapping(q)
    path = output/'same-call-material.json'
    path.write_text(json.dumps(mapping.material(q), default=_json_value, allow_nan=False))
    material = json.loads(path.read_text())
    vertices = [[Fraction(x) for x in row] for row in material['variable_rational_vertices']]
    variable = integrate(vertices, material['variable_faces'], request.material_density)
    fixed = json.loads((frozen_path.parent/'fixed.binary64.json').read_text())
    # Saved fixed vertices already use normalized world coordinates; rotate and
    # translate them to the registered target frame independently of mapping.
    axes = np.asarray(request.target_axes)
    target = (np.asarray(fixed['vertices'])-np.asarray(request.target_origin)) @ axes.T
    fixed_moments = integrate([[Fraction(float(x)) for x in row] for row in target], fixed['faces'], 1.)
    independent_error = float(np.max(np.abs(variable+fixed_moments-mapped['raw_moments'])))
    baseline = json.loads((prior/'point-000-material.json').read_text())
    if not np.array_equal(baseline['coefficients'], q):
        raise ValueError('COMPARISON_STATE_MISMATCH')
    old = baseline['evaluation']
    comparison_m = float(np.max(np.abs(np.asarray(old['total_raw_moments'])-mapped['raw_moments'])))
    comparison_j = float(np.max(np.abs(np.asarray(old['coefficient_jacobian'])-mapped['moment_jacobian'])))
    report = dict(independent_error=independent_error,
        cgal_moment_difference=comparison_m, cgal_jacobian_difference=comparison_j,
        independent_passed=independent_error<=1e-10,
        numerical_backend_comparison_passed=comparison_m<=1e-10 and comparison_j<=1e-6,
        elapsed_seconds=time.monotonic()-start, case_admission=False,
        derivative_admission=False, exact_entity_comparison='NOT_PERFORMED',
        scope='ONE_REGISTERED_FIXED_STATE_NOT_OPTIMIZATION')
    (output/'report.json').write_text(json.dumps(report))
    if not report['independent_passed'] or not report['numerical_backend_comparison_passed']:
        raise ValueError('FIXED_STATE_COMPARISON_FAILED')


if __name__=='__main__':
    main(Path(sys.argv[1]))
