"""Four identical boundary audits at the actual stalled Heart state.

No Boolean, source reconstruction or Newton. Recreates the original fixed
proof in the exact physical-to-target frame and verifies its saved identity.
"""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import _make_it_stand as native

HERE = Path(__file__).resolve().parent


def main(output, mode='phases'):
    output.mkdir(parents=True, exist_ok=False)
    run = HERE/'heart-guarded-formal-seed-0-job-01/output/run'
    config = json.loads((run/'configuration.json').read_text())
    records = [json.loads(line) for line in (run/'events.jsonl').read_text().splitlines()]
    evaluation = 2 if mode == 'classify' else 45
    reference = next(row for row in records if row.get('kind') == 'MAP_SAVED' and row['evaluation'] == evaluation)
    material_path = run/reference['evaluation_file']
    if hashlib.sha256(material_path.read_bytes()).hexdigest() != reference['sha256']:
        raise ValueError('SAVED_MATERIAL_CHANGED')
    material = json.loads(material_path.read_text())
    frozen_path = HERE.parent/'faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json'
    frozen = json.loads(frozen_path.read_text())
    input_paths = {Path(path).name: Path(path) for path in frozen['input_sources']}
    outer = json.loads(input_paths['outer.binary64.json'].read_text())
    design = json.loads(input_paths['design.binary64.json'].read_text())
    vertices = np.asarray(outer['vertices'])
    low, high = vertices.min(axis=0), vertices.max(axis=0)
    spans = high-low
    origin = [Fraction(float(x)) for x in low+spans*.5]
    length = Fraction(float(max(spans)))
    request = config['request']
    if float(length) != request['original_units_per_normalized_unit']:
        raise ValueError('NORMALIZATION_LENGTH_CHANGED')
    target = [Fraction(float(x)) for x in request['target_origin']]
    axes = [[Fraction(float(x)) for x in row] for row in request['target_axes']]
    rows = []
    for vertex in design['vertices']:
        model = [(Fraction(float(vertex[i]))-origin[i])/length-target[i] for i in range(3)]
        exact = [sum(model[i]*axes[i][j] for i in range(3)) for j in range(3)]
        rows.append([f'{value.numerator}/{value.denominator}' for value in exact])
    start = time.monotonic()
    proof = native.verify_fixed_rational_boundary(rows, np.asarray(design['faces'], dtype=np.int32))
    preparation = time.monotonic()-start
    if proof.reference_identity != config['fixed_boundary_proof_identity']:
        raise ValueError('ORIGINAL_FIXED_PROOF_IDENTITY_CHANGED')
    start = time.monotonic()
    from make_it_stand.verification.spin_problem import _canonical_rational_vertices, _mesh_input_array
    exact = _canonical_rational_vertices(material['variable_rational_vertices'])
    faces = _mesh_input_array(material['variable_faces'], indices=True)
    conversion = time.monotonic()-start
    if mode == 'classify':
        from make_it_stand.verification.spin_problem import rational_solver_candidate_snapshot
        from make_it_stand.verification.faster_material import _raw_moments
        start = time.monotonic()
        baseline, baseline_mode = rational_solver_candidate_snapshot(exact, faces, fixed_boundary_proof=proof)
        baseline_seconds = time.monotonic()-start
        facts = baseline.exact_region().facts
        (output/'baseline-facts.json').write_text(json.dumps(dict(mode=baseline_mode,
            seconds=baseline_seconds, facts=facts), allow_nan=False))
        samples = []
        for index in range(4):
            start = time.monotonic()
            result = native.audit_exact_rational_material_boundary(exact, faces, proof)
            elapsed = time.monotonic()-start
            (output/f'classification-{index}.json').write_text(json.dumps(result, allow_nan=False))
            checks = {key: result[key] == facts[key] for key in
                ('parents', 'depths', 'original_signs', 'material_components', 'bounded_void_components')}
            checks['actual_boundary_digest'] = result['boundary']['input_digest'] == facts['input_digest']
            if not result['completed'] or not result['classification_valid'] or not all(checks.values()):
                raise ValueError('EXACT_POINT_LOCATION_DISAGREES_WITH_NEF:'+json.dumps(checks))
            samples.append(dict(sample=index, wall_seconds=elapsed, checks=checks, performance=result['performance']))
        moments = _raw_moments(baseline, request['material_density'])
        error = float(np.max(np.abs(np.asarray(moments)-material['variable_raw_moments'])))
        if error > 1e-10:
            raise ValueError('INDEPENDENT_BASELINE_MOMENT_MISMATCH')
        report = dict(schema='heart-initial-signed-boundary-classification-comparison-v1', evaluation=evaluation,
            material_sha256=reference['sha256'], baseline_seconds=baseline_seconds,
            samples=samples, independent_moment_error=error, passed=True,
            production_profile_admission=False, derivative_certification=False, case_complete=False,
            scope='One actual registered initial boundary compared with complete exact Nef classification')
        (output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False))
        print(json.dumps(dict(baseline_seconds=baseline_seconds, candidate_seconds=[s['wall_seconds'] for s in samples],
                             independent_moment_error=error)), flush=True)
        return
    if mode == 'compare':
        from make_it_stand.verification.spin_problem import rational_solver_candidate_snapshot
        from make_it_stand.verification.faster_material import _raw_moments
        start = time.monotonic()
        baseline, baseline_mode = rational_solver_candidate_snapshot(exact, faces, fixed_boundary_proof=proof)
        baseline_seconds = time.monotonic()-start
        (output/'baseline-facts.json').write_text(json.dumps(dict(mode=baseline_mode,
            seconds=baseline_seconds, facts=baseline.exact_region().facts), allow_nan=False))
        samples = []
        for index in range(4):
            start = time.monotonic()
            candidate, candidate_mode = rational_solver_candidate_snapshot(exact, faces,
                fixed_boundary_proof=proof, allow_separated_roots=True)
            elapsed = time.monotonic()-start
            equal = (candidate.rational_vertices == baseline.rational_vertices and
                     all(a.tobytes() == b.tobytes() for a, b in zip(candidate.arrays(), baseline.arrays())))
            hierarchy = (candidate.separation_proof['parents'] == baseline.exact_region().facts['parents'])
            if not equal or not hierarchy or candidate_mode != 'EXACT_RATIONAL_SEPARATED_OUTWARD_BOUNDARIES':
                raise ValueError('OPTIMIZED_BOUNDARY_OR_REGION_CLASSIFICATION_CHANGED')
            samples.append(dict(sample=index, wall_seconds=elapsed, boundary_equal=equal,
                                exact_region_hierarchy_equal=hierarchy,
                                separation_proof=candidate.separation_proof))
            (output/f'comparison-sample-{index}.json').write_text(json.dumps(samples[-1], allow_nan=False))
        start = time.monotonic()
        moments = _raw_moments(candidate, request['material_density'])
        independent_seconds = time.monotonic()-start
        moment_error = float(np.max(np.abs(np.asarray(moments)-material['variable_raw_moments'])))
        if moment_error > 1e-10:
            raise ValueError('INDEPENDENT_MATERIAL_MOMENTS_CHANGED')
        report = dict(schema='heart-multiple-root-snapshot-comparison-v1',
            material_sha256=reference['sha256'], fixed_proof_identity=proof.reference_identity,
            baseline_seconds=baseline_seconds, samples=samples,
            moment_error=moment_error, independent_integration_seconds=independent_seconds,
            actual_matrix_or_jacobian_changed=False, optimization_invocations=0,
            formal_profile_admission=False, case_complete=False)
        (output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False))
        print(json.dumps(dict(baseline_seconds=baseline_seconds,
            candidate_seconds=[x['wall_seconds'] for x in samples], moment_error=moment_error)), flush=True)
        return
    if mode != 'phases':
        raise ValueError('UNREGISTERED_BOUNDARY_PROFILE_MODE')
    samples = []
    baseline = None
    for index in range(4):
        start = time.monotonic()
        facts = native.audit_rational_boundary_with_fixed_proof(exact, faces, 'AnyNonzero', proof, True)
        seconds = time.monotonic()-start
        (output/f'raw-boundary-facts-{index}.json').write_text(json.dumps(facts, allow_nan=False))
        if not facts['boundary_valid']:
            raise ValueError('REGISTERED_BOUNDARY_INVALID')
        science = {key: value for key, value in facts.items() if key != 'performance'}
        if baseline is None:
            baseline = science
        elif science != baseline:
            raise ValueError('REPEATED_BOUNDARY_AUDIT_CHANGED_FACTS')
        samples.append(dict(sample=index, process_cold=index == 0, wall_seconds=seconds,
                            performance=facts['performance']))
        # Persist after every actual return, including if a later sample fails.
        (output/f'sample-{index}.json').write_text(json.dumps(samples[-1], allow_nan=False))
    report = dict(schema='heart-saved-boundary-phase-profile-v1', evaluation=45,
        material_sha256=reference['sha256'], fixed_proof_identity=proof.reference_identity,
        native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
        fixed_preparation_seconds=preparation, python_canonicalization_seconds=conversion,
        samples=samples, scientific_facts=baseline, optimization_invocations=0,
        production_snapshot_route=('EXACT_RATIONAL_SINGLE_OUTWARD_BOUNDARY'
            if tuple(baseline['component_volume_signs']) == (1,) else 'EXACT_NEF_SIGNED_REGION'),
        derivative_certification=False, case_complete=False)
    (output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False))
    print(json.dumps(dict(seconds=[s['wall_seconds'] for s in samples],
                         last_phases=samples[-1]['performance'])), flush=True)


if __name__ == '__main__':
    main(Path(sys.argv[1]), sys.argv[2] if len(sys.argv) > 2 else 'phases')
