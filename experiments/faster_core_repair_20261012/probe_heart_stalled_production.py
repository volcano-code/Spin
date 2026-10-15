"""Registered 81-point production probe at Heart's actual accepted state 45."""
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from make_it_stand import spin
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator, _json_value
from make_it_stand.verification.cherchi_run_journal import RunJournal
from make_it_stand.verification.spin_problem import rational_solver_candidate_snapshot
from independent_tetra_moments import integrate
from probe_formal_cherchi_derivatives import weights

HERE = Path(__file__).resolve().parent


def main(output):
    run = HERE/'heart-guarded-formal-seed-0-job-01/output/run'
    config = json.loads((run/'configuration.json').read_text())
    events = [json.loads(row) for row in (run/'events.jsonl').read_text().splitlines()]
    ref = next(row for row in events if row.get('kind') == 'MAP_SAVED' and row['evaluation'] == 45)
    saved_bytes = (run/ref['evaluation_file']).read_bytes()
    if hashlib.sha256(saved_bytes).hexdigest() != ref['sha256']:
        raise ValueError('SAVED_STALLED_MATERIAL_CHANGED')
    saved = json.loads(saved_bytes)
    base = np.asarray(saved['quadric_coefficients'])
    samples, lookup, layouts = [base.copy()], {tuple(base): 0}, []
    for direction in range(10):
        for h in (3e-6, 1.5e-6, .75e-6):
            indices, positions = [], []
            for shift in (-h, -h/2, 0., h/2, h):
                q = base.copy(); q[direction] += shift
                key = tuple(q)
                if key not in lookup:
                    lookup[key] = len(samples); samples.append(q)
                indices.append(lookup[key]); positions.append(float(q[direction]-base[direction]))
            layouts.append(dict(direction=direction, scale=h, indices=indices, positions=positions))
    if len(samples) != 81:
        raise ValueError('REGISTERED_POINT_COUNT_CHANGED')
    frozen = json.loads((HERE.parent/'faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json').read_text())
    paths = {Path(path).name: Path(path) for path in frozen['input_sources']}
    for path, sha in frozen['input_sources'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
            raise ValueError('FROZEN_INPUT_CHANGED')
    outer, design, fixed = [json.loads(paths[name+'.binary64.json'].read_text())
                            for name in ('outer', 'design', 'fixed')]
    v = np.asarray(outer['vertices']); low, high = v.min(0), v.max(0)
    normalization = dict(origin=(low+(high-low)*.5).tolist(), L=float(max(high-low)))
    request_data = dict(config['request'])
    request_data['settings'] = spin.SpinSolverConfig(**request_data['settings'])
    request = spin.SpinItFasterRequest(**request_data)
    settings = spin.SpinEvaluatorConfig(**config['evaluator'])
    output.mkdir(parents=True, exist_ok=False)
    durations = Counter()
    with RunJournal(output, durability_profile='EVALUATION_CHECKPOINT_V1', json_value=_json_value) as journal:
        def stage(name, operation):
            start = time.monotonic_ns()
            journal.emit(dict(kind='STAGE_STARTED', stage=name, monotonic_ns=start))
            try:
                value = operation()
            except Exception as error:
                journal.emit(dict(kind='STAGE_FAILED', stage=name,
                    elapsed_ns=time.monotonic_ns()-start, error=f'{type(error).__name__}: {error}'))
                raise
            elapsed = time.monotonic_ns()-start
            durations[name] += elapsed
            journal.emit(dict(kind='STAGE_RETURNED', stage=name, elapsed_ns=elapsed))
            return value
        journal.publish('manifest.json', dict(schema='heart-stalled-production-81-v1', evaluation=45,
            original_material_sha256=ref['sha256'], configuration_sha256=hashlib.sha256(
                (run/'configuration.json').read_bytes()).hexdigest(), samples=[q.tolist() for q in samples],
            layouts=layouts, request=config['request'], evaluator=config['evaluator'],
            optimization_invocations=0, case_complete=False, branch_certification=False))
        producer = stage('prepare_same_input_producer', lambda: InProcessCherchiEvaluator(
            request, settings, normalization=normalization, design_vertices_physical=design['vertices'],
            fixed_boundary_reuse=True, rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1',
            stage_observer=stage))
        if producer.fixed_boundary_proof.reference_identity != config['fixed_boundary_proof_identity']:
            raise ValueError('FIXED_BOUNDARY_PROOF_CHANGED')
        exact_fixed = tuple(tuple(Fraction(float(x)) for x in row) for row in fixed['vertices'])
        target_fixed = producer._target_view(exact_fixed)[1]
        independent_fixed = stage('independent_fixed_integral', lambda:
            integrate(target_fixed, fixed['faces'], request.material_density))
        if np.max(np.abs(independent_fixed-producer.fixed)) > 1e-10:
            raise ValueError('FIXED_MATERIAL_MOMENTS_CHANGED')
        moments, identities, J = [], [], None
        for index, q in enumerate(samples):
            mapped = stage('production_geometry_evaluation', lambda: producer(q))
            material = producer.material(q)
            if index == 0:
                for field in ('variable_rational_vertices', 'variable_faces', 'moment_jacobian',
                              'quadric_source_identity', 'geometry_hash', 'topology_hash'):
                    if material[field] != saved[field]:
                        raise ValueError('STALLED_BASE_MAPPING_CHANGED:'+field)
                J = np.asarray(mapped['moment_jacobian'])
            journal.publish(f'point-{index:03d}-material.json', material)
            reloaded = json.loads((output/f'point-{index:03d}-material.json').read_text())
            snapshot, mode = stage('exact_boundary_classification', lambda:
                rational_solver_candidate_snapshot(reloaded['variable_rational_vertices'],
                    reloaded['variable_faces'], fixed_boundary_proof=producer.fixed_boundary_proof,
                    allow_signed_classification=True))
            independent = stage('independent_saved_material_integral', lambda:
                integrate([[Fraction(x) for x in row] for row in reloaded['variable_rational_vertices']],
                    reloaded['variable_faces'], request.material_density)+independent_fixed)
            error = float(np.max(np.abs(independent-np.asarray(mapped['raw_moments']))))
            if error > 1e-10:
                raise ValueError('INDEPENDENT_MATERIAL_INTEGRAL_CHANGED')
            moments.append(independent)
            identities.append(mapped['topology_hash'])
            journal.emit(dict(kind='POINT_RETURNED', index=index, topology_hash=identities[-1],
                independent_error=error, boundary_mode=mode))
            print('POINT_RETURNED', index, error, flush=True)
        checks = []
        for layout in layouts:
            numerical = weights(layout['positions'])@np.asarray([moments[i] for i in layout['indices']])
            error = float(np.max(np.abs(numerical-J[:, layout['direction']])))
            checks.append(dict(**layout, error=error, passed=error <= 1e-6))
        passed = all(check['passed'] for check in checks)
        journal.publish('report.json', dict(schema='heart-stalled-production-81-v1', numerical_passed=passed,
            checks=checks, point_count=len(samples), topology_labels=len(set(identities)),
            branch_certification=False, case_complete=False, optimization_invocations=0,
            stage_ns=dict(durations),
            scope='Actual accepted state45; finite numerical sample, not a neighborhood or branch certificate'))
        if not passed:
            raise ValueError('STALLED_PRODUCTION_FD_FAILED_NO_STEP_RETUNING')


if __name__ == '__main__':
    main(Path(sys.argv[1]))
