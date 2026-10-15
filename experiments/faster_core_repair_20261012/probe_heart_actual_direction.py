"""Actual Cherchi materials on preregistered saved Heart geodesics, no Newton."""
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
from audit_heart_stalled_prefix import coordinate, symbolic_model
from audit_heart_trial_prediction import residual
from independent_tetra_moments import integrate

HERE = Path(__file__).resolve().parent


def main(output):
    registration_path = HERE/('heart-actual-direction-probe-registration-05.json' if '--local-suffix' in sys.argv else
                             'heart-actual-direction-probe-registration-04.json'
                             if '--backtrack-suffix' in sys.argv else 'heart-actual-direction-probe-registration-03.json')
    registration = json.loads(registration_path.read_text())
    run = HERE/'heart-guarded-formal-seed-0-job-01/output/run'
    for name in ('events', 'configuration'):
        path = run/(name+'.jsonl' if name == 'events' else name+'.json')
        if hashlib.sha256(path.read_bytes()).hexdigest() != registration[name+'_sha256']:
            raise ValueError('REGISTERED_PREFIX_CHANGED')
    config = json.loads((run/'configuration.json').read_text())
    frozen = json.loads((HERE.parent/'faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json').read_text())
    files = {Path(p).name: Path(p) for p in frozen['input_sources']}
    for path, sha in frozen['input_sources'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
            raise ValueError('FROZEN_INPUT_CHANGED')
    outer, design = [json.loads(files[n+'.binary64.json'].read_text()) for n in ('outer', 'design')]
    v = np.asarray(outer['vertices']); low, high = v.min(0), v.max(0)
    normalization = dict(origin=(low+(high-low)*.5).tolist(), L=float(max(high-low)))
    data = dict(config['request']); data['settings'] = spin.SpinSolverConfig(**data['settings'])
    request = spin.SpinItFasterRequest(**data)
    settings = spin.SpinEvaluatorConfig(**config['evaluator'])
    output.mkdir(parents=True, exist_ok=False)
    with RunJournal(output, durability_profile='EVALUATION_CHECKPOINT_V1', json_value=_json_value) as journal:
        def stage(name, operation):
            start = time.monotonic_ns()
            journal.emit(dict(kind='STAGE_STARTED', stage=name, monotonic_ns=start))
            try:
                value = operation()
            except Exception as error:
                journal.emit(dict(kind='STAGE_FAILED', stage=name, error=str(error)))
                raise
            journal.emit(dict(kind='STAGE_RETURNED', stage=name, elapsed_ns=time.monotonic_ns()-start))
            return value
        journal.publish('registration.json', registration)
        producer = stage('prepare_producer', lambda: InProcessCherchiEvaluator(
            request, settings, normalization=normalization,
            design_vertices_physical=design['vertices'], fixed_boundary_reuse=True,
            rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1', stage_observer=stage,
            integral_profile=('CERTIFIED_MPFR256_GMP_FALLBACK_V1' if '--certified-integral' in sys.argv else 'EXACT_GMP_V1')))
        if producer.fixed_boundary_proof.reference_identity != config['fixed_boundary_proof_identity']:
            raise ValueError('FIXED_BOUNDARY_IDENTITY_CHANGED')
        model = symbolic_model(request.center_of_mass_weight)
        benchmark_args = [x for x in sys.argv if x.startswith('--benchmark-state=')]
        if benchmark_args:
            index = int(benchmark_args[0].split('=')[1])
            if index not in (2,38,45,46): raise ValueError('UNREGISTERED_BENCHMARK_STATE')
            path = run/f'evaluation-{index:06d}.json'
            reference = next(r for r in json.loads('['+','.join((run/'events.jsonl').read_text().splitlines())+']')
                             if r.get('kind') == 'MAP_SAVED' and r['evaluation'] == index)
            if hashlib.sha256(path.read_bytes()).hexdigest() != reference['sha256']:
                raise ValueError('BENCHMARK_MATERIAL_CHANGED')
            saved = json.loads(path.read_text()); q = saved['quadric_coefficients']; samples=[]
            for repeat in range(4):
                started=time.monotonic()
                mapped=stage('benchmark_evaluation',lambda:producer(q))
                elapsed=time.monotonic()-started; material=producer.material(q)
                for field in ('variable_rational_vertices','variable_faces','moment_jacobian',
                              'quadric_source_identity','geometry_hash','topology_hash','raw_moments'):
                    if material[field] != saved[field]:raise ValueError('BENCHMARK_MAPPING_CHANGED:'+field)
                if repeat == 0:
                    journal.publish('material.json',material)
                    independent=stage('independent_variable_integral',lambda:integrate(
                        [[Fraction(x) for x in row] for row in material['variable_rational_vertices']],
                        material['variable_faces'],request.material_density))
                    error=float(np.max(np.abs(independent+producer.fixed-np.asarray(mapped['raw_moments']))))
                    if error>1e-10:raise ValueError('BENCHMARK_INDEPENDENT_MOMENTS_CHANGED')
                samples.append(dict(repeat=repeat,elapsed_seconds=elapsed))
                print('BENCHMARK_POINT_RETURNED',index,repeat,elapsed,flush=True)
            journal.publish('report.json',dict(schema='heart-fixed-state-arithmetic-benchmark-v1',
                evaluation=index,samples=samples,independent_variable_error=error,
                arithmetic_identity=producer.arithmetic_identity,material_sha256=reference['sha256'],
                full_saved_mapping_equal=True,optimization_invocations=0,case_complete=False))
            return
        events = [json.loads(line) for line in (run/'events.jsonl').read_text().splitlines()]
        base_events = {r['event']['evaluation_index']: r['event'] for r in events
                       if r.get('kind') == 'NEWTON_EVENT' and r['event'].get('type') == 'evaluation'
                       and r['event'].get('evaluation_payload')}
        bases, results, base_branches = {}, [], {}
        def branch():
            source = producer.last['source']
            return json.loads(json.dumps(dict(surface_type=source['surface_type'],
                cap_face_count=source['cap_face_count'], interior_is_positive=source['interior_is_positive'],
                operation=producer.last['operation'], source_normal_edges=source['source_normal_edges']),
                default=_json_value, allow_nan=False))
        for index in sorted({p['base_evaluation'] for p in registration['probes']}):
            path = run/f'evaluation-{index:06d}.json'
            sha = next(p['base_material_sha256'] for p in registration['probes'] if p['base_evaluation']==index)
            if hashlib.sha256(path.read_bytes()).hexdigest() != sha:
                raise ValueError('BASE_MATERIAL_CHANGED')
            saved = json.loads(path.read_text()); q = saved['quadric_coefficients']
            mapped = stage('replay_base', lambda: producer(q))
            material = producer.material(q)
            for field in ('variable_rational_vertices', 'variable_faces', 'moment_jacobian',
                          'quadric_source_identity', 'geometry_hash', 'topology_hash'):
                if material[field] != saved[field]:
                    raise ValueError('BASE_MAPPING_CHANGED:'+field)
            raw = np.asarray(mapped['raw_moments'])
            if not np.allclose(raw, base_events[index]['evaluation_payload']['solver_raw_moments'],
                               rtol=1e-12, atol=1e-14):
                raise ValueError('SOLVER_MOMENT_COORDINATES_DIFFER')
            bases[index] = (coordinate(q), raw, np.asarray(mapped['moment_jacobian']), np.asarray(q))
            base_branches[index] = branch()
        for i, probe in enumerate(registration['probes']):
            q = probe['coefficients']; p = coordinate(q)
            mapped = stage('direction_material', lambda: producer(q))
            material = producer.material(q)
            if probe['existing_trial']:
                path = run/f"evaluation-{probe['original_trial']:06d}.json"
                if hashlib.sha256(path.read_bytes()).hexdigest() != probe['original_trial_material_sha256']:
                    raise ValueError('ORIGINAL_TRIAL_CHECKPOINT_CHANGED')
                original = json.loads(path.read_text())
                for field in ('quadric_coefficients', 'variable_rational_vertices', 'variable_faces',
                              'moment_jacobian', 'quadric_source_identity', 'geometry_hash', 'topology_hash'):
                    if json.loads(json.dumps(material[field], default=_json_value, allow_nan=False)) != original[field]:
                        raise ValueError('ORIGINAL_TRIAL_MAPPING_CHANGED:'+field)
            journal.publish(f'point-{i:03d}-material.json', material)
            independent_variable = stage('independent_variable_integral', lambda: integrate(
                [[Fraction(x) for x in row] for row in material['variable_rational_vertices']],
                material['variable_faces'], request.material_density))
            actual = np.asarray(mapped['raw_moments'])
            integral_error = float(np.max(np.abs(independent_variable+producer.fixed-actual)))
            if integral_error > 1e-10:
                raise ValueError('INDEPENDENT_VARIABLE_MOMENTS_MISMATCH')
            bp, br, J, bq = bases[probe['base_evaluation']]
            change, prediction = actual-br, J@(np.asarray(q)-bq)
            error = np.linalg.norm(change-prediction)
            row = dict(probe=probe, moment_linearization_error=float(error),
                relative_linearization_error=float(error/max(np.linalg.norm(change), np.linalg.norm(prediction), 1e-300)),
                actual_residual=residual(model, actual, p),
                base_residual=residual(model, br, bp), independent_variable_error=integral_error,
                base_source_branch=base_branches[probe['base_evaluation']], source_branch=branch(),
                derivative_defined=mapped['derivative_defined'],
                derivative_identity=mapped['derivative_identity'], geometry_hash=mapped['geometry_hash'])
            results.append(row); journal.emit(dict(kind='DIRECTION_POINT_RETURNED', **row))
            print('DIRECTION_POINT_RETURNED', i, row['actual_residual'], flush=True)
        journal.publish('report.json', dict(schema='heart-actual-direction-probe-v1', results=results,
            optimization_invocations=0, derivative_admission=False, case_complete=False,
            scope='Registered finite secants and independent variable integrals; fixed moments reused from matched producer'))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
