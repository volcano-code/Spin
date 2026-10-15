"""Replay one recorded failure; never initialize or restart an optimizer."""
import dataclasses
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
CASE = Path(__file__).parent / 'asymmetric-candidate-v1'
OLD_BUILD = ROOT / 'build/heart-exact-input-20261011-v2'
BUILD = ROOT / 'build/heart-fixed-facet-20261011-v1'


def save(path, value):
    with path.open('x') as out:
        json.dump(value, out, indent=2, allow_nan=False)
        out.flush()
        os.fsync(out.fileno())


def worker(mode, job):
    if mode == 'build':
        BUILD.mkdir(exist_ok=False)
        commands = json.loads((OLD_BUILD / 'commands.json').read_text())
        commands = [[arg.replace(str(OLD_BUILD), str(BUILD)) for arg in row]
                    for row in commands]
        save(BUILD / 'commands.json', commands)
        for command in commands:
            subprocess.run(command, check=True)
        save(job / 'build.json', {'status': 'DIAGNOSTIC_OVERLAY_BUILT',
             'behavior_change': 'EXACT_STATIONARY_FACET_MEMBERSHIP', 'commands': commands})
        return
    from make_it_stand.spin import (SpinItFasterRequest, SpinSolverConfig,
        SpinEvaluatorConfig, evaluate_spin_it_faster_fullrank_quadric)
    from make_it_stand.verification.spin_common import digest
    import _make_it_stand
    frozen = json.loads((CASE / 'frozen-requests-v3.json').read_text())
    expected = frozen.pop('frozen_payload_digest')
    assert digest(frozen) == expected
    for path, want in frozen['input_sources'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == want, path
    run = json.loads((CASE / 'solve-seed-0-v3/solver-result.json').read_text())
    event_index, event = next((i, e) for i, e in enumerate(run['events'])
        if 'birth facet membership' in e.get('detail', ''))
    payload = frozen['requests'][0]['request']
    request = SpinItFasterRequest(**{**payload,
        'settings': SpinSolverConfig(**payload['settings'])})
    save(job / 'replay-input.json', {'frozen_digest': expected,
        'recorded_event_index': event_index, 'recorded_event': event,
        'native': _make_it_stand.__file__, 'purpose': 'DIAGNOSTIC_ONLY',
        'optimizer_invoked': False})
    value = evaluate_spin_it_faster_fullrank_quadric(request,
        event['quadric_coefficients'],
        evaluator_settings=SpinEvaluatorConfig(**frozen['evaluator_settings']),
        diagnostic_stage_observer=lambda stage: print(stage, flush=True))
    result = dataclasses.asdict(value)
    save(job / 'evaluation.json', result)
    save(job / 'conclusion.json', {'derivative_defined': value.derivative_defined,
        'diagnostic': value.derivative_diagnostic,
        'scientific_acceptance': 'NOT_GRANTED', 'optimizer_invoked': False})
    print(value.derivative_diagnostic, flush=True)
    if mode == 'difference':
        # A single registered direction is a local diagnostic, not the complete
        # 81-point production protocol or a material/optimizer acceptance gate.
        import numpy as np
        q = np.asarray(event['quadric_coefficients'], dtype=float)
        q /= np.linalg.norm(q)
        direction = np.eye(10)[4] - q[4]*q
        direction /= np.linalg.norm(direction)
        scales = (1.e-4, 5.e-5, 2.5e-5)
        save(job / 'difference-protocol.json', {'direction': direction.tolist(),
            'scales': scales, 'scope': 'ONE_DIRECTION_FLOAT_GEOMETRY_DIAGNOSTIC',
            'not_complete_production_derivative_acceptance': True})
        assert value.derivative_defined, value.derivative_diagnostic

        def integrate(evaluation):
            v = np.asarray(evaluation.material_vertices)
            f = np.asarray(evaluation.material_faces)
            a,b,c = (v[f[:,i]] for i in range(3))
            volume = np.einsum('ij,ij->i', a, np.cross(b,c))/6.
            s = a+b+c
            first = np.sum(volume[:,None]*s/4., axis=0)
            second = np.einsum('i,ijk->jk', volume,
                (np.einsum('ij,ik->ijk',s,s)+np.einsum('ij,ik->ijk',a,a)+
                 np.einsum('ij,ik->ijk',b,b)+np.einsum('ij,ik->ijk',c,c))/20.)
            return request.material_density*np.array([volume.sum(),*first,
                *np.diag(second),second[1,2],second[0,2],second[0,1]])

        jacobian = np.asarray(value.coefficient_jacobian)
        checks = []
        for index,h in enumerate(scales):
            moments = []
            positions = []
            for sign in (-1,1):
                point = q+sign*h*direction
                point /= np.linalg.norm(point)
                sample = evaluate_spin_it_faster_fullrank_quadric(request, point,
                    evaluator_settings=SpinEvaluatorConfig(**frozen['evaluator_settings']))
                save(job / f'difference-{index}-{sign}.json', dataclasses.asdict(sample))
                moments.append(integrate(sample)); positions.append(point)
                print(f'difference scale={h} sign={sign} defined={sample.derivative_defined}', flush=True)
            measured = (moments[1]-moments[0])/(2*h)
            predicted = jacobian@((positions[1]-positions[0])/(2*h))
            errors = np.abs(measured-predicted)/np.maximum(1.,np.abs(predicted))
            checks.append({'scale':h,'measured':measured.tolist(),
                'predicted':predicted.tolist(),'errors':errors.tolist(),
                'maximum_error':float(errors.max())})
        save(job / 'difference-report.json', {'checks':checks,
            'all_errors_at_most_1e_minus_6': all(x['maximum_error']<=1e-6 for x in checks),
            'full_derivative_acceptance': 'NOT_GRANTED'})


def main():
    if sys.argv[1] == 'worker':
        worker(sys.argv[2], Path(sys.argv[3]))
        return
    from make_it_stand.verification.source_job import run_source_job
    mode = sys.argv[1]
    env = dict(PYTHONPATH=str(ROOT / 'python') + ':' + str(BUILD))
    identities = {'script': __file__,
        'source': str(ROOT / 'src/spin/PolyhedralMomentDifferential.cpp'),
        'old_commands': str(OLD_BUILD / 'commands.json'),
        'frozen': str(CASE / 'frozen-requests-v3.json'),
        'recorded_run': str(CASE / 'solve-seed-0-v3/solver-result.json')}
    if mode != 'build':
        identities['native'] = str(BUILD / '_make_it_stand.cpython-313-x86_64-linux-gnu.so')
    receipt = run_source_job(CASE / ('fixed-facet-' + mode + '-v1'),
        ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
         'OMP_NUM_THREADS=1', 'OPENBLAS_NUM_THREADS=1', 'MKL_NUM_THREADS=1',
         sys.executable, '-B', __file__, 'worker', mode, '{job}'],
        identities, environment=env)
    print(receipt, flush=True)


if __name__ == '__main__':
    main()
