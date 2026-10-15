"""Registered seed-only Noodles CGAL run, using existing finite-sample evidence.

This composes input, derivative, source-chain and integration evidence. It does
not turn their finite sample coverage into a global neighborhood proof, nor
does numerical convergence authorize final material, OBJ or dynamics.
"""
import hashlib
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
NATIVE = HERE / 'cgal-proof-capture-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
PATHS = {
    'frozen': 'noodles-balanced-request-freeze-job-01/output/request.json',
    'manifest': 'noodles-algebraic-production-moments-job-01/output/manifest.json',
    'numerical': 'noodles-algebraic-production-moments-job-01/output/report.json',
    'source': 'noodles-interval-algebraic-source-derivative-job-01/output/report.json',
    'source_manifest': 'noodles-interval-algebraic-source-derivative-job-01/output/manifest.json',
    'branch': 'noodles-native-tape-certification-job-01/output/report.json',
    'closure': 'small-cgal-proof-closure-job-01/output/report.json',
    'input': 'noodles-balanced-shared-export-job-02/output/manifest.json',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_evidence(here):
    paths = {key: here / path for key, path in PATHS.items()}
    evidence = {key: json.loads(path.read_text()) for key, path in paths.items()}
    for path, expected in evidence['frozen']['input_sources'].items():
        if sha(path) != expected:
            raise ValueError('FROZEN_INPUT_IDENTITY_CHANGED: ' + path)
    evidence.update(request_sha256=sha(paths['frozen']), manifest_sha256=sha(paths['manifest']),
                    source_sha256=sha(paths['source']), source_manifest_sha256=sha(paths['source_manifest']),
                    native_sha256=sha(NATIVE), evidence_sha256={k: sha(p) for k, p in paths.items()})
    return evidence


def check_evidence(e, seed):
    if seed not in (0, 1):
        raise ValueError('REGISTERED_SEED_REQUIRED')
    if (e['manifest']['request_sha256'] != e['request_sha256'] or
            e['branch']['manifest_sha256'] != e['manifest_sha256'] or
            e['manifest']['source_gate_sha256'] != e['source_sha256'] or
            e['manifest']['source_manifest_sha256'] != e['source_manifest_sha256']):
        raise ValueError('EVIDENCE_IDENTITY_MISMATCH')
    if e['branch']['native_sha256'] != e['native_sha256']:
        raise ValueError('CERTIFIED_BUILD_MISMATCH')
    if (e['branch']['complete_frozen_sample_certification'] is not True or
            e['branch']['passed_indices'] != list(range(81))):
        raise ValueError('COMPLETE_FINITE81_CERTIFICATE_REQUIRED')
    if (e['numerical']['numerical_checks_passed'] is not True or
            e['numerical']['independent_integration_passed'] is not True or
            e['numerical']['point_count'] != 81 or len(e['numerical']['checks']) != 30 or
            not all(x['passed'] for x in e['numerical']['checks']) or
            e['source']['source_checks_passed'] is not True):
        raise ValueError('PRODUCTION_NUMERICAL_GATES_REQUIRED')
    if e['closure']['closure_passed'] is not True or e['input']['derived_input_partition_passed'] is not True:
        raise ValueError('INPUT_AND_SAME_CALL_CLOSURE_REQUIRED')
    rows = [r for r in e['frozen']['requests'] if r['settings']['deterministic_seed'] == seed]
    if len(rows) != 1:
        raise ValueError('UNIQUE_REGISTERED_REQUEST_REQUIRED')
    row = rows[0]
    required = dict(maximum_iterations=60, maximum_restarts=8,
                    maximum_function_evaluations_per_restart=60, equation_tolerance=1e-12)
    if any(row['settings'][k] != v for k, v in required.items()):
        raise ValueError('FROZEN_FULL_SOLVE_BUDGET_MISMATCH')
    evaluator = dict(e['frozen']['evaluator_settings'], algebraic_finite_map=True)
    if evaluator != e['manifest']['settings']:
        raise ValueError('EVALUATOR_IDENTITY_MISMATCH')
    if e['frozen']['alpha'] != .5 or e['frozen']['t_max'] != .5:
        raise ValueError('FROZEN_NEWTON_PROTOCOL_MISMATCH')
    return dict(request=row, evaluator=evaluator, evidence_sha256=e['evidence_sha256'],
                native_sha256=e['native_sha256'], initialization='FROZEN_REQUEST_SEED_ONLY',
                finite_sample_branch_certification=True, global_neighborhood_certification=False,
                status='ADMITTED_REGISTERED_CORE_RUN_NOT_FINAL_CASE_SUCCESS')


def main(output, seed):
    from make_it_stand import spin
    sys.path.insert(0, str(HERE.parent / 'faster_heart_contact_20261011'))
    from rounded_solver_journal import run_durable_fullrank, save
    start = time.monotonic()
    config = check_evidence(load_evidence(HERE), seed)
    if sha(spin._native_module().__file__) != config['native_sha256']:
        raise ValueError('LOADED_BUILD_MISMATCH')
    output.mkdir(parents=True, exist_ok=False)
    save(output / 'configuration.json', config)
    payload = dict(config['request'])
    payload['settings'] = spin.SpinSolverConfig(**payload['settings'])
    request = spin.SpinItFasterRequest(**payload)
    evaluator = spin.SpinEvaluatorConfig(**config['evaluator'])
    setup_seconds = time.monotonic() - start
    result = run_durable_fullrank(request, evaluator, output / 'run', capture_material_proof=True)
    updates = sum(float(t.get('accepted_geodesic_step', 0)) > 0 for t in result.get('trace', []))
    report = dict(status=result['status'], seed=seed, function_evaluations=result['function_evaluations'],
                  raw_equation_residual=result.get('raw_equation_residual'), accepted_updates=updates,
                  setup_seconds=setup_seconds, solve_seconds=time.monotonic()-start-setup_seconds,
                  final_material_acceptance='PENDING', repeat_acceptance='PENDING',
                  obj_acceptance='PENDING', dynamics_acceptance='PENDING', complete_case=False)
    save(output / 'report.json', report)
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main(Path(sys.argv[1]), int(sys.argv[2]))
