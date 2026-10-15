"""Bind existing finite Heart evidence to a branch-protected formal attempt.

Preserves the failed original five-point check and the absence of neighborhood
certification. No source construction or Boolean is performed here.
"""
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def assess_finite_evidence(*, preparation, cap_numeric, cap_classification,
        small_numeric, small_classification, reference, cap_events, small_events):
    def require(value, message):
        if not value:
            raise ValueError(message)
    require(preparation.get('passed') is True and
            preparation.get('request_comparison', {}).get('equivalent') is True and
            len(preparation.get('rows', [])) == 2 and all(
                row['material_geometry_equal'] and row['jacobian_error'] <= 1e-10 and
                row['moment_error'] <= 1e-10 for row in preparation['rows']),
            'FORMAL_PREPARATION_NOT_VALID')
    require(cap_numeric.get('point_count') == 81 and cap_numeric.get('numerical_passed') is True and
            cap_numeric['independent_error'] <= 1e-10 and all(
                row['passed'] and row['error'] <= 1e-6 for row in cap_numeric['checks']),
            'CAP_DERIVATIVE_EVIDENCE_INCOMPLETE')
    require(small_numeric.get('point_count') == 81 and small_numeric.get('numerical_passed') is False and
            small_numeric['independent_error'] <= 1e-10, 'SMALL_FACE_ORIGINAL_GATE_CHANGED')
    failed = [row for row in small_numeric['checks'] if not row['passed']]
    require(bool(failed) and all(row['error'] <= 1e-6 for row in small_numeric['checks'] if row['passed']),
            'INCONSISTENT_ORIGINAL_FD_CHECKS')
    require(reference.get('status') == 'SUPPORTED' and reference.get('wrong_jacobian_status') == 'INCONSISTENT' and
            math.isfinite(float(reference['maximum_derivative_error'])) and
            float(reference['maximum_derivative_error']) <= 1e-6 and
            reference.get('old_five_point_gate_passed') is False and
            reference.get('neighborhood_certification') is False and
            all(row['direction'] == reference['direction'] for row in failed),
            'FAILED_DIRECTION_NOT_COVERED_BY_POINTWISE_REFERENCE')
    def classify(report, events):
        rows = report.get('rows', [])
        require(report.get('finite_classification_complete') is True and report.get('point_count') == 81 and
                len(rows) == 81 and sorted(row['index'] for row in rows) == list(range(81)),
                'INCOMPLETE_FINITE_CLASSIFICATION')
        observed = {row['index']: row for row in events if row.get('status') == 'RETURNED'}
        require(len(observed) == 81 and sorted(observed) == list(range(81)), 'MISSING_PRODUCTION_SAMPLE')
        allowed = {'IDENTICAL_OR_RENUMBERED', 'EXACT_COPLANAR_RETRIANGULATION', 'REAL_CONSTRUCTION_BRANCH_CHANGE'}
        for row in rows:
            require(row.get('status') == 'FINITE_SAMPLE_CLASSIFIED' and row.get('own_construction_replayed') is True and
                    row.get('classification') in allowed, 'UNSUPPORTED_FINITE_SAMPLE')
            if row['classification'] == 'REAL_CONSTRUCTION_BRANCH_CHANGE':
                require(observed[row['index']]['topology_hash'] != observed[0]['topology_hash'],
                        'REAL_BRANCH_CHANGE_NOT_DETECTED')
        return dict(Counter(row['classification'] for row in rows))
    return dict(finite_evidence_ready=True, cap_classifications=classify(cap_classification, cap_events),
        small_classifications=classify(small_classification, small_events),
        failed_original_checks=[dict(direction=row['direction'], scale=row['scale'], error=row['error']) for row in failed],
        old_five_point_gate_passed=False, neighborhood_certification=False, case_complete=False,
        scope='Finite evidence for guarded solve attempts; no neighborhood smoothness or case completion')


def main(output, *, certified=False):
    from make_it_stand import spin
    load = lambda path: json.loads(path.read_text())
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    jobs = dict(preparation='heart-formal-mapping-prepare-job-01',
        cap_numeric='heart-cherchi-heart-first-cap-failure-derivative-job-01',
        cap_classification='heart-cherchi-heart-first-cap-failure-connection-job-02',
        small_numeric='heart-cherchi-heart-first-small-face-rejection-derivative-job-01',
        small_classification='heart-small-face-finite-classification-job-01',
        reference='heart-projective-independent-reference-job-01')
    bindings = {key: HERE/name/'output/report.json' for key, name in jobs.items()}
    data = {key: load(path) for key, path in bindings.items()}
    for key, name in jobs.items():
        result_path = HERE/name/'result.json'
        config_path = HERE/name/'config.json'
        result, config = load(result_path), load(config_path)
        expected_status = ('PROCESS_FAILED', 1) if key == 'small_numeric' else ('PROCESS_EXITED', 0)
        if (result['status'], result['returncode']) != expected_status or result['unit'] != config['unit']:
            raise ValueError('ORIGINAL_EVIDENCE_RESOURCE_STATUS_CHANGED')
        bindings[key+'_resource_result'] = result_path
        bindings[key+'_resource_config'] = config_path
    for label in ('cap', 'small'):
        path = HERE/jobs[label+'_numeric']/'output/events.jsonl'
        bindings[label+'_events'] = path
        data[label+'_events'] = [json.loads(line) for line in path.read_text().splitlines()]
        for row in data[label+'_classification']['rows']:
            path = HERE/jobs[label+'_numeric']/f"output/point-{row['index']:03d}-material.json"
            if sha(path) != row['material_sha256']:
                raise ValueError('FINITE_MATERIAL_ATTACHMENT_CHANGED')
            bindings[f"{label}_material_{row['index']:03d}"] = path
    assessment = assess_finite_evidence(**data)
    frozen = HERE.parent/'faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json'
    clearance = HERE/'heart-derived-clearance-transfer-01.json'
    bindings.update(frozen_request=frozen, clearance=clearance)
    if (data['preparation']['request_sha256'] != sha(frozen) or
            load(clearance)['request_sha256'] != sha(frozen) or not load(clearance)['at_least_required']):
        raise ValueError('FROZEN_HEART_CLEARANCE_NOT_BOUND')
    for index, path in enumerate(load(frozen)['input_sources']):
        if sha(Path(path)) != load(frozen)['input_sources'][path]:
            raise ValueError('FROZEN_HEART_INPUT_CHANGED')
        bindings['input_'+str(index)] = Path(path)
    test = ROOT/'tests/python/test_cherchi_branch_relinearization.py'
    native = Path(spin._native_module().__file__).resolve()
    expected_native = HERE/('formal-native-core-build-job-12/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so' if certified else
                           'formal-native-core-build-job-08/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so')
    if sha(native) != sha(expected_native):
        raise ValueError('UNREGISTERED_RUNTIME_NATIVE')
    runtime = subprocess.run([sys.executable, '-m', 'pytest', '-q', str(test)], cwd=ROOT,
        capture_output=True, text=True, timeout=30)
    if runtime.returncode != 0:
        raise ValueError('RUNTIME_BRANCH_REGRESSION_FAILED:'+runtime.stdout+runtime.stderr)
    bindings.update(runtime_test=test, native=native, builder=Path(__file__),
        spin=ROOT/'python/make_it_stand/spin.py',
        adapter=ROOT/'python/make_it_stand/verification/cherchi_inprocess.py',
        material=ROOT/'python/make_it_stand/verification/faster_material.py',
        coordinate_contract=HERE/'case_coordinate_contract.py')
    if certified:
        arithmetic=HERE/'formal-rational-arithmetic-build-job-05/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
        bindings['certified_arithmetic']=arithmetic
        for state in (2,38,45,46):
            job=HERE/f'heart-arithmetic-05-certified-integral-state-{state}-job-01'
            sample=load(job/'output/report.json');result=load(job/'result.json')
            if (result['status']!='PROCESS_EXITED' or result['returncode']!=0 or
                sample['full_saved_mapping_equal'] is not True or len(sample['samples'])!=4 or
                sample['arithmetic_identity']['module_sha256']!=sha(arithmetic)):
                raise ValueError('CERTIFIED_MAPPING_NOT_VERIFIED')
            for name in ('config.json','result.json','output/report.json'):
                bindings[f'certified_state_{state}_{name}']=job/name
        first=HERE/'certified-integral-small-closure-job-02'
        repeat=HERE/'certified-integral-small-closure-repeat-01'
        a,b=[load(p/'output/report.json') for p in (first,repeat)]
        if not all(a[k]==b[k] for k in ('material_digest','scientific_event_digest','raw_equation_residual','accepted_updates')):
            raise ValueError('CERTIFIED_CLOSURE_REPEAT_CHANGED')
        for label,p in (('first',first),('repeat',repeat)):
            r=load(p/'result.json');v=load(p/'output/report.json')
            if r['status']!='PROCESS_EXITED' or r['returncode']!=0 or not v['material_closed_loop_passed']:
                raise ValueError('CERTIFIED_CLOSURE_NOT_COMPLETE')
            for name in ('config.json','result.json','output/report.json'):
                bindings[f'certified_closure_{label}_{name}']=p/name
        facets=HERE/'heart-exact-facet-deltas-03.json'
        if load(facets)['finite_support_classification_complete'] is not True:
            raise ValueError('STALLED_FINITE_CLASSIFICATION_INCOMPLETE')
        bindings['stalled_finite_classification']=facets
    output.mkdir(parents=True, exist_ok=False)
    (output/'runtime-test-result.json').write_text(json.dumps(dict(returncode=runtime.returncode,
        stdout=runtime.stdout, stderr=runtime.stderr, native_sha256=sha(native))))
    report = dict(assessment, schema='heart-guarded-solve-admission-v1', guarded_solve_admission=True,
        runtime_branch_regression_passed=True, input_digest=data['preparation']['input_digest'],
        request_sha256=sha(frozen), source_calls=0, boolean_calls=0, optimization_invocations=0,
        runtime_policy='FRESH_CANDIDATE_JACOBIAN_RESIDUAL_DECREASE_RELINEARIZE_ELSE_BACKTRACK',
        bindings={key: dict(path=str(path.resolve()), sha256=sha(path)) for key, path in bindings.items()})
    if certified:
        report.update(integral_profile='CERTIFIED_MPFR256_GMP_FALLBACK_V1',
                      production_precision_scope='CORRECTLY_ROUNDED_INTEGRAL_EXACT_GEOMETRY_GMP_DERIVATIVE')
    (output/'report.json').write_text(json.dumps(report, sort_keys=True, allow_nan=False))
    print(json.dumps({key: value for key, value in report.items() if key != 'bindings'}), flush=True)


if __name__ == '__main__':
    main(Path(sys.argv[1]),certified='--certified-integral' in sys.argv)
