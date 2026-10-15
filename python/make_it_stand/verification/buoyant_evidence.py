"""Re-read bounded buoyant execution evidence before admitting method support.

Local provenance consistency, not authenticated logs, compiler attestation,
author-case reproduction, or manufacture approval. No solver is called here.
The first supported execution protocol is the saved single-case A→B worker.
"""
from pathlib import Path

from .buoyant_application import buoyant_case_plan
from .buoyant_material import audit_saved_buoyant_case
from .run_evidence import PaperRunEvidenceV2, check_freshness, _record
from .source_checkpoint_store import _fingerprints
from .spin_common import parse, require, digest

SCOPE = 'SAVED_BUOYANT_EXECUTION_AND_MATERIAL_METHOD_V1'
OPERATION = 'READ_SAVED_BUOYANT_EXECUTION'
JOB_FILES = ('config.json', 'kernel.json', 'result.json', 'complete.json')


def _read(path):
    return parse(Path(path).read_text(encoding='utf-8'))


def assert_independent_buoyant_runs(records):
    require(len(records) == 2, 'two independent native runs are required')
    first, second = records
    require(first['execution_sha256'] != second['execution_sha256'] and
            (first['boot_id'], first['unit']) != (second['boot_id'], second['unit']),
            'reused original evidence is the same native run, not an independent repeat')


def _job_receipt(root, request, observation, execution, application):
    config, kernel, result, complete = [_read(root/name) for name in JOB_FILES]
    require(config['schema'] == 'source-job-v1' and
            result['schema'] == 'source-job-result-v1', 'unsupported source-job receipt')
    require(complete == {'result_sha256': _record(root/'result.json')['sha256']},
            'source-job completion digest mismatch')
    require(result['status'] == 'PROCESS_EXITED' and result['returncode'] == 0 and
            result['unit'] == config['unit'] and kernel['boot_id'] == config['boot_id'],
            'original execution did not complete in the recorded source job')
    limits = config['limits']
    require(0 < limits['maximum_wall_ns'] <= 3600_000_000_000 and
            0 <= result['elapsed_ns'] <= limits['maximum_wall_ns'] and
            0 < kernel['memory_max'] == limits['maximum_memory_bytes'] <= 24*1024**3 and
            kernel['swap_max'] == 0 and
            0 <= result['accounted_disk_bytes'] <= limits['maximum_disk_bytes'] and
            limits['minimum_start_free_bytes'] >= 100*1024**3 and
            limits['minimum_free_bytes'] >= 20*1024**3,
            'source-job resource envelope is invalid')
    service = result['service']
    require(service['Result'] == 'success' and service['ExecMainStatus'] == '0' and
            int(service['MemoryMax']) == kernel['memory_max'] and
            int(service['MemorySwapMax']) == 0 and int(service['MemorySwapPeak']) == 0 and
            0 <= int(service['MemoryPeak']) <= kernel['memory_max'] and
            service['ControlGroup'] == kernel['cgroup'] and
            kernel['cgroup'].endswith('/'+config['unit']), 'source-job resource receipt mismatch')
    require(config['identity'] == _fingerprints(config['identity_files']),
            'STALE_RESULT: original source-job effective inputs changed')
    pinned = {row['path']: row for row in config['identity'].values()}
    # Bind effective request files and the actually observed loaded native binary.
    for rows in observation['files'].values():
        for row in rows:
            require(pinned.get(row['path']) == row, 'source job did not bind execution input')
    for row in observation['runtime']['native_modules'].values():
        require(pinned.get(row['path']) == row, 'source job did not bind loaded native module')
    argv = config['argv']
    require(len(argv) >= 4 and argv[-3:] == ['worker', request['asset_files'][0], str(application.parent)] and
            argv[-4] in request['runner_files'] and argv[-4] in pinned and
            execution == application.parent/'evidence.json' and
            application.parent == root/'outputs', 'source-job worker/output binding mismatch')
    return {'execution_sha256': _record(execution)['sha256'],
            'boot_id': config['boot_id'], 'unit': config['unit'],
            'elapsed_ns': result['elapsed_ns'], 'memory_peak_bytes': int(service['MemoryPeak']),
            'swap_peak_bytes': int(service['MemorySwapPeak'])}


def _validate_saved_ledger(inputs, cases, recorded_events):
    """Bind request/run labels and every recorded stage, not caller pass flags."""
    require(len(inputs) == len(cases) == 1, 'supported native worker requires exactly one case')
    plan = buoyant_case_plan(inputs[0], 1, 1, [0])
    case = cases[0]
    require(all(case.get(key) == item for key, item in plan.items() if key != 'runs'),
            'saved application differs from frozen request')
    require(len(case['runs']) == 1 and all(case['runs'][0].get(k) == v for k,v in plan['runs'][0].items()),
            'saved application lost an original run or initialization label')
    require(case.get('asset_class') == inputs[0].get('asset_class'), 'asset classification changed')
    events = [{'iteration': i, 'status': 'ok' if s['status'] == 'CONVERGED' else 'failed',
               'data': {'stage': s['stage'], 'trace': s['trace'],
                        'material_digest': s['artifact_scientific_digest']}}
              for run in case['runs'] for i,s in enumerate(run['stages'])]
    require(events == recorded_events, 'original execution events differ from saved ledger')
    return events


def audit_buoyant_execution(*, application, bundle_root, execution_evidence, source_job_root):
    execution = Path(execution_evidence)
    require(execution.is_file(), 'original execution evidence is required')
    execution = execution.resolve(strict=True)
    application = Path(application).resolve(strict=True)
    bundle = Path(bundle_root).resolve(strict=True)
    job = Path(source_job_root).resolve(strict=True)
    evidence = PaperRunEvidenceV2.from_json(execution.read_text())
    freshness = check_freshness(evidence)
    require(freshness.status == 'CURRENT', 'STALE_RESULT: original execution '+str(freshness))
    value = evidence.as_dict()
    request = value['begin']['request']
    require(request['parameters'] == {'operation': 'ACTUAL_SEALED_HULL_BUOYANT_A_THEN_B',
                                      'starts': 1, 'max_restarts': 1} and
            request['seed_protocol'] == {'seed': 0, 'protocol': 'LOCKED_EXPLICIT_PLANES_NOT_RANDOM_TABLE_STATISTICS'} and
            request['native_modules'] == ['_make_it_stand'] and request['pure_python'] is False and
            len(request['asset_files']) == 1, 'unsupported original buoyant execution protocol')
    provenance = _job_receipt(job, request, value['begin']['observation'], execution, application)
    require(bundle == application.parent/'materials', 'original material root mismatch')
    inputs = _read(request['asset_files'][0])['cases']
    cases = _read(application)['cases']
    events = _validate_saved_ledger(inputs, cases, value['scientific']['events'])
    case = cases[0]
    expected_outputs = {application, *[p.resolve(strict=True) for p in bundle.rglob('*') if p.is_file()]}
    outputs = {Path(row['path']).resolve(strict=True) for row in value['scientific']['outputs']}
    require(outputs == expected_outputs, 'original execution output inventory differs from saved materials')
    report = audit_saved_buoyant_case(case, bundle_root=bundle)
    return {'schema': 'saved-buoyant-execution-audit-v1', 'scope': SCOPE,
            'provenance': provenance, 'cases': [{'id': case['id'], **report}],
            'input_sha256': digest(inputs), 'application_sha256': digest({'cases': cases}),
            'event_sha256': digest(events), 'build_source_attestation': 'NOT_VERIFIED',
            'verification_effect': 'SUPPORTING_ONLY_NO_CATALOGUE_REQUIREMENT_VERIFIED'}


def scientific_signature(report):
    """Process identity and resource measurements are not scientific outputs."""
    return {key: value for key,value in report.items() if key != 'provenance'}
