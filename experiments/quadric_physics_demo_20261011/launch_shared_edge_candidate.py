"""Bounded serial exact-predicate experiment, preserving the frozen Heart build."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from launch_native_geometry_profile import paused_dispatcher

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def main(output=None):
    sys.path.insert(0,str(ROOT/'python'))
    from make_it_stand.verification.source_job import run_source_job
    output=Path(output) if output is not None else HERE/'three-dimensional-faster-performance-v1/edge-predicate-candidate-v1'
    output.mkdir(exist_ok=False)
    with (output/'events.jsonl').open('x') as stream:
        def event(status,**values):
            row=dict(status=status,**values)
            stream.write(json.dumps(row)+'\n');stream.flush();os.fsync(stream.fileno())
            print(json.dumps(row),flush=True)
        with paused_dispatcher(1766941,'rounded_production_sequence.py'):
            queue=ROOT/'experiments/faster_heart_contact_20261011/rounded-clearance-v3-production-sequence-01/events.jsonl'
            rows=[json.loads(row) for row in queue.read_text().splitlines()]
            active=Path(next(row['job'] for row in reversed(rows) if row['status']=='POINT_STARTED'))
            event('WAIT_CURRENT_HEART_POINT',job=str(active))
            deadline=time.monotonic()+3660
            while not (active/'complete.json').exists():
                if time.monotonic()>deadline:raise RuntimeError('ACTIVE_POINT_NOT_TERMINAL')
                time.sleep(2)
            worker=HERE/'shared_edge_candidate_worker.py'
            inputs={'launcher':Path(__file__),'scheduler':HERE/'launch_native_geometry_profile.py',
                'worker':worker,'predicate':ROOT/'include/verification/ExactSharedEdgeIntersection.h',
                'audit':ROOT/'src/verification/ExactBoundaryAudit.cpp',
                'audit_header':ROOT/'include/verification/ExactBoundaryAudit.h',
                'rational_header':ROOT/'include/verification/ExactRationalCoordinates.h',
                'sha_header':ROOT/'include/verification/Sha256.h',
                'oracle':HERE/'shared_edge_predicate_test.cpp','regression':ROOT/'tests/exact_boundary_audit_regression.cpp',
                'benchmark':HERE/'exact_boundary_benchmark.cpp','compiler':Path('/usr/bin/c++')}
            def job(name,operation,extra,identity):
                event('JOB_STARTED',name=name,operation=operation)
                result=run_source_job(output/name,
                    ['/home/lin/anaconda3/bin/python3','-B',str(worker),operation,'--output','{job}/output',*extra],
                    identity)
                event('JOB_TERMINAL',name=name,status_code=result.status,returncode=result.returncode)
                if result.status!='PROCESS_EXITED' or result.returncode!=0:
                    raise RuntimeError('CANDIDATE_JOB_FAILED:'+name)
            job('build-job-01','build',[],inputs)
            built=output/'build-job-01/output'
            boundary=HERE/'noodles-dense-v5-domain-stage-04/output/design.rational.json'
            for name in ('candidate','baseline'):
                job(name+'-audit-job-01','benchmark',['--source',str(boundary),'--binary',str(built/(name+'-benchmark'))],
                    dict(inputs,source=boundary,binary=built/(name+'-benchmark'),build_receipt=built/'build.json'))
            candidate=json.loads((output/'candidate-audit-job-01/output/audit.json').read_text())
            baseline=json.loads((output/'baseline-audit-job-01/output/audit.json').read_text())
            match=candidate['facts']==baseline['facts']
            report=dict(scope='EXACT_BOUNDARY_CANDIDATE_NOT_FULL_MATERIAL_OR_OPTIMIZATION_ADMISSION',
                source_sha256=hashlib.sha256(boundary.read_bytes()).hexdigest(),same_facts=match,
                baseline_seconds=baseline['audit_seconds'],candidate_seconds=candidate['audit_seconds'],
                speedup=baseline['audit_seconds']/candidate['audit_seconds'],
                geometry_admission=False,optimization_admission=False)
            (output/'comparison.json').write_text(json.dumps(report,indent=2))
            event('BOUNDARY_COMPARISON_FINISHED',**report)
            if not match:raise RuntimeError('EXACT_BOUNDARY_FACTS_CHANGED')
        event('HEART_DISPATCHER_RESUMED')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--directory',type=Path)
    main(parser.parse_args().directory)
