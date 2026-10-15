#!/usr/bin/env python3
"""Check saved buoyant execution, job receipt and material without optimizing."""
import argparse
from pathlib import Path
import time

import numpy as np
from make_it_stand import spin
from make_it_stand.verification.buoyant_evidence import audit_buoyant_execution, SCOPE, OPERATION, JOB_FILES
from make_it_stand.verification.run_evidence import begin_run, finish_run, EXPLICIT_EFFECTIVE_SOURCE_POLICY
from make_it_stand.verification.spin_common import canonical


def main():
    parser = argparse.ArgumentParser()
    for flag in ('application', 'bundle-root', 'execution-evidence', 'source-job-root', 'output'):
        parser.add_argument('--'+flag, type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[3]
    paths = {key: getattr(args,key).resolve(strict=True) for key in
             ('application','bundle_root','execution_evidence','source_job_root')}
    assets = {paths['application'], paths['execution_evidence'],
              *[paths['source_job_root']/name for name in JOB_FILES],
              *[p.resolve(strict=True) for p in paths['bundle_root'].rglob('*') if p.is_file()]}
    args.output.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    parameters = {'scope': SCOPE, 'operation': OPERATION, **{k:str(v) for k,v in paths.items()}}
    begin = begin_run(repository_root=root, runner_files=[Path(__file__)],
        checker_files=[*sorted((root/'python/make_it_stand/verification').glob('*.py')),
                       root/'python/make_it_stand/spin.py'], asset_files=sorted(assets),
        build_configuration_files=[root/'CMakeLists.txt'], parameters=parameters,
        seed_protocol={'seed':0, 'protocol':'NO_OPTIMIZATION_ORIGINAL_EXECUTION_RECHECK'},
        dependency_roots={'numpy':Path(np.__file__).resolve().parent},
        native_modules=[spin._native_module()], pure_python=False,
        source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY)
    (args.output/'begin.json').write_text(begin.to_json()+'\n')
    report = audit_buoyant_execution(**paths)
    scientific = args.output/'scientific.json'
    scientific.write_text(canonical(report)+'\n')
    events = [{'iteration':i,'status':'ok' if case['qualified_runs'] else 'failed',
               'data':case} for i,case in enumerate(report['cases'])]
    evidence = finish_run(begin, events=events, output_files=[scientific],
                         measurements={'elapsed_seconds':time.perf_counter()-start})
    (args.output/'evidence.json').write_text(evidence.to_json()+'\n')
    print(canonical({'scope':SCOPE,'qualified_runs':sum(c['qualified_runs'] for c in report['cases'])}))


if __name__ == '__main__':
    main()
