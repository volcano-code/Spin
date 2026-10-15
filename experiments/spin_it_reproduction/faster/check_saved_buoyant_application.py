#!/usr/bin/env python3
"""Recheck saved labelled buoyant geometry; no optimization or paper promotion."""
import argparse
import json
from pathlib import Path
import time

import numpy as np

from make_it_stand import spin
from make_it_stand.verification.buoyant_material import audit_saved_buoyant_case
from make_it_stand.verification.run_evidence import begin_run, finish_run, EXPLICIT_EFFECTIVE_SOURCE_POLICY


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--application', type=Path, required=True)
    parser.add_argument('--bundle-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[3]
    application = args.application.resolve(strict=True)
    bundles = args.bundle_root.resolve(strict=True)
    payload = json.loads(application.read_text())
    cases = payload['cases']
    if not cases or len({c['id'] for c in cases}) != len(cases):
        raise ValueError('missing or duplicate buoyant cases')
    assets = {application}
    for case in cases:
        for run in case['runs']:
            for stage in run['stages']:
                if stage.get('relative_bundle') is None:
                    continue
                relative = Path(stage['relative_bundle'])
                path = (bundles/relative).resolve(strict=True)
                if relative.is_absolute() or bundles not in path.parents:
                    raise ValueError('bundle path escapes evidence root')
                assets.update(p.resolve(strict=True) for p in path.rglob('*') if p.is_file())
    args.output.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    begin = begin_run(repository_root=root, runner_files=[Path(__file__)],
        checker_files=[*sorted((root/'python/make_it_stand/verification').glob('*.py')),
                       root/'python/make_it_stand/spin.py'],
        asset_files=sorted(assets), build_configuration_files=[root/'CMakeLists.txt'],
        parameters={'scope':'SAVED_BUOYANT_MATERIAL_METHOD_ONLY',
                    'operation':'READ_SAVED_BUOYANT_APPLICATION', 'bundle_root':str(bundles)},
        seed_protocol={'seed':0, 'protocol':'NO_OPTIMIZATION_SAVED_MATERIAL_RELOAD'},
        dependency_roots={'numpy':Path(np.__file__).resolve().parent},
        native_modules=[spin._native_module()], pure_python=False,
        source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY)
    (args.output/'begin.json').write_text(begin.to_json()+'\n')
    reports = [{'id':case['id'], **audit_saved_buoyant_case(case, bundle_root=bundles)} for case in cases]
    report = {'schema':'saved-buoyant-application-audit-v1', 'cases':reports,
              'scope':'MATERIAL_AND_INTERNAL_LEDGER_CONSISTENCY_NOT_SOURCE_EXECUTION_OR_PAPER_ATTESTATION'}
    scientific = args.output/'scientific.json'
    scientific.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n')
    events = [{'iteration':i, 'status':'ok' if row['qualified'] else 'failed',
               'data':{'case':case['id'], **row}}
              for case in reports for i,row in enumerate(case['runs'])]
    evidence = finish_run(begin, events=events, output_files=[scientific],
                         measurements={'elapsed_seconds':time.perf_counter()-start})
    (args.output/'evidence.json').write_text(evidence.to_json()+'\n')
    print(json.dumps({'cases':len(cases), 'qualified_runs':sum(c['qualified_runs'] for c in reports)}))


if __name__ == '__main__':
    main()
