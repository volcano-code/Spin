"""Recompute run freshness, repeat identities and the historical trace prefix.

Record comparison only; this program cannot grant paper or artifact acceptance.
"""
import argparse
import json
from pathlib import Path
from teapot import STRATEGIES,check_files,file_hashes,write


def audit(base, historical):
    rows=[];sources=file_hashes([Path(__file__),Path(__file__).with_name('teapot.py')])
    observed=['before','after','input_beta','beta','output_beta','initial_phi',
              'final_phi','evaluation','trace','final_kkt']
    for strategy in STRATEGIES:
        records=[]
        old_path=historical/strategy/'native-run.json'
        old=json.loads(old_path.read_text())['native_run']['solve']
        sources.update(file_hashes([old_path]))
        for run in (2,3):
            folder=base/f'{strategy}-run{run}'
            files=json.loads((folder/'outputs.json').read_text())
            check_files(files);sources.update(files)
            sources.update(file_hashes([folder/'outputs.json']))
            current=json.loads((folder/'source.json').read_text())['files']
            check_files(current);sources.update(current)
            cfg=json.loads((folder/'config.json').read_text())
            ref=json.loads(Path(cfg['reference_config']).read_text())
            expected={**ref['solver'],'termination_policy':'PAPER_POST_TRANSITION_V1'}
            if (cfg['solver_settings']!=expected or cfg['source_settings']!=ref['geometry'] or
                    cfg['boundary_settings']!=ref['boundary'] or
                    cfg['minimum_outer_wall']!=ref['minimum_wall']):
                raise ValueError('UNDECLARED_NUMERIC_TREATMENT')
            record=json.loads((folder/'native-run.json').read_text())
            if record['settings']['solver_settings']!=expected:
                raise ValueError('ACTUAL_NATIVE_SETTINGS_MISMATCH')
            summary=json.loads((folder/'summary.json').read_text());solve=record['solve']
            prefix_equal=len(solve['passes'])<=len(old['passes']) and all(
                p[k]==old['passes'][i][k] for i,p in enumerate(solve['passes']) for k in observed)
            records.append(dict(run=run,summary=summary,solve=solve,
                                historical_scientific_prefix_equal=prefix_equal))
        a,b=records
        rows.append(dict(strategy=strategy,status=a['solve']['status'],
            passes=len(a['solve']['passes']),
            normalized_binary_residual=a['solve']['binary_check']['normalized_constraint_residual'],
            independent_moment_error=a['solve']['binary_check']['normalized_moment_error'],
            terminal_relaxation_kkt_applies=a['solve']['terminal_relaxation_kkt_applies'],
            repeat_full_solve_equal=a['solve']==b['solve'],
            repeat_trace_hash_equal=a['solve']['trace_hash']==b['solve']['trace_hash'],
            repeat_output_hash_equal=a['solve']['output_hash']==b['solve']['output_hash'],
            historical_scientific_prefix_equal=all(r['historical_scientific_prefix_equal'] for r in records),
            elapsed_seconds=[r['summary']['elapsed_seconds'] for r in records],
            peak_rss_kib=[r['summary']['max_rss_kib'] for r in records],
            output_hash=a['solve']['output_hash'],trace_hash=a['solve']['trace_hash']))
    first=base/'inertia-run1/source.json'
    stale=[]
    if first.is_file():
        for name,h in json.loads(first.read_text())['files'].items():
            if not Path(name).is_file() or file_hashes([name])[str(Path(name).resolve())]!=h:
                stale.append(name)
    check_files(sources)
    return dict(scope='EXECUTION_RECORD_COMPARISON_NOT_PAPER_ACCEPTANCE',
        project_status='PARTIAL_BLOCKED',original_paper_case=False,
        current_runs_source_freshness='CURRENT_AT_AUDIT',cases=rows,
        earlier_inertia_run1=dict(status='STALE_RESULT' if stale else 'CURRENT_AT_AUDIT',
                                 mismatched_files=stale),files=sources)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--runs',type=Path,required=True)
    p.add_argument('--historical',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    report=audit(a.runs,a.historical);write(a.output,report)
    for row in report['cases']:
        print(row['strategy'],row['status'],'passes',row['passes'],
              'binary residual',row['normalized_binary_residual'],
              'repeat',row['repeat_full_solve_equal'],
              'old prefix',row['historical_scientific_prefix_equal'])
