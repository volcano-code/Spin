"""Read-only derivation from native run records; never grants paper acceptance."""
import argparse
import hashlib
import json
from pathlib import Path
from resolution_rounding_audit import rounding_diagnostics
from teapot import STRATEGIES,write


def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def verify_required_outputs(folder, hashes):
    required={str((folder/name).resolve()) for name in
              ('native-run.json','config.json','summary.json','source.json')}
    if not required.issubset(hashes):
        raise ValueError('MISSING_EVIDENCE_COVERAGE: '+str(sorted(required-set(hashes))))


def verify_native_settings(actual, declared):
    for key in ('source_settings','solver_settings','boundary_settings','density','minimum_outer_wall'):
        if key not in actual or key not in declared or actual[key]!=declared[key]:
            raise ValueError('NATIVE_SETTINGS_MISMATCH: '+key)


def pass_metrics(passes):
    last=passes[-1] if passes else {}
    kkt=last.get('final_kkt') or {}
    return dict(largest_solved_partition=max((len(v['before']) for v in passes),default=None),
        maximum_solved_leaf_level=max((v['level'] for step in passes for v in step['before']),default=None),
        final_partition_level_counts=last.get('output_level_counts'),
        final_fractional_level_counts=last.get('output_fractional_level_counts'),
        last_kkt=kkt.get('full'),last_equality=kkt.get('equality_primal'),
        last_solve_status=last.get('status'),
        last_solve_iterations=len(last['trace']) if 'trace' in last else None)


def repeat_comparison(first, second):
    available=all(v.get('native_record_sha256') for v in (first,second))
    return dict(repeat_comparison_status='AVAILABLE' if available else 'UNAVAILABLE',
        repeat_native_record_bytes_equal=(first['native_record_sha256']==second['native_record_sha256'])
        if available else None)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    rows=[];files={str(Path(__file__).resolve()):sha(__file__)}
    helper=Path(__file__).with_name('resolution_rounding_audit.py').resolve();files[str(helper)]=sha(helper)
    # Explicit expected identities keep a resource-skipped repeat in the ledger.
    matrix=[(s,g,o,r,'DEPTH_ONLY') for g,o in [(6,6),(7,6),(7,7)]
            for s in STRATEGIES for r in (1,2)]
    matrix += [('inertia',6,6,r,'ITERATION_960') for r in (1,2)]
    for strategy,g,o,repeat,group in matrix:
        name=f'{strategy}-g{g}-o{o}'+('-i960' if group=='ITERATION_960' else '')+f'-r{repeat}'
        folder=a.root/name
        if not (folder/'summary.json').is_file():
            rows.append(dict(name=name,group=group,status='NOT_COMPLETED',repeat=repeat,
                             artifacts_present=folder.exists()))
            continue
        manifests=[folder/'outputs.json',folder/'source.json']
        for manifest in manifests:
            hashes=json.loads(manifest.read_text())
            if manifest.name=='source.json':hashes=hashes['files']
            else:verify_required_outputs(folder,hashes)
            for name_,expected in hashes.items():
                if not Path(name_).is_file() or sha(name_)!=expected:
                    raise ValueError('STALE_RESULT: '+name_)
                files[name_]=expected
            files[str(manifest.resolve())]=sha(manifest)
        summary=json.loads((folder/'summary.json').read_text())
        cfg=json.loads((folder/'config.json').read_text())
        native=json.loads((folder/'native-run.json').read_text());s=native['solve']
        verify_native_settings(native['settings'],cfg)
        if s['status']!=summary['status'] or s['output_hash']!=summary['output_hash']:
            raise ValueError('NATIVE_SUMMARY_MISMATCH')
        if (cfg['source_settings']['maximum_depth']!=g or
            cfg['solver_settings']['maximum_optimization_level']!=o or
            cfg['solver_settings']['objective']!=STRATEGIES[strategy] or
            cfg['solver_settings']['slqp']['maximum_iterations']!=(960 if group=='ITERATION_960' else 240)):
            raise ValueError('EXPERIMENT_IDENTITY_MISMATCH')
        row=dict(name=name,group=group,strategy=strategy,geometry_depth=g,optimization_max=o,
            repeat=repeat,status=s['status'],diagnostic=s['diagnostic'],passes=len(s['passes']),
            actual_geometry_depth=native['source']['stats']['maximum_visited_level'],
            **pass_metrics(s['passes']),
            source_hash=native['source']['source_hash'],
            source_partition_hash=native['source']['partition_hash'],
            unresolved_boundary_volume=native['source']['unresolved_boundary_volume'],
            conservative_carveable_volume=native['source']['carveable_moments'][0],
            binary_check=s['binary_check'],rounding=rounding_diagnostics(native),
            elapsed_seconds=summary['elapsed_seconds'],max_rss_kib=summary['max_rss_kib'],
            native_record_sha256=sha(folder/'native-run.json'),trace_hash=s['trace_hash'],output_hash=s['output_hash'])
        if repeat==2:
            first=next(v for v in rows if v['name']==name[:-1]+'1')
            row.update(repeat_comparison(first,row))
        rows.append(row)
        print(name,s['status'],'max_solved',row['maximum_solved_leaf_level'],
              'binary',None if s['binary_check'] is None else s['binary_check']['normalized_constraint_residual'],flush=True)
        del native,s
    for file in [a.root/'PLAN.md',a.root/'ITERATION_PILOT.md']:
        files[str(file.resolve())]=sha(file)
    write(a.output,dict(scope='EXPERIMENT_AND_TRACE_DIAGNOSTICS_NOT_PAPER_ACCEPTANCE',
        project_status='PARTIAL_BLOCKED',original_paper_case=False,expected_depth_runs=18,
        expected_separate_budget_runs=2,cases=rows,files=files))
