"""Summarize genuine failed attempts; never upgrade relaxed moments to material."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from teapot import STRATEGIES,write,file_hashes,check_files
from admission import check_case_context


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--campaign',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(exist_ok=False)
    rows=[];sources=file_hashes([Path(__file__),Path(__file__).with_name('admission.py')]);common=None
    fig,axes=plt.subplots(2,3,figsize=(13,7),layout='constrained')
    for i,name in enumerate(STRATEGIES):
        folder=a.campaign/name
        record=json.loads((folder/'native-run.json').read_text());summary=json.loads((folder/'summary.json').read_text())
        config=json.loads((folder/'config.json').read_text());upstream=json.loads((folder/'source.json').read_text())['files']
        check_files(upstream);sources.update(upstream)
        sources.update(file_hashes([folder/n for n in ['native-run.json','summary.json','config.json','source.json']]))
        if common is None:common=record['problem_digest']
        check_case_context(record,config,summary,name,common)
        solve=record['native_run']['solve'];passes=[];seen={};repeats=[]
        for step in solve['passes']:
            beta=np.array(step['output_beta']);cells=step['after'];eps=config['solver']['fractional_epsilon']
            level=config['solver']['maximum_optimization_level']
            fraction=[c['level'] for b,c in zip(beta,cells) if eps<=b<=1-eps]
            after_at_max=all(x==level for x in fraction)
            key=step['after_hash']
            if key in seen:repeats.append(dict(previous=seen[key],current=step['index'],partition_hash=key))
            seen[key]=step['index']
            moment=step['evaluation']['final_moments']
            passes.append(dict(index=step['index'],variables=len(step['beta']),status=step['status'],
                kkt=step['final_kkt']['full'],equality=step['final_kkt']['equality_primal'],
                relaxed_mass_g=moment[0]*18.75,relaxed_center_height_mm=moment[3]/moment[0]*25,
                levels=step['level_counts'],after_fractional_levels=fraction,
                all_post_transition_fractional_at_maximum=after_at_max,
                before_hash=step['before_hash'],after_hash=key))
        row=dict(strategy=name,status=solve['status'],sealed=summary['sealed'],diagnostic=solve['diagnostic'],
            elapsed_seconds=summary['elapsed_seconds'],max_rss_kib=summary['max_rss_kib'],passes=passes,
            repeated_partition_hashes=repeats,binary_check=solve['binary_check'],
            actual_native_objective=record['solver_settings']['objective'])
        rows.append(row)
        titles=['Inertia only','Inertia + lowering','Inertia + mass/height']
        x=[s['index'] for s in passes]
        axes[0,i].plot(x,[s['variables'] for s in passes],'o-',color='#527fa3')
        axes[0,i].set(title=titles[i],xlabel='Adaptive pass (zero based)',ylabel='Optimization variables')
        axes[0,i].grid(alpha=.2)
        axes[1,i].semilogy(x,[s['kkt'] for s in passes],'o-',color='#b87039')
        axes[1,i].axhline(1e-7,color='#777777',linestyle='--',linewidth=1)
        axes[1,i].set(xlabel='Adaptive pass',ylabel='Last continuous KKT residual')
        axes[1,i].grid(alpha=.2)
        axes[1,i].text(.03,.97,solve['status'].replace('_','\n'),transform=axes[1,i].transAxes,
            va='top',fontsize=9,color='#a32828',bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
    fig.suptitle('Reconstructed Utah Teapot | Three actual native attempts\nContinuous convergence does not imply valid binary material',fontsize=15)
    fig.savefig(a.output/'optimization-diagnostics.png',dpi=150);plt.close(fig)
    check_files(sources)
    write(a.output/'source.json',dict(files=sources))
    write(a.output/'report.json',dict(status='NO_OPTIMIZED_MATERIAL_ADMITTED' if not any(r['sealed'] for r in rows) else 'PARTIAL',
        problem_digest=common,source_unchanged=True,original_paper_case=False,cases=rows,
        geometry_depth=6,maximum_optimization_level=5,paper_level9=False,
        termination_note='Post-transition fractional-level condition is reported for diagnosis only. No paper-pass, changed solver status, or fabricated binary output is granted.'))
    print('DIAGNOSTICS WRITTEN',flush=True)
