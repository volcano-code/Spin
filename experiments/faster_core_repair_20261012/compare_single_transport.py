"""Same registered states, old/new independent processes, full numerical comparison."""
import json,os,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]


def capture(path):
    from dataclasses import asdict
    import numpy as np
    from make_it_stand import spin
    frozen=json.loads((HERE/'same-input-job-03/output/frozen-inputs.json').read_text())
    payload=dict(frozen['request']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
    request=spin.SpinItFasterRequest(**payload)
    cfg=dict(frozen['settings']);cfg['exact_binary64_material_boundary']=True
    settings=spin.SpinEvaluatorConfig(**cfg);rows=[]
    for name,coefficients in frozen['states'].items():
        q=np.asarray(coefficients,dtype=float);q/=np.linalg.norm(q)
        for repetition in range(4):
            events=[];start=time.monotonic()
            result=spin.evaluate_spin_it_faster_fullrank_quadric(request,q,evaluator_settings=settings,
                diagnostic_stage_observer=events.append)
            rows.append(dict(state=name,repetition=repetition,seconds=time.monotonic()-start,
                events=events,result=asdict(result)))
    with path.open('x') as stream:json.dump(rows,stream,allow_nan=False)


def compare(output):
    import numpy as np
    output.mkdir(parents=True,exist_ok=False)
    for role,module in (('before','component-selector-build-job-01'),('after','single-transport-build-job-01')):
        env=dict(os.environ,PYTHONPATH=str(ROOT/'python')+':'+str(HERE/module/'output'),
            LD_PRELOAD='/usr/lib/x86_64-linux-gnu/libstdc++.so.6')
        subprocess.run([sys.executable,'-B',__file__,'--capture',str(output/(role+'.json'))],env=env,check=True)
    before=json.loads((output/'before.json').read_text());after=json.loads((output/'after.json').read_text());rows=[]
    for a,b in zip(before,after):
        same=all(a['result'][key]==b['result'][key] for key in (
            'variable_raw_moments','total_raw_moments','geometry_hash','material_vertices','material_faces',
            'derivative_defined','derivative_diagnostic','coefficient_jacobian'))
        count=lambda r:sum(x in ('material-transport-begin','canonical-material-transport-begin') for x in r['events'])
        rows.append(dict(state=a['state'],repetition=a['repetition'],numerical_material_equal=same,
            before_transport_calls=count(a),after_transport_calls=count(b),
            before_seconds=a['seconds'],after_seconds=b['seconds']))
    passed=len(rows)==12 and all(r['numerical_material_equal'] and r['before_transport_calls']==2 and r['after_transport_calls']==1 for r in rows)
    (output/'report.json').write_text(json.dumps(dict(passed=passed,rows=rows,core_admission=False,
        scope='Exact-material identical-transport overlay only; scientific identities change with build')))
    if not passed:raise SystemExit(1)


if __name__=='__main__':
    if sys.argv[1]=='--capture':capture(Path(sys.argv[2]))
    else:compare(Path(sys.argv[1]))
