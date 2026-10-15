"""Nine registered source points: existing algebraic finite-map contrast only."""
import json,sys,os,time,hashlib
from pathlib import Path
import numpy as np
import _make_it_stand as native
from probe_inprocess_derivatives import weights
from run_small_inprocess_solver import safe
HERE=Path(__file__).resolve().parent
saved=HERE/'noodles-interval-source-derivative-job-01/output'
config=json.loads((HERE/'noodles-bound-seed0-production-probe-job-01/output/run/mapping-config.json').read_text())
r=config['request'];e=config['evaluator'];design=(np.asarray(r['design_vertices'])-r['target_origin'])@np.asarray(r['target_axes'])
s=native.SpinTriangulationSettings()
for k in ('icosphere_subdivisions','maximum_icosphere_subdivisions','maximum_radius_doublings','retain_fixed_normal_provenance','truncation_radius','paper_projective_hyperboloid','algebraic_finite_map'):setattr(s,k,e[k])
s.algebraic_finite_map=True
layouts=[x for x in json.loads((saved/'manifest.json').read_text())['layouts'] if x['direction'] in (1,9) and x['scale']==7.5e-7]
ids=sorted({i for x in layouts for i in x['indices']})
if len(ids)!=9:raise ValueError('REGISTERED_POINT_COUNT_MISMATCH')
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False);records={};start=time.monotonic()
with (out/'events.jsonl').open('x') as journal:
 for i in ids:
    q=json.loads((saved/f'point-{i:03d}.json').read_text())['coefficients']
    source=native.evaluate_spin_it_faster_source_motion(np.asarray(q),design,s)
    with (out/f'point-{i:03d}.json').open('x') as f:json.dump(safe(dict(coefficients=q,source=source)),f)
    if not source['defined'] or not source['side_selection_defined']:raise ValueError(source['diagnostic'])
    records[i]=source;journal.write(json.dumps(dict(index=i,seconds=time.monotonic()-start))+'\n');journal.flush();os.fsync(journal.fileno())
b=records[0];j=np.asarray(b['vertex_velocity']).reshape(-1,3,10);checks=[]
for x in layouts:
 n=np.tensordot(weights(x['positions']),np.asarray([records[i]['vertices'] for i in x['indices']]),axes=1)
 error=np.linalg.norm(n-j[:,:,x['direction']],axis=1)/np.maximum(1,np.linalg.norm(j[:,:,x['direction']],axis=1))
 checks.append(dict(direction=x['direction'],scale=x['scale'],maximum_error=float(max(error)),worst_vertex=int(np.argmax(error)),passed=bool(np.all(error<=1e-6))))
report=dict(scope='NINE_SOURCE_POINT_EXISTING_ALGEBRAIC_FINITE_MAP_DIAGNOSIS',checks=checks,point_count=len(ids),native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),seconds=time.monotonic()-start,production_admission=False,optimization_admission=False)
with (out/'report.json').open('x') as f:json.dump(report,f,indent=2)
print(json.dumps(report),flush=True)
