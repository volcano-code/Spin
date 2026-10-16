import json,sys,time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
import _make_it_stand as n
from make_it_stand.verification.rational_text import fraction_from_rational_text
p=json.loads((HERE/'inertia-boolean65-capture-job-07/output/capture/boolean-evaluation65.json').read_text())
r=json.loads((HERE/'cherchi-traversal-diagnostic-run-job-05/output/report.json').read_text())['result']
def decimal(x):
    a,b=x.split('/');return str(int(a,0))+'/'+str(int(b,0))
q=[[decimal(x) for x in row] for row in r['exact_vertices']]
out=HERE/'boolean65-candidate-verification-01.json';report={}
try:
    t=time.monotonic();audit=dict(n.audit_exact_rational_triangle_boundary(q,np.array(r['faces'],dtype=np.int32),'AnyNonzero'))
    report['boundary']={k:audit[k] for k in ('boundary_valid','intersecting_pair_count','degenerate_faces','component_volume_signs')}
    if not audit['boundary_valid']:raise ValueError('CANDIDATE_EXACT_BOUNDARY_INVALID')
    f=np.array(r['faces'],dtype=np.int32)
    # Cherchi may emit the opposite global component orientation. Apply only
    # the established closed-component convention, recording every reversal.
    flips=[]
    for i,component in enumerate(audit['face_components']):
        if audit['component_volume_signs'][component]<0:f[i]=f[i,::-1];flips.append(i)
    a=n.make_exact_material_region(np.array(p['design_vertices']),np.array(p['design_faces'],dtype=np.int32))
    b=n.make_exact_material_region(np.array(p['source_vertices']),np.array(p['source_faces'],dtype=np.int32))
    c=n.make_exact_material_region_from_rationals(q,f)
    operation='difference' if p['operation'] in ('subtraction','difference') else p['operation']
    report['relation']=dict(n.relate_exact_boolean_operation_to_material(a,b,c,operation));report['orientation_reversed_faces']=flips
    report['seconds']=time.monotonic()-t
except Exception as ex:report['error']=repr(ex)
out.write_text(json.dumps(report));print(json.dumps(report),flush=True)
