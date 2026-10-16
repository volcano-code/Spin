"""New equivalent interior source mapping; original source failure remains."""
import json,sys
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
import _cherchi_core_boolean as bridge
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parents[1]/'python'))
from make_it_stand.verification.cherchi_inprocess import exact_fraction,construction_constraint_key
from make_it_stand.verification.cherchi_motion_differential import recover_construction_vertex_motion
from probe_inprocess_derivatives import weights
folder=HERE/'heart-radial-source-derivatives-job-01/output'
config=json.loads((HERE/'heart-certified-formal-seed-1-job-01/output/run/configuration.json').read_text());r=config['request']
design=(np.array(r['design_vertices'])-r['target_origin'])@np.array(r['target_axes']).T
margin=float(np.linalg.norm(np.ptp(design,axis=0))*.1);lo=design.min(axis=0)-margin;hi=design.max(axis=0)+margin
bv=np.array([[lo[0],lo[1],lo[2]],[hi[0],lo[1],lo[2]],[hi[0],hi[1],lo[2]],[lo[0],hi[1],lo[2]],[lo[0],lo[1],hi[2]],[hi[0],lo[1],hi[2]],[hi[0],hi[1],hi[2]],[lo[0],hi[1],hi[2]]])
bf=[[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],[1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]]
base=np.load(folder/'point-000.npz');samples=[]
for path in sorted(folder.glob('point-*.npz')):
    p=np.load(path);delta=p['coefficients']-base['coefficients']
    if np.count_nonzero(delta)==0 or (np.count_nonzero(delta)==1 and delta[5]!=0):samples.append((path,p))
out=HERE/'heart-box-clip-source-probe-01';out.mkdir(exist_ok=False);records=[]
try:
    for path,s in samples:
        sv=s['vertices'];packet=bridge.boolean(bv.tolist(),bf,sv.tolist(),s['faces'].tolist(),'intersection')
        lookup={}
        for obj,vs in enumerate([bv,sv]):
            for i,v in enumerate(vs):lookup.setdefault(tuple(Q(float(x)) for x in v),[]).append(dict(input_object=obj,input_vertex=i))
        provenance=[]
        for item in packet['construction_records']:
            operands=[];positions=[]
            for row in item['operand_coordinates_exact']:
                point=tuple(exact_fraction(x) for x in row);matches=lookup.get(point,[])
                if len(matches)!=1:raise ValueError('CLIPPED_SOURCE_PROVENANCE_NOT_UNIQUE')
                operands.append(dict(source_vertices=matches,status='UNIQUE'));positions.append(list(map(float,point)))
            provenance.append(dict(vertex_index=item['vertex_index'],point_type=item['point_type'],parameter_t=item['parameter_t'],operands=operands,operand_coordinates=positions))
        keys=[construction_constraint_key(x) for x in provenance]
        if len(set(map(str,keys)))!=len(keys):raise ValueError('CLIPPED_SOURCE_FEATURE_COLLISION')
        order=sorted(range(len(keys)),key=lambda i:str(keys[i]));points=np.array(packet['vertices'])[order]
        velocity=recover_construction_vertex_motion(provenance,bv,sv,s['velocity'].reshape(-1,3,10),np.array(packet['vertices']))[order]
        records.append(dict(offset=float(s['coefficients'][5]-base['coefficients'][5]),keys=[keys[i] for i in order],points=points,velocity=velocity))
        (out/(path.stem+'.json')).write_text(json.dumps(packet));print(path.name,len(points),flush=True)
    baseline=next(x for x in records if x['offset']==0);checks=[]
    assert all(x['keys']==baseline['keys'] for x in records),'CLIPPED_FEATURE_BRANCH_CHANGED'
    for h in [3e-6,1.5e-6,.75e-6]:
        selected=[min(records,key=lambda x:abs(x['offset']-s)) for s in [-h,-h/2,0,h/2,h]]
        w=weights([x['offset'] for x in selected]);fd=np.tensordot(w,np.array([x['points'] for x in selected]),axes=1);j=baseline['velocity'][:,:,5]
        error=np.linalg.norm(fd-j,axis=1)/np.maximum(1,np.linalg.norm(j,axis=1));checks.append(dict(scale=h,error=float(error.max())))
    report=dict(checks=checks,passed=all(x['error']<=1e-6 for x in checks),domain_box_margin=margin,scope='NEW_CLIPPED_MAPPING_DIRECTION5_ONLY',case_admission=False)
except Exception as ex:report=dict(error=repr(ex),case_admission=False)
(out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
