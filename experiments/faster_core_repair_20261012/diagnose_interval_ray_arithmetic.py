"""Frozen vertex 93 only: independent arithmetic; no new map/Boolean calls."""
import hashlib,json,subprocess,sys
from pathlib import Path
import mpmath as mp
import numpy as np
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
source=HERE/'inspect_interval_normals.cpp';exe=output/'inspect-normals'
template=json.loads((ROOT/'build/heart-fixed-facet-20261011-v1/commands.json').read_text())[0]
eigen=[x for x in template if 'eigen' in x.lower() and x.startswith('-I')]
command=['/usr/bin/c++','-O2','-std=c++17','-I'+str(ROOT/'src'),*eigen,str(source),'-o',str(exe)]
subprocess.run(command,check=True,capture_output=True)
normals=json.loads(subprocess.check_output([str(exe)],text=True))
saved=HERE/'noodles-interval-source-derivative-job-01/output'
manifest=json.loads((saved/'manifest.json').read_text())
records=[json.loads((saved/f'point-{i:03d}.json').read_text()) for i in range(81)]
def exact_float(x):
    a,b=float(x).as_integer_ratio();return mp.mpf(a)/b
def replay(q, digits):
    with mp.workdps(digits):
        q=list(map(exact_float,q));a=mp.matrix([[q[4],q[9],q[8]],[q[9],q[5],q[7]],[q[8],q[7],q[6]]]);B=a**-1
        b=mp.matrix(q[1:4]);n=mp.matrix(list(map(exact_float,normals[0])));d=mp.matrix(list(map(exact_float,normals[1])))-n
        dot=lambda x,y:(x.T*y)[0]
        aa=dot(d,B*d);bb=2*dot(n,B*d);cc=dot(n,B*n)
        roots=sorted([(-bb-mp.sqrt(bb*bb-4*aa*cc))/(2*aa),(-bb+mp.sqrt(bb*bb-4*aa*cc))/(2*aa)])
        middle=(roots[0]+roots[1])/2;finite_normal=n+middle*d
        origin=-B*b/2+B*finite_normal/mp.sqrt(4*dot(finite_normal,B*finite_normal)/(dot(b,B*b)-4*q[0]))
        direction=B*(n+roots[0]*d)
        radius=exact_float(records[0]['source']['truncation_radius'])
        t=(-dot(origin,direction)+mp.sqrt(dot(origin,direction)**2-dot(direction,direction)*(dot(origin,origin)-radius**2)))/dot(direction,direction)
        return list(origin),list(origin+t*direction),roots
rows=[]
for layout in manifest['layouts']:
    if layout['direction'] not in (1,9):continue
    ids=layout['indices'];positions=list(map(exact_float,layout['positions']))
    with mp.workdps(100):
        # Solve exact stored-position polynomial reproduction, degrees 0..4.
        matrix=mp.matrix([[x**k for x in positions] for k in range(5)])
        weight=mp.lu_solve(matrix,mp.matrix([0,1,0,0,0]))
        values=[replay(records[i]['coefficients'],100) for i in ids]
        numeric=mp.matrix([sum(weight[j]*values[j][1][k] for j in range(5)) for k in range(3)])
        analytic=mp.matrix([exact_float(x) for x in np.asarray(records[0]['source']['vertex_velocity']).reshape(-1,3,10)[93,:,layout['direction']]])
        err=mp.norm(numeric-analytic)/max(1,mp.norm(analytic))
        rows.append(dict(direction=layout['direction'],scale=layout['scale'],high_precision_replay_vs_production_jacobian=float(err),numerical=[float(x) for x in numeric]))
base=replay(records[0]['coefficients'],100); repeat=replay(records[0]['coefficients'],180)
report=dict(scope='FIXED_INTERVAL_RAY_ARITHMETIC_DIAGNOSIS_NOT_PRODUCTION_GATE',vertex=93,anchor=92,normal_edge=[40,135],root_ordinal=0,
    normals=normals,source_manifest_sha256=hashlib.sha256((saved/'manifest.json').read_bytes()).hexdigest(),
    normal_generator_sha256=hashlib.sha256((ROOT/'src/spin/detail/FixedNormalSphere.h').read_bytes()).hexdigest(),commands=[command],
    anchor_difference=[float(base[0][k]-exact_float(records[0]['source']['vertices'][92][k])) for k in range(3)],
    ray_difference=[float(base[1][k]-exact_float(records[0]['source']['vertices'][93][k])) for k in range(3)],
    arithmetic_repeat_difference=max(float(abs(base[1][k]-repeat[1][k])) for k in range(3)),checks=rows,
    new_mapping_evaluations=0,production_admission=False,optimization_admission=False)
with (output/'report.json').open('x') as stream:json.dump(report,stream,indent=2)
print(json.dumps(report),flush=True)
