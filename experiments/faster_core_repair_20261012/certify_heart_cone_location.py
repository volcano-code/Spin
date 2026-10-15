"""Exact interval apex enclosure and stable ray parity on production input."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp
from audit_heart_quadric_transition import matrix


def iv(x):
    return x if isinstance(x, tuple) else (F(x), F(x))


def add(a,b):
    a,b=iv(a),iv(b);return a[0]+b[0],a[1]+b[1]


def neg(a):
    a=iv(a);return -a[1],-a[0]


def sub(a,b): return add(a,neg(b))


def mul(a,b):
    a,b=iv(a),iv(b);p=[x*y for x in a for y in b];return min(p),max(p)


def div(a,b):
    b=iv(b)
    if b[0]<=0<=b[1]: raise ValueError('INTERVAL_DIVISOR_CONTAINS_ZERO')
    return mul(a,(1/b[1],1/b[0]))


def dot(a,b):
    s=iv(0)
    for x,y in zip(a,b):s=add(s,mul(x,y))
    return s


def cross(a,b):
    return [sub(mul(a[1],b[2]),mul(a[2],b[1])),
            sub(mul(a[2],b[0]),mul(a[0],b[2])),
            sub(mul(a[0],b[1]),mul(a[1],b[0]))]


def polynomial(expr,t,interval):
    result=iv(0)
    for c in sp.Poly(expr,t).all_coeffs():result=add(mul(result,interval),F(str(c)))
    return result


def parity(vertices,faces,point):
    vertices=[[F(x) for x in v] for v in vertices]
    d=[F(1),F(17),F(257)];count=0
    for face in faces:
        a,b,c=[vertices[i] for i in face]
        if any(point[k][0]>max(a[k],b[k],c[k]) for k in range(3)):continue
        e1=[b[k]-a[k] for k in range(3)];e2=[c[k]-a[k] for k in range(3)]
        h=cross(d,e2);det=dot(e1,h)
        offset=[sub(point[k],a[k]) for k in range(3)]
        if det==(0,0):
            plane=dot(cross(e1,e2),offset)
            if plane[0]<=0<=plane[1]:raise ValueError('RAY_COPLANAR_INTERVAL')
            continue
        u=div(dot(offset,h),det)
        if u[1]<0 or u[0]>1:continue
        q=cross(offset,e1);v=div(dot(d,q),det)
        uv=add(u,v)
        if v[1]<0 or uv[0]>1:continue
        distance=div(dot(e2,q),det)
        if distance[1]<0:continue
        if distance[0]>0 and u[0]>0 and v[0]>0 and uv[1]<1:count+=1
        else:raise ValueError('RAY_BOUNDARY_INTERVAL_AMBIGUOUS')
    return dict(ray=[str(x) for x in d], crossings=count, inside=bool(count%2),
                whole_enclosure_same_parity=True)


def main(output):
    here=Path(__file__).resolve().parent;run=here/'heart-guarded-formal-seed-0-job-01/output/run'
    proof_path=here/'heart-exact-quadric-transition-01.json';proof=json.loads(proof_path.read_text())
    config_path=run/'configuration.json';request=json.loads(config_path.read_text())['request']
    qs=[]
    for state in proof['states']:
        path=run/f"evaluation-{state['evaluation']:06d}.json"
        if hashlib.sha256(path.read_bytes()).hexdigest()!=state['material_sha256']:raise ValueError('MATERIAL_CHANGED')
        qs.append([sp.Rational(x) for x in json.loads(path.read_text())['quadric_coefficients']])
    t=sp.Symbol('t');A,b,_=matrix([(1-t)*a+t*z for a,z in zip(*qs)])
    root=proof['degeneracy_parameter_intervals'][0];interval=(F(root[0]),F(root[1]))
    den=polynomial(A.det(),t,interval);num=-A.adjugate()*b/2
    center=[div(polynomial(x,t,interval),den) for x in num]
    point=[add(dot(center,request['target_axes'][k]),request['target_origin'][k]) for k in range(3)]
    results={name:parity(request[name+'_vertices'],request[name+'_faces'],point) for name in ('design','outer')}
    result=dict(schema='heart-cone-exact-location-v1',
        proof_sha256=hashlib.sha256(proof_path.read_bytes()).hexdigest(),
        configuration_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest(),
        center_target_interval=[[str(x) for x in v] for v in center], queries=results,
        scope='Production binary64 request entity, exact rational operations; input validity inherited from frozen admission',
        apex_excluded_from_design=not results['design']['inside'],
        stagnation_root_cause_certified=False,case_complete=False)
    with output.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False)
    print(json.dumps(results),flush=True)


if __name__=='__main__':main(Path(sys.argv[1]))
