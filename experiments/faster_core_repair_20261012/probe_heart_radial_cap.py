"""Offline radial cap experiment on the saved rejected interface."""
import json
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
import _make_it_stand as native


def ears(ids, vertices, axis):
    keep=[a for a in range(3) if a!=int(np.argmax(abs(axis)))]
    xy={i:[Q(float(vertices[i,a])) for a in keep] for i in ids}
    def turn(a,b,c):
        x,y,z=xy[a],xy[b],xy[c]
        return (y[0]-x[0])*(z[1]-x[1])-(y[1]-x[1])*(z[0]-x[0])
    area=sum(xy[a][0]*xy[b][1]-xy[a][1]*xy[b][0] for a,b in zip(ids,ids[1:]+ids[:1]))
    if not area:raise ValueError('ZERO_PROJECTED_AREA')
    sign=1 if area>0 else -1
    remaining=list(ids);faces=[]
    while len(remaining)>3:
        for k,b in enumerate(remaining):
            a,c=remaining[k-1],remaining[(k+1)%len(remaining)]
            if sign*turn(a,b,c)<=0:continue
            if any(sign*turn(a,b,x)>=0 and sign*turn(b,c,x)>=0 and sign*turn(c,a,x)>=0
                   for x in remaining if x not in (a,b,c)):continue
            faces.append([c,b,a]);remaining.pop(k);break
        else:raise ValueError('NO_CERTIFIED_EAR')
    faces.append(remaining[::-1]);return faces


def main():
    root=Path(__file__).resolve().parent
    error=json.loads((root/'heart-source-failure-replay-job-03/stdout.log').read_text())['error']
    fields=dict(a.split('=',1) for a in error.split(';')[1:] if '=' in a)
    def rows(k,t):return [[t(x) for x in a.split(',') if x] for a in fields[k].split('|') if a]
    vertices=np.array(rows('witness_vertices',float));oldfaces=np.array(rows('witness_faces',int))
    loops=rows('witness_loops',int);count=int(fields['interface_face_count'])
    faces=oldfaces[:count].tolist();offset=int(oldfaces[:count].max())+1
    q=np.array(json.loads((root/'heart-seed1-source-failure-registration-01.json').read_text())['event']['quadric_coefficients'])
    a=np.array([[q[4],q[9],q[8]],[q[9],q[5],q[7]],[q[8],q[7],q[6]]])
    center=-np.linalg.solve(a,q[1:4])/2;_,axes=np.linalg.eigh(a)
    request=json.loads((root/'heart-certified-formal-seed-1-job-01/output/run/configuration.json').read_text())['request']
    design=(np.asarray(request['design_vertices'])-np.asarray(request['target_origin']))@np.asarray(request['target_axes']).T
    def exact_dot(point,axis):return sum(Q(float(x))*Q(float(y)) for x,y in zip(point,axis))
    margins=[]
    for loop in loops:
        points=vertices[loop]-center;axis=axes[:,2]
        if np.mean(points@axis)<0:axis=-axis
        if min(points@axis)<=0:raise ValueError('RADIAL_PROJECTION_NOT_SEPARATED')
        upper=list(range(offset,offset+len(loop)));offset+=len(loop)
        vertices[upper]=center+points*(2*float(fields['truncation_radius'])/(points@axis))[:,None]
        support=max(exact_dot(point,axis) for point in design)
        margin=min(exact_dot(vertices[i],axis)-support for i in loop+upper)
        if margin<=0:raise ValueError('CAP_NOT_SEPARATED_FROM_DESIGN')
        margins.append(float(margin))
        for k,b in enumerate(loop):
            nxt=(k+1)%len(loop)
            faces.extend([[loop[nxt],b,upper[k]],[loop[nxt],upper[k],upper[nxt]]])
        faces.extend(ears(upper,vertices,axis))
    report=native.audit_exact_triangle_boundary(vertices,np.array(faces,dtype=np.int32),'OutwardPositive')
    print(json.dumps({k:v for k,v in report.items() if k in ('boundary_valid','closed','degenerate_faces','intersecting_pairs','intersection_error')}))
    print(json.dumps(dict(exact_cap_halfspace_margins=margins,derivative_verified=False,production_adopted=False)))


if __name__=='__main__':main()
