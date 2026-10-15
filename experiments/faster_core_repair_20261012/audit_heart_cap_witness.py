"""Offline cap witness analysis; exact segment tests on represented projection."""
from fractions import Fraction as Q
import json
from pathlib import Path
import sys
import numpy as np


def turn(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def strict_crossings(points):
    exact = [[Q(float(x)) for x in row] for row in points]
    result = []
    n = len(exact)
    for i in range(n):
        for j in range(i+1, n):
            if j in (i, (i+1)%n) or (j+1)%n == i:
                continue
            a,b,c,d = exact[i],exact[(i+1)%n],exact[j],exact[(j+1)%n]
            if turn(a,b,c)*turn(a,b,d)<0 and turn(c,d,a)*turn(c,d,b)<0:
                result.append([i,j])
    return result


def audit(path):
    data=json.loads(Path(path).read_text())
    fields=dict(x.split('=',1) for x in data['error'].split(';')[1:] if '=' in x)
    def rows(name, cast):
        return [[cast(x) for x in row.split(',') if x] for row in fields[name].split('|') if row]
    vertices=np.array(rows('witness_vertices',float))
    faces=np.array(rows('witness_faces',int))
    loops=rows('witness_loops',int)
    interface_faces=int(fields['interface_face_count'])
    offset=int(faces[:interface_faces].max())+1
    reports=[]
    for loop in loops:
        original=vertices[loop]
        upper=vertices[offset:offset+len(loop)];offset+=len(loop)
        area=np.cross(original,np.roll(original,-1,axis=0)).sum(0)
        axis=area/np.linalg.norm(area)
        if axis@original.sum(0)<0:axis=-axis
        omitted=int(np.argmax(abs(axis)));keep=[i for i in range(3) if i!=omitted]
        reports.append(dict(vertex_count=len(loop),axis=axis.tolist(),
            projected_upper_strict_crossings=strict_crossings(upper[:,keep]),
            minimum_extrusion_offset=float(np.min((upper-original)@axis))))
    return dict(schema='heart-cap-witness-analysis-v1',loops=reports,
        scope='Exact proper crossings of stored binary64 projected cap edges; contacts and 3D shell intersections require separate audit')


if __name__=='__main__':
    print(json.dumps(audit(sys.argv[1]),indent=2))
