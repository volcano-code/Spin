"""Read-only exact edge/triangle witnesses in frozen binary64 source parts.

This is a targeted diagnostic, not a complete self-intersection admission gate.
No native Boolean, geometry repair, optimization or simulation is invoked.
"""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
import trimesh


def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def dot(a,b):return sum(x*y for x,y in zip(a,b))


def segment_triangle(p,q,a,b,c):
    """Exact transversal strict triangle/strict segment crossing, no tolerances."""
    e1,e2,d=sub(b,a),sub(c,a),sub(q,p)
    h=cross(d,e2);det=dot(e1,h)
    if not det:return None
    s=sub(p,a);u=dot(s,h)/det;tvec=cross(s,e1)
    v=dot(d,tvec)/det;t=dot(e2,tvec)/det
    if 0<t<1 and 0<u and 0<v and u+v<1:
        return dict(segment_parameter=str(t),barycentric=[str(1-u-v),str(u),str(v)],
            point=[str(x+t*y) for x,y in zip(p,d)])
    return None


def diagnose(obj):
    obj=Path(obj);mesh=trimesh.load(obj,force='mesh',process=False)
    vertices=np.asarray(mesh.vertices);faces=np.asarray(mesh.faces)
    triangles=vertices[faces];lo=triangles.min(axis=1);hi=triangles.max(axis=1)
    exact={};witnesses=[];candidates=0
    def row(i):
        if i not in exact:exact[i]=tuple(Q(float(x)) for x in vertices[i])
        return exact[i]
    for i,face in enumerate(faces):
        possible=np.flatnonzero(np.all(hi[i]>=lo,axis=1)&np.all(hi>=lo[i],axis=1))
        for j in possible:
            if j<=i or set(face)&set(faces[j]):continue
            candidates+=1
            for edge_face,triangle_face in ((face,faces[j]),(faces[j],face)):
                a,b,c=map(row,triangle_face)
                for p,q in zip(edge_face,np.roll(edge_face,-1)):
                    hit=segment_triangle(row(p),row(q),a,b,c)
                    if hit:
                        witnesses.append(dict(faces=[int(i),int(j)],
                            edge_vertices=[int(p),int(q)],triangle_vertices=list(map(int,triangle_face)),**hit))
    return dict(schema='frozen-noodles-tube-diagnostic-v1',source=str(obj.resolve()),
        source_sha256=hashlib.sha256(obj.read_bytes()).hexdigest(),vertices=len(vertices),faces=len(faces),
        disjoint_face_aabb_candidates=candidates,strict_exact_edge_triangle_witnesses=witnesses,
        status='SELF_INTERSECTION_WITNESSED' if witnesses else 'UNKNOWN_NO_STRICT_WITNESS',
        complete_admission=False,geometry_modified=False,native_invocations=0)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--obj',type=Path,required=True)
    p.add_argument('--output',type=Path)
    a=p.parse_args();record=diagnose(a.obj)
    if a.output:
        with a.output.open('x') as f:json.dump(record,f,indent=2)
        print(json.dumps(dict(status=record['status'],witness_count=len(record['strict_exact_edge_triangle_witnesses']),output=str(a.output))))
    else:print(json.dumps(record,indent=2))
