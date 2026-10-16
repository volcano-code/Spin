"""Certify outer = whole + cavity using identical chains and exact nesting rays."""
import json,hashlib
from pathlib import Path
import numpy as np
import gmpy2 as g
from collections import Counter
HERE=Path(__file__).resolve().parent;out=HERE/'noodles-terminal-interior-patch-01'
p=json.loads((out/'candidate.json').read_text());v=p['vertices'];f=p['faces']
a=json.loads((out/'report.json').read_text())['audit']
assert a['boundary_valid'] and a['component_volume_signs']==[1,-1,-1]
components=[[t for t,c in zip(f,a['face_components']) if c==k] for k in range(3)]
source=HERE/'noodles-terminal-recovery-job-03/output/material'
ov=np.load(source/'mesh_0000_vertices.npy').tolist();of=np.load(source/'mesh_0000_faces.npy').tolist()
def key(vertices,t):
    points=[tuple(vertices[i]) for i in t]
    return min(tuple(points[i:]+points[:i]) for i in range(3))
assert Counter(key(ov,t) for t in of)==Counter(key(v,t) for t in components[0]),'OUTER_IDENTITY_CHANGED'
q=[[g.mpq(float(x)) for x in row] for row in v]
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def inside(point,faces):
    for direction in [(1,2,3),(2,3,5),(3,5,7)]:
        hits=0;ambiguous=False
        for tri in faces:
            x,y,z=[q[i] for i in tri];e1=sub(y,x);e2=sub(z,x);h=cross(direction,e2);det=dot(e1,h)
            if det==0:continue
            s=sub(point,x);u=dot(s,h)/det
            if u<0 or u>1:continue
            tvec=cross(s,e1);w=dot(direction,tvec)/det
            if w<0 or u+w>1:continue
            t=dot(e2,tvec)/det
            if t<=0:continue
            if u==0 or w==0 or u+w==1:ambiguous=True;break
            hits+=1
        if not ambiguous:return dict(inside=bool(hits%2),direction=direction,positive_crossings=hits)
    raise ValueError('NO_NONDEGENERATE_EXACT_RAY')
nesting={}
for c in [1,2]:
    point=q[components[c][0][0]]
    for target in [0,3-c]:
        witness=inside(point,components[target]);nesting[f'{c}_in_{target}']=witness
        assert witness['inside']==(target==0),'INVALID_REGION_NESTING'
cavity=[list(reversed(t)) for c in [1,2] for t in components[c]]
used=sorted({i for t in cavity for i in t});mapping={old:new for new,old in enumerate(used)}
path=out/'cavity.obj'
with path.open('x') as stream:
    for i in used:stream.write('v '+' '.join(format(x,'.17g') for x in v[i])+'\n')
    for t in cavity:stream.write('f '+' '.join(str(mapping[i]+1) for i in t)+'\n')
report=dict(outer_identity='EXACT_SAME_ORIENTED_TRIANGLES',global_boundary_valid=True,shared_cavity_interface='IDENTICAL_OPPOSITELY_ORIENTED_TRIANGLES',nesting=nesting,whole_cavity_outer_partition_passed=True,fixed_retained_partition='PENDING',qualified_full_material_package=False,cavity_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
(out/'whole-partition-check.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
