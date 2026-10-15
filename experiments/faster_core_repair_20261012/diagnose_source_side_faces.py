"""Read-only exact binary64 face witnesses; does not repair or admit geometry."""
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys


def main(source,output):
    data=json.loads(source.read_text());s=data['snapshots'][0]
    q=list(map(Q,data['normalized_coefficients']))
    vertices=[list(map(Q,p)) for p in s['source_vertices']]
    A=((q[4],q[9],q[8]),(q[9],q[5],q[7]),(q[8],q[7],q[6]))
    def dot(a,b):return sum(x*y for x,y in zip(a,b))
    def mat(v):return [dot(row,v) for row in A]
    def value(v):return q[0]+dot(q[1:4],v)+dot(v,mat(v))
    rows=[]
    for index,face in enumerate(s['source_faces']):
        if s['source_face_is_cap'][index]:continue
        a,b,c=[vertices[i] for i in face]
        u=[y-x for x,y in zip(a,b)];v=[y-x for x,y in zip(a,c)]
        n=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]]
        center=[(a[k]+b[k]+c[k])/3 for k in range(3)]
        alignment=dot(n,[2*x+y for x,y in zip(mat(center),q[1:4])])
        rows.append(dict(face=index,vertices=face,alignment_exact=str(alignment),
            alignment_sign=1 if alignment>0 else -1 if alignment<0 else 0,
            vertex_quadric_values=[float(value(p)) for p in (a,b,c)],
            projective_ray_provenance=[s['source_projective_rays'][i] for i in face]))
    report=dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        coefficients=data['normalized_coefficients'],rows=rows,
        positive_faces=[r['face'] for r in rows if r['alignment_sign']>0],
        scope='Exact represented face/gradient signs; not inside/outside certification',
        geometry_modified=False,core_admission=False)
    with output.open('x') as stream:json.dump(report,stream)


if __name__=='__main__':main(Path(sys.argv[1]),Path(sys.argv[2]))
