"""Inspect captured tapes only: no Boolean, optimization, or mesh repair."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

def cross(a,b):
    return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]

def incidence(faces):
    edges=Counter();oriented=Counter();links=defaultdict(list)
    for a,b,c in faces:
        for v,x,y in ((a,b,c),(b,c,a),(c,a,b)):
            edges[tuple(sorted((v,x)))]+=1
            oriented[(v,x)]+=1
            links[v].append((x,y))
    bad=[]
    for vertex, pairs in links.items():
        graph=defaultdict(list)
        for a,b in pairs:graph[a].append(b);graph[b].append(a)
        remaining=set(graph);components=0
        while remaining:
            stack=[remaining.pop()];components+=1
            while stack:
                for other in graph[stack.pop()]:
                    if other in remaining:remaining.remove(other);stack.append(other)
        if components!=1 or any(len(v)!=2 for v in graph.values()):
            bad.append(dict(vertex=vertex,link_components=components,
                            bad_link_degrees={str(k):len(v) for k,v in graph.items() if len(v)!=2}))
    return dict(vertices_used=len(links),edges=len(edges),faces=len(faces),
                edges_not_two=sum(n!=2 for n in edges.values()),
                inconsistent_oriented_edges=sum(oriented[(a,b)]!=1 or oriented[(b,a)]!=1 for a,b in edges),
                bad_vertex_links=bad)

def material_report(snapshot):
    vertices=np.asarray(snapshot['material_vertices']);faces=np.asarray(snapshot['material_faces'])
    exact=snapshot['exact_material_vertices']
    area=np.cross(vertices[faces[:,1]]-vertices[faces[:,0]],vertices[faces[:,2]]-vertices[faces[:,0]])
    float_zero=np.flatnonzero(np.all(area==0,axis=1)).tolist()
    float_groups=defaultdict(list);exact_groups=defaultdict(list)
    for i,row in enumerate(vertices):float_groups[tuple(row)].append(i)
    for i,row in enumerate(exact):exact_groups[tuple(row)].append(i)
    aliases=[v for v in float_groups.values() if len(v)>1]
    duplicates=[v for v in exact_groups.values() if len(v)>1]
    witnesses=[]
    for f in float_zero:
        ids=faces[f].tolist();p=[[Q(x) for x in exact[i]] for i in ids]
        exact_area=cross([p[1][k]-p[0][k] for k in range(3)],[p[2][k]-p[0][k] for k in range(3)])
        represented=[[Q(float(x)) for x in vertices[i]] for i in ids]
        represented_area=cross([represented[1][k]-represented[0][k] for k in range(3)],
                               [represented[2][k]-represented[0][k] for k in range(3)])
        witnesses.append(dict(face=f,vertices=ids,birth_face=snapshot['material_birth_faces'][f],
                              exact_area_zero=all(x==0 for x in exact_area),
                              exact_double_area=[str(x) for x in exact_area],
                              binary64_exact_area_zero=all(x==0 for x in represented_area),
                              binary64_exact_double_area=[str(x) for x in represented_area]))
    return dict(stage=snapshot['stage'],incidence=incidence(faces.tolist()),
                float_aliases=aliases,exact_duplicate_vertices=duplicates,
                float_zero_area_faces_exact_witnesses=witnesses)

def source_report(packet):
    source=packet['snapshots'][0];vertices=np.asarray(source['source_vertices']);faces=source['source_faces']
    q=[Q(float(x)) for x in packet['normalized_coefficients']]
    rows=((4,9,8),(9,5,7),(8,7,6));positive=[];negative=[];zero=[]
    witnesses=[]
    for f,ids in enumerate(faces):
        if source['source_face_is_cap'][f]:continue
        p=[[Q(float(x)) for x in vertices[i]] for i in ids]
        area=cross([p[1][k]-p[0][k] for k in range(3)],[p[2][k]-p[0][k] for k in range(3)])
        dots=[]
        for v in p:
            gradient=[q[k+1]+2*sum(q[rows[k][j]]*v[j] for j in range(3)) for k in range(3)]
            dots.append(sum(area[k]*gradient[k] for k in range(3)))
        centroid_dot=sum(dots)/3
        (positive if centroid_dot>0 else negative if centroid_dot<0 else zero).append(f)
        if centroid_dot>0:
            witnesses.append(dict(face=f,vertices=ids,gradient_dots=[str(v) for v in dots],
                                  rays=[source['source_projective_rays'][i] for i in ids],
                                  normal_edges=[source['source_normal_edges'][i] for i in ids],
                                  fixed_normal_ids=[source['source_normal_vertices'][i] for i in ids]))
    return dict(incidence=incidence(faces),positive_faces=positive,negative_faces=negative,zero_faces=zero,
                positive_face_witnesses=witnesses,
                note='Signs are exact for the recorded binary64 source/coefficients; no majority orientation decision.')

if __name__=='__main__':
    directory=Path(sys.argv[1]);output=Path(sys.argv[2])
    inputs={name:directory/(name+'.json') for name in ('source_side','material_boundary')}
    source=json.loads(inputs['source_side'].read_text());material=json.loads(inputs['material_boundary'].read_text())
    report=dict(schema='same-call-failure-tape-analysis-v1',production_admission=False,
                source_sha256={k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in inputs.items()},
                source=source_report(source),
                materials=[material_report(s) for s in material['snapshots'] if s['stage'].startswith('boolean-')])
    with output.open('x') as stream:json.dump(report,stream,indent=2)
    print(json.dumps(dict(source_positive_faces=report['source']['positive_faces'],
                         materials=[dict(stage=m['stage'],incidence=m['incidence'],
                                         aliases=len(m['float_aliases']),exact_duplicates=len(m['exact_duplicate_vertices']),
                                         float_zero_faces=len(m['float_zero_area_faces_exact_witnesses'])) for m in report['materials']]),indent=2))
