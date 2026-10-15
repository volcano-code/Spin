"""Read-only coordinate/erasure witness diagnosis; no Boolean or mesh repair."""
import json
from pathlib import Path
import gmpy2
import numpy as np
HERE=Path(__file__).resolve().parent
source=HERE/'noodles-balanced-full-input-job-01/output/outer-rational.json'
packet=json.loads(source.read_text());scale=gmpy2.mpq(1,40)
exact=[tuple(gmpy2.mpq(x)*scale for x in row) for row in packet['vertices']]
vertices=np.asarray([[float(x) for x in row] for row in exact]);faces=np.asarray(packet['faces'])
unique,inverse=np.unique(vertices,axis=0,return_inverse=True);mapped=inverse[faces]
collapsed=np.any(np.diff(np.sort(mapped,axis=1),axis=1)==0,axis=1)
groups={};removed=set(map(int,np.flatnonzero(collapsed)));pairs=[]
for index in np.flatnonzero(~collapsed):
    triangle=mapped[index];key=tuple(sorted(map(int,triangle)))
    parity=sum(int(triangle[a]>triangle[b]) for a,b in ((0,1),(0,2),(1,2)))%2
    groups.setdefault(key,[]).append((int(index),parity))
for occurrences in groups.values():
    if len(occurrences)==2 and occurrences[0][1]!=occurrences[1][1]:
        pair=[v[0] for v in occurrences];pairs.append(pair);removed.update(pair)
    elif len(occurrences)>1:raise ValueError('NONCANCELLING_DUPLICATE_FACE')
keep=[i for i in range(len(faces)) if i not in removed]
edges=set()
for face in mapped[keep]:
    edges.update(tuple(sorted((int(face[i]),int(face[(i+1)%3])))) for i in range(3))
used=set(map(int,mapped[keep].ravel()))
spans=[max(p[i] for p in exact)-min(p[i] for p in exact) for i in range(3)]
squared=sum(x*x for x in spans)
with gmpy2.context(precision=256,round=gmpy2.RoundDown):length=float(gmpy2.sqrt(squared))
while gmpy2.mpq(length)**2>squared:length=float(np.nextafter(length,0))
rows=[]
for index in sorted(removed):
    candidates=sorted(set(map(int,mapped[index]))&used)
    points=[exact[int(v)] for v in faces[index]]
    vertex_bound=min(max(sum(abs(x-gmpy2.mpq(float(y))) for x,y in zip(p,unique[v]))
        for p in points) for v in candidates)
    segments=[]
    for a in candidates:
        for b in candidates:
            if a>=b or (a,b) not in edges:continue
            start=tuple(gmpy2.mpq(float(x)) for x in unique[a]);end=tuple(gmpy2.mpq(float(x)) for x in unique[b])
            direction=tuple(y-x for x,y in zip(start,end));den=sum(x*x for x in direction)
            bound=gmpy2.mpq(0)
            for point in points:
                t=min(gmpy2.mpq(1),max(gmpy2.mpq(0),sum((x-y)*d for x,y,d in zip(point,start,direction))/den))
                bound=max(bound,sum(abs(x-(a+t*d)) for x,a,d in zip(point,start,direction)))
            segments.append(dict(edge=[a,b],distance_l1_upper_bound=float(bound)))
    rows.append(dict(face=index,collapsed=bool(collapsed[index]),vertex_bound_m=float(vertex_bound),
        surviving_edge_witnesses=segments))
worst=sorted(rows,key=lambda x:x['vertex_bound_m'],reverse=True)[:4]
print(json.dumps(dict(role='REPRESENTATION_WITNESS_DIAGNOSIS_ONLY',removed_faces=len(removed),
    collapsed_faces=int(np.count_nonzero(collapsed)),opposite_pairs=len(pairs),length_m=length,
    maximum_boundary_error_allowed_m=length/1e6,worst_faces=worst,admission=False),indent=2))
