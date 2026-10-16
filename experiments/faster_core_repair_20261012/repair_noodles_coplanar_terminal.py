"""Apply existing exact patch equivalence rules; audit every candidate."""
import json,sys
from fractions import Fraction
from pathlib import Path
from collections import Counter
import numpy as np
import _make_it_stand as native
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification._coplanar_boundary import propose_coplanar_diagonal_repairs
source=HERE/'noodles-terminal-recovery-job-03/output/material'
out=HERE/'noodles-terminal-coplanar-03';out.mkdir(exist_ok=False)
exact=[[Fraction(x) for x in row] for row in json.loads((source/'mesh_0004_rational_vertices.json').read_text())]
v=[[float(x) for x in row] for row in exact]
faces=json.loads((HERE/'noodles-terminal-coplanar-02/proposal-4.json').read_text())['faces']
seen=set()
for step in range(50):
    import hashlib
    key=hashlib.sha256(np.asarray(faces,dtype=np.int32).tobytes()).hexdigest()
    if key in seen:
        print('CONNECTIVITY_CYCLE',step,flush=True);break
    seen.add(key)
    audit=dict(native.audit_exact_triangle_boundary(np.asarray(v),np.asarray(faces,dtype=np.int32),'AnyNonzero'))
    (out/f'audit-{step}.json').write_text(json.dumps(audit,default=lambda x:x.tolist()))
    if audit['boundary_valid']:
        with (out/'complete-candidate.obj').open('x') as stream:
            for row in v:stream.write('v '+' '.join(format(x,'.17g') for x in row)+'\n')
            for row in faces:stream.write('f '+' '.join(str(x+1) for x in row)+'\n')
        print('BOUNDARY_VALID',step,flush=True);break
    pairs=list(audit['intersecting_pairs']);adj={}
    for i,f in enumerate(faces):
        for a,b in zip(f,f[1:]+f[:1]):adj.setdefault(tuple(sorted((a,b))),[]).append(i)
    affected={i for pair in pairs for i in pair}|set(audit['degenerate_faces'])
    for i in affected:
        f=faces[i]
        for a,b in zip(f,f[1:]+f[:1]):
            pairs.extend([i,j] for j in adj[tuple(sorted((a,b)))] if i!=j)
    proposal=propose_coplanar_diagonal_repairs(exact,v,faces,pairs,face_regions=['material']*len(faces),allow_positive_display_patch=True)
    (out/f'proposal-{step}.json').write_text(json.dumps(proposal))
    print(step,'intersections',audit['intersecting_pair_count'],'patches',len(proposal['patches']),dict(Counter(x['reason'] for x in proposal['skipped_pairs'])),flush=True)
    if not proposal['patches']:break
    faces=proposal['faces']
