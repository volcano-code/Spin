"""Accept exact patch flips only after strict global intersection decrease."""
import json,sys
from pathlib import Path
from fractions import Fraction
import numpy as np
import _make_it_stand as native
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification._coplanar_boundary import propose_coplanar_diagonal_repairs
s=HERE/'noodles-terminal-recovery-job-03/output/material'
e=[[Fraction(x) for x in r] for r in json.loads((s/'mesh_0004_rational_vertices.json').read_text())];v=[[float(x) for x in r] for r in e]
f=json.loads((HERE/'noodles-terminal-coplanar-03/proposal-2.json').read_text())['faces']
out=HERE/'noodles-terminal-greedy-01';out.mkdir(exist_ok=False)
def audit(f):return dict(native.audit_exact_triangle_boundary(np.asarray(v),np.asarray(f,dtype=np.int32),'AnyNonzero'))
for step in range(30):
    a=audit(f);count=a['intersecting_pair_count']
    (out/f'accepted-{step}.json').write_text(json.dumps(dict(faces=f,audit=a),default=lambda x:x.tolist()))
    print('ACCEPTED',step,count,flush=True)
    if a['boundary_valid']:break
    adj={}
    for i,t in enumerate(f):
        for x,y in zip(t,t[1:]+t[:1]):adj.setdefault(tuple(sorted((x,y))),[]).append(i)
    hints=set()
    for i in {i for p in a['intersecting_pairs'] for i in p}:
        for x,y in zip(f[i],f[i][1:]+f[i][:1]):
            for j in adj[tuple(sorted((x,y)))]:
                if i!=j:hints.add(tuple(sorted((i,j))))
    changed=False
    for pair in sorted(hints):
        p=propose_coplanar_diagonal_repairs(e,v,f,[list(pair)],face_regions=['material']*len(f),allow_positive_display_patch=True)
        if not p['patches']:continue
        b=audit(p['faces'])
        if b['intersecting_pair_count']<count and not b['degenerate_faces']:
            f=p['faces'];changed=True;break
    if not changed:
        print('NO_SINGLE_EXACT_PATCH_DESCENT',flush=True);break
