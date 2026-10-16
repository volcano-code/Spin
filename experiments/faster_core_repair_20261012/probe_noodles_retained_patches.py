"""Exact planar fan candidates for retained boundary; no admission by construction."""
import json,sys
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
import _make_it_stand as native
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'python'))
from make_it_stand.verification._coplanar_boundary import _normal,_dot,_sub,_boundary
source=HERE/'noodles-terminal-recovery-job-03/output/material'
e=[tuple(Q(x) for x in r) for r in json.loads((source/'mesh_0008_rational_vertices.json').read_text())]
f=np.load(source/'mesh_0008_faces.npy').tolist()
pairs=json.loads((HERE/'noodles-remaining-parts-audit-01.json').read_text())['retained']['intersecting_pairs']
adj={}
for i,t in enumerate(f):
    for a,b in zip(t,t[1:]+t[:1]):adj.setdefault(tuple(sorted((a,b))),[]).append(i)
covered=set();removed=set();added=[];proofs=[];rejected=[]
for seed in sorted({i for p in pairs for i in p}):
    if seed in covered:continue
    n=_normal(e,f[seed]);origin=e[f[seed][0]];seen={seed};todo=[seed]
    while todo:
        i=todo.pop()
        for a,b in zip(f[i],f[i][1:]+f[i][:1]):
            for j in adj[tuple(sorted((a,b)))]:
                if j not in seen and all(_dot(n,_sub(e[k],origin))==0 for k in f[j]):seen.add(j);todo.append(j)
    covered.update(seen);old=[f[i] for i in sorted(seen)];boundary=_boundary(old)
    edges=[(a,b) if s==1 else (b,a) for (a,b),s in boundary.items()]
    if not edges or not all(abs(s)==1 for s in boundary.values()) or len({a for a,b in edges})!=len(edges) or len({b for a,b in edges})!=len(edges):
        rejected.append(dict(seed=seed,reason='NOT_SIMPLE_DISK'));continue
    bv=sorted({i for edge in edges for i in edge});center=tuple(sum(e[i][j] for i in bv)/len(bv) for j in range(3))
    k=len(e);e.append(center);new=[(a,b,k) for a,b in edges]
    if not all(_dot(n,_normal(e,t))>0 for t in old+new):
        e.pop();rejected.append(dict(seed=seed,reason='CENTER_NOT_STRICT_KERNEL'));continue
    assert _boundary(old)==_boundary(new)
    removed.update(seen);added.extend(new);proofs.append(dict(seed=seed,old_indices=sorted(seen),center=[str(x) for x in center]))
nf=[t for i,t in enumerate(f) if i not in removed]+added
v=np.array([[float(x) for x in r] for r in e]);a=dict(native.audit_exact_triangle_boundary(v,np.array(nf,dtype=np.int32),'AnyNonzero'))
out=HERE/'noodles-retained-interior-patches-01';out.mkdir(exist_ok=False)
(out/'candidate.json').write_text(json.dumps(dict(vertices=v.tolist(),faces=nf,rational_vertices=[[str(x) for x in r] for r in e])))
(out/'report.json').write_text(json.dumps(dict(audit=a,proofs=proofs,rejected=rejected,joint_partition_certified=False),default=lambda x:x.tolist()))
print(json.dumps(dict(patches=len(proofs),rejected=len(rejected),valid=a['boundary_valid'],pairs=a['intersecting_pairs'])),flush=True)
