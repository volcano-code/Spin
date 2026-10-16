"""Exact oriented planar disk replacement with an authenticated interior vertex."""
import json,sys
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
import _make_it_stand as native
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'python'))
from make_it_stand.verification._coplanar_boundary import _normal,_dot,_sub,_boundary
e=[[Q(x) for x in p] for p in json.loads((HERE/'noodles-terminal-recovery-job-03/output/material/mesh_0004_rational_vertices.json').read_text())]
f=json.loads((HERE/'noodles-terminal-greedy-01/accepted-14.json').read_text())['faces']
components=json.loads((HERE/'noodles-terminal-planar-components-01.json').read_text())
out=HERE/'noodles-terminal-interior-patch-01';out.mkdir(exist_ok=False)
proofs=[]
for component in [components[0],components[2]]:
    ids=component['indices']; old=[f[i] for i in ids]; boundary=_boundary(old)
    edges=[(a,b) if s==1 else (b,a) for (a,b),s in boundary.items()]
    assert all(abs(s)==1 for s in boundary.values())
    assert len({a for a,b in edges})==len(edges)==len({b for a,b in edges})
    bv=sorted({a for edge in edges for a in edge})
    center=tuple(sum(e[i][j] for i in bv)/len(bv) for j in range(3))
    n=_normal(e,old[0]); origin=e[old[0][0]]
    assert _dot(n,_sub(center,origin))==0
    k=len(e);e.append(center);new=[(a,b,k) for a,b in edges]
    assert _boundary(old)==_boundary(new)
    assert all(_dot(n,_normal(e,t))>0 for t in old+new),'CENTER_OUTSIDE_KERNEL'
    # Positive oriented fan and identical simple disk boundary certify exact surface equality.
    proofs.append(dict(old_indices=ids,old_faces=old,new_faces=new,center=[str(x) for x in center],relation='EXACT_ORIENTED_PLANAR_DISK_EQUALITY'))
    removed=set(ids);f=[t for i,t in enumerate(f) if i not in removed]+[list(t) for t in new]
    # Original face indices of the next component shift after replacement.
    for other in components:
        if other is component:continue
        other['indices']=[i-sum(j<i for j in removed) for i in other['indices']]
v=np.array([[float(x) for x in p] for p in e]);a=dict(native.audit_exact_triangle_boundary(v,np.array(f,dtype=np.int32),'AnyNonzero'))
(out/'report.json').write_text(json.dumps(dict(audit=a,proofs=proofs,qualified_delivery=False),default=lambda x:x.tolist()))
(out/'candidate.json').write_text(json.dumps(dict(vertices=v.tolist(),faces=f,rational_vertices=[[str(x) for x in p] for p in e])))
print(a,flush=True)
