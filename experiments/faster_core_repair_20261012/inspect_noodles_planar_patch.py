"""Extract maximal exact coplanar components at the remaining display failures."""
import json,sys
from pathlib import Path
from fractions import Fraction
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'python'))
from make_it_stand.verification._coplanar_boundary import _normal,_dot,_sub,_boundary
e=[[Fraction(x) for x in p] for p in json.loads((HERE/'noodles-terminal-recovery-job-03/output/material/mesh_0004_rational_vertices.json').read_text())]
f=json.loads((HERE/'noodles-terminal-greedy-01/accepted-14.json').read_text())['faces']
adj={}
for i,t in enumerate(f):
    for a,b in zip(t,t[1:]+t[:1]):adj.setdefault(tuple(sorted((a,b))),[]).append(i)
components=[]
for seed in [87559,87560,91123,91124]:
    n=_normal(e,f[seed]);origin=e[f[seed][0]]; seen={seed}; todo=[seed]
    while todo:
        i=todo.pop()
        for a,b in zip(f[i],f[i][1:]+f[i][:1]):
            for j in adj[tuple(sorted((a,b)))]:
                if j not in seen and all(_dot(n,_sub(e[k],origin))==0 for k in f[j]):
                    seen.add(j);todo.append(j)
    boundary=_boundary([f[i] for i in seen]);bv={i for edge in boundary for i in edge};allv={k for i in seen for k in f[i]}
    entry=dict(seed=seed,indices=sorted(seen),faces=[f[i] for i in sorted(seen)],boundary=[[a,b,s] for (a,b),s in boundary.items()],interior_vertices=sorted(allv-bv),vertices={str(k):[str(x) for x in e[k]] for k in allv})
    components.append(entry)
    print(seed,'faces',len(seen),'boundary',len(boundary),'interior',sorted(allv-bv),flush=True)
(HERE/'noodles-terminal-planar-components-01.json').write_text(json.dumps(components))
