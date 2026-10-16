"""Compare exact oriented boundary chains in every supporting plane."""
import json,hashlib
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np
import gmpy2 as g
HERE=Path(__file__).resolve().parent;out=HERE/'noodles-terminal-interior-patch-01';source=HERE/'noodles-terminal-recovery-job-03/output/material'
p=json.loads((out/'candidate.json').read_text())
old=json.loads((source/'mesh_0004_rational_vertices.json').read_text());of=np.load(source/'mesh_0004_faces.npy').tolist()
def chains(rows,faces):
    v=[tuple(g.mpq(x) for x in r) for r in rows];groups=defaultdict(Counter)
    for t in faces:
        a,b,c=(v[i] for i in t);u=[b[i]-a[i] for i in range(3)];w=[c[i]-a[i] for i in range(3)]
        n=(u[1]*w[2]-u[2]*w[1],u[2]*w[0]-u[0]*w[2],u[0]*w[1]-u[1]*w[0]);lead=next(x for x in n if x)
        plane=tuple(x/lead for x in n)+(-sum(n[i]*a[i] for i in range(3))/lead,)
        for x,y in [(a,b),(b,c),(c,a)]:
            edge=tuple(sorted((x,y)));groups[plane][edge]+=1 if x<y else -1
    return {plane:{edge:s for edge,s in edges.items() if s} for plane,edges in groups.items()}
a=chains(old,of);b=chains(p['rational_vertices'],p['faces']);assert a==b,'EXACT_ORIENTED_SURFACE_CHAIN_CHANGED'
report=dict(passed=True,relation='IDENTICAL_FINITE_ORIENTED_PLANAR_SURFACE_CHAINS',supporting_planes=len(a),authority_faces=len(of),derived_exact_faces=len(p['faces']),source_sha256=hashlib.sha256((source/'mesh_0004_rational_vertices.json').read_bytes()).hexdigest(),candidate_sha256=hashlib.sha256((out/'candidate.json').read_bytes()).hexdigest())
(out/'exact-surface-chain-check.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
