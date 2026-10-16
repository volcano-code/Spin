"""Search two exact-equivalent flips as one larger patch, auditing globally."""
import json,sys,hashlib
from pathlib import Path
from fractions import Fraction
import numpy as np
import _make_it_stand as native
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification._coplanar_boundary import propose_coplanar_diagonal_repairs
s=HERE/'noodles-terminal-recovery-job-03/output/material'
e=[[Fraction(x) for x in row] for row in json.loads((s/'mesh_0004_rational_vertices.json').read_text())];v=[[float(x) for x in row] for row in e]
initial=json.loads((HERE/'noodles-terminal-greedy-01/accepted-14.json').read_text())['faces']
out=HERE/'noodles-terminal-patch-search-01';out.mkdir(exist_ok=False)
def hints(f,focus):
    adj={}
    for i,t in enumerate(f):
        for a,b in zip(t,t[1:]+t[:1]):adj.setdefault(tuple(sorted((a,b))),[]).append(i)
    pairs=set()
    for _ in range(2):
        expanded=set(focus)
        for i in focus:
            for a,b in zip(f[i],f[i][1:]+f[i][:1]):
                ids=adj[tuple(sorted((a,b)))];expanded.update(ids)
                if len(ids)==2:pairs.add(tuple(sorted(ids)))
        focus=expanded
    return sorted(pairs)
def audit(f):return dict(native.audit_exact_triangle_boundary(np.asarray(v),np.asarray(f,dtype=np.int32),'AnyNonzero'))
seen=set();count=0
for pair in hints(initial,{87559,87560,91123}):
    p=propose_coplanar_diagonal_repairs(e,v,initial,[list(pair)],face_regions=['material']*len(initial),allow_positive_display_patch=True)
    if not p['patches']:continue
    first=p['faces']
    for second in hints(first,set(pair)|{87559,87560,91123}):
        q=propose_coplanar_diagonal_repairs(e,v,first,[list(second)],face_regions=['material']*len(first),allow_positive_display_patch=True)
        if not q['patches']:continue
        key=hashlib.sha256(np.asarray(q['faces'],dtype=np.int32).tobytes()).hexdigest()
        if key in seen:continue
        seen.add(key);a=audit(q['faces']);count+=1
        print(count,pair,second,a['intersecting_pair_count'],flush=True)
        if a['boundary_valid']:
            (out/'passed.json').write_text(json.dumps(dict(faces=q['faces'],audit=a,patches=p['patches']+q['patches']),default=lambda x:x.tolist()))
            print('FULL_BOUNDARY_PASSED',flush=True);sys.exit(0)
(out/'report.json').write_text(json.dumps(dict(status='NO_TWO_FLIP_SOLUTION',tested=count,qualified_delivery=False)))
