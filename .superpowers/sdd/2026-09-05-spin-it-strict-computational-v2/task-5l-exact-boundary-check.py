"""Exact Fraction reference from represented coordinates; no project imports."""
import itertools
import json
import math
import sys
from fractions import Fraction as F
from pathlib import Path

actual=json.loads(Path(sys.argv[1]).read_text())
assert actual["scope"]=="SYNTHETIC_BOUNDARY_ONLY"
keys=[tuple(k) for k in actual["keys"]]
assert keys==list(itertools.product(range(2,7),repeat=3))
index={key:i for i,key in enumerate(keys)}
nodes=[[F(x) for x in row] for row in actual["nodes"]]
handles=actual["handles"]
expected_handles={index[k] for k in keys if all(x%2==0 for x in k) and any(x in (2,6) for x in k)}
assert len(handles)==len(set(handles))==26 and set(handles)==expected_handles
handle_index={keys[n]:h for h,n in enumerate(handles)}
rows={r["node"]:r["entries"] for r in actual["boundary"]}
assert len(rows)==len(actual["boundary"])==98
assert set(rows)=={i for i,k in enumerate(keys) if any(x in (2,6) for x in k)}
maximum=0.; exact_entries=0; shared_edge_count=0
for n,entries in rows.items():
    key=keys[n]
    axis=next(a for a in range(3) if key[a] in (2,6))
    axes=[a for a in range(3) if a!=axis]
    low=list(key)
    for a in axes: low[a]=min((key[a]//2)*2,4)
    ratios=[]
    for a in axes:
        high=low.copy();high[a]+=2
        ratios.append((nodes[n][a]-nodes[index[tuple(low)]][a])/(nodes[index[tuple(high)]][a]-nodes[index[tuple(low)]][a]))
    u,v=ratios
    assert 0<=u<=1 and 0<=v<=1
    expected=[F(0)]*26
    for uv,value in {(0,0):1-max(u,v),(2,2):min(u,v),(2,0):max(u-v,0),(0,2):max(v-u,0)}.items():
        corner=low.copy()
        for a,offset in zip(axes,uv):corner[a]+=offset
        expected[handle_index[tuple(corner)]]+=value
    got=[0.]*26
    assert len(entries)==len({h for h,w in entries})
    for h,w in entries:
        assert 0<=h<26 and math.isfinite(w) and w>0
        got[h]=w
    maximum=max(maximum,max(abs(g-float(e)) for g,e in zip(got,expected)))
    exact_entries+=26
    assert abs(sum(got)-1)<1e-14
    if sum(k in (2,6) for k in key)==2 and any(k%2 for k in key):
        assert len(entries)==2, "shared edge must use only its actual endpoints"
        shared_edge_count+=1
assert maximum<=1e-15,maximum
assert shared_edge_count==24
if len(sys.argv)>2:
    assert Path(sys.argv[1]).read_bytes()==Path(sys.argv[2]).read_bytes(),"repeat differs"
print(json.dumps(dict(scope="EXACT_REFERENCE_VS_NUMERICAL_SYNTHETIC_BOUNDARY",rule=actual["rule"],
                     rows=len(rows),entries_compared=exact_entries,shared_edge_nodes=shared_edge_count,
                     maximum_absolute_error=maximum,repeat_checked=len(sys.argv)>2),sort_keys=True))
