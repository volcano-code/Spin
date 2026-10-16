"""Compare the same terminal boundary before and after coordinate conversion."""
import json,time
from pathlib import Path
import numpy as np
import _make_it_stand as native
HERE=Path(__file__).resolve().parent
source=HERE/'noodles-terminal-recovery-job-03/output/material'
r=json.loads((source/'mesh_0004_rational_vertices.json').read_text())
f=np.load(source/'mesh_0004_faces.npy')
start=time.monotonic()
result=dict(native.audit_exact_rational_triangle_boundary(r,f,'AnyNonzero'))
out=HERE/'noodles-terminal-exact-boundary-comparison-01.json'
with out.open('x') as stream:json.dump(dict(seconds=time.monotonic()-start,result=result),stream,default=lambda x:x.tolist())
print({k:result.get(k) for k in ('boundary_valid','self_intersection_completed','self_intersection_free','intersecting_pair_count')},flush=True)
