"""Independent OBJ reload, exact coordinate bounds, and exact tetrahedral moments."""
import json,hashlib,sys
from pathlib import Path
import numpy as np
import gmpy2 as g
import _make_it_stand as native
import _cherchi_rational_moments as arithmetic
HERE=Path(__file__).resolve().parent
source=HERE/'noodles-terminal-recovery-job-03/output/material'
out=HERE/'noodles-terminal-interior-patch-01'
p=json.loads((out/'candidate.json').read_text());v=np.array(p['vertices']);f=np.array(p['faces'],dtype=np.int32)
path=out/'complete.obj'
with path.open('x') as stream:
    stream.write('# metres; optimized terminal; derived whole boundary; partition acceptance pending\n')
    for row in v:stream.write('v '+' '.join(format(x,'.17g') for x in row)+'\n')
    for row in f:stream.write('f '+' '.join(str(int(x)+1) for x in row)+'\n')
rv=[];rf=[]
for line in path.read_text().splitlines():
    if line.startswith('v '):rv.append([float(x) for x in line.split()[1:]])
    elif line.startswith('f '):rf.append([int(x)-1 for x in line.split()[1:]])
assert np.array_equal(v,rv) and np.array_equal(f,rf)
audit=dict(native.audit_exact_triangle_boundary(np.array(rv),np.array(rf,dtype=np.int32),'AnyNonzero'))
exact=[[g.mpq(x) for x in row] for row in p['rational_vertices']]
display=[[g.mpq(float(x)) for x in row] for row in rv]
delta=max(sum(abs(a-b) for a,b in zip(x,y)) for x,y in zip(exact,display))
span=sum(max(row[j] for row in exact)-min(row[j] for row in exact) for j in range(3))
swept=len(f)*(delta*span**2+12*delta**2*span+6*delta**3)
original=np.load(source/'mesh_0000_vertices.npy');L=float(np.linalg.norm(np.ptp(original,axis=0)))
assert delta<=g.mpq(L)/10**6 and swept<=g.mpq(L)**3/10**6
def tokens(rows):return [[str(x.numerator)+'/'+str(x.denominator) for x in row] for row in rows]
authority=json.loads((source/'mesh_0004_rational_vertices.json').read_text())
af=np.load(source/'mesh_0004_faces.npy')
am=np.array(arithmetic.integrate_moments(authority,af,1.0))
em=np.array(arithmetic.integrate_moments(tokens(exact),f,1.0))
dm=np.array(arithmetic.integrate_moments(tokens(display),f,1.0))
scales=np.array([L**3]+[L**4]*3+[L**5]*6)
error=float(np.max(abs(dm-am)/scales));exact_error=float(np.max(abs(em-am)/scales))
assert error<=1e-6 and exact_error<=1e-10
report=dict(boundary_valid=audit['boundary_valid'],intersecting_pair_count=audit['intersecting_pair_count'],obj_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),roundtrip_identical=True,length_scale=L,boundary_distance_l1_bound_exact=str(delta),symmetric_difference_volume_bound_exact=str(swept),moment_scaled_error=error,exact_retriangulation_moment_scaled_error=exact_error,authority_moments=am.tolist(),derived_moments=dm.tolist(),derived_partition='PENDING',qualified_material_delivery=False)
(out/'independent-output-check.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
