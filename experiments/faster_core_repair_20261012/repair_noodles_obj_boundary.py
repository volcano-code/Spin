"""Try the registered bounded alias conversion on the saved exact terminal."""
import json,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'experiments/quadric_physics_demo_20261011'))
from core_boolean_mesh import bounded_rational_core_mesh
source=HERE/'noodles-terminal-recovery-job-03/output/material'
out=HERE/'noodles-terminal-bounded-obj-02';out.mkdir(exist_ok=False)
v=np.load(source/'mesh_0000_vertices.npy')
length=float(np.nextafter(np.linalg.norm(np.ptp(v,axis=0)),0))
r=json.loads((source/'mesh_0004_rational_vertices.json').read_text())
f=np.load(source/'mesh_0004_faces.npy')
mesh,proof=bounded_rational_core_mesh(r,f,length_scale=length,boundary_witness_policy='SURVIVING_EDGE_V2')
(out/'representation.json').write_text(json.dumps(proof))
mesh.export(out/'complete-candidate.obj')
import _make_it_stand as native
audit=dict(native.audit_exact_triangle_boundary(mesh.vertices,np.asarray(mesh.faces,dtype=np.int32),'AnyNonzero'))
(out/'audit.json').write_text(json.dumps(audit,default=lambda x:x.tolist()))
print(json.dumps(dict(boundary_valid=audit['boundary_valid'],intersecting_pair_count=audit.get('intersecting_pair_count'),removed_pairs=proof['removed_opposite_face_pairs'],qualified_delivery=False)))
