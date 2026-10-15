"""Derive and independently check a mesh only after the rational boundary passes."""
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
import trimesh
import _make_it_stand as native
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
demo=ROOT/'experiments/quadric_physics_demo_20261011'
sys.path.insert(0,str(demo))
from core_boolean_mesh import bounded_rational_core_mesh
from independent_tetra_moments import integrate
audit_path=HERE/'noodles-repaired-rational-audit-job-01/output/report.json'
if not json.loads(audit_path.read_text())['boundary_valid']:raise ValueError('RATIONAL_BOUNDARY_NOT_ADMITTED')
source=HERE/'noodles-repaired-rational-suffix-job-01/output/domain-081-rational.json'
packet=json.loads(source.read_text());case=demo/'noodles-dense-v5'
manifest_path=case/'manifest.json';manifest=json.loads(manifest_path.read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
parts=[]
for name in ['outer-ellipsoid']+manifest['construction']['fixed_parts']:
    entry=manifest['parts'][name];path=case/entry['model_obj']
    if sha(path)!=entry['sha256']['model']:raise ValueError('FROZEN_PART_CHANGED')
    parts.append(trimesh.load(path,force='mesh',process=False))
lo=np.min([m.bounds[0] for m in parts],axis=0);hi=np.max([m.bounds[1] for m in parts],axis=0)
length=float(np.nextafter(np.linalg.norm(hi-lo),0))
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
mesh,representation=bounded_rational_core_mesh(packet['vertices'],packet['faces'],length_scale=length)
def save(name,value):
    with (output/name).open('x') as stream:json.dump(value,stream,default=lambda x:x.tolist(),allow_nan=False)
save('representation.json',representation)
save('domain-081.json',dict(schema='noodles-core-binary64-domain-v1',profile='CGAL_CORE_RATIONAL_FROM71_V1',
    vertices=mesh.vertices.tolist(),faces=mesh.faces.tolist(),
    excluded_parts=manifest['construction']['fixed_parts'],manifest_sha256=sha(manifest_path),
    authority_sha256=sha(source),geometry_admission=False,optimization_admission=False))
mesh.export(output/'design-candidate.obj')
boundary=native.audit_exact_triangle_boundary(mesh.vertices,mesh.faces,'AnyNonzero')
save('boundary-audit.json',dict(boundary))
report=dict(authority_sha256=sha(source),audit_sha256=sha(audit_path),
    derived_boundary_valid=boundary['boundary_valid'],seconds=time.monotonic()-started,
    geometry_admission=False,optimization_admission=False,
    pending=['actual outer union','actual fixed material','whole input partition','capacity check'])
if boundary['boundary_valid']:
    independent=integrate(mesh.vertices,mesh.faces,1.)
    reference=np.asarray(native.integrate_raw_mass_moments(mesh.vertices,mesh.faces,1.))
    report.update(design_moments=independent.tolist(),moment_max_error=float(np.max(np.abs(independent-reference))))
else:report['blocker']='DERIVED_BOUNDARY_REJECTED'
save('report.json',report);print(json.dumps(report),flush=True)
