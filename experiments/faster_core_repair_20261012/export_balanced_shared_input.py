"""Same-boundary SI input exports; not optimized material or collision admission."""
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import gmpy2
import numpy as np
import trimesh
import _make_it_stand as native
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
DEMO=ROOT/'experiments/quadric_physics_demo_20261011'
HEART=ROOT/'experiments/faster_heart_contact_20261011'
sys.path[:0]=[str(DEMO),str(HEART)]
from core_boolean_mesh import bounded_rational_core_mesh
from shared_input_boundary import assemble_fixed_boundary
from rounded_input import exact_moments,write_mesh
prior=HERE/'noodles-balanced-full-input-job-01/output'
moments=HERE/'noodles-balanced-input-moments-job-01/output'
manifest_path=DEMO/'noodles-dense-v5/manifest.json'
manifest=json.loads(manifest_path.read_text());spec=manifest['spec']
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
record=json.loads((moments/'report.json').read_text())
if record['manifest_sha256']!=sha(manifest_path):raise ValueError('FROZEN_MANIFEST_CHANGED')
if record['capacity_status']=='EXCLUDED_BY_FIRST_MOMENT_BOUND':raise ValueError('INPUT_CAPACITY_EXCLUDED')
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
policy=sys.argv[2] if len(sys.argv)==3 else 'SURVIVING_VERTEX_V1'
def save(name,value):
    with (output/name).open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.flush();os.fsync(stream.fileno())
def event(stage,**details):
    with (output/'events.jsonl').open('a') as stream:
        stream.write(json.dumps(dict(stage=stage,seconds=time.monotonic()-started,**details))+'\n')
        stream.flush();os.fsync(stream.fileno())
try:
    packets={}
    for name in ('outer','design','fixed'):
        source=prior/(name+'-rational.json')
        if sha(source)!=record['material_sha256'][name]:raise ValueError('MOMENTS_SOURCE_CHANGED:'+name)
        packets[name]=json.loads(source.read_text())
    scale=gmpy2.mpq(str(spec['model_unit_m']))
    exact_points=[tuple(gmpy2.mpq(x)*scale for x in row) for row in packets['outer']['vertices']]
    spans=[max(p[i] for p in exact_points)-min(p[i] for p in exact_points) for i in range(3)]
    squared=sum(x*x for x in spans)
    with gmpy2.context(precision=256,round=gmpy2.RoundDown):lower=float(gmpy2.sqrt(squared))
    while gmpy2.mpq(lower)**2>squared:lower=float(np.nextafter(lower,0))
    if lower<=0:raise ValueError('INVALID_FROZEN_LENGTH_SCALE')
    meshes={};representations={}
    for name in ('outer','design'):
        event('BOUNDED_SHARED_SOURCE_BEGIN',material=name,boundary_witness_policy=policy)
        rows=[[str((gmpy2.mpq(x)*scale).numerator)+'/'+str((gmpy2.mpq(x)*scale).denominator)
               for x in row] for row in packets[name]['vertices']]
        meshes[name],representations[name]=bounded_rational_core_mesh(rows,packets[name]['faces'],
            length_scale=lower,boundary_witness_policy=policy)
        save(name+'-representation.json',representations[name])
    # Every design vertex/triangle is exactly reused with the opposite orientation.
    # Containment is independently checked below; this assembly alone grants nothing.
    meshes['fixed']=assemble_fixed_boundary(meshes['outer'],meshes['design'])
    limit=gmpy2.mpq(lower)**3/gmpy2.mpq(10**6)
    swept=sum(gmpy2.mpq(representations[name]['swept_volume_upper_bound_exact']) for name in ('outer','design'))
    if swept>limit:raise ValueError('ASSEMBLED_FIXED_VOLUME_BOUND_EXCEEDED')
    save('shared-assembly-correspondence.json',dict(
        profile='SAME_DERIVED_OUTER_AND_DESIGN_BOUNDARY_FIXED_V1',
        outer_vertices=len(meshes['outer'].vertices),design_vertices=len(meshes['design'].vertices),
        design_in_fixed_vertex_offset=len(meshes['outer'].vertices),design_faces_reversed=True,
        length_scale_lower_bound_m=lower,authority_bbox_diagonal_squared_m2_exact=str(squared),
        fixed_swept_volume_upper_bound_exact=str(swept),independent_material_admission=False))
    regions={};files={};errors={};powers=[3]+[4]*3+[5]*6
    for name in ('outer','design','fixed'):
        event('OBJ_WRITE_RELOAD_BEGIN',material=name)
        path=output/(name+'-si-candidate.obj');write_mesh(path,meshes[name])
        loaded=trimesh.load(path,force='mesh',process=False)
        if not np.array_equal(loaded.vertices,meshes[name].vertices) or not np.array_equal(loaded.faces,meshes[name].faces):
            raise ValueError('OBJ_RELOAD_CHANGED_SHARED_BOUNDARY:'+name)
        event('DERIVED_EXACT_REGION_BEGIN',material=name)
        regions[name]=native.make_exact_material_region(loaded.vertices,loaded.faces.astype(np.int32))
        save(name+'-derived-facts.json',dict(regions[name].facts))
        if not regions[name].facts['completed'] or not regions[name].facts['valid']:
            raise ValueError('DERIVED_REGION_REJECTED:'+name)
        event('DERIVED_INDEPENDENT_MOMENTS_BEGIN',material=name)
        values=exact_moments(loaded.vertices,loaded.faces)
        reference=[gmpy2.mpq(x) for x in json.loads((moments/(name+'-moments-exact.json')).read_text())]
        errors[name]=[float(abs(v-r*scale**p)/gmpy2.mpq(lower)**p) for v,r,p in zip(values,reference,powers)]
        if max(errors[name])>1e-6:raise ValueError('DERIVED_MOMENTS_CHANGED:'+name)
        files[name]=dict(name=path.name,sha256=sha(path),role='INPUT_NOT_OPTIMIZED',length_unit='m')
    event('DERIVED_INDEPENDENT_PARTITION_BEGIN')
    partition=native.check_exact_material_partition(regions['outer'],[regions['fixed'],regions['design']])
    save('derived-partition.json',dict(partition))
    if not partition['completed'] or not partition['complete_partition']:
        raise ValueError('DERIVED_SHARED_PARTITION_REJECTED')
    save('manifest.json',dict(schema='balanced-noodles-shared-input-obj-v1',files=files,
        authority_material_sha256=record['material_sha256'],input_moments_sha256=sha(moments/'report.json'),
        model_manifest_sha256=sha(manifest_path),length_unit='m',density_kg_m3=spec['density_kg_m3'],
        boundary_witness_policy=policy,
        maximum_normalized_moment_errors=errors,derived_input_partition_passed=True,
        capacity_status=record['capacity_status'],feasible_witness=False,
        solver_admission=False,collision_admission=False,optimization_invocations=0,
        printing_qualification='NOT_VERIFIED',native_sha256=sha(Path(native.__file__))))
    event('SHARED_DERIVED_INPUT_COMPLETE_OPTIMIZATION_PENDING')
except Exception as exc:
    save('failure.json',dict(type=type(exc).__name__,message=str(exc),seconds=time.monotonic()-started,
        derived_input_admission=False,optimization_admission=False,collision_admission=False))
    event('SHARED_DERIVED_INPUT_REJECTED',reason=str(exc));raise
