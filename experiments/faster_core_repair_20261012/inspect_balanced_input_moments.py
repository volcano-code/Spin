"""Independent input integration/capacity after all three boundaries pass."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import gmpy2
import trimesh
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
HEART=ROOT/'experiments/faster_heart_contact_20261011'
DEMO=ROOT/'experiments/quadric_physics_demo_20261011'
sys.path[:0]=[str(HEART),str(DEMO),str(ROOT/'experiments/heart_stem_contact')]
from rounded_input import exact_moments
from saved_input_capacity import capacity_from_fixed_moments
from core_boolean_mesh import binary64_rational_token,manifest_design_seed
from noodles_admission import boundary_components
from contact_matrix import physical_properties
case=DEMO/'noodles-dense-v5';manifest_path=case/'manifest.json'
manifest=json.loads(manifest_path.read_text());spec=manifest['spec']
if spec['target_axis']!=[0.,0.,1.]:raise ValueError('UNSUPPORTED_CAPACITY_TARGET_AXIS')
if spec['density_kg_m3']<=0 or spec['model_unit_m']<=0:raise ValueError('INVALID_DENSITY_OR_UNIT')
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
prior=HERE/'noodles-balanced-full-input-job-01/output';packets={};identities={}
started=time.monotonic()
def save(name,value):
    with (output/name).open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.flush();os.fsync(stream.fileno())
def event(stage,**details):
    with (output/'events.jsonl').open('a') as stream:
        stream.write(json.dumps(dict(stage=stage,elapsed_seconds=time.monotonic()-started,**details))+'\n')
        stream.flush();os.fsync(stream.fileno())
for name in ('design','outer','fixed'):
    source=prior/(name+'-rational.json');audit=HERE/f'noodles-balanced-{name}-audit-job-01/output'
    if not json.loads((audit/'report.json').read_text())['boundary_valid']:
        raise ValueError('BOUNDARY_NOT_VALID:'+name)
    if sha(source)!=json.loads((audit/'receipt.json').read_text())['input_sha256']:
        raise ValueError('BOUNDARY_AUDIT_INPUT_CHANGED:'+name)
    packets[name]=json.loads(source.read_text());identities[name]=sha(source)
if boundary_components(packets['outer']['faces'])!=1:raise ValueError('OUTER_NOT_SINGLE_CONNECTED')
event('THREE_SAVED_BOUNDARIES_BOUND',material_sha256=identities)
# Independent necessary subset assumption: original inner lies inside original outer.
seed=manifest_design_seed(manifest);operands=[]
for name in (seed,'outer-ellipsoid'):
    entry=manifest['parts'][name];path=case/entry['model_obj']
    if sha(path)!=entry['sha256']['model']:raise ValueError('FROZEN_SOURCE_CHANGED')
    mesh=trimesh.load(path,force='mesh',process=False)
    operands.append(([[binary64_rational_token(v) for v in row] for row in mesh.vertices],mesh.faces))
text=output/'inner-outer-operands.txt'
with text.open('x') as stream:
    for vertices,faces in operands:
        stream.write(f'{len(vertices)} {len(faces)}\n')
        for row in vertices:stream.write(' '.join(row)+'\n')
        for row in faces:stream.write(' '.join(map(str,row))+'\n')
exe=HERE/'exact-input-union-build-job-01/output/exact-input-boolean'
with (output/'membership-stderr.log').open('x') as errors:
    subprocess.run([str(exe),str(text),str(output/'inner-minus-outer.json'),'difference'],stderr=errors,check=True)
if json.loads((output/'inner-minus-outer.json').read_text())['faces']:
    raise ValueError('ORIGINAL_INNER_NOT_CONTAINED_IN_OUTER')
event('ORIGINAL_INNER_SUBSET_OUTER_KERNEL_CHECK_PASSED')
values={};physics={}
for name,packet in packets.items():
    event('EXACT_MOMENTS_BEGIN',material=name)
    values[name]=exact_moments(packet['vertices'],packet['faces'])
    if values[name][0]<=0:raise ValueError('NONPOSITIVE_INPUT_MASS:'+name)
    save(name+'-moments-exact.json',list(map(str,values[name])))
    p=physical_properties([float(x) for x in values[name]],
        mass_unit_kg=spec['density_kg_m3']*spec['model_unit_m']**3,length_unit_m=spec['model_unit_m'])
    physics[name]=dict(mass_kg=float(p[0]),center_m=p[1].tolist(),inertia_center_kg_m2=p[2].tolist(),
        principal_moments_kg_m2=p[3].tolist(),principal_axes=p[4].tolist())
    event('EXACT_MOMENTS_COMPLETE',material=name,mass_kg=physics[name]['mass_kg'])
delta=[a-b-c for a,b,c in zip(values['outer'],values['fixed'],values['design'])]
if any(x!=0 for x in delta):raise ValueError('EXACT_PARTITION_MOMENTS_MISMATCH')
event('EXACT_TEN_MOMENT_PARTITION_PASSED')
event('LATERAL_CAPACITY_BEGIN')
capacity=capacity_from_fixed_moments(values['fixed'],packets['design']['vertices'],
    packets['design']['faces'],spec['target_origin_model'])
save('capacity.json',capacity)
save('physical-properties.json',dict(materials=physics,material_sha256=identities,
    length_unit_m=spec['model_unit_m'],density_kg_m3=spec['density_kg_m3'],
    scope='INPUT_NOT_OPTIMIZED',collision_admission=False))
save('report.json',dict(schema='balanced-noodles-independent-input-moments-v1',
    material_sha256=identities,manifest_sha256=sha(manifest_path),
    exact_partition_moment_delta=list(map(str,delta)),outer_connected=True,
    original_inner_subset_outer_kernel_check=True,capacity_status=capacity['status'],
    feasible_witness=False,independent_partition_region_certificate=False,
    derived_output_admission=False,optimization_admission=False,
    gmpy2_version=gmpy2.version(),gmp_version=gmpy2.mp_version(),
    seconds=time.monotonic()-started))
event('INPUT_MOMENTS_CAPACITY_COMPLETE',capacity_status=capacity['status'],optimization_admission=False)
