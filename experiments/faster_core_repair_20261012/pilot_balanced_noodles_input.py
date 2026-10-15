"""Eight frozen parts: exact balanced/sequential same-set pilot, NOT full input."""
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
import trimesh
from balanced_exact_input import balanced_union
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
DEMO=ROOT/'experiments/quadric_physics_demo_20261011'
sys.path.insert(0,str(DEMO))
from core_boolean_mesh import binary64_rational_token,manifest_design_seed
case=DEMO/'noodles-dense-v5';manifest_path=case/'manifest.json'
manifest=json.loads(manifest_path.read_text());seed=manifest_design_seed(manifest)
names=manifest['construction']['fixed_parts'][:8]
if len(names)!=8 or len(set(names))!=8:raise ValueError('INVALID_FIXED_PILOT_PARTS')
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
executable=HERE/'exact-input-union-build-job-01/output/exact-input-boolean'
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
events=output/'events.jsonl';counter=0
def event(data):
    with events.open('a') as stream:
        stream.write(json.dumps(data,sort_keys=True,allow_nan=False)+'\n');stream.flush()
def load(name):
    entry=manifest['parts'][name];path=case/entry['model_obj']
    if sha(path)!=entry['sha256']['model']:raise ValueError('FROZEN_PART_CHANGED:'+name)
    mesh=trimesh.load(path,force='mesh',process=False)
    return dict(vertices=[[binary64_rational_token(v) for v in row] for row in mesh.vertices],
                faces=mesh.faces.tolist())
def operation(a,b,kind):
    global counter
    counter+=1;input_path=output/f'operands-{counter:03d}.txt';result=output/f'result-{counter:03d}.json'
    with input_path.open('x') as stream:
        for packet in (a,b):
            stream.write(f"{len(packet['vertices'])} {len(packet['faces'])}\n")
            for row in packet['vertices']:
                # Canonical native grammar only; integer/rational value unchanged.
                values=[Q(v) for v in row]
                stream.write(' '.join(str(v.numerator)+'/'+str(v.denominator) for v in values)+'\n')
            for row in packet['faces']:stream.write(' '.join(map(str,row))+'\n')
    event(dict(event='boolean_begin',index=counter,operation=kind,
        operand_sha256=sha(input_path),input_face_counts=[len(a['faces']),len(b['faces'])]))
    started=time.monotonic()
    with (output/f'stderr-{counter:03d}.log').open('x') as errors:
        subprocess.run([str(executable),str(input_path),str(result),kind],stderr=errors,check=True)
    elapsed=time.monotonic()-started;packet=json.loads(result.read_text())
    event(dict(event='boolean_end',index=counter,operation=kind,seconds=elapsed,
        result_sha256=sha(result),output_faces=len(packet['faces'])))
    return packet
event(dict(event='pilot_begin',role='EIGHT_PART_INPUT_CONSTRUCTION_DIAGNOSIS_ONLY',
    manifest_sha256=sha(manifest_path),seed=seed,parts=names,
    executable_sha256=sha(executable),complete_model=False))
parts=[load(name) for name in names];inner=load(seed);outer=load('outer-ellipsoid')
started=time.monotonic()
fixed_union=balanced_union(parts,lambda a,b:operation(a,b,'union'))
design=operation(inner,fixed_union,'difference')
balanced_seconds=time.monotonic()-started
started=time.monotonic();sequential=inner
for part in parts:sequential=operation(sequential,part,'difference')
sequential_seconds=time.monotonic()-started
left=operation(design,sequential,'difference');right=operation(sequential,design,'difference')
if left['faces'] or right['faces']:raise ValueError('BALANCED_AND_SEQUENTIAL_SETS_DIFFER')
# The actual fixed shell is assembled outer minus design, not fixed_union.
assembled=operation(outer,fixed_union,'union')
fixed=operation(assembled,design,'difference')
intersection=operation(fixed,design,'intersection')
if intersection['faces']:raise ValueError('FIXED_DESIGN_POSITIVE_VOLUME_OVERLAP')
for name,packet in [('pilot-design',design),('pilot-outer',assembled),('pilot-fixed',fixed)]:
    with (output/(name+'.json')).open('x') as stream:json.dump(packet,stream,allow_nan=False)
report=dict(role='EIGHT_PART_INPUT_CONSTRUCTION_DIAGNOSIS_ONLY',parts=names,seed=seed,
    manifest_sha256=sha(manifest_path),balanced_seconds=balanced_seconds,
    sequential_seconds=sequential_seconds,empty_bidirectional_set_difference=True,
    exact_kernel_partition_overlap_empty=True,complete_model=False,
    independent_boundary_audit=False,geometry_admission=False,optimization_admission=False)
with (output/'report.json').open('x') as stream:json.dump(report,stream,indent=2)
event(dict(event='pilot_complete',**report));print(json.dumps(report),flush=True)
