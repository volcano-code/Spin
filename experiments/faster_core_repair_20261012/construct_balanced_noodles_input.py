"""Frozen full input, same-set balanced exact union; no solver admission."""
from fractions import Fraction as Q
import hashlib
import json
import os
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
names=manifest['construction']['fixed_parts']
if len(names)!=81 or len(set(names))!=81:raise ValueError('FROZEN_81_PARTS_REQUIRED')
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
executable=HERE/'exact-input-union-build-job-01/output/exact-input-boolean'
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
events=output/'events.jsonl';counter=0
def event(data):
    with events.open('a') as stream:
        stream.write(json.dumps(data,sort_keys=True,allow_nan=False)+'\n')
        stream.flush();os.fsync(stream.fileno())
def load(name):
    entry=manifest['parts'][name];path=case/entry['model_obj']
    if sha(path)!=entry['sha256']['model']:raise ValueError('FROZEN_PART_CHANGED:'+name)
    mesh=trimesh.load(path,force='mesh',process=False)
    return dict(vertices=[[binary64_rational_token(v) for v in row] for row in mesh.vertices],
                faces=mesh.faces.tolist(),leaf_parts=[name])
def operation(a,b,kind):
    global counter
    counter+=1;input_path=output/f'operands-{counter:03d}.txt';result=output/f'result-{counter:03d}.json'
    with input_path.open('x') as stream:
        for packet in (a,b):
            stream.write(f"{len(packet['vertices'])} {len(packet['faces'])}\n")
            for row in packet['vertices']:
                values=[Q(v) for v in row]
                stream.write(' '.join(str(v.numerator)+'/'+str(v.denominator) for v in values)+'\n')
            for row in packet['faces']:stream.write(' '.join(map(str,row))+'\n')
    ancestry=[a['leaf_parts'],b['leaf_parts']]
    event(dict(event='boolean_begin',index=counter,operation=kind,leaf_parts=ancestry,
        operand_sha256=sha(input_path),input_face_counts=[len(a['faces']),len(b['faces'])]))
    started=time.monotonic()
    with (output/f'stderr-{counter:03d}.log').open('x') as errors:
        subprocess.run([str(executable),str(input_path),str(result),kind],stderr=errors,check=True)
    packet=json.loads(result.read_text());packet['leaf_parts']=a['leaf_parts']+b['leaf_parts']
    event(dict(event='boolean_end',index=counter,operation=kind,seconds=time.monotonic()-started,
        result_sha256=sha(result),output_faces=len(packet['faces']),leaf_parts=packet['leaf_parts']))
    return packet
pilot=HERE/'noodles-balanced-eight-part-pilot-job-01/output/report.json'
audit=HERE/'noodles-balanced-eight-part-audit-job-01/output/report.json'
if not json.loads(pilot.read_text())['empty_bidirectional_set_difference'] or not json.loads(audit.read_text())['boundary_valid']:
    raise ValueError('SMALL_RECIPE_GATE_NOT_PASSED')
event(dict(event='construction_begin',profile='CGAL_BALANCED_EXACT_FROZEN_INPUT_V1',
    manifest_sha256=sha(manifest_path),seed=seed,parts=names,
    executable_sha256=sha(executable),pilot_sha256=sha(pilot),audit_sha256=sha(audit),
    geometry_admission=False,optimization_admission=False))
started=time.monotonic()
parts=[load(name) for name in names];inner=load(seed);outer=load('outer-ellipsoid')
fixed_union=balanced_union(parts,lambda a,b:operation(a,b,'union'))
if fixed_union['leaf_parts']!=names:raise ValueError('PART_COVERAGE_CHANGED')
design=operation(inner,fixed_union,'difference')
assembled=operation(outer,fixed_union,'union')
fixed=operation(assembled,design,'difference')
for name,packet in [('fixed-parts-union',fixed_union),('design',design),('outer',assembled),('fixed',fixed)]:
    with (output/(name+'-rational.json')).open('x') as stream:
        json.dump(packet,stream,allow_nan=False);stream.flush();os.fsync(stream.fileno())
report=dict(profile='CGAL_BALANCED_EXACT_FROZEN_INPUT_V1',seed=seed,parts=names,
    actual_boolean_count=counter,construction_seconds=time.monotonic()-started,
    complete_81_part_construction=True,independent_boundary_audit=False,
    geometry_admission=False,optimization_admission=False,
    output_sha256={name:sha(output/(name+'-rational.json')) for name in ['design','outer','fixed','fixed-parts-union']})
with (output/'report.json').open('x') as stream:json.dump(report,stream,indent=2)
event(dict(event='construction_complete',**report));print(json.dumps(report),flush=True)
