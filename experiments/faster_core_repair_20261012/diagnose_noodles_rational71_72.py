"""Registered two-step test: remove only the intermediate step71 float feedback."""
import hashlib
import json
from pathlib import Path
import sys
import time
import gmpy2
import numpy as np
import trimesh
import _make_it_stand as native
from noodles_vertex_link_witness import link_components

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
demo=ROOT/'experiments/quadric_physics_demo_20261011';case=demo/'noodles-dense-v5'
prefix_path=demo/'noodles-dense-v5-core-preprocessing-v1/domain-job-01/output/domain-070.json'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest=json.loads((case/'manifest.json').read_text());prefix=json.loads(prefix_path.read_text())
if prefix['excluded_parts']!=manifest['construction']['fixed_parts'][:70] or prefix['manifest_sha256']!=sha(case/'manifest.json'):
    raise ValueError('FROZEN_70_PREFIX_MISMATCH')
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
def save(name,data):
    with (output/name).open('x') as stream:json.dump(data,stream,default=lambda x:x.tolist(),allow_nan=False)
target_path=HERE/'noodles-rational-authority-diagnosis-job-01/output/vertex-link-witness.json'
target=np.array([float(gmpy2.mpq(x)) for x in json.loads(target_path.read_text())['witnesses'][0]['exact_coordinate']])
packet=None;started=time.monotonic();rows=[]
for step in (71,72):
    name=manifest['construction']['fixed_parts'][step-1];entry=manifest['parts'][name]
    path=case/entry['model_obj']
    if sha(path)!=entry['sha256']['model']:raise ValueError('FIXED_PART_CHANGED')
    part=trimesh.load(path,force='mesh',process=False)
    if step==71:
        packet=native.libigl_exact_boolean_construction(np.asarray(prefix['vertices']),
            np.asarray(prefix['faces'],dtype=np.int32),part.vertices,part.faces,'difference')
    else:
        cutter=[[str(gmpy2.mpq(float(x))) for x in row] for row in part.vertices]
        packet=native.libigl_exact_boolean_construction_from_rationals(packet['vertices'],
            np.asarray(packet['faces'],dtype=np.int32),cutter,part.faces,'difference')
    save(f'original-{step:03d}-rational.json',dict(packet))
    if not packet['succeeded']:raise ValueError('ACTUAL_BOOLEAN_FAILED')
    vertices=np.array([[float(gmpy2.mpq(x)) for x in row] for row in packet['vertices']])
    matches=np.flatnonzero(np.all(vertices==target,axis=1)).tolist()
    stars=[dict(vertex=i,exact_coordinate=packet['vertices'][i],**link_components(packet['faces'],i)) for i in matches]
    rows.append(dict(step=step,part=name,part_sha256=sha(path),seconds=time.monotonic()-started,stars=stars))
save('report.json',dict(schema='rational71-72-feedback-diagnosis-v1',rows=rows,
    prefix_sha256=sha(prefix_path),native_sha256=sha(Path(native.__file__)),
    hypothesis='The intermediate step71 binary64 feedback creates the step72 bad vertex link.',
    geometry_admission=False,optimization_admission=False,
    limitation='Only local target star checked; not a full geometry certificate.'))
print(json.dumps([dict(step=row['step'],stars=[dict(vertex=s['vertex'],single_cycle=s['is_single_cycle'],
    components=len(s['components'])) for s in row['stars']]) for row in rows]),flush=True)
