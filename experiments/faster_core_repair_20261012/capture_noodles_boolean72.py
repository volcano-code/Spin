"""One actual preserved-prefix Boolean, with original rational local-link evidence."""
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

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]
demo=ROOT/'experiments/quadric_physics_demo_20261011'
case=demo/'noodles-dense-v5'
predecessor=demo/'noodles-dense-v5-core-preprocessing-v1/domain-job-01/output/domain-071.json'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest=json.loads((case/'manifest.json').read_text())
prefix=json.loads(predecessor.read_text())
if len(prefix['excluded_parts'])!=71 or prefix['excluded_parts']!=manifest['construction']['fixed_parts'][:71] or prefix['manifest_sha256']!=sha(case/'manifest.json'):
    raise ValueError('FROZEN_BOOLEAN72_CONTEXT_MISMATCH')
name=manifest['construction']['fixed_parts'][71]
entry=manifest['parts'][name]; part_path=case/entry['model_obj']
if sha(part_path)!=entry['sha256']['model']:raise ValueError('FIXED_PART_CHANGED')
witness=HERE/'noodles-rational-authority-diagnosis-job-01/output/vertex-link-witness.json'
exact_target=json.loads(witness.read_text())['witnesses'][0]['exact_coordinate']
target=np.array([float(gmpy2.mpq(x)) for x in exact_target])
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
part=trimesh.load(part_path,force='mesh',process=False)
start=time.monotonic()
packet=native.libigl_exact_boolean_construction(np.asarray(prefix['vertices']),
    np.asarray(prefix['faces'],dtype=np.int32),part.vertices,part.faces,'difference')
def save(name,data):
    with (output/name).open('x') as stream:
        json.dump(data,stream,default=lambda x:x.tolist(),allow_nan=False)
save('original-rational-output.json',dict(packet))
if not packet['succeeded']:raise ValueError('ACTUAL_BOOLEAN72_FAILED')
vertices=np.array([[float(gmpy2.mpq(x)) for x in row] for row in packet['vertices']])
matches=np.flatnonzero(np.all(vertices==target,axis=1)).tolist()
stars=[dict(vertex=i,exact_coordinate=packet['vertices'][i],
    same_exact_coordinate=all(gmpy2.mpq(a)==gmpy2.mpq(b) for a,b in zip(packet['vertices'][i],exact_target)),
    **link_components(packet['faces'],i)) for i in matches]
save('local-link-report.json',dict(schema='original-rational-boolean72-link-v1',
    predecessor_sha256=sha(predecessor),part=name,part_sha256=sha(part_path),
    native_sha256=sha(Path(native.__file__)),witness_sha256=sha(witness),seconds=time.monotonic()-start,
    stars=stars,geometry_admission=False,optimization_admission=False,
    limitation='Local link evidence only; no complete boundary or region certification.'))
print(json.dumps(dict(part=name,stars=[dict(vertex=s['vertex'],single_cycle=s['is_single_cycle'],
    components=len(s['components']),same_exact_coordinate=s['same_exact_coordinate']) for s in stars])),flush=True)
