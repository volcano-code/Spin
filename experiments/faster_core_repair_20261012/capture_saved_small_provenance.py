"""Only two necessary provenance replays, binding each to a saved exact boundary."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from make_it_stand import spin
from cherchi_core_map import CherchiCoreMap
from run_small_inprocess_solver import safe
HERE=Path(__file__).resolve().parent
run=HERE/'small-top-inprocess-solver-sequence-01/seed-1-job-01/output'
inputs=HERE/'small-top-derivative-probe-job-02/output'
config=json.loads((run/'configuration.json').read_text())
payload=dict(config['request']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
mapping=CherchiCoreMap(spin.SpinItFasterRequest(**payload),spin.SpinEvaluatorConfig(**config['evaluator']))
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False);rows=[]
def triangles(points,faces):
    result=[]
    for face in faces:
        tri=tuple(points[int(i)] for i in face)
        result.append(min(tri,tri[1:]+tri[:1],tri[2:]+tri[:2]))
    return sorted(result)
for index in (0,2):
    path=inputs/f'point-{index:03d}-material.json';saved=json.loads(path.read_text())
    evaluation=mapping(saved['actual_coefficients']);last=mapping.last
    old=tuple(tuple(Fraction(x) for x in row) for row in saved['vertices'])
    if triangles(old,saved['faces'])!=triangles(last['exact_vertices'],last['faces']):
        raise ValueError('REPLAY_NOT_SAME_EXACT_ORIENTED_TRIANGLE_BOUNDARY')
    if len(set(old))!=len(old):raise ValueError('SAVED_COORDINATE_CORRESPONDENCE_AMBIGUOUS')
    old_indices={p:i for i,p in enumerate(old)}
    index_map=[old_indices[p] for p in last['exact_vertices']]
    record=dict(index=index,saved_geometry_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        identical_exact_oriented_triangle_boundary=True,new_vertex_to_saved_vertex=index_map,
        original_construction_records=last['original_construction_records'],
        resolved_construction_records=last['resolved_construction_records'],
        source=last['source'],operation=last['operation'],evaluation=evaluation,
        full_branch_certification=False,optimization_admission=False)
    with (output/f'point-{index:03d}.json').open('x') as stream:json.dump(safe(record),stream,allow_nan=False)
    rows.append(dict(index=index,saved_vertex_count=len(old),records=len(record['original_construction_records'])))
with (output/'report.json').open('x') as stream:
    json.dump(dict(points=rows,material_boolean_calls=2,native_sha256=hashlib.sha256(Path(mapping.native.__file__).read_bytes()).hexdigest(),
        adapter_sha256=hashlib.sha256((HERE/'cherchi_core_map.py').read_bytes()).hexdigest(),
        bridge_sha256=hashlib.sha256(Path(mapping.bridge.__file__).read_bytes()).hexdigest(),
        purpose='Necessary missing source-feature evidence only; original81 numerical batch not rerun.',
        full_branch_certification=False,optimization_admission=False),stream,indent=2)
print(json.dumps(rows),flush=True)
