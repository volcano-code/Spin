"""Compare both index patterns at EACH fixed geometry, never different entities."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from cherchi_core_map import construction_constraint_key
from exact_planar_patches import patch_boundary_signature,coalesce_collinear_cycle,oriented_planar_chain_signature
HERE=Path(__file__).resolve().parent;capture=HERE/'saved-small-provenance-capture-01'
inputs=HERE/'small-top-derivative-probe-job-02/output'
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
data=[]
for index in (0,2):
    p=inputs/f'point-{index:03d}-material.json';e=capture/f'point-{index:03d}.json'
    saved=json.loads(p.read_text());evidence=json.loads(e.read_text())
    if evidence['saved_geometry_sha256']!=hashlib.sha256(p.read_bytes()).hexdigest():raise ValueError('SAVED_SOURCE_MISMATCH')
    keys=[None]*len(saved['vertices'])
    for record in evidence['resolved_construction_records']:
        keys[evidence['new_vertex_to_saved_vertex'][record['vertex_index']]]=construction_constraint_key(record)
    if None in keys or len(set(keys))!=len(keys):raise ValueError('NONBIJECTIVE_SOURCE_FEATURES')
    data.append((index,saved,keys))
if set(data[0][2])!=set(data[1][2]):raise ValueError('CONSTRUCTION_FEATURE_SET_CHANGED')
canonical={key:i for i,key in enumerate(sorted(data[0][2]))}
patterns=[[[canonical[keys[v]] for v in face] for face in saved['faces']] for _,saved,keys in data]
rows=[]
for index,saved,keys in data:
    points=[None]*len(keys)
    for i,key in enumerate(keys):points[canonical[key]]=saved['vertices'][i]
    signatures=[];chain_signatures=[];null_counts=[]
    row=dict(index=index,geometry_unchanged=True)
    for faces in patterns:
        chain=oriented_planar_chain_signature(points,faces)
        chain_signatures.append(chain['signature'])
        null_counts.append(len(chain['exact_null_triangle_indices']))
    row.update(exact_oriented_surface_chain_equal=chain_signatures[0]==chain_signatures[1],
               exact_null_triangle_counts_in_proof_view=null_counts)
    try:
        for faces in patterns:
            cert=patch_boundary_signature(points,faces)
            signature=[]
            for patch in cert['patches']:
                loops=[]
                for cycle in patch:
                    loop=tuple(coalesce_collinear_cycle(points,cycle))
                    loops.append(min(loop[i:]+loop[:i] for i in range(len(loop))))
                signature.append(tuple(sorted(loops)))
            signatures.append(sorted(signature))
        row.update(exact_coplanar_fixed_geometry_boundary_equal=signatures[0]==signatures[1],
                   patch_counts=[len(s) for s in signatures])
    except ValueError as exc:
        row.update(exact_coplanar_fixed_geometry_boundary_equal=False,exception=str(exc))
    rows.append(row)
report=dict(rows=rows,source_feature_correspondence_bijective=True,new_boolean_calls=0,
    finite_pair_retriangulation_certified=all(r['exact_oriented_surface_chain_equal'] for r in rows),
    full81_branch_certification=False,optimization_admission=False,
    scope='Point000/002 only; each comparison uses one identical rational vertex geometry.')
with (output/'report.json').open('x') as stream:json.dump(report,stream,indent=2)
print(json.dumps(report),flush=True)
