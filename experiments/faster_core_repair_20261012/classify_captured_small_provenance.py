"""Exact finite-sample retriangulation evidence from same-operation source recipes."""
import hashlib
import json
from pathlib import Path
import sys
from exact_planar_patches import patch_boundary_signature,coalesce_collinear_cycle
from cherchi_core_map import construction_constraint_key
HERE=Path(__file__).resolve().parent
captured=HERE/'saved-small-provenance-capture-01'
inputs=HERE/'small-top-derivative-probe-job-02/output'
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False);rows=[]
for index in (0,2):
    geometry_path=inputs/f'point-{index:03d}-material.json'
    evidence_path=captured/f'point-{index:03d}.json'
    saved=json.loads(geometry_path.read_text());evidence=json.loads(evidence_path.read_text())
    if evidence['saved_geometry_sha256']!=hashlib.sha256(geometry_path.read_bytes()).hexdigest() or not evidence['identical_exact_oriented_triangle_boundary']:
        raise ValueError('SOURCE_EVIDENCE_NOT_BOUND_TO_SAVED_GEOMETRY')
    labels=[None]*len(saved['vertices'])
    for record in evidence['resolved_construction_records']:
        operands=[]
        for operand in record['operands']:
            if operand['status']!='UNIQUE' or len(operand['source_vertices'])!=1:
                raise ValueError('AMBIGUOUS_SOURCE_RECIPE')
            candidate=operand['source_vertices'][0]
            operands.append((candidate['input_object'],candidate['input_vertex']))
        old_index=evidence['new_vertex_to_saved_vertex'][record['vertex_index']]
        labels[old_index]=json.dumps(construction_constraint_key(record),separators=(',',':'))
    if None in labels or len(set(labels))!=len(labels):raise ValueError('NONBIJECTIVE_SOURCE_RECIPE_KEYS')
    certificates=patch_boundary_signature(saved['vertices'],saved['faces'])
    boundaries=[];removed=[]
    for patch in certificates['patches']:
        loops=[]
        for cycle in patch:
            core=coalesce_collinear_cycle(saved['vertices'],cycle)
            removed.extend(v for v in cycle if v not in core)
            tokens=tuple(labels[v] for v in core)
            loops.append(min(tokens[i:]+tokens[:i] for i in range(len(tokens))))
        boundaries.append(tuple(sorted(loops)))
    row=dict(index=index,geometry_sha256=evidence['saved_geometry_sha256'],
        evidence_sha256=hashlib.sha256(evidence_path.read_bytes()).hexdigest(),
        unique_same_call_source_keys=True,patch_count=len(boundaries),
        exact_collinear_removed_in_certificate_only=sorted(set(removed)),
        exact_coplanar_removed_edge_count=len(certificates['removed_edges']),
        oriented_patch_boundaries=sorted(boundaries))
    rows.append(row)
    with (output/f'point-{index:03d}.json').open('x') as stream:json.dump(row,stream)
report=dict(points=[r['index'] for r in rows],new_boolean_calls=0,
    exact_feature_patch_boundaries_equal=rows[0]['oriented_patch_boundaries']==rows[1]['oriented_patch_boundaries'],
    full_branch_certification=False,optimization_admission=False,
    scope='Two frozen finite samples; exact witness coalescing, not a complete81 or interval certificate.')
with (output/'report.json').open('x') as stream:json.dump(report,stream,indent=2)
print(json.dumps(report),flush=True)
