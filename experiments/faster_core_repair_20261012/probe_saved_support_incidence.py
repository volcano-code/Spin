"""Two preserved samples, exact source-feature incidence; no material Boolean."""
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import sys
import numpy as np
from make_it_stand import spin
from exact_planar_patches import cross,subtract,dot,point_on_triangle,patch_boundary_signature,coalesce_collinear_cycle

HERE=Path(__file__).resolve().parent;inputs=HERE/'small-top-derivative-probe-job-02/output'
run=HERE/'small-top-inprocess-solver-sequence-01/seed-1-job-01/output'
config=json.loads((run/'configuration.json').read_text());req=config['request'];settings=config['evaluator']
design=(np.array(req['design_vertices'])-req['target_origin'])@np.array(req['target_axes'])
native=spin._native_module();tri=native.SpinTriangulationSettings()
for name in ('icosphere_subdivisions','maximum_icosphere_subdivisions','maximum_radius_doublings',
    'retain_fixed_normal_provenance','truncation_radius','paper_projective_hyperboloid','algebraic_finite_map'):
    setattr(tri,name,settings[name])
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
rows=[]
for index in (0,2):
    path=inputs/f'point-{index:03d}-material.json';saved=json.loads(path.read_text())
    source=native.evaluate_spin_it_faster_source_motion(saved['actual_coefficients'],design,tri)
    branch=dict(surface_type=source['surface_type'],radius=source['truncation_radius'],cap_count=source['cap_face_count'],
        faces=np.asarray(source['faces']).tolist(),normal_vertices=np.asarray(source['source_normal_vertices']).tolist(),
        normal_edges=np.asarray(source['source_normal_edges']).tolist(),
        projective_ray_ids=np.asarray(source['source_projective_rays']).tolist())
    if branch!=saved['source_branch']:raise ValueError('FROZEN_SOURCE_BRANCH_NOT_REPRODUCED')
    planes=[]
    for operand,(vertices,faces) in enumerate(((design,req['design_faces']),(source['vertices'],source['faces']))):
        vertices=[tuple(Fraction(float(x)) for x in row) for row in vertices];groups={}
        for face_index,face in enumerate(faces):
            triangle=tuple(vertices[int(v)] for v in face);a,b,c=triangle
            n=cross(subtract(b,a),subtract(c,a))
            if n==(0,0,0):raise ValueError('EXACT_ZERO_SOURCE_FACE')
            pivot=next(x for x in n if x!=0);d=dot(n,a)
            key=tuple(x/pivot for x in n)+(d/pivot,)
            if key not in groups:groups[key]=dict(identity=(operand,face_index),normal=n,d=d,triangles=[])
            groups[key]['triangles'].append(triangle)
        planes.extend(groups.values())
    signatures=[];ranks=[]
    for coordinates in saved['vertices']:
        p=tuple(Fraction(x) for x in coordinates);supports=[]
        for plane in planes:
            if dot(plane['normal'],p)==plane['d'] and any(point_on_triangle(p,t) for t in plane['triangles']):
                supports.append(plane)
        signatures.append(tuple(sorted(s['identity'] for s in supports)))
        ranks.append(any(dot(cross(a['normal'],b['normal']),c['normal'])!=0
            for a,b,c in itertools.combinations(supports,3)))
    unique=len(set(signatures))==len(signatures);correspondence=unique and all(ranks)
    row=dict(index=index,input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        source_branch_reproduced=True,unique_support_keys=unique,all_rank3=all(ranks),support_keys=signatures,
        unsupported_vertices=[i for i,s in enumerate(signatures) if not s],
        correspondence_certified=correspondence)
    patches=patch_boundary_signature(saved['vertices'],saved['faces'])['patches']
    reduced=[[coalesce_collinear_cycle(saved['vertices'],cycle) for cycle in patch] for patch in patches]
    retained=sorted({v for patch in reduced for cycle in patch for v in cycle})
    retained_keys=[signatures[v] for v in retained]
    core_valid=len(set(retained_keys))==len(retained_keys) and all(ranks[v] for v in retained)
    row.update(retained_corner_count=len(retained),exact_redundant_vertex_count=len(signatures)-len(retained),
               correspondence_certified=core_valid)
    if core_valid:
        def oriented_signature(cycle):
            loop=tuple(signatures[v] for v in cycle)
            return min(loop[i:]+loop[:i] for i in range(len(loop)))
        row['coalesced_patch_boundaries']=sorted(tuple(sorted(oriented_signature(cycle) for cycle in patch)) for patch in reduced)
    rows.append(row)
    with (output/f'point-{index:03d}.json').open('x') as stream:json.dump(row,stream)
report=dict(points=[r['index'] for r in rows],material_boolean_calls=0,source_replay_calls=2,
    native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
    matched_patch_boundaries=all(r.get('coalesced_patch_boundaries')==rows[0].get('coalesced_patch_boundaries') for r in rows)
        and all(r['correspondence_certified'] for r in rows),
    full_branch_certification=False,optimization_admission=False,
    scope='Exact support incidence on two saved finite samples; not interval/neighborhood certification.')
with (output/'report.json').open('x') as stream:json.dump(report,stream,indent=2)
print(json.dumps(report),flush=True)
