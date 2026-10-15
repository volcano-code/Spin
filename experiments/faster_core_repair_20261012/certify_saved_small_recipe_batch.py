"""Offline exact source-recipe and surface-chain checks; zero material Boolean calls."""
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
import numpy as np
from make_it_stand import spin
from cherchi_core_map import replay_lpi_constraint,construction_constraint_key
from exact_planar_patches import oriented_planar_chain_signature
HERE=Path(__file__).resolve().parent;inputs=HERE/'small-top-derivative-probe-job-02/output'
capture=HERE/'saved-small-provenance-capture-01/point-000.json'
evidence=json.loads(capture.read_text());base=json.loads((inputs/'point-000-material.json').read_text())
config=json.loads((HERE/'small-top-inprocess-solver-sequence-01/seed-1-job-01/output/configuration.json').read_text())
req=config['request'];settings=config['evaluator'];native=spin._native_module()
receipt=json.loads((capture.parent/'report.json').read_text())
if hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest()!=receipt['native_sha256']:
    raise ValueError('SOURCE_REPLAY_BUILD_MISMATCH')
if hashlib.sha256((inputs/'point-000-material.json').read_bytes()).hexdigest()!=evidence['saved_geometry_sha256']:
    raise ValueError('BASE_CAPTURE_MISMATCH')
design=(np.asarray(req['design_vertices'])-req['target_origin'])@np.asarray(req['target_axes'])
tri=native.SpinTriangulationSettings()
for name in ('icosphere_subdivisions','maximum_icosphere_subdivisions','maximum_radius_doublings',
    'retain_fixed_normal_provenance','truncation_radius','paper_projective_hyperboloid','algebraic_finite_map'):
    setattr(tri,name,settings[name])
records=sorted(evidence['resolved_construction_records'],key=lambda r:r['vertex_index'])
if len(set(construction_constraint_key(r) for r in records))!=len(records):raise ValueError('AMBIGUOUS_BASE_FEATURES')
inverse={old:new for new,old in enumerate(evidence['new_vertex_to_saved_vertex'])}
base_pattern=[[inverse[v] for v in face] for face in base['faces']]
manifest=json.loads((inputs/'manifest.json').read_text());output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
passed=[]
with (output/'events.jsonl').open('x') as events:
    for index,q in enumerate(manifest['samples']):
        path=inputs/f'point-{index:03d}-material.json';saved=json.loads(path.read_text())
        try:
            if saved['actual_coefficients']!=q:raise ValueError('SAVED_PARAMETER_MISMATCH')
            source=native.evaluate_spin_it_faster_source_motion(q,design,tri)
            branch=dict(surface_type=source['surface_type'],radius=source['truncation_radius'],cap_count=source['cap_face_count'],
                faces=np.asarray(source['faces']).tolist(),normal_vertices=np.asarray(source['source_normal_vertices']).tolist(),
                normal_edges=np.asarray(source['source_normal_edges']).tolist(),projective_ray_ids=np.asarray(source['source_projective_rays']).tolist())
            if branch!=saved['source_branch']:raise ValueError('SOURCE_BRANCH_REPLAY_MISMATCH')
            operands=[[tuple(Fraction(float(x)) for x in row) for row in vertices] for vertices in (design,source['vertices'])]
            points=[]
            for r in records:
                if r['point_type']=='EXPLICIT3D':
                    vertex=r['operands'][0]['source_vertices'][0];p=operands[vertex['input_object']][vertex['input_vertex']]
                else:p=replay_lpi_constraint(r,operands)
                points.append(p)
            old=[tuple(Fraction(x) for x in row) for row in saved['vertices']]
            if len(set(points))!=len(points) or set(points)!=set(old):raise ValueError('SOURCE_RECIPE_NOT_BIJECTIVE_TO_SAVED_MATERIAL')
            correspondence={p:i for i,p in enumerate(points)}
            current=[[correspondence[old[v]] for v in face] for face in saved['faces']]
            a=oriented_planar_chain_signature(points,base_pattern);b=oriented_planar_chain_signature(points,current)
            if a['signature']!=b['signature']:raise ValueError('FIXED_GEOMETRY_ORIENTED_SURFACE_CHAIN_DIFFERS')
            row=dict(index=index,status='OFFLINE_RECIPE_CHAIN_CERTIFIED',material_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                exact_null_counts_in_proof_view=[len(a['exact_null_triangle_indices']),len(b['exact_null_triangle_indices'])],
                reference_recipe_to_saved_vertex=[old.index(p) for p in points])
            passed.append(index)
        except Exception as exc:
            events.write(json.dumps(dict(index=index,status='OFFLINE_CERTIFICATION_FAILED',exception_type=type(exc).__name__,exception=str(exc)))+'\n')
            events.flush();os.fsync(events.fileno());raise
        events.write(json.dumps(row)+'\n');events.flush();os.fsync(events.fileno())
    report=dict(point_count=len(passed),material_boolean_calls=0,source_replay_calls=len(passed),
        complete_frozen_sample_source_recipe_chain_certification=len(passed)==manifest['point_count']==81,
        interval_neighborhood_certification=False,optimization_admission=False,
        manifest_sha256=hashlib.sha256((inputs/'manifest.json').read_bytes()).hexdigest(),
        base_capture_sha256=hashlib.sha256(capture.read_bytes()).hexdigest(),
        scope='All saved finite81 samples only; no geometry repair, new derivative sweep or complete solver admission.')
    with (output/'report.json').open('x') as stream:json.dump(report,stream,indent=2)
    print(json.dumps(report),flush=True)
