"""Finite saved samples: exact construction replay and oriented chain equality.

No Boolean or source reconstruction. Does not certify a parameter neighborhood.
"""
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from cherchi_core_map import replay_lpi_constraint
from make_it_stand.verification.cherchi_inprocess import construction_constraint_key,exact_fraction
from exact_planar_patches import oriented_planar_chain_signature

HERE=Path(__file__).resolve().parent
inputs=HERE/'formal-cherchi-derivative-job-01/output'
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
manifest=json.loads((inputs/'manifest.json').read_text())
base=json.loads((inputs/'point-000-material.json').read_text())
material=base['same_call_material']
records=sorted(material['material_proof']['resolved_construction_records'],key=lambda r:r['vertex_index'])
if [r['vertex_index'] for r in records]!=list(range(len(records))):raise ValueError('NONCONTIGUOUS_RECORDS')
if len(set(construction_constraint_key(r) for r in records))!=len(records):raise ValueError('AMBIGUOUS_FEATURES')
design=json.loads((HERE/'small-top-inprocess-solver-sequence-01/seed-1-job-01/output/configuration.json').read_text())['request']['design_vertices']
passed=[]
with (output/'events.jsonl').open('x') as events:
    for index,q in enumerate(manifest['samples']):
        path=inputs/f'point-{index:03d}-material.json'
        saved=json.loads(path.read_text());m=saved['same_call_material']
        try:
            if saved['actual_coefficients']!=q or m['quadric_coefficients']!=q:raise ValueError('PARAMETER_MISMATCH')
            if saved['source_branch']!=base['source_branch']:raise ValueError('SOURCE_BRANCH_CHANGED')
            operands=[[tuple(Fraction(float(x)) for x in row) for row in vertices]
                      for vertices in (design,m['quadric_source_vertices'])]
            points=[]
            for r in records:
                if r['point_type']=='EXPLICIT3D':
                    v=r['operands'][0]['source_vertices'][0]
                    p=operands[v['input_object']][v['input_vertex']]
                else:p=replay_lpi_constraint(r,operands)
                points.append(p)
            actual=[tuple(exact_fraction(x) for x in row) for row in m['material_proof']['exact_vertices']]
            if len(set(points))!=len(points) or set(points)!=set(actual):raise ValueError('RECIPE_NOT_BIJECTIVE')
            correspondence={p:i for i,p in enumerate(points)}
            current=[[correspondence[actual[v]] for v in face] for face in m['variable_faces']]
            reference=material['variable_faces']
            a=oriented_planar_chain_signature(points,reference)
            b=oriented_planar_chain_signature(points,current)
            if a['signature']!=b['signature']:raise ValueError('ORIENTED_SURFACE_CHAIN_DIFFERS')
            row=dict(index=index,status='FINITE_SAMPLE_CERTIFIED',
                classification='IDENTICAL_OR_RENUMBERED' if current==reference else 'EXACT_COPLANAR_RETRIANGULATION',
                material_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            passed.append(row)
        except Exception as exc:
            events.write(json.dumps(dict(index=index,status='FAILED',error=f'{type(exc).__name__}: {exc}'))+'\n')
            events.flush();os.fsync(events.fileno());raise
        events.write(json.dumps(row)+'\n');events.flush();os.fsync(events.fileno())
with (output/'report.json').open('x') as stream:
    json.dump(dict(point_count=len(passed),finite_sample_certification=len(passed)==81,
        neighborhood_certification=False,boolean_calls=0,source_calls=0,
        manifest_sha256=hashlib.sha256((inputs/'manifest.json').read_bytes()).hexdigest()),stream)
