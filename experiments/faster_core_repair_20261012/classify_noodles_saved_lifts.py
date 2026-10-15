"""Triage all existing samples using exact recipe lifts, not missing native tapes.

Float-view correspondence cannot certify the original unrecorded exact Boolean
vertices. This pilot only decides whether provenance capture has a coherent
fixed recipe to test, without rerunning any numerical difference or Boolean.
"""
import hashlib
import json
import os
from fractions import Fraction
from pathlib import Path
import sys
import time
import numpy as np
from birth_face_recipes import replay_vertex_recipe
from exact_planar_patches import oriented_planar_chain_signature

HERE = Path(__file__).resolve().parent
out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=False)
saved = HERE/'noodles-algebraic-production-moments-job-01/output'
sources = HERE/'noodles-interval-algebraic-source-derivative-job-01/output'
manifest = json.loads((saved/'manifest.json').read_text())
req = json.loads((HERE/'noodles-balanced-request-freeze-job-01/output/request.json').read_text())['requests'][0]
design = (np.asarray(req['design_vertices'])-req['target_origin']) @ np.asarray(req['target_axes'])
design = [tuple(Fraction(float(x)) for x in p) for p in design]
base = json.loads((HERE/'noodles-same-call-proof-capture-job-01/output/point-000.json').read_text())
recipes = json.loads((HERE/'noodles-offline-proof-pair-job-01/output/point-000-recipes.json').read_text())['recipes']
start=time.monotonic(); passed=[]
report=dict(scope='EXACT_RECIPE_LIFT_AND_BINARY64_VIEW_TRIAGE_ONLY', material_boolean_calls=0,
    original_exact_material_certification=False, optimization_admission=False)
with (out/'events.jsonl').open('x') as log:
    try:
        for index in range(81):
            actual=json.loads((saved/f'point-{index:03d}-material.json').read_text())
            source=json.loads((sources/f'point-{index:03d}.json').read_text())
            if actual['coefficients'] != manifest['samples'][index] or source['coefficients'] != actual['coefficients'] or source['source']['scientific_hash'] != actual['evaluation']['source_scientific_hash']:
                raise ValueError('SOURCE_MATERIAL_SAMPLE_IDENTITY_MISMATCH')
            if source['source']['faces'] != base['proof']['source_faces']:
                raise ValueError('SOURCE_CONNECTIVITY_CHANGED')
            operands=[design,[tuple(Fraction(float(x)) for x in p) for p in source['source']['vertices']]]
            points=[replay_vertex_recipe(r,operands) for r in recipes]
            rounded=[tuple(float(x) for x in p) for p in points]
            old=[tuple(p) for p in actual['evaluation']['material_vertices']]
            if len(set(rounded)) != len(points) or len(set(old)) != len(old) or set(rounded) != set(old):
                raise ValueError('EXACT_LIFT_BINARY64_VIEW_NOT_BIJECTIVE')
            lookup={p:i for i,p in enumerate(rounded)}
            current=[[lookup[old[v]] for v in face] for face in actual['evaluation']['material_faces']]
            a=oriented_planar_chain_signature(points,base['evaluation']['material_faces'])
            b=oriented_planar_chain_signature(points,current)
            if a['signature'] != b['signature']:
                raise ValueError('EXACT_LIFT_ORIENTED_CHAIN_DIFFERS')
            passed.append(index)
            log.write(json.dumps(dict(index=index,status='LIFT_TRIAGE_PASSED_NOT_NATIVE_TAPE_CERTIFIED',seconds=time.monotonic()-start))+'\n')
            log.flush();os.fsync(log.fileno())
    except Exception as exc:
        report.update(failure_index=index,exception_type=type(exc).__name__,exception=str(exc))
    report.update(passed_indices=passed,all_lifts_passed=len(passed)==81,seconds=time.monotonic()-start,
        manifest_sha256=hashlib.sha256((saved/'manifest.json').read_bytes()).hexdigest())
    with (out/'report.json').open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps(report),flush=True)
if not report['all_lifts_passed']:sys.exit(1)
