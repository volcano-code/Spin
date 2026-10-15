"""Capture ONLY missing exact tapes and certify the existing finite81 samples.

No finite-difference weights/integration sweep, changed parameters or new
optimization. Original numerical records must match field-for-field.
"""
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import numpy as np
from make_it_stand import spin
from birth_face_recipes import replay_vertex_recipe
from exact_planar_patches import oriented_planar_chain_signature
from run_small_inprocess_solver import safe

HERE=Path(__file__).resolve().parent
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
saved=HERE/'noodles-algebraic-production-moments-job-01/output'
sources=HERE/'noodles-interval-algebraic-source-derivative-job-01/output'
captures=HERE/'noodles-same-call-proof-capture-job-01/output'
triage=json.loads((HERE/'noodles-saved-lift-triage-job-01/output/report.json').read_text())
if not triage['all_lifts_passed']:raise ValueError('ALL_SAVED_LIFT_TRIAGE_REQUIRED')
manifest=json.loads((saved/'manifest.json').read_text())
frozen=json.loads((HERE/'noodles-balanced-request-freeze-job-01/output/request.json').read_text())
row=dict(next(x for x in frozen['requests'] if x['settings']['deterministic_seed']==0))
row['settings']=spin.SpinSolverConfig(**row['settings']);request=spin.SpinItFasterRequest(**row)
native,nreq=spin._native_request(request)
settings=spin._native_evaluator_settings(native,spin.SpinEvaluatorConfig(**manifest['settings']))
base=json.loads((captures/'point-000.json').read_text())
recipes=json.loads((HERE/'noodles-offline-proof-pair-job-01/output/point-000-recipes.json').read_text())['recipes']
design=(np.asarray(request.design_vertices)-request.target_origin)@np.asarray(request.target_axes)
design=[tuple(Fraction(float(x)) for x in p) for p in design]
base_source=json.loads((sources/'point-000.json').read_text())['source']
branch_keys=('faces','face_is_cap','source_normal_vertices','source_normal_edges','source_projective_rays',
             'surface_type','truncation_radius','interior_is_positive')
start=time.monotonic();passed=[];calls=0
report=dict(scope='SAME_PROFILE_FINITE81_EXACT_RECIPE_ORIENTED_CHAIN_CERTIFICATE',
    interval_neighborhood_certification=False,optimization_admission=False,
    manifest_sha256=hashlib.sha256((saved/'manifest.json').read_bytes()).hexdigest(),
    native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest())
with (out/'events.jsonl').open('x') as journal:
    def event(data):journal.write(json.dumps(data)+'\n');journal.flush();os.fsync(journal.fileno())
    try:
        for index in range(81):
            event(dict(index=index,status='STARTED',seconds=time.monotonic()-start))
            old_path=saved/f'point-{index:03d}-material.json';old=json.loads(old_path.read_text())
            source=json.loads((sources/f'point-{index:03d}.json').read_text())
            if old['coefficients'] != manifest['samples'][index] or source['coefficients'] != old['coefficients'] or source['source']['scientific_hash'] != old['evaluation']['source_scientific_hash']:
                raise ValueError('SAVED_PARAMETER_SOURCE_IDENTITY_MISMATCH')
            if any(source['source'][k] != base_source[k] for k in branch_keys):
                raise ValueError('SOURCE_CONSTRUCTION_BRANCH_CHANGED')
            if index in (0,1):
                captured=json.loads((captures/f'point-{index:03d}.json').read_text())
                if captured['saved_sha256'] != hashlib.sha256(old_path.read_bytes()).hexdigest():
                    raise ValueError('CAPTURE_REFERENCE_CHANGED')
                proof=captured['proof'];actual=captured['evaluation']
            else:
                calls+=1
                result=native.evaluate_spin_it_faster_fullrank_quadric(nreq,np.asarray(old['coefficients']),settings,capture_material_proof=True)
                proof=safe(result.pop('material_proof_attachment'));actual=safe(result)
                with (out/f'point-{index:03d}-tape.json').open('x') as f:
                    json.dump(dict(index=index,coefficients=old['coefficients'],proof=proof,evaluation=actual,
                        saved_sha256=hashlib.sha256(old_path.read_bytes()).hexdigest()),f,allow_nan=False)
                    f.flush();os.fsync(f.fileno())
            if actual != old['evaluation']:raise ValueError('CAPTURE_CHANGED_ORIGINAL_SCIENTIFIC_FIELDS')
            if proof['geometry_hash'] != actual['geometry_hash'] or proof['source_faces'] != source['source']['faces']:
                raise ValueError('SAME_CALL_PROOF_IDENTITY_MISMATCH')
            operands=[design,[tuple(Fraction(float(x)) for x in p) for p in source['source']['vertices']]]
            points=[replay_vertex_recipe(r,operands) for r in recipes]
            exact=[tuple(Fraction(x) for x in p) for p in proof['exact_material_vertices']]
            if len(set(points)) != len(points) or len(set(exact)) != len(exact) or set(points) != set(exact):
                raise ValueError('ACTUAL_EXACT_MATERIAL_NOT_RECIPE_BIJECTION')
            lookup={p:i for i,p in enumerate(points)}
            current=[[lookup[exact[v]] for v in face] for face in actual['material_faces']]
            a=oriented_planar_chain_signature(points,base['evaluation']['material_faces'])
            b=oriented_planar_chain_signature(points,current)
            if a['signature'] != b['signature']:raise ValueError('ACTUAL_EXACT_ORIENTED_CHAIN_DIFFERS')
            passed.append(index)
            event(dict(index=index,status='ACTUAL_EXACT_SAMPLE_CERTIFIED',seconds=time.monotonic()-start,
                actual_geometry_hash=actual['geometry_hash'],same_original_scientific_fields=True,
                proof_only_null_triangles=[len(a['exact_null_triangle_indices']),len(b['exact_null_triangle_indices'])]))
    except Exception as exc:
        report.update(failure_index=index,exception_type=type(exc).__name__,exception=str(exc))
        event(dict(index=index,status='CERTIFICATION_FAILED',exception=str(exc)))
    report.update(passed_indices=passed,material_boolean_calls=calls,complete_frozen_sample_certification=len(passed)==81,
                  seconds=time.monotonic()-start)
    with (out/'report.json').open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps(report),flush=True)
if not report['complete_frozen_sample_certification']:sys.exit(1)
