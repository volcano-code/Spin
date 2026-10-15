"""Two necessary provenance captures, not another derivative sweep or solve."""
import hashlib,json,os,sys,time
from pathlib import Path
import numpy as np
from make_it_stand import spin
from run_small_inprocess_solver import safe
HERE=Path(__file__).resolve().parent;saved=HERE/'noodles-algebraic-production-moments-job-01/output'
manifest=json.loads((saved/'manifest.json').read_text());frozen=json.loads((HERE/'noodles-balanced-request-freeze-job-01/output/request.json').read_text())
row=dict(next(x for x in frozen['requests'] if x['settings']['deterministic_seed']==0));row['settings']=spin.SpinSolverConfig(**row['settings'])
native,request=spin._native_request(spin.SpinItFasterRequest(**row));settings=spin._native_evaluator_settings(native,spin.SpinEvaluatorConfig(**manifest['settings']))
base=json.loads((saved/'point-000-material.json').read_text())['evaluation']['material_faces']
changed=next(i for i in range(1,81) if json.loads((saved/f'point-{i:03d}-material.json').read_text())['evaluation']['material_faces']!=base)
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False);start=time.monotonic();captures=[]
with (out/'events.jsonl').open('x') as journal:
 def event(row):journal.write(json.dumps(row)+'\n');journal.flush();os.fsync(journal.fileno())
 for index in (0,changed):
    path=saved/f'point-{index:03d}-material.json';old=json.loads(path.read_text());event(dict(index=index,status='CAPTURE_STARTED'))
    result=native.evaluate_spin_it_faster_fullrank_quadric(request,np.asarray(old['coefficients']),settings,capture_material_proof=True)
    proof=safe(result.pop('material_proof_attachment'));actual=safe(result)
    if actual!=old['evaluation']:raise ValueError('BINDING_ONLY_CAPTURE_CHANGED_SCIENTIFIC_RESULT')
    if len(proof['exact_material_vertices'])!=len(actual['material_vertices']):raise ValueError('MISSING_SAME_CALL_EXACT_COORDINATES')
    with (out/f'point-{index:03d}.json').open('x') as f:json.dump(dict(index=index,coefficients=old['coefficients'],evaluation=actual,proof=proof,saved_sha256=hashlib.sha256(path.read_bytes()).hexdigest()),f,allow_nan=False)
    captures.append(index);event(dict(index=index,status='SAME_CALL_PROOF_CAPTURED',seconds=time.monotonic()-start))
 with (out/'report.json').open('x') as f:json.dump(dict(captured_indices=captures,material_boolean_calls=2,all_original_scientific_fields_equal=True,native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),manifest_sha256=hashlib.sha256((saved/'manifest.json').read_bytes()).hexdigest(),scope='Two registered differing-connectivity points only; no whole branch certificate',full_material_branch_certification=False,optimization_admission=False),f,indent=2)
