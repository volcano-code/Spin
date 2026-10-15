"""Actual saved formal run evidence, with one corruption per negative test."""
import importlib
import importlib.util
import json
from pathlib import Path
import shutil
import hashlib
import os
import pytest

ROOT=Path(__file__).resolve().parents[2]
RUN=Path(os.environ.get('MAKE_IT_STAND_FORMAL_EVIDENCE_RUN',
    str(ROOT/'experiments/faster_core_repair_20261012/formal-cherchi-small-closure-job-04/output/run')))

def api():
    name='make_it_stand.verification.cherchi_run_evidence'
    assert importlib.util.find_spec(name) is not None, 'formal evidence verifier missing'
    return importlib.import_module(name)

def test_actual_completed_run_is_valid():
    result=api().verify_cherchi_run_evidence(RUN)
    # Independently count the public solver trace, not the verifier's event
    # tape. Seed0 has 24 updates; the historical seed1 fixture has 17.
    native=json.loads((RUN/'result.json').read_text())['identity']['native_run']
    expected_updates=sum(row['accepted_geodesic_step']>0 for row in native['trace'])
    assert expected_updates>0
    assert result['record_integrity_valid'] and result['accepted_updates']==expected_updates
    assert not result['core_admission']

@pytest.mark.parametrize('fault',['missing_log','tampered_material','wrong_backend','wrong_source','missing_update','wrong_terminal','stale_jacobian'])
def test_corrupted_run_is_rejected(tmp_path,fault):
    module=api();p=tmp_path/'run';shutil.copytree(RUN,p)
    events=[json.loads(x) for x in (p/'events.jsonl').read_text().splitlines()]
    if fault=='missing_log':(p/'events.jsonl').unlink()
    elif fault=='missing_update':
        events=[e for e in events if not(e['kind']=='NEWTON_EVENT' and e['event']['type']=='evaluation' and (e['event'].get('step') or 0)>0)]
        (p/'events.jsonl').write_text(''.join(json.dumps(e)+'\n' for e in events))
    elif fault=='wrong_terminal':
        r=json.loads((p/'terminal-material-reference.json').read_text());r['evaluation_file']='evaluation-000001.json'
        (p/'terminal-material-reference.json').write_text(json.dumps(r))
    elif fault=='wrong_backend':
        c=json.loads((p/'configuration.json').read_text());c['boolean_backend_identity']='cgal'
        (p/'configuration.json').write_text(json.dumps(c))
    else:
        f=p/'evaluation-000018.json';m=json.loads(f.read_text())
        if fault=='tampered_material':m['variable_vertices'][0][0]+=0.01
        elif fault=='wrong_source':m['quadric_source_identity']='wrong'
        else:m['moment_jacobian'][0][0]+=0.001
        f.write_text(json.dumps(m))
    with pytest.raises(ValueError):module.verify_cherchi_run_evidence(p)

def test_scientific_projection_only_ignores_named_costs():
    module=api()
    a={'event':{'objective':2.,'evaluator_wall_seconds':1.,'moment_jacobian':[[1.]],'backend':'cherchi'},'monotonic_ns':10}
    b=json.loads(json.dumps(a));b['monotonic_ns']=20;b['event']['evaluator_wall_seconds']=9.
    assert module.scientific_digest(a)==module.scientific_digest(b)
    for key,value in [('moment_jacobian',[[2.]]),('backend','cgal'),('objective',3.)]:
        b=json.loads(json.dumps(a));b['event'][key]=value
        assert module.scientific_digest(a)!=module.scientific_digest(b)

def test_rehashed_stale_jacobian_still_disagrees_with_solver_payload(tmp_path):
    module=api();p=tmp_path/'run';shutil.copytree(RUN,p)
    name='evaluation-000018.json';f=p/name;m=json.loads(f.read_text())
    m['moment_jacobian'][0][0]+=0.001;f.write_text(json.dumps(m))
    sha=hashlib.sha256(f.read_bytes()).hexdigest()
    events=[json.loads(x) for x in (p/'events.jsonl').read_text().splitlines()]
    for e in events:
        if e.get('evaluation_file')==name:e['sha256']=sha
    (p/'events.jsonl').write_text(''.join(json.dumps(e)+'\n' for e in events))
    terminal=json.loads((p/'terminal-material-reference.json').read_text());terminal['sha256']=sha
    (p/'terminal-material-reference.json').write_text(json.dumps(terminal))
    with pytest.raises(ValueError,match='JACOBIAN_ATTACHMENT_MISMATCH'):
        module.verify_cherchi_run_evidence(p)


def test_claimed_native_arithmetic_must_match_saved_materials(tmp_path):
    module=api();p=tmp_path/'run';shutil.copytree(RUN,p)
    config=json.loads((p/'configuration.json').read_text())
    config['rational_arithmetic_backend']='NATIVE_GMP_SURFACE_V1'
    config['rational_arithmetic_identity']=dict(backend='NATIVE_GMP_SURFACE_V1',
        algorithm='DIV_SURFACE_V1_REYNOLDS_V1_GMP_BINARY64_NEAREST_EVEN_V1',module_sha256='a'*64)
    (p/'configuration.json').write_text(json.dumps(config))
    with pytest.raises(ValueError,match='ARITHMETIC_ATTACHMENT_MISMATCH'):
        module.verify_cherchi_run_evidence(p)


def test_claimed_trial_boundary_profile_must_match_actual_solver_record(tmp_path):
    module=api();p=tmp_path/'run';shutil.copytree(RUN,p)
    config=json.loads((p/'configuration.json').read_text())
    actual=config.get('solver_boundary_profile','EXACT_NEF_REGION_V1')
    config['solver_boundary_profile']={'EXACT_NEF_REGION_V1':'EXACT_RATIONAL_BOUNDARY_V1',
        'EXACT_RATIONAL_BOUNDARY_V1':'EXACT_NEF_REGION_V1'}[actual]
    (p/'configuration.json').write_text(json.dumps(config))
    with pytest.raises(ValueError,match='SOLVER_BOUNDARY_PROFILE_MISMATCH'):
        module.verify_cherchi_run_evidence(p)


def test_claimed_durability_profile_must_match_actual_run(tmp_path):
    module=api();p=tmp_path/'run';shutil.copytree(RUN,p)
    config=json.loads((p/'configuration.json').read_text())
    actual=config.get('durability_profile','PER_EVENT_FSYNC_V1')
    config['durability_profile']={'PER_EVENT_FSYNC_V1':'EVALUATION_CHECKPOINT_V1',
        'EVALUATION_CHECKPOINT_V1':'PER_EVENT_FSYNC_V1'}[actual]
    (p/'configuration.json').write_text(json.dumps(config))
    with pytest.raises(ValueError,match='RUN_DURABILITY_PROFILE_MISMATCH'):
        module.verify_cherchi_run_evidence(p)


def test_input_reuse_profile_cannot_be_added_by_changing_only_configuration(tmp_path):
    module=api();p=tmp_path/'run';shutil.copytree(RUN,p)
    config=json.loads((p/'configuration.json').read_text())
    config['input_validation_profile']='CONSTRUCTOR_PROOF_REUSE_V1'
    (p/'configuration.json').write_text(json.dumps(config))
    with pytest.raises(ValueError,match='INPUT_VALIDATION_PROFILE_MISMATCH'):
        module.verify_cherchi_run_evidence(p)


def test_input_reuse_profile_must_bind_the_actual_problem_identity(tmp_path):
    module=api();p=tmp_path/'run';shutil.copytree(RUN,p)
    config=json.loads((p/'configuration.json').read_text())
    config['input_validation_profile']='CONSTRUCTOR_PROOF_REUSE_V1'
    config['input_validation_problem_digest']='0'*64
    result=json.loads((p/'result.json').read_text())
    result['identity']['input_validation_profile']='CONSTRUCTOR_PROOF_REUSE_V1'
    (p/'configuration.json').write_text(json.dumps(config))
    (p/'result.json').write_text(json.dumps(result))
    with pytest.raises(ValueError,match='INPUT_VALIDATION_PROBLEM_IDENTITY_MISMATCH'):
        module.verify_cherchi_run_evidence(p)


def test_terminal_strategy_cannot_be_changed_only_in_configuration(tmp_path):
    module=api();p=tmp_path/'run';shutil.copytree(RUN,p)
    config=json.loads((p/'configuration.json').read_text())
    config['terminal_region_profile']='OWNED_INTERIOR_DIFFERENCE_V1'
    (p/'configuration.json').write_text(json.dumps(config))
    with pytest.raises(ValueError,match='TERMINAL_REGION_PROFILE_MISMATCH'):
        module.verify_cherchi_run_evidence(p)
