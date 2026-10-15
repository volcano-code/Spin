"""Read-only formal Cherchi log integrity, not an independent physics auditor."""
import hashlib
import json
from pathlib import Path
import re
import math
import numpy as np

TIMING_FIELDS=frozenset(('monotonic_ns','elapsed_ns','evaluator_wall_seconds',
                         'residual_evaluation_seconds','derivative_evaluation_seconds'))

def scientific_projection(value):
    if isinstance(value,dict):return {k:scientific_projection(v) for k,v in value.items() if k not in TIMING_FIELDS}
    if isinstance(value,list):return [scientific_projection(v) for v in value]
    return value

def _digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

def scientific_digest(value):return _digest(scientific_projection(value))


def _check_jacobian_attachment(material, event):
    payload = event.get('evaluation_payload', {})
    J = np.asarray(material['moment_jacobian'], dtype=float)
    solver_J = np.asarray(payload.get('solver_moment_jacobian'), dtype=float)
    s = np.asarray(payload.get('solver_coordinates'), dtype=float)
    if J.shape != (10, 10) or solver_J.shape != (10, 10) or s.shape != (10,):
        raise ValueError('JACOBIAN_ATTACHMENT_MISMATCH')
    if event.get('solver_coordinate_identity') != 'literal_monomial_unit_sphere_v1':
        raise ValueError('UNSUPPORTED_SOLVER_COORDINATES')
    inverse_scale = np.array([1.]*7+[.5]*3)
    p = s*inverse_scale; length = np.linalg.norm(p)
    if not math.isfinite(length) or length <= 0:
        raise ValueError('INVALID_SOLVER_COORDINATES')
    p = p/length
    expected = J@((np.eye(10)-np.outer(p, p))*inverse_scale[None, :]/length)
    if (not np.isfinite(expected).all() or not np.isfinite(solver_J).all() or
            np.max(np.abs(expected-solver_J)) > 1e-10):
        raise ValueError('JACOBIAN_ATTACHMENT_MISMATCH')


def load_cherchi_terminal_prefix(directory):
    """Read a durable numerical terminal for independent postprocessing.

    Supports historical interrupted journals without a native-result file.
    Returns the actual saved terminal event, never a synthesized solver run.
    Record integrity is not independent material or solver certification.
    """
    root = Path(directory)
    def require(condition, code):
        if not condition:
            raise ValueError(code)
    def read(name):
        path = root/name
        require(path.is_file() and not path.is_symlink(), 'MISSING_OR_LINKED_RUN_FILE:'+name)
        return json.loads(path.read_text(), parse_constant=lambda x:
                         (_ for _ in ()).throw(ValueError('NONFINITE_JSON')))
    config = read('configuration.json')
    require(config.get('schema') == 'cherchi-inprocess-formal-run-v1', 'UNSUPPORTED_RUN_CONFIGURATION')
    log = root/'events.jsonl'
    require(log.is_file() and not log.is_symlink(), 'MISSING_RUN_LOG')
    events = [json.loads(line, parse_constant=lambda x:
                        (_ for _ in ()).throw(ValueError('NONFINITE_JSON')))
              for line in log.read_text().splitlines()]
    newton = [row['event'] for row in events if row.get('kind') == 'NEWTON_EVENT']
    require(bool(newton) and newton[-1].get('type') == 'attempt_terminal' and
            newton[-1].get('reason') == 'converged', 'MISSING_NUMERICAL_CONVERGENCE')
    accepted = [row for row in newton if row.get('type') == 'evaluation' and
                row.get('outcome') == 'accepted' and (row.get('step') or 0.) > 0.]
    require(bool(accepted), 'NO_ACCEPTED_UPDATE')
    terminal = accepted[-1]
    residual = terminal.get('raw_equation_residual'); mu = terminal.get('multiplier_mu')
    require(residual is not None and math.isfinite(residual) and 0 <= residual <= 1e-12 and
            mu is not None and math.isfinite(mu) and mu > 0, 'NUMERICAL_TERMINAL_NOT_QUALIFIED')
    require(terminal.get('evaluation_index') == newton[-1].get('evaluation_index') and
            terminal.get('restart') == newton[-1].get('restart'), 'TERMINAL_EVENT_MISMATCH')
    mapped = [row for row in events if row.get('kind') in ('MAP_SAVED', 'MAP_REJECTED')]
    require([row['evaluation'] for row in mapped] == list(range(1, len(mapped)+1)), 'MISSING_MAP_EVENT')
    selected = None
    for record in mapped:
        if record['kind'] != 'MAP_SAVED':
            continue
        name = record['evaluation_file']
        require(name == f"evaluation-{record['evaluation']:06d}.json", 'INVALID_MATERIAL_REFERENCE')
        path = root/name
        require(path.is_file() and not path.is_symlink(), 'MISSING_OR_LINKED_RUN_FILE:'+name)
        require(hashlib.sha256(path.read_bytes()).hexdigest() == record['sha256'], 'MATERIAL_FILE_DIGEST_MISMATCH')
        if record['evaluation'] == terminal['evaluation_index']:
            selected = record
    require(selected is not None, 'TERMINAL_NOT_SAVED')
    material = read(selected['evaluation_file'])
    payload = terminal.get('evaluation_payload', {})
    require(material['quadric_coefficients'] == selected['coefficients'] == terminal['quadric_coefficients'],
            'TERMINAL_COEFFICIENT_MISMATCH')
    require(material['geometry_hash'] == selected['geometry_hash'] == terminal['material_geometry_hash'] ==
            payload.get('geometry_hash'), 'TERMINAL_GEOMETRY_MISMATCH')
    for field in ('boolean_backend_identity', 'evaluation_map_identity'):
        require(material[field] == config[field] == payload.get(field), 'MATERIAL_CONTEXT_MISMATCH')
    require(material.get('rational_arithmetic_identity') == config.get('rational_arithmetic_identity'),
            'ARITHMETIC_ATTACHMENT_MISMATCH')
    require(material.get('derivative_defined') is True, 'DERIVATIVE_UNAVAILABLE')
    require(material['quadric_source_identity'] == _digest(dict(vertices=material['quadric_source_vertices'],
            faces=material['quadric_source_faces'], coefficients=material['quadric_coefficients'])),
            'SOURCE_IDENTITY_MISMATCH')
    require(material['geometry_hash'] == _digest(dict(vertices=material['material_proof']['exact_vertices'],
            faces=material['variable_faces'])), 'ACTUAL_GEOMETRY_DIGEST_MISMATCH')
    observed = payload.get('solver_raw_moments', [])
    require(len(observed) == 10 and len(material['raw_moments']) == 10 and
            all(math.isfinite(a) and math.isfinite(b) and abs(a-b) <= 1e-10
                for a, b in zip(observed, material['raw_moments'])) and
            payload.get('variable_raw_moments') == material['variable_raw_moments'],
            'MOMENT_ATTACHMENT_MISMATCH')
    _check_jacobian_attachment(material, terminal)
    if (root/'terminal-material-reference.json').exists():
        reference = read('terminal-material-reference.json')
        require(all(reference.get(key) == selected[key] for key in ('evaluation_file', 'sha256', 'geometry_hash')) and
                reference.get('coefficients') == terminal['quadric_coefficients'], 'TERMINAL_REFERENCE_MISMATCH')
    if (root/'native-result.json').exists():
        native = read('native-result.json')
        require(native.get('configuration_sha256') == hashlib.sha256((root/'configuration.json').read_bytes()).hexdigest(),
                'NATIVE_CONFIGURATION_MISMATCH')
        require(native['native_run'].get('status') == 'CONVERGED' and
                native['native_run'].get('quadric_coefficients') == terminal['quadric_coefficients'] and
                scientific_projection(native['native_run']['events']) == scientific_projection(newton),
                'NATIVE_EVENT_LOG_MISMATCH')
    return dict(configuration=config, material=material, terminal_event=terminal,
        terminal_reference=dict(coefficients=terminal['quadric_coefficients'], **{
            key: selected[key] for key in ('evaluation_file', 'sha256', 'geometry_hash')}),
        accepted_updates=len(accepted), record_integrity_valid=True,
        event_log_sha256=hashlib.sha256(log.read_bytes()).hexdigest(),
        configuration_sha256=hashlib.sha256((root/'configuration.json').read_bytes()).hexdigest(),
        material_audit_complete=False, run_complete=False,
        scope='Numerical terminal checkpoint only; independent postprocessing required')

def verify_cherchi_run_evidence(directory):
    """Require a completed converged run with real updates and bound materials.

    Checks accidental corruption/inconsistent attachments, not authenticity
    against an adversary able to rewrite all evidence and hashes together.
    """
    root=Path(directory)
    def require(condition,code):
        if not condition:raise ValueError(code)
    def read(name):
        p=root/name
        require(p.is_file() and not p.is_symlink(),'MISSING_OR_LINKED_RUN_FILE:'+name)
        try:return json.loads(p.read_text(),parse_constant=lambda x:(_ for _ in ()).throw(ValueError('NONFINITE_JSON')))
        except (OSError,json.JSONDecodeError) as exc:raise ValueError('INVALID_RUN_FILE:'+name) from exc
    config=read('configuration.json');complete=read('complete.json');result=read('result.json')
    terminal_profile=config.get('terminal_region_profile','INDEXED_BOUNDARY_REBUILD_V1')
    require(terminal_profile in ('INDEXED_BOUNDARY_REBUILD_V1','OWNED_INTERIOR_DIFFERENCE_V1') and
            terminal_profile==result['identity'].get('terminal_region_profile','INDEXED_BOUNDARY_REBUILD_V1'),
            'TERMINAL_REGION_PROFILE_MISMATCH')
    input_profile=config.get('input_validation_profile','RECHECK_INPUT_V1')
    require(input_profile in ('RECHECK_INPUT_V1','CONSTRUCTOR_PROOF_REUSE_V1') and
            input_profile==result['identity'].get('input_validation_profile','RECHECK_INPUT_V1'),
            'INPUT_VALIDATION_PROFILE_MISMATCH')
    if input_profile=='CONSTRUCTOR_PROOF_REUSE_V1':
        input_digest=config.get('input_validation_problem_digest')
        require(bool(re.fullmatch(r'[0-9a-f]{64}',str(input_digest))) and
                input_digest==result['identity'].get('problem_digest'),
                'INPUT_VALIDATION_PROBLEM_IDENTITY_MISMATCH')
    durability=config.get('durability_profile','PER_EVENT_FSYNC_V1')
    require(durability in ('PER_EVENT_FSYNC_V1','EVALUATION_CHECKPOINT_V1') and
            durability==result.get('durability_profile','PER_EVENT_FSYNC_V1'),
            'RUN_DURABILITY_PROFILE_MISMATCH')
    boundary_profile=config.get('solver_boundary_profile','EXACT_NEF_REGION_V1')
    require(boundary_profile in ('EXACT_NEF_REGION_V1','EXACT_RATIONAL_BOUNDARY_V1','EXACT_RATIONAL_FIXED_SUBSET_V1','EXACT_RATIONAL_FIXED_SEPARATED_ROOTS_V1','EXACT_RATIONAL_FIXED_CLASSIFIED_V1') and
            boundary_profile==result['identity'].get('solver_boundary_profile','EXACT_NEF_REGION_V1'),
            'SOLVER_BOUNDARY_PROFILE_MISMATCH')
    if boundary_profile in ('EXACT_RATIONAL_FIXED_SUBSET_V1','EXACT_RATIONAL_FIXED_SEPARATED_ROOTS_V1','EXACT_RATIONAL_FIXED_CLASSIFIED_V1'):
        proof_identity=config.get('fixed_boundary_proof_identity')
        require(bool(re.fullmatch(r'[0-9a-f]{64}',str(proof_identity))) and
                proof_identity==result['identity'].get('fixed_boundary_proof_identity'),
                'FIXED_BOUNDARY_PROOF_IDENTITY_MISMATCH')
    validation_profile=config.get('validation_profile','INDEPENDENT_EACH_EVALUATION_V1')
    require(validation_profile in ('INDEPENDENT_EACH_EVALUATION_V1','LAYERED_SAME_CALL_V1') and
            validation_profile==result['identity'].get('validation_profile','INDEPENDENT_EACH_EVALUATION_V1'),
            'MATERIAL_VALIDATION_PROFILE_MISMATCH')
    if validation_profile=='LAYERED_SAME_CALL_V1':
        require(bool(re.fullmatch(r'[0-9a-f]{64}', str(result['identity'].get(
            'fixed_material_preparation_identity','')))), 'MISSING_FIXED_MATERIAL_PREPARATION')
    backend=config['boolean_backend_identity'];map_identity=config['evaluation_map_identity']
    arithmetic_backend=config.get('rational_arithmetic_backend','PYTHON_FRACTION_V1')
    arithmetic_identity=config.get('rational_arithmetic_identity')
    require(arithmetic_backend in ('PYTHON_FRACTION_V1','NATIVE_GMP_SURFACE_V1'),'ARITHMETIC_BACKEND_MISMATCH')
    if arithmetic_backend=='NATIVE_GMP_SURFACE_V1':
        require(isinstance(arithmetic_identity,dict) and
            arithmetic_identity.get('backend')==arithmetic_backend and
            arithmetic_identity.get('algorithm')=='DIV_SURFACE_V1_REYNOLDS_V1_GMP_BINARY64_NEAREST_EVEN_V1' and
            bool(re.fullmatch(r'[0-9a-f]{64}',str(arithmetic_identity.get('module_sha256','')))),
            'ARITHMETIC_IDENTITY_MISMATCH')
    else:require(arithmetic_identity is None,'ARITHMETIC_IDENTITY_MISMATCH')
    require(bool(re.fullmatch(r'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:[0-9a-f]{64}',backend)),
            'BACKEND_IDENTITY_MISMATCH')
    require(complete.get('complete') is True and complete['status']==result['status']=='CONVERGED','RUN_NOT_CONVERGED')
    event_path=root/'events.jsonl'
    require(event_path.is_file() and not event_path.is_symlink(),'MISSING_RUN_LOG')
    try:events=[json.loads(line) for line in event_path.read_text().splitlines()]
    except (OSError,json.JSONDecodeError) as exc:raise ValueError('INVALID_RUN_LOG') from exc
    require(bool(events) and events[-1].get('kind')=='RUN_TERMINATED','INCOMPLETE_RUN_LOG')
    require(events[-1].get('status')=='CONVERGED','LOG_STATUS_MISMATCH')
    saved=[e for e in events if e['kind']=='MAP_SAVED'];rejected=[e for e in events if e['kind']=='MAP_REJECTED']
    indices=[e['evaluation'] for e in events if e['kind'] in ('MAP_SAVED','MAP_REJECTED')]
    require(indices==list(range(1,complete['evaluation_count']+1)),'MISSING_MAP_EVENT')
    require(events[-1]['evaluation_count']==complete['evaluation_count'],'EVALUATION_COUNT_MISMATCH')
    newton=[e['event'] for e in events if e['kind']=='NEWTON_EVENT']
    native=result['identity']['native_run']
    require(scientific_projection(newton)==scientific_projection(native['events']),'NATIVE_EVENT_LOG_MISMATCH')
    updates=[e for e in newton if e['type']=='evaluation' and e['outcome']=='accepted' and (e.get('step') or 0)>0]
    require(bool(updates),'NO_ACCEPTED_UPDATE')
    materials={}
    for record in saved:
        name=record['evaluation_file']
        require(name==f"evaluation-{record['evaluation']:06d}.json",'INVALID_MATERIAL_REFERENCE')
        material=read(name)
        require(material.get('rational_arithmetic_identity')==arithmetic_identity,'ARITHMETIC_ATTACHMENT_MISMATCH')
        require(hashlib.sha256((root/name).read_bytes()).hexdigest()==record['sha256'],'MATERIAL_FILE_DIGEST_MISMATCH')
        require(material['quadric_coefficients']==record['coefficients'],'MATERIAL_COEFFICIENT_MISMATCH')
        require(material['geometry_hash']==record['geometry_hash'],'MATERIAL_GEOMETRY_MISMATCH')
        require(material['boolean_backend_identity']==backend and material['evaluation_map_identity']==map_identity,'MATERIAL_CONTEXT_MISMATCH')
        require(material.get('derivative_defined') is True,'DERIVATIVE_UNAVAILABLE')
        require(material['quadric_source_identity']==_digest(dict(vertices=material['quadric_source_vertices'],
            faces=material['quadric_source_faces'],coefficients=material['quadric_coefficients'])),'SOURCE_IDENTITY_MISMATCH')
        require(material['geometry_hash']==_digest(dict(vertices=material['material_proof']['exact_vertices'],
            faces=material['variable_faces'])),'ACTUAL_GEOMETRY_DIGEST_MISMATCH')
        matches=[e for e in newton if e['type']=='evaluation' and e.get('evaluation_index')==record['evaluation']]
        require(len(matches)==1,'MATERIAL_EVALUATION_EVENT_MISSING')
        event=matches[0];payload=event.get('evaluation_payload',{})
        require(event['quadric_coefficients']==material['quadric_coefficients'] and
                payload.get('geometry_hash')==material['geometry_hash'] and
                payload.get('boolean_backend_identity')==backend and
                payload.get('evaluation_map_identity')==map_identity,'EVALUATION_ATTACHMENT_MISMATCH')
        # The sealing bridge independently reintegrates the fixed shell; its
        # target-frame sum can differ by roundoff from the evaluator's sum.
        observed=payload.get('solver_raw_moments',[])
        require(len(observed)==10 and len(material['raw_moments'])==10 and
                all(math.isfinite(a) and math.isfinite(b) and abs(a-b)<=1e-10
                    for a,b in zip(observed,material['raw_moments'])) and
                payload.get('variable_raw_moments')==material['variable_raw_moments'],'MOMENT_ATTACHMENT_MISMATCH')
        J=np.asarray(material['moment_jacobian'],dtype=float)
        solver_J=np.asarray(payload.get('solver_moment_jacobian'),dtype=float)
        s=np.asarray(payload.get('solver_coordinates'),dtype=float)
        require(J.shape==solver_J.shape==(10,10) and s.shape==(10,), 'JACOBIAN_ATTACHMENT_MISMATCH')
        require(event['solver_coordinate_identity']=='literal_monomial_unit_sphere_v1', 'UNSUPPORTED_SOLVER_COORDINATES')
        inverse_scale=np.array([1.]*7+[.5]*3)
        p=s*inverse_scale;length=np.linalg.norm(p)
        require(math.isfinite(length) and length>0,'INVALID_SOLVER_COORDINATES')
        p=p/length
        geometry_J=(np.eye(10)-np.outer(p,p))*inverse_scale[None,:]/length
        expected=J@geometry_J
        require(np.isfinite(expected).all() and np.isfinite(solver_J).all() and
                np.max(np.abs(expected-solver_J))<=1e-10,'JACOBIAN_ATTACHMENT_MISMATCH')
        materials[name]=material
    terminal=read('terminal-material-reference.json');name=terminal['evaluation_file']
    require(name in materials,'TERMINAL_NOT_SAVED')
    material=materials[name]
    require(terminal['sha256']==hashlib.sha256((root/name).read_bytes()).hexdigest(),'TERMINAL_DIGEST_MISMATCH')
    require(terminal['coefficients']==native['quadric_coefficients']==material['quadric_coefficients'],'TERMINAL_COEFFICIENT_MISMATCH')
    require(updates[-1]['quadric_coefficients']==terminal['coefficients'],'TERMINAL_NOT_LAST_ACCEPTED_UPDATE')
    require(native['raw_equation_residual']<=1e-12,'RAW_RESIDUAL_NOT_CONVERGED')
    return dict(record_integrity_valid=True,accepted_updates=len(updates),saved_materials=len(saved),
        rejected_evaluations=len(rejected),scientific_event_digest=scientific_digest(newton),
        raw_log_digest=hashlib.sha256(event_path.read_bytes()).hexdigest(),
        core_admission=False,scope='Run record integrity only; independent material/derivative admission separate')
