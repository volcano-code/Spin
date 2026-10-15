"""Frozen 81 actual CGAL maps; independent integration of saved boundaries."""
from dataclasses import asdict
import hashlib,json,os,sys,time
from pathlib import Path
import numpy as np
from make_it_stand import spin
from independent_tetra_moments import integrate
from probe_inprocess_derivatives import weights
from run_small_inprocess_solver import safe
HERE=Path(__file__).resolve().parent
frozen_path=HERE/'noodles-balanced-request-freeze-job-01/output/request.json'
source_manifest_path=HERE/'noodles-interval-algebraic-source-derivative-job-01/output/manifest.json'
source_gate_path=source_manifest_path.parent/'report.json'
frozen=json.loads(frozen_path.read_text());manifest=json.loads(source_manifest_path.read_text());gate=json.loads(source_gate_path.read_text())
if not gate['source_checks_passed'] or gate['point_count']!=81:raise ValueError('SOURCE_GATE_REQUIRED')
for path,digest in frozen['input_sources'].items():
    if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=digest:raise ValueError('FROZEN_INPUT_CHANGED:'+path)
payload=dict(next(x for x in frozen['requests'] if x['settings']['deterministic_seed']==0))
payload['settings']=spin.SpinSolverConfig(**payload['settings']);request=spin.SpinItFasterRequest(**payload)
config=dict(frozen['evaluator_settings']);config['algebraic_finite_map']=True;settings=spin.SpinEvaluatorConfig(**config)
samples=[json.loads((source_manifest_path.parent/f'point-{i:03d}.json').read_text())['coefficients'] for i in range(81)]
if len({tuple(q) for q in samples})!=81:raise ValueError('UNIQUE_PRODUCTION_POINT_COUNT_MISMATCH')
fixed_path=HERE/'noodles-balanced-request-freeze-job-01/output/fixed.binary64.json'
fixed=json.loads(fixed_path.read_text());rotation=np.asarray(request.target_axes).T;origin=np.asarray(request.target_origin)
fixed_vertices=(np.asarray(fixed['vertices'])-origin)@rotation.T
independent_fixed=integrate(fixed_vertices,fixed['faces'],1.0)
m=request.fixed_shell_moments;mass=m[0];first=np.asarray(m[1:4])
second=np.array([[m[4],m[9],m[8]],[m[9],m[5],m[7]],[m[8],m[7],m[6]]])
first_target=rotation@(first-mass*origin)
second_target=rotation@(second-np.outer(origin,first)-np.outer(first,origin)+mass*np.outer(origin,origin))@rotation.T
expected=np.array([mass,*first_target,second_target[0,0],second_target[1,1],second_target[2,2],second_target[1,2],second_target[0,2],second_target[0,1]])
fixed_error=float(np.max(np.abs(expected-independent_fixed)))
if fixed_error>1e-10:raise ValueError('INDEPENDENT_FIXED_MATERIAL_MISMATCH')
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
with (out/'manifest.json').open('x') as stream:
    json.dump(dict(scope='ACTUAL_CGAL_PRODUCTION_MOMENT_JACOBIAN_CHECK',explicit_profile='PAPER_INTERVAL_ALGEBRAIC_FINITE_CORE_DIAGNOSIS',point_count=81,layouts=manifest['layouts'],samples=samples,settings=config,request_sha256=hashlib.sha256(frozen_path.read_bytes()).hexdigest(),source_manifest_sha256=hashlib.sha256(source_manifest_path.read_bytes()).hexdigest(),source_gate_sha256=hashlib.sha256(source_gate_path.read_bytes()).hexdigest(),native_sha256=hashlib.sha256(Path(spin._native_module().__file__).read_bytes()).hexdigest(),independent_fixed_error=fixed_error,optimization_admission=False),stream)
moments=[];integration_errors=[];jac=None;started=time.monotonic()
with (out/'events.jsonl').open('x') as journal:
    def event(row):journal.write(json.dumps(row,allow_nan=False)+'\n');journal.flush();os.fsync(journal.fileno())
    for i,q in enumerate(samples):
        event(dict(index=i,status='PRODUCTION_POINT_STARTED',seconds=time.monotonic()-started))
        try:
            result=spin.evaluate_spin_it_faster_fullrank_quadric(request,np.asarray(q),evaluator_settings=settings,
                diagnostic_stage_observer=lambda stage:event(dict(index=i,status='STAGE',stage=stage,seconds=time.monotonic()-started)))
            if not result.derivative_defined:raise ValueError('PRODUCTION_DERIVATIVE_UNDEFINED:'+result.derivative_diagnostic)
            path=out/f'point-{i:03d}-material.json'
            with path.open('x') as stream:json.dump(safe(dict(coefficients=q,evaluation=asdict(result))),stream,allow_nan=False)
            saved=json.loads(path.read_text())['evaluation']
            independent=integrate(saved['material_vertices'],saved['material_faces'],request.material_density)+independent_fixed
            error=float(np.max(np.abs(independent-np.asarray(saved['total_raw_moments']))))
            integration_errors.append(error);moments.append(independent)
            if i==0:jac=np.asarray(result.coefficient_jacobian)
            event(dict(status='PRODUCTION_POINT_RETURNED',index=i,seconds=time.monotonic()-started,independent_raw_moments=independent.tolist(),independent_integration_error=error,geometry_hash=result.geometry_hash))
        except Exception as exc:
            event(dict(status='PRODUCTION_POINT_FAILED',index=i,exception_type=type(exc).__name__,exception=str(exc)));raise
checks=[]
for layout in manifest['layouts']:
    numerical=weights(layout['positions'])@np.asarray([moments[i] for i in layout['indices']])
    errors=np.abs(numerical-jac[:,layout['direction']])
    checks.append(dict(**layout,errors=errors.tolist(),maximum_error=float(max(errors)),passed=bool(np.all(errors<=1e-6))))
report=dict(scope='ACTUAL_CGAL_PRODUCTION_MOMENT_JACOBIAN_CHECK',point_count=81,checks=checks,numerical_checks_passed=all(x['passed'] for x in checks),independent_integration_max_error=max(integration_errors),independent_integration_passed=max(integration_errors)<=1e-6,full_material_branch_certification=False,optimization_admission=False,seconds=time.monotonic()-started)
with (out/'report.json').open('x') as stream:json.dump(report,stream,indent=2)
print(json.dumps(dict(numerical_checks_passed=report['numerical_checks_passed'],maximum_error=max(x['maximum_error'] for x in checks))),flush=True)
