"""New registered 81-point actual CGAL map; saved geometry integrated independently."""
import argparse
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import numpy as np
from make_it_stand import spin
from independent_tetra_moments import integrate
from probe_inprocess_derivatives import weights
from run_small_inprocess_solver import safe

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]


def main(output,state_name,material_policy='legacy'):
    heart=ROOT/'experiments/faster_heart_contact_20261011'
    sys.path.insert(0,str(heart))
    from core_faster_heart import core_settings
    frozen_path=heart/'rounded-clearance-v3-faster-frozen-v2/request.json'
    frozen=json.loads(frozen_path.read_text())
    for path,digest in frozen['input_sources'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=digest:raise ValueError('FROZEN_INPUT_CHANGED')
    payload=dict(frozen['requests'][0]);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
    request=spin.SpinItFasterRequest(**payload)
    config=core_settings(frozen['evaluator_settings']);config['paper_projective_hyperboloid']=True
    if material_policy == 'exact':
        config['exact_binary64_material_boundary']=True
    settings=spin.SpinEvaluatorConfig(**config)
    state=json.loads((HERE/'frozen-heart-states.json').read_text())['records'][state_name]
    base=np.asarray(state['coefficients']);samples=[base.copy()];lookup={tuple(base):0};layouts=[]
    for direction in range(10):
        for H in (3e-6,1.5e-6,.75e-6):
            ids=[];positions=[]
            for shift in (-H,-H/2,0,H/2,H):
                q=base.copy();q[direction]+=shift;key=tuple(q)
                if key not in lookup:lookup[key]=len(samples);samples.append(q)
                ids.append(lookup[key]);positions.append(float(q[direction]-base[direction]))
            layouts.append(dict(direction=direction,scale=H,indices=ids,positions=positions))
    if len(samples)!=81:raise ValueError('UNIQUE_PRODUCTION_POINT_COUNT_MISMATCH')
    fixed_path=next(Path(path) for path in frozen['input_sources'] if Path(path).name=='fixed.binary64.json')
    fixed=json.loads(fixed_path.read_text());rotation=np.asarray(request.target_axes).T
    fixed_vertices=(np.asarray(fixed['vertices'])-request.target_origin)@rotation.T
    independent_fixed=integrate(fixed_vertices,fixed['faces'],1.0)
    # The request is already normalized with fixed density 1; bind and check it.
    mass=request.fixed_shell_moments[0];first=np.asarray(request.fixed_shell_moments[1:4]);origin=np.asarray(request.target_origin)
    second=np.array(((request.fixed_shell_moments[4],request.fixed_shell_moments[9],request.fixed_shell_moments[8]),
                     (request.fixed_shell_moments[9],request.fixed_shell_moments[5],request.fixed_shell_moments[7]),
                     (request.fixed_shell_moments[8],request.fixed_shell_moments[7],request.fixed_shell_moments[6])))
    first_target=rotation@(first-mass*origin)
    second_target=rotation@(second-np.outer(origin,first)-np.outer(first,origin)+mass*np.outer(origin,origin))@rotation.T
    expected_fixed=np.array((mass,*first_target,second_target[0,0],second_target[1,1],second_target[2,2],
                            second_target[1,2],second_target[0,2],second_target[0,1]))
    fixed_error=float(np.max(np.abs(independent_fixed-expected_fixed)))
    if fixed_error>1e-10:raise ValueError('INDEPENDENT_FIXED_MATERIAL_MISMATCH')
    output.mkdir(parents=True,exist_ok=False)
    with (output/'manifest.json').open('x') as stream:
        json.dump(dict(scope='ACTUAL_CGAL_PRODUCTION_MOMENT_JACOBIAN_CHECK',state=state_name,
            base_step=3e-6,point_count=81,layouts=layouts,samples=[q.tolist() for q in samples],
            settings=config,request_sha256=hashlib.sha256(frozen_path.read_bytes()).hexdigest(),
            native_sha256=hashlib.sha256(Path(spin._native_module().__file__).read_bytes()).hexdigest(),
            independent_fixed_error=fixed_error,independent_fixed_density=1.0),stream)
    moments=[];integration_errors=[];jac=None;started=time.monotonic()
    with (output/'events.jsonl').open('x') as events:
        for i,q in enumerate(samples):
            try:
                result=spin.evaluate_spin_it_faster_fullrank_quadric(request,q,evaluator_settings=settings)
                if not result.derivative_defined:raise ValueError('PRODUCTION_DERIVATIVE_UNDEFINED:'+result.derivative_diagnostic)
                record=dict(coefficients=q.tolist(),evaluation=asdict(result))
                path=output/f'point-{i:03d}-material.json'
                with path.open('x') as stream:json.dump(safe(record),stream,allow_nan=False)
                # Re-read persisted geometry; never integrate an optimizer-provided moment vector.
                saved=json.loads(path.read_text())['evaluation']
                independent=integrate(saved['material_vertices'],saved['material_faces'],request.material_density)+independent_fixed
                integration_errors.append(float(np.max(np.abs(independent-np.asarray(result.total_raw_moments)))))
                moments.append(independent)
                if i==0:jac=np.asarray(result.coefficient_jacobian)
                row=dict(status='PRODUCTION_POINT_RETURNED',index=i,seconds=time.monotonic()-started,
                    coefficients=q.tolist(),independent_raw_moments=independent.tolist(),
                    independent_integration_error=integration_errors[-1],geometry_hash=result.geometry_hash)
            except Exception as exc:
                row=dict(status='PRODUCTION_POINT_FAILED',index=i,exception_type=type(exc).__name__,exception=str(exc))
                events.write(json.dumps(row)+'\n');events.flush();os.fsync(events.fileno());raise
            events.write(json.dumps(row)+'\n');events.flush();os.fsync(events.fileno())
    checks=[]
    for layout in layouts:
        numerical=weights(layout['positions'])@np.asarray([moments[i] for i in layout['indices']])
        errors=np.abs(numerical-jac[:,layout['direction']])
        checks.append(dict(**layout,errors=errors.tolist(),maximum_error=float(np.max(errors)),
                           passed=bool(np.all(errors<=1e-6))))
    with (output/'report.json').open('x') as stream:
        json.dump(dict(state=state_name,point_count=81,checks=checks,
            numerical_checks_passed=all(x['passed'] for x in checks),
            independent_integration_max_error=max(integration_errors),
            independent_integration_passed=max(integration_errors)<=1e-6,
            full_material_branch_certification=False,optimization_admission=False),stream,indent=2)
    print(json.dumps(dict(state=state_name,numerical_checks_passed=all(x['passed'] for x in checks),
        maximum_error=max(x['maximum_error'] for x in checks))),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--state',choices=('heart-first-cap-failure','heart-first-small-face-rejection'),required=True)
    p.add_argument('--material-policy',choices=('legacy','exact'),default='legacy')
    a=p.parse_args();main(a.output,a.state,a.material_policy)
