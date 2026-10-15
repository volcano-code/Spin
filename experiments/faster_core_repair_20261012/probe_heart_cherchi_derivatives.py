"""Registered Heart Cherchi finite differences; no optimization admission."""
from fractions import Fraction
import json
import os
from pathlib import Path
import sys
import time
import numpy as np
from make_it_stand import spin
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator,_json_value
from independent_tetra_moments import integrate
from probe_formal_cherchi_derivatives import weights
HERE=Path(__file__).resolve().parent

def main(output,state):
    if state not in ('heart-first-cap-failure','heart-first-small-face-rejection'):
        raise ValueError('UNREGISTERED_HEART_STATE')
    prior=HERE/'heart-cherchi-coordinate-map-job-01/output'
    if not json.loads((prior/'report.json').read_text())['passed']:raise ValueError('COORDINATE_GATE_FAILED')
    config=json.loads((prior/'registration.json').read_text())
    payload=dict(config['normalized']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
    request=spin.SpinItFasterRequest(**payload);settings=spin.SpinEvaluatorConfig(**config['settings'])
    base=np.asarray(json.loads((prior/f'{state}-normalized.json').read_text())['quadric_coefficients'])
    samples=[base.copy()];lookup={tuple(base):0};layouts=[]
    for direction in range(10):
        for h in (3e-6,1.5e-6,.75e-6):
            indices=[];positions=[]
            for shift in (-h,-h/2,0.,h/2,h):
                q=base.copy();q[direction]+=shift;key=tuple(q)
                if key not in lookup:lookup[key]=len(samples);samples.append(q)
                indices.append(lookup[key]);positions.append(float(q[direction]-base[direction]))
            layouts.append(dict(direction=direction,scale=h,indices=indices,positions=positions))
    if len(samples)!=81:raise ValueError('POINT_COUNT_MISMATCH')
    output.mkdir(parents=True,exist_ok=False)
    (output/'manifest.json').write_text(json.dumps(dict(state=state,point_count=81,
        samples=[q.tolist() for q in samples],layouts=layouts,coordinate_frame='NORMALIZED_TARGET',
        request=config['normalized'],settings=config['settings'],case_admission=False)))
    mapping=InProcessCherchiEvaluator(request,settings,
        normalization=dict(origin=config['origin'],L=config['length']),
        design_vertices_physical=config['original']['design_vertices'],
        rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1')
    frozen=json.loads((HERE.parents[1]/'experiments/faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json').read_text())
    fixed_path=next(Path(p) for p in frozen['input_sources'] if Path(p).name=='fixed.binary64.json')
    fixed=json.loads(fixed_path.read_text())
    exact_fixed=tuple(tuple(Fraction(float(x)) for x in row) for row in fixed['vertices'])
    target_fixed=mapping._target_view(exact_fixed)[1]
    independent_fixed=integrate(target_fixed,fixed['faces'],1.)
    if np.max(np.abs(independent_fixed-mapping.fixed))>1e-10:raise ValueError('FIXED_MOMENTS_MISMATCH')
    moments=[];errors=[];jac=None
    with (output/'events.jsonl').open('x') as events:
        for i,q in enumerate(samples):
            start=time.monotonic()
            try:
                mapped=mapping(q);material=mapping.material(q)
                path=output/f'point-{i:03d}-material.json'
                path.write_text(json.dumps(material,default=_json_value,allow_nan=False))
                saved=json.loads(path.read_text())
                exact=[[Fraction(x) for x in row] for row in saved['variable_rational_vertices']]
                independent=integrate(exact,saved['variable_faces'],request.material_density)+independent_fixed
                error=float(np.max(np.abs(independent-np.asarray(mapped['raw_moments']))))
                if error>1e-10:raise ValueError('INDEPENDENT_INTEGRATION_MISMATCH')
                errors.append(error);moments.append(independent)
                if i==0:jac=np.asarray(mapped['moment_jacobian'])
                row=dict(index=i,status='RETURNED',seconds=time.monotonic()-start,
                    topology_hash=mapped['topology_hash'],independent_moments=independent.tolist(),integration_error=error)
            except Exception as exc:
                events.write(json.dumps(dict(index=i,status='FAILED',error=f'{type(exc).__name__}: {exc}'))+'\n')
                events.flush();os.fsync(events.fileno());raise
            events.write(json.dumps(row)+'\n');events.flush();os.fsync(events.fileno())
    checks=[]
    for layout in layouts:
        numerical=weights(layout['positions'])@np.asarray([moments[i] for i in layout['indices']])
        error=float(np.max(np.abs(numerical-jac[:,layout['direction']])))
        checks.append(dict(**layout,error=error,passed=error<=1e-6))
    report=dict(state=state,checks=checks,numerical_passed=all(c['passed'] for c in checks),
        independent_error=max(errors),point_count=len(samples),branch_certification=False,case_admission=False)
    (output/'report.json').write_text(json.dumps(report))
    if not report['numerical_passed']:raise ValueError('PRODUCTION_FD_FAILED_NO_STEP_RETUNING')

if __name__=='__main__':main(Path(sys.argv[1]),sys.argv[2])
