"""81 source points, fixed failed seed state; not a material Jacobian gate."""
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import numpy as np
import _make_it_stand as native
from probe_inprocess_derivatives import weights
from run_small_inprocess_solver import safe
HERE=Path(__file__).resolve().parent
run=HERE/'noodles-bound-seed0-production-probe-job-01/output/run'
config_path=run/'mapping-config.json';config=json.loads(config_path.read_text())
request=config['request'];evaluator=config['evaluator']
first=json.loads((run/'mapping-events.jsonl').read_text().splitlines()[0]);base=np.asarray(first['coefficients'])
design=(np.asarray(request['design_vertices'])-request['target_origin'])@np.asarray(request['target_axes'])
settings=native.SpinTriangulationSettings()
for name in ('icosphere_subdivisions','maximum_icosphere_subdivisions','maximum_radius_doublings',
    'retain_fixed_normal_provenance','truncation_radius','paper_projective_hyperboloid','algebraic_finite_map'):
    setattr(settings,name,evaluator[name])
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
samples=[base.copy()];index={tuple(base):0};layouts=[];scales=[3e-6,1.5e-6,.75e-6]
for direction in range(10):
    for H in scales:
        ids=[];positions=[]
        for shift in (-H,-H/2,0,H/2,H):
            q=base.copy();q[direction]+=shift;key=tuple(q)
            if key not in index:index[key]=len(samples);samples.append(q)
            ids.append(index[key]);positions.append(float(q[direction]-base[direction]))
        layouts.append(dict(direction=direction,scale=H,indices=ids,positions=positions))
if len(samples)!=81:raise ValueError('SOURCE_POINT_COUNT_MISMATCH')
with (output/'manifest.json').open('x') as stream:
    json.dump(dict(scope='SOURCE_COORDINATE_DERIVATIVE_ONLY_NOT_MATERIAL_MOMENTS',
        base=base.tolist(),scales=scales,layouts=layouts,
        native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
        input_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest()),stream)
records=[];started=time.monotonic()
with (output/'events.jsonl').open('x') as stream:
    def event(row):stream.write(json.dumps(safe(row),allow_nan=False)+'\n');stream.flush();os.fsync(stream.fileno())
    for i,q in enumerate(samples):
        try:
            result=native.evaluate_spin_it_faster_source_motion(q,design,settings)
            if not result['defined'] or not result['side_selection_defined']:raise ValueError(result['diagnostic'])
            with (output/f'point-{i:03d}.json').open('x') as saved:
                json.dump(safe(dict(coefficients=q,source=result)),saved,allow_nan=False)
            records.append(result);event(dict(index=i,status='SOURCE_RETURNED',seconds=time.monotonic()-started))
        except Exception as exc:event(dict(index=i,status='SOURCE_FAILED',exception=str(exc)));raise
centre=records[0];jac=np.asarray(centre['vertex_velocity']).reshape(len(centre['vertices']),3,10)
topology=True
for result in records:
    topology&=result['truncation_radius']==centre['truncation_radius']
    topology&=result['source_projective_rays']==centre['source_projective_rays']
    topology&=result['source_normal_vertices']==centre['source_normal_vertices']
    topology&=np.array_equal(result['faces'],centre['faces'])
checks=[]
for layout in layouts:
    numerical=np.tensordot(weights(layout['positions']),np.asarray([records[i]['vertices'] for i in layout['indices']]),axes=1)
    error=np.linalg.norm(numerical-jac[:,:,layout['direction']],axis=1)
    denominator=np.maximum(1.,np.linalg.norm(jac[:,:,layout['direction']],axis=1))
    checks.append(dict(**layout,maximum_relative_error=float(np.max(error/denominator)),passed=bool(np.all(error/denominator<=1e-6))))
audit=dict(native.audit_exact_triangle_boundary(centre['vertices'],centre['faces'],'OutwardPositive'))
report=dict(scope='SOURCE_COORDINATE_DERIVATIVE_ONLY_NOT_MATERIAL_MOMENTS',point_count=81,
    source_boundary_audit=audit,discrete_branch_stable=topology,checks=checks,
    source_checks_passed=topology and audit['boundary_valid'] and all(row['passed'] for row in checks),
    production_moment_admission=False,optimization_admission=False,seconds=time.monotonic()-started)
with (output/'report.json').open('x') as stream:json.dump(report,stream,indent=2)
print(json.dumps(dict(source_checks_passed=report['source_checks_passed'],branch_stable=topology,
    maximum_relative_error=max(row['maximum_relative_error'] for row in checks),optimization_admission=False)),flush=True)
