"""New explicit algebraic-finite profile, same 81 source-only positions."""
import json,sys,os,time,hashlib
from pathlib import Path
import numpy as np
import _make_it_stand as native
from probe_inprocess_derivatives import weights
from run_small_inprocess_solver import safe
HERE=Path(__file__).resolve().parent;saved=HERE/'noodles-interval-source-derivative-job-01/output'
config=json.loads((HERE/'noodles-bound-seed0-production-probe-job-01/output/run/mapping-config.json').read_text())
r=config['request'];e=config['evaluator'];design=(np.asarray(r['design_vertices'])-r['target_origin'])@np.asarray(r['target_axes'])
s=native.SpinTriangulationSettings()
for k in ('icosphere_subdivisions','maximum_icosphere_subdivisions','maximum_radius_doublings','retain_fixed_normal_provenance','truncation_radius','paper_projective_hyperboloid','algebraic_finite_map'):setattr(s,k,e[k])
s.algebraic_finite_map=True
manifest=json.loads((saved/'manifest.json').read_text());out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
with (out/'manifest.json').open('x') as f:json.dump(dict(**manifest,explicit_profile='PAPER_INTERVAL_ALGEBRAIC_FINITE_SOURCE_DIAGNOSIS',algebraic_finite_map=True,parent_manifest_sha256=hashlib.sha256((saved/'manifest.json').read_bytes()).hexdigest()),f)
records=[];start=time.monotonic()
with (out/'events.jsonl').open('x') as journal:
 for i in range(81):
    try:
        q=json.loads((saved/f'point-{i:03d}.json').read_text())['coefficients']
        source=native.evaluate_spin_it_faster_source_motion(np.asarray(q),design,s)
        with (out/f'point-{i:03d}.json').open('x') as f:json.dump(safe(dict(coefficients=q,source=source)),f,allow_nan=False)
        if not source['defined'] or not source['side_selection_defined']:raise ValueError(source['diagnostic'])
        records.append(source);status='SOURCE_RETURNED'
    except Exception as exc:
        journal.write(json.dumps(dict(index=i,status='SOURCE_FAILED',exception=str(exc)))+'\n');journal.flush();os.fsync(journal.fileno());raise
    journal.write(json.dumps(dict(index=i,status=status,seconds=time.monotonic()-start))+'\n');journal.flush();os.fsync(journal.fileno())
b=records[0];j=np.asarray(b['vertex_velocity']).reshape(-1,3,10);checks=[]
stable=all(x['truncation_radius']==b['truncation_radius'] and x['source_projective_rays']==b['source_projective_rays'] and x['source_normal_vertices']==b['source_normal_vertices'] and np.array_equal(x['faces'],b['faces']) for x in records)
for x in manifest['layouts']:
 n=np.tensordot(weights(x['positions']),np.asarray([records[i]['vertices'] for i in x['indices']]),axes=1)
 error=np.linalg.norm(n-j[:,:,x['direction']],axis=1)/np.maximum(1,np.linalg.norm(j[:,:,x['direction']],axis=1))
 checks.append(dict(**x,maximum_relative_error=float(max(error)),worst_vertex=int(np.argmax(error)),passed=bool(np.all(error<=1e-6))))
audit=dict(native.audit_exact_triangle_boundary(b['vertices'],b['faces'],'OutwardPositive'))
report=dict(scope='SOURCE_COORDINATE_DERIVATIVE_ONLY_NOT_MATERIAL_MOMENTS',explicit_profile='PAPER_INTERVAL_ALGEBRAIC_FINITE_SOURCE_DIAGNOSIS',checks=checks,point_count=81,discrete_branch_stable=stable,source_boundary_audit=audit,source_checks_passed=stable and audit['boundary_valid'] and all(x['passed'] for x in checks),native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),seconds=time.monotonic()-start,production_moment_admission=False,optimization_admission=False)
with (out/'report.json').open('x') as f:json.dump(report,f,indent=2)
print(json.dumps(dict(source_checks_passed=report['source_checks_passed'],maximum_relative_error=max(x['maximum_relative_error'] for x in checks))),flush=True)
