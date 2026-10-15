"""Input-defined 81-point Faster map registration and one actual point job."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'python'))


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def register(request_file,out,native_file=None):
    import numpy as np
    from make_it_stand.verification.surface_derivatives import quadric_five_point_layout
    frozen=json.loads(Path(request_file).read_text())
    for path,digest in frozen['input_sources'].items():
        if sha(path)!=digest:raise ValueError('FROZEN_REQUEST_SOURCE_CHANGED:'+path)
    vertices=np.asarray(frozen['requests'][0]['design_vertices'])
    low,high=vertices.min(axis=0),vertices.max(axis=0)
    center=(low+high)/2
    radii=.7*(high-low)/2
    diagonal=np.array([1.,1.1,1.2])/radii**2
    base=np.r_[np.dot(diagonal,center**2)-1,-2*diagonal*center,diagonal,[0.,0.,0.]]
    steps=.001*np.maximum(1.,abs(base))
    trials,rows=quadric_five_point_layout(base,steps)
    native=Path(native_file).resolve() if native_file else ROOT/'build/spin2014-objectives-quadric-20261007/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    sources={str(Path(request_file).resolve()):sha(request_file),str(native):sha(native),
        str(Path(__file__).resolve()):sha(__file__),
        str(ROOT/'python/make_it_stand/spin.py'):sha(ROOT/'python/make_it_stand/spin.py'),
        str(Path(__file__).with_name('rounded_input.py')):sha(Path(__file__).with_name('rounded_input.py')),
        str(ROOT/'python/make_it_stand/verification/surface_derivatives.py'):sha(ROOT/'python/make_it_stand/verification/surface_derivatives.py')}
    payload=dict(schema='rounded-heart-faster-81-point-registration-v1',request_path=str(Path(request_file).resolve()),native_path=str(native),
        sources=sources,base=base.tolist(),steps=steps.tolist(),
        points=[list(x) for x in trials],rows=[dict(column=c,factor=f,indices=list(ids)) for c,f,ids in rows],
        initial_state_rule='INPUT_DESIGN_BBOX_70_PERCENT_RADII_DISTINCT_AXIS_FACTORS_1_1.1_1.2',
        step_rule='0.001_TIMES_MAX_1_ABS_RAW_COEFFICIENT',
        scientific_tolerance=1e-6,production_queries=81,derivative_admission=False,optimizer_invocations=0)
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    (out/'registration.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')


def point(registration,index,out):
    import gmpy2
    import numpy as np
    from make_it_stand.spin import SpinItFasterRequest,SpinSolverConfig,SpinEvaluatorConfig,evaluate_spin_it_faster_fullrank_quadric
    from make_it_stand import spin
    from rounded_input import exact_moments
    from make_it_stand.verification.spin_problem import MeshSnapshot
    registration=Path(registration);out=Path(out)
    registered=json.loads(registration.read_text())
    if not 0<=index<81:raise ValueError('POINT_INDEX_OUTSIDE_REGISTRATION')
    for path,digest in registered['sources'].items():
        if sha(path)!=digest:raise ValueError('REGISTERED_PRODUCTION_SOURCE_CHANGED:'+path)
    if Path(spin._native_module().__file__).resolve()!=Path(registered['native_path']).resolve():
        raise ValueError('LOADED_NATIVE_DIFFERS_FROM_REGISTRATION')
    frozen=json.loads(Path(registered['request_path']).read_text())
    for path,digest in frozen['input_sources'].items():
        if sha(path)!=digest:raise ValueError('GEOMETRY_BINDING_CHANGED:'+path)
    payload=dict(frozen['requests'][0]);payload['settings']=SpinSolverConfig(**payload['settings'])
    request=SpinItFasterRequest(**payload)
    settings=SpinEvaluatorConfig(**frozen['evaluator_settings'])
    out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    with (out/'events.jsonl').open('x') as stream:
        import os
        def stage(name):
            row=dict(stage=name,index=index,elapsed_seconds=time.monotonic()-started)
            stream.write(json.dumps(row)+'\n');stream.flush();os.fsync(stream.fileno());print(json.dumps(row),flush=True)
        stage('PRODUCTION_EVALUATION_STARTED')
        evaluated=evaluate_spin_it_faster_fullrank_quadric(request,registered['points'][index],evaluator_settings=settings,diagnostic_stage_observer=stage)
        packet=asdict(evaluated)
        (out/'actual-evaluation.json').write_text(json.dumps(packet,allow_nan=False)+'\n')
        authority=evaluated.rational_authority
        stage('INDEPENDENT_RATIONAL_MATERIAL_INTEGRATION')
        variable=exact_moments(authority['material_vertices'],evaluated.material_faces)
        total=[m*gmpy2.mpq(request.material_density)+gmpy2.mpq(x) for m,x in zip(variable,request.fixed_shell_moments)]
        stage('INDEPENDENT_MATERIAL_ENTITY')
        facts=MeshSnapshot.from_rationals(authority['material_vertices'],evaluated.material_faces).exact_region().facts
        errors=[float(abs(a-gmpy2.mpq(x))/max(gmpy2.mpq(1),abs(a))) for a,x in zip(total,evaluated.total_raw_moments)]
        record=dict(schema='rounded-heart-faster-production-point-v1',index=index,
            registration_sha256=sha(registration),raw_coefficients=registered['points'][index],
            normalized_coefficients=list(evaluated.normalized_coefficients),
            independent_total_moments=[str(x) for x in total],moment_errors=errors,
            derivative_defined=evaluated.derivative_defined,derivative_diagnostic=evaluated.derivative_diagnostic,
            branch_identity=evaluated.branch_identity,material_facts=facts,
            evaluation_sha256=sha(out/'actual-evaluation.json'),
            elapsed_seconds=time.monotonic()-started,derivative_admission=False,optimizer_invocations=0)
        (out/'point.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
        stage('POINT_COMPLETE_NO_BATCH_VERDICT')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('operation',choices=['register','point','launch'])
    parser.add_argument('--request',type=Path);parser.add_argument('--registration',type=Path)
    parser.add_argument('--index',type=int,default=0);parser.add_argument('--output',type=Path);parser.add_argument('--job',type=Path)
    parser.add_argument('--native',type=Path)
    args=parser.parse_args()
    if args.operation=='register':register(args.request,args.output,args.native);return
    if args.operation=='point':point(args.registration,args.index,args.output);return
    from make_it_stand.verification.source_job import run_source_job
    registration=args.registration.resolve();data=json.loads(registration.read_text())
    frozen=json.loads(Path(data['request_path']).read_text())
    files={**{str(i):Path(p) for i,p in enumerate(data['sources'])},
           **{'input.'+str(i):Path(p) for i,p in enumerate(frozen['input_sources'])},'registration':registration}
    native=Path(data['native_path']).parent
    result=run_source_job(args.job.resolve(),['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
        '/home/lin/anaconda3/bin/python3','-B',str(Path(__file__).resolve()),'point',
        '--registration',str(registration),'--index',str(args.index),'--output','{job}/output'],files,
        environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native)})
    print(json.dumps(dict(status=result.status,returncode=result.returncode)),flush=True)
    raise SystemExit(0 if result.status=='PROCESS_EXITED' and result.returncode==0 else 1)


if __name__=='__main__':main()
