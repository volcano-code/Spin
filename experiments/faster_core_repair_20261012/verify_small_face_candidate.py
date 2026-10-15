"""One registered production replay plus audit of already captured source."""
import hashlib,json,sys,time
from pathlib import Path
import numpy as np
from make_it_stand import spin
HERE=Path(__file__).resolve().parent
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
registration=json.loads((HERE/'noodles-first-core-failures-registration-01/registration.json').read_text())
raw=Path(registration['mapping_config']).read_bytes()
if hashlib.sha256(raw).hexdigest()!=registration['mapping_config_sha256']:raise ValueError('STALE_CONFIG')
config=json.loads(raw);payload=dict(config['request']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
native,request=spin._native_request(spin.SpinItFasterRequest(**payload))
settings=spin._native_evaluator_settings(native,spin.SpinEvaluatorConfig(**config['evaluator']))
def convert(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)
def save(name,value):
    with (out/name).open('x') as stream:json.dump(value,stream,default=convert,allow_nan=False)
source=json.loads((HERE/'core-registered-failure-capture-job-01/output/source_side.json').read_text())['snapshots'][0]
audit=native.audit_exact_triangle_boundary(np.asarray(source['source_vertices']),np.asarray(source['source_faces']), 'OutwardPositive')
save('source-side-boundary-audit.json',audit)
start=time.monotonic()
def stage(name):
    with (out/'stages.jsonl').open('a') as stream:
        stream.write(json.dumps(dict(stage=name,monotonic_seconds=time.monotonic()-start))+'\n')
settings.diagnostic_stage_observer=stage
q=np.asarray(registration['failures']['material_boundary']['record']['coefficients'])
report=dict(schema='small-face-fixed-state-verification-v1',production_admission=False,optimization_invocations=0,
            native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest())
try:
    result=native.evaluate_spin_it_faster_fullrank_quadric(request,q,settings,True)
    save('actual-evaluation.json',result)
    report.update(evaluation_completed=True,derivative_defined=result['derivative_defined'],
                  diagnostic=result['derivative_diagnostic'],material_faces=len(result['material_faces']),
                  geometry_hash=result['geometry_hash'])
except Exception as error:
    report.update(evaluation_completed=False,error=str(error),error_type=type(error).__name__)
report.update(wall_seconds=time.monotonic()-start,source_boundary_valid=audit['boundary_valid'])
save('report.json',report)
print(json.dumps(report),flush=True)
