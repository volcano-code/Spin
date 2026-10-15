"""Replay only the two registered states through the actual production evaluator."""
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from make_it_stand import spin

HERE=Path(__file__).resolve().parent
registration_path=HERE/'noodles-first-core-failures-registration-01/registration.json'
registration=json.loads(registration_path.read_text())
config_path=Path(registration['mapping_config'])
if hashlib.sha256(config_path.read_bytes()).hexdigest()!=registration['mapping_config_sha256']:
    raise ValueError('REGISTERED_CONFIG_CHANGED')
config=json.loads(config_path.read_text())
payload=dict(config['request']);payload['settings']=spin.SpinSolverConfig(**payload['settings'])
native,request=spin._native_request(spin.SpinItFasterRequest(**payload))
settings=spin._native_evaluator_settings(native,spin.SpinEvaluatorConfig(**config['evaluator']))
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)

def serializable(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)

summaries=[]
for name in ('source_side','material_boundary'):
    record=registration['failures'][name]['record']
    start=time.monotonic()
    result=native.capture_spin_it_faster_fullrank_geometry(request,np.asarray(record['coefficients']),settings)
    result.update(registration_sha256=hashlib.sha256(registration_path.read_bytes()).hexdigest(),
                  original_evaluation_index=record['index'],original_error=record['exception_message'],
                  native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
                  wall_seconds=time.monotonic()-start,optimization_invocations=0)
    path=out/(name+'.json')
    with path.open('x') as stream:json.dump(result,stream,default=serializable,allow_nan=False)
    summary=dict(name=name,evaluation_completed=result['evaluation_completed'],error=result['error'],
                 original_error_matches=result['error']==record['exception_message'],
                 stages=[s['stage'] for s in result['snapshots']],seconds=result['wall_seconds'],
                 sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    summaries.append(summary)
    print(json.dumps(summary),flush=True)
    if result['evaluation_completed'] or not summary['original_error_matches']:
        raise ValueError('REGISTERED_FAILURE_NOT_REPRODUCED:'+name)
with (out/'summary.json').open('x') as stream:
    json.dump(dict(schema='registered-same-call-failure-capture-v1',production_admission=False,
                   optimization_invocations=0,states=summaries),stream,indent=2)
