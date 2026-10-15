"""Diagnose the failed registered stencil without new evaluations or admission."""
import hashlib
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
p=HERE/'formal-cherchi-derivative-job-01/output'
r=json.loads((p/'report.json').read_text())
rows=[json.loads(s) for s in (p/'events.jsonl').read_text().splitlines()]
m=json.loads((p/'manifest.json').read_text())
base=m['base'][0]
indices=sorted({i for layout in m['layouts'] if layout['direction']==0 for i in layout['indices']},
               key=lambda i:m['samples'][i][0])
nodes=np.array([m['samples'][i][0]-base for i in indices]);scale=max(abs(nodes))
A=np.vstack([(nodes/scale)**k for k in range(len(nodes))])
rhs=np.zeros(len(nodes));rhs[1]=1
unit=np.linalg.solve(A,rhs)
if max(abs(A@unit-rhs))>1e-10:raise ValueError('DIAGNOSTIC_WEIGHTS_INVALID')
estimate=(unit/scale)@np.array([rows[i]['raw_moments'] for i in indices])
observed=[x['maximum_error'] for x in r['checks'] if x['direction']==0]
record=dict(role='SAVED_DATA_DIAGNOSIS_ONLY',original_registered_gate_passed=r['passed'],
    original_five_point_errors=observed,error_reduction_ratios=[observed[i]/observed[i+1] for i in range(2)],
    nine_node_diagnostic_max_error=float(max(abs(estimate-np.array(r['production_jacobian'])[:,0]))),
    new_evaluations=0,grants_derivative_admission=False,
    interpretation='Consistent with fourth-order stencil truncation; does not override original three-scale failure.',
    input_sha256={name:hashlib.sha256((p/name).read_bytes()).hexdigest() for name in ('manifest.json','events.jsonl','report.json')})
with (HERE/'formal-cherchi-fd-diagnosis-01.json').open('x') as stream:json.dump(record,stream,indent=2)
print(json.dumps(record))
