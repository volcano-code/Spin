"""Independent high-precision coordinate-chain diagnostic; does not waive gate."""
import json
from pathlib import Path
import numpy as np
import mpmath as mp
HERE=Path(__file__).resolve().parent
root=HERE/'spin2014-inertia-actual-source-guard-seed1-job-04/output/run'
events=[json.loads(x).get('event',{}) for x in (root/'events.jsonl').read_text().splitlines()]
event=next(x for x in events if x.get('type')=='evaluation' and x.get('evaluation_index')==38)
material=json.loads((root/'evaluation-000038.json').read_text());payload=event['evaluation_payload']
mp.mp.dps=100
J=mp.matrix(material['moment_jacobian']);s=mp.matrix(payload['solver_coordinates']);D=mp.diag([1]*7+[mp.mpf('.5')]*3)
p=D*s;length=mp.norm(p);p/=length;exact=J*(mp.eye(10)-p*p.T)*D/length
solver=np.array(payload['solver_moment_jacobian']);j=np.array(material['moment_jacobian']);q=np.array(payload['solver_coordinates'])*np.array([1]*7+[.5]*3);l=np.linalg.norm(q);q/=l
numpy=j@((np.eye(10)-np.outer(q,q))*np.array([1]*7+[.5]*3)[None,:]/l)
native_error=max(abs(mp.mpf(float(solver[a,b]))-exact[a,b]) for a in range(10) for b in range(10))
numpy_error=max(abs(mp.mpf(float(numpy[a,b]))-exact[a,b]) for a in range(10) for b in range(10))
r=dict(evaluation=38,precision_digits=100,max_native_magnitude=float(abs(solver).max()),native_reference_absolute_error=str(native_error),numpy_reference_absolute_error=str(numpy_error),native_numpy_absolute_difference=float(abs(solver-numpy).max()),attachment_absolute_threshold=1e-10,attachment_passed=False,action='Preserve rejection; require verified native arithmetic replay or revised explicitly approved attachment contract. No derivative threshold changed.')
(HERE/'jacobian-attachment38-reference-01.json').write_text(json.dumps(r,indent=2));print(json.dumps(r),flush=True)
