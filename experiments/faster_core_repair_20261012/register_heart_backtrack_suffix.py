"""Preregister the next four untested original half-step trials at state45."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from audit_heart_stalled_prefix import coordinate


def main(output):
    here=Path(__file__).resolve().parent
    prior_path=here/'heart-actual-direction-probe-registration-03.json'
    prior=json.loads(prior_path.read_text())
    trial=next(p for p in prior['probes'] if p['base_evaluation']==45 and p['existing_trial'])
    run=here/'heart-guarded-formal-seed-0-job-01/output/run'
    base_path=run/'evaluation-000045.json'
    if hashlib.sha256(base_path.read_bytes()).hexdigest()!=trial['base_material_sha256']:
        raise ValueError('BASE_MATERIAL_CHANGED')
    p=coordinate(json.loads(base_path.read_text())['quadric_coefficients'])
    t=coordinate(trial['coefficients']);v=t-(p@t)*p;v/=np.linalg.norm(v)
    probes=[]
    steps = ((.00390625,.001953125,.000244140625,.000030517578125,.000003814697265625)
             if '--local-suffix' in sys.argv else (.0625,.03125,.015625,.0078125))
    for step in steps:
        target=p*np.cos(step)+v*np.sin(step);q=target.copy();q[7:]*=.5;q/=np.linalg.norm(q)
        if np.linalg.norm(coordinate(q)-target)>1e-12:raise ValueError('GEODESIC_CHANGED')
        probes.append(dict(trial,step=step,fraction=step/.5,coefficients=q.tolist(),existing_trial=False))
    registration=dict(prior,probes=probes,
        prior_registration_sha256=hashlib.sha256(prior_path.read_bytes()).hexdigest(),
        scope='Preregistered paper half-step subsequence on saved state45 direction; diagnostic only, not a Newton restart')
    with output.open('x') as stream:json.dump(registration,stream,indent=2,allow_nan=False)
    print('REGISTERED_STEPS',[p['step'] for p in probes])


if __name__=='__main__':main(Path(sys.argv[1]))
