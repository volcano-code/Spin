"""Register half/quarter probes on saved actual geodesics; never initialize Newton."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from audit_heart_stalled_prefix import coordinate


def register(source):
    payload = (source/'events.jsonl').read_bytes()
    records = [json.loads(line) for line in payload.splitlines()]
    maps = {r['evaluation']: r for r in records if r.get('kind') == 'MAP_SAVED'}
    for r in maps.values():
        if hashlib.sha256((source/r['evaluation_file']).read_bytes()).hexdigest() != r['sha256']:
            raise ValueError('MATERIAL_CHECKPOINT_CHANGED')
    events = [r['event'] for r in records if r.get('kind') == 'NEWTON_EVENT'
              and r['event'].get('type') == 'evaluation' and r['event'].get('evaluation_payload')]
    probes = []
    for base_index in (38, 45):
        base = next(e for e in events if e['evaluation_index'] == base_index)
        if base['outcome'] != 'accepted':
            raise ValueError('BASE_NOT_ACTUALLY_ACCEPTED')
        trial = next(e for e in events if e['evaluation_index'] > base_index
                     and e['restart'] == base['restart'] and (e.get('step') or 0) > 0
                     and (base_index == 45 or e['outcome'] == 'accepted'))
        p, t = coordinate(base['quadric_coefficients']), coordinate(trial['quadric_coefficients'])
        cosine = float(p@t)
        tangent = t-cosine*p
        sine = np.linalg.norm(tangent)
        arc = float(np.arctan2(sine, cosine))
        if sine == 0 or abs(arc-trial['step']) > 1e-10:
            raise ValueError('SAVED_TRIAL_GEODESIC_MISMATCH')
        direction = tangent/sine
        for fraction in (1., .5, .25):
            step = fraction*trial['step']
            target = p*np.cos(step)+direction*np.sin(step)
            q = target.copy(); q[7:] *= .5
            # Literal coefficients live on their own unit sphere; scaling
            # preserves the material and recovers the dual coordinate above.
            q /= np.linalg.norm(q)
            if fraction == 1.:
                # Keep literal binary64 operands of the recorded trial. A
                # normalize/denormalize roundtrip can move them by one ulp.
                q = np.asarray(trial['quadric_coefficients'], dtype=float).copy()
            if np.linalg.norm(coordinate(q)-target) > 1e-12:
                raise ValueError('COORDINATE_ROUNDTRIP_FAILED')
            probes.append(dict(base_evaluation=base_index, original_trial=trial['evaluation_index'],
                               base_material_sha256=maps[base_index]['sha256'],
                               original_trial_material_sha256=maps[trial['evaluation_index']]['sha256'],
                               fraction=fraction, step=float(step), coefficients=q.tolist(),
                               existing_trial=fraction == 1., independent_geometry_evaluated=False))
    return dict(schema='heart-actual-direction-probe-registration-v1',
                events_sha256=hashlib.sha256(payload).hexdigest(),
                configuration_sha256=hashlib.sha256((source/'configuration.json').read_bytes()).hexdigest(),
                probes=probes, optimization_invocations=0, derivative_admission=False,
                scope='Actual saved trial geodesic, not a replacement for original three-scale FD protocol')


if __name__ == '__main__':
    here = Path(__file__).resolve().parent
    result = register(here/'heart-guarded-formal-seed-0-job-01/output/run')
    with Path(sys.argv[1]).open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
    print(json.dumps([dict(base=p['base_evaluation'], step=p['step']) for p in result['probes']]))
