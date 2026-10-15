"""Compare saved actual trial moments against the preceding accepted Jacobian."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from audit_heart_stalled_prefix import coordinate, symbolic_model


def residual(model, r, p):
    _, df, _, g, dg, _ = model
    grad, G = np.asarray(df(*r)).ravel(), np.asarray(dg(*r))
    multipliers = np.linalg.lstsq(np.column_stack((-G.T, p)), -grad, rcond=None)[0]
    value = np.r_[grad + multipliers[-1]*p-G.T@multipliers[:-1], np.asarray(g(*r)).ravel()]
    return float(np.linalg.norm(value))


def audit(source):
    config_bytes = (source/'configuration.json').read_bytes()
    request = json.loads(config_bytes)['request']
    if request['spinning_top_objective'] != 'INERTIA_AND_LOWERING':
        raise ValueError('UNSUPPORTED_OBJECTIVE')
    event_bytes = (source/'events.jsonl').read_bytes()
    records = [json.loads(line) for line in event_bytes.splitlines()]
    maps = {row['evaluation']: row for row in records if row.get('kind') == 'MAP_SAVED'}
    for row in maps.values():
        if hashlib.sha256((source/row['evaluation_file']).read_bytes()).hexdigest() != row['sha256']:
            raise ValueError('CHECKPOINT_CHANGED')
    model = symbolic_model(request['center_of_mass_weight'])
    bases, trials = {}, []
    for record in records:
        if record.get('kind') != 'NEWTON_EVENT':
            continue
        event = record['event']
        payload = event.get('evaluation_payload')
        if event['type'] != 'evaluation' or not payload:
            continue
        restart = event['restart']
        p, r = coordinate(event['quadric_coefficients']), np.asarray(payload['solver_raw_moments'])
        base = bases.get(restart)
        if base is not None and (event.get('step') or 0) > 0:
            b, bp, br = base
            dr, dp = r-br, p-bp
            predicted_delta = np.asarray(b['evaluation_payload']['solver_moment_jacobian']) @ dp
            error = np.linalg.norm(dr-predicted_delta)
            predicted_r = br+predicted_delta
            trials.append(dict(base_evaluation=b['evaluation_index'], evaluation=event['evaluation_index'],
                step=event['step'], coordinate_chord=float(np.linalg.norm(dp)),
                moment_change=float(np.linalg.norm(dr)), predicted_change=float(np.linalg.norm(predicted_delta)),
                linearization_error=float(error), relative_linearization_error=float(error/max(
                    np.linalg.norm(dr), np.linalg.norm(predicted_delta), np.finfo(float).tiny)),
                actual_raw_residual=residual(model,r,p),
                predicted_raw_residual=(residual(model,predicted_r,p) if predicted_r[0]>0 else None),
                base_raw_residual=b['raw_equation_residual'], outcome=event['outcome'], reason=event['reason'],
                same_recorded_topology=(b['material_topology_hash']==event['material_topology_hash'])))
        if event['outcome'] == 'accepted':
            bases[restart] = event,p,r
    return dict(schema='heart-saved-trial-linear-prediction-v1', source=str(source),
        configuration_sha256=hashlib.sha256(config_bytes).hexdigest(),
        events_sha256=hashlib.sha256(event_bytes).hexdigest(), trials=trials,
        geometry_derivative_certification=False, convergence_certificate=False,
        scope='Saved secant comparison and independent residual prediction; topology labels are not branch certificates')


if __name__ == '__main__':
    here = Path(__file__).resolve().parent
    report = audit(here/'heart-guarded-formal-seed-0-job-01/output/run')
    with Path(sys.argv[1]).open('x') as stream:
        json.dump(report,stream,indent=2,allow_nan=False)
    print(json.dumps(report['trials'][-8:],indent=2))
