"""Read-only scientific audit of saved Heart Newton systems (no geometry calls)."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np


def audit(run):
    path = Path(run) / 'events.jsonl'
    events = [r['event'] for line in path.open()
              if (r := json.loads(line)).get('kind') == 'NEWTON_EVENT']
    attempts = []
    for restart in sorted({e['restart'] for e in events}):
        rows = [e for e in events if e['restart'] == restart]
        accepted = [e for e in rows if e['type'] == 'evaluation'
                    and e.get('outcome') == 'accepted']
        attempts.append(dict(restart=restart,
            accepted_updates=sum(e.get('step', 0) > 0 for e in accepted),
            best_accepted_raw=min((e['raw_equation_residual'] for e in accepted), default=None),
            terminal=[e['reason'] for e in rows if e['type'] == 'attempt_terminal']))
    systems = []
    residual_components = []
    for event in events:
        if event['type'] != 'newton_linear_solve':
            continue
        rhs = np.asarray(event['newton_system_payload']['right_hand_side'])
        residual_components.append(dict(
            restart=event['restart'], evaluation_index=event['evaluation_index'],
            stationarity_norm=float(np.linalg.norm(rhs[:10])),
            constraint_residual=(-rhs[10:-1]).tolist(),
            constraint_norm=float(np.linalg.norm(rhs[10:-1])),
            raw_residual=event['raw_equation_residual'],
            multiplier_mu=event['multiplier_mu']))
    with mp.workdps(90):
        for index in (38, 45, 56):
            found = [e for e in events if e['type'] == 'newton_linear_solve'
                     and e['evaluation_index'] == index]
            if not found:
                continue
            e = found[0]
            v = e['newton_system_payload']
            a, b = np.array(v['matrix']), np.array(v['right_hand_side'])
            r, c = np.array(v['row_scales']), np.array(v['column_scales'])
            scaled = a / r[:, None] / c[None, :]
            double = np.linalg.solve(scaled, b / r) / c
            exact_input_a = mp.matrix([[mp.mpf(float(z)) for z in row] for row in a])
            exact_input_b = mp.matrix([mp.mpf(float(z)) for z in b])
            reference = np.array([float(z) for z in mp.lu_solve(exact_input_a, exact_input_b)])
            q = np.array(e['quadric_coefficients']); q[7:] *= 2; q /= np.linalg.norm(q)
            tangent = reference[:10] - q * np.dot(q, reference[:10])
            trial = next(z for z in events if z['type'] == 'evaluation'
                         and z['evaluation_index'] > index and z['restart'] == e['restart'])
            qt = np.array(trial['quadric_coefficients']); qt[7:] *= 2; qt /= np.linalg.norm(qt)
            t = trial['step']
            prediction = q * np.cos(t) + tangent / np.linalg.norm(tangent) * np.sin(t)
            systems.append(dict(evaluation_index=index,
                condition_original=float(np.linalg.cond(a)),
                condition_equilibrated=float(np.linalg.cond(scaled)),
                double_reference_relative_error=float(np.linalg.norm(double-reference)/np.linalg.norm(reference)),
                tangent_norm=float(np.linalg.norm(tangent)),
                next_trial_index=trial['evaluation_index'],
                next_trial_max_coordinate_error=float(np.max(np.abs(prediction-qt))),
                native_backward_error=v['linear_backward_error']))
    return dict(schema='heart-saved-newton-audit-v1',
                events_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                attempts=attempts, systems=systems, residual_components=residual_components,
                reference_decimal_digits=90,
                scope='Saved binary64 systems and next actual sphere trial only; no new derivative or convergence certification')


if __name__ == '__main__':
    import sys
    print(json.dumps(audit(sys.argv[1]), sort_keys=True, indent=2))
