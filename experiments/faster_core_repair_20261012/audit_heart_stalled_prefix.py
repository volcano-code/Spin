"""Independent algebra and timing audit of the stopped Heart production prefix.

Reads the original checkpointed moments/Jacobians; invokes no geometry or
optimizer. This does not certify those Jacobians by finite differences.
"""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent


def symbolic_model(weight):
    r = sp.symbols('m x y z xx yy zz yz xz xy')
    m, x, y, z, xx, yy, zz, yz, xz, xy = r
    a, b, c = yy+zz-(y*y+z*z)/m, xx+zz-(x*x+z*z)/m, xx+yy-(x*x+y*y)/m
    d = -xy+x*y/m
    f = (a*a+b*b+2*d*d)/(c*c)+sp.Rational(str(weight))*z/m
    g = sp.Matrix([x, y, -xz+x*z/m, -yz+y*z/m])
    return tuple(sp.lambdify(r, v, 'numpy', cse=True) for v in
                 (f, sp.Matrix([sp.diff(f, t) for t in r]), sp.hessian(f, r),
                  g, g.jacobian(r), [sp.hessian(v, r) for v in g]))


def coordinate(q):
    p = np.asarray(q, dtype=float).copy()
    p[7:] *= 2
    return p/np.linalg.norm(p)


def relative_error(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.max(np.abs(a-b))/max(1., np.max(np.abs(a)), np.max(np.abs(b))))


def independent_state(event, model):
    f, df, hf, g, dg, hg = model
    payload = event['evaluation_payload']
    r = np.asarray(payload['solver_raw_moments'])
    p = coordinate(event['quadric_coefficients'])
    J = np.asarray(payload['solver_moment_jacobian'])
    grad, G = np.asarray(df(*r)).ravel(), np.asarray(dg(*r))
    constraints = np.asarray(g(*r)).ravel()
    multipliers = np.linalg.lstsq(np.column_stack((-G.T, p)), -grad, rcond=None)[0]
    lam, mu = multipliers[:-1], multipliers[-1]
    residual = np.r_[grad+mu*p-G.T@lam, constraints]
    H = np.asarray(hf(*r))-np.einsum('i,ijk->jk', lam, np.asarray(hg(*r)))
    A = np.block([[H@J+mu*np.eye(10), -G.T, p[:, None]],
                  [G@J, np.zeros((4, 4)), np.zeros((4, 1))],
                  [2*p[None, :], np.zeros((1, 4)), np.zeros((1, 1))]])
    rhs = np.r_[-residual, 0.]
    return dict(p=p, r=r, J=J, matrix=A, rhs=rhs, mu=float(mu),
                residual=float(np.linalg.norm(residual)), objective=float(f(*r)),
                constraints=constraints)


def audit(directory):
    config = json.loads((directory/'configuration.json').read_text())
    request = config['request']
    if request['mode'] != 'SPINNING_TOP' or request['spinning_top_objective'] != 'INERTIA_AND_LOWERING':
        raise ValueError('unsupported objective for this independent model')
    records = [json.loads(line) for line in (directory/'events.jsonl').read_text().splitlines()]
    maps = {row['evaluation']: row for row in records if row.get('kind') == 'MAP_SAVED'}
    for record in maps.values():
        actual = hashlib.sha256((directory/record['evaluation_file']).read_bytes()).hexdigest()
        if actual != record['sha256']:
            raise ValueError('checkpoint digest mismatch')
    events = [row['event'] for row in records if row.get('kind') == 'NEWTON_EVENT']
    evaluations = [e for e in events if e['type'] == 'evaluation' and e.get('evaluation_payload')]
    accepted = [e for e in evaluations if e['outcome'] == 'accepted']
    updates = [e for e in accepted if (e.get('step') or 0) > 0]
    last = updates[-1]
    successor = next(e for e in evaluations if e['evaluation_index'] > last['evaluation_index'])
    selected = [accepted[0], updates[0], last, successor]
    model = symbolic_model(request['center_of_mass_weight'])
    states = {e['evaluation_index']: independent_state(e, model) for e in evaluations}
    errors = []
    for e in evaluations:
        s = states[e['evaluation_index']]
        errors.append(dict(evaluation=e['evaluation_index'],
            coordinate=relative_error(s['p'], e['evaluation_payload']['solver_coordinates']),
            objective=relative_error(s['objective'], e['objective']),
            multiplier=relative_error(s['mu'], e['multiplier_mu']),
            raw_residual=relative_error(s['residual'], e['raw_equation_residual'])))
    matrices = []
    for event in events:
        payload = event.get('newton_system_payload')
        if not payload:
            continue
        index = event['evaluation_index']
        s = states[index]
        native_matrix = np.asarray(payload['matrix'])
        native_rhs = np.asarray(payload['right_hand_side'])
        row_scales = np.max(np.abs(s['matrix']), axis=1)
        equilibrated = s['matrix']/row_scales[:, None]
        equilibrated /= np.max(np.abs(equilibrated), axis=0)[None, :]
        solution = np.linalg.solve(s['matrix'], s['rhs'])
        tangent = solution[:10]-s['p']*np.dot(s['p'], solution[:10])
        nxt = next((e for e in evaluations if e['restart'] == event['restart']
                    and e['evaluation_index'] > index), None)
        row = dict(evaluation=index, matrix_relative_error=relative_error(s['matrix'], native_matrix),
                   rhs_relative_error=relative_error(s['rhs'], native_rhs),
                   condition=float(np.linalg.cond(s['matrix'])),
                   reversible_equilibrated_condition=float(np.linalg.cond(equilibrated)),
                   linear_backward_error=float(np.linalg.norm(s['matrix']@solution-s['rhs'], np.inf)/
                       (np.linalg.norm(s['matrix'], np.inf)*np.linalg.norm(solution, np.inf)+np.linalg.norm(s['rhs'], np.inf))))
        if nxt is not None:
            step = nxt['step']
            predicted = np.cos(step)*s['p']+np.sin(step)*tangent/np.linalg.norm(tangent)
            predicted /= np.linalg.norm(predicted)
            row['next_evaluation'] = nxt['evaluation_index']
            row['geodesic_error'] = relative_error(predicted, coordinate(nxt['quadric_coefficients']))
        matrices.append(row)
    timing, pending = {}, []
    for row in records:
        if row.get('kind') == 'STAGE_STARTED':
            pending.append([row['stage'], row['monotonic_ns'], 0])
        elif row.get('kind') == 'STAGE_RETURNED':
            if not pending or pending[-1][0] != row['stage']:
                raise ValueError('unmatched stage nesting')
            name, start, children = pending.pop()
            elapsed = row['elapsed_ns']
            if elapsed < children:
                raise ValueError('negative exclusive stage time')
            item = timing.setdefault(name, dict(calls=0, inclusive_seconds=0., exclusive_seconds=0.))
            item['calls'] += 1
            item['inclusive_seconds'] += elapsed/1e9
            item['exclusive_seconds'] += (elapsed-children)/1e9
            if pending:
                pending[-1][2] += elapsed
    maxima = {key: max(e[key] for e in errors) for key in ('coordinate', 'objective', 'multiplier', 'raw_residual')}
    matrix_ok = all(e['matrix_relative_error'] <= 1e-9 and e['rhs_relative_error'] <= 1e-9
                    and e.get('geodesic_error', 0) <= 1e-6 for e in matrices)
    return dict(schema='heart-stalled-prefix-independent-algebra-v1',
        directory=str(directory), configuration_sha256=hashlib.sha256((directory/'configuration.json').read_bytes()).hexdigest(),
        verified_material_checkpoints=len(maps), valid_evaluations=len(evaluations), accepted_updates=len(updates),
        converged=False, geometry_derivative_certification=False, case_complete=False,
        model_relative_errors=maxima, model_consistent=all(v <= 1e-9 for v in maxima.values()),
        native_matrices_consistent=matrix_ok, matrices=matrices, evaluations=errors,
        diagnostic_states=[dict(evaluation=e['evaluation_index'], coefficients=e['quadric_coefficients'],
            checkpoint=maps[e['evaluation_index']], raw_residual=e['raw_equation_residual'],
            mu=e['multiplier_mu'], step=e['step'], outcome=e['outcome'], reason=e['reason']) for e in selected],
        accepted_trajectory=[dict(evaluation=e['evaluation_index'], raw=e['raw_equation_residual'], step=e['step']) for e in accepted],
        timing=timing, unfinished_stages=[item[0] for item in pending],
        scope='Independent Eq.14 algebra and saved native matrix/next valid trial action; not a geometry derivative or convergence certificate')


if __name__ == '__main__':
    source = HERE/'heart-guarded-formal-seed-0-job-01/output/run'
    destination = Path(sys.argv[1])
    report = audit(source)
    with destination.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps({k: report[k] for k in ('valid_evaluations', 'accepted_updates', 'model_relative_errors',
        'model_consistent', 'native_matrices_consistent', 'unfinished_stages')}))
