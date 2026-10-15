"""Exact binary64-coefficient algebra along a saved accepted projective arc."""
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


def matrix(q):
    A = sp.Matrix([[q[4], q[9]/2, q[8]/2],
                   [q[9]/2, q[5], q[7]/2], [q[8]/2, q[7]/2, q[6]]])
    b = sp.Matrix(q[1:4])
    H = A.row_join(b/2).col_join((b.T/2).row_join(sp.Matrix([[q[0]]])))
    return A, b, H


def audit(run):
    events_bytes = (run/'events.jsonl').read_bytes()
    records = [json.loads(line) for line in events_bytes.splitlines()]
    maps = {r['evaluation']: r for r in records if r.get('kind') == 'MAP_SAVED'}
    states, coefficients = [], []
    for index in (38, 45):
        reference = maps[index]; path = run/reference['evaluation_file']
        if hashlib.sha256(path.read_bytes()).hexdigest() != reference['sha256']:
            raise ValueError('SAVED_MATERIAL_CHANGED')
        event = next(r['event'] for r in records if r.get('kind') == 'NEWTON_EVENT'
                     and r['event'].get('type') == 'evaluation'
                     and r['event']['evaluation_index'] == index)
        if event['outcome'] != 'accepted':
            raise ValueError('NOT_ACCEPTED')
        q = [sp.Rational(x) for x in json.loads(path.read_text())['quadric_coefficients']]
        A, b, H = matrix(q)
        center = -A.inv()*b/2
        level = H.det()/A.det()
        states.append(dict(evaluation=index, material_sha256=reference['sha256'],
                           center_exact=[str(x) for x in center], center_level_exact=str(level),
                           center_level_float=float(level)))
        coefficients.append(q)
    t = sp.Symbol('t')
    q = [(1-t)*a+t*b for a,b in zip(*coefficients)]
    A, _, H = matrix(q)
    det_A, det_H = sp.Poly(A.det(), t), sp.Poly(H.det(), t)
    a_roots = int(det_A.count_roots(0, 1))
    h_roots = int(det_H.count_roots(0, 1))
    intervals = [[str(lo),str(hi),multiplicity] for (lo,hi),multiplicity in
                 det_H.intervals(eps=sp.Rational(1, 10**12)) if hi>=0 and lo<=1]
    return dict(schema='heart-exact-quadric-transition-v1', states=states,
        events_sha256=hashlib.sha256(events_bytes).hexdigest(),
        A_determinant_coefficients=[str(x) for x in det_A.all_coeffs()],
        augmented_determinant_coefficients=[str(x) for x in det_H.all_coeffs()],
        A_singular_roots_on_segment=a_roots, quadric_degenerate_roots_on_segment=h_roots,
        degeneracy_parameter_intervals=intervals,
        interpretation='Positive projective interpolation traces the same minor spherical arc; coefficient scale does not change zero locus',
        domain_intersection_certified=False, nondifferentiability_certified=False,
        stagnation_root_cause_certified=False, case_complete=False)


if __name__ == '__main__':
    here = Path(__file__).resolve().parent
    report = audit(here/'heart-guarded-formal-seed-0-job-01/output/run')
    with Path(sys.argv[1]).open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps({k: v for k, v in report.items() if k not in
                     ('states','A_determinant_coefficients','augmented_determinant_coefficients')}))
