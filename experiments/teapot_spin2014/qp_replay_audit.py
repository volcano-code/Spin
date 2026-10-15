"""Independent direct-BFGS certificate check of a fixed native QP replay.

This is not nonlinear convergence, binary validation or geometry acceptance.
Engineering tolerances are the unchanged 1e-9 QP and 1e-10 box/dual gates.
"""
import numpy as np
from direction_replay import hessian_action


def check_replay(value, *, bounded):
    t = value['terminal']; q = t['finite_box']
    if q.get('checked') is not True or q.get('status') != 'BOX_QP_FIRST_ORDER_CHECKED':
        raise ValueError('QP_NOT_CHECKED')
    x = np.asarray(t['state'], dtype=float); n = len(x)
    g = np.asarray(t['gradient'], dtype=float)
    c = np.asarray(t['constraints'], dtype=float)
    a = np.asarray(t['constraint_jacobian'], dtype=float)
    d = np.asarray(q['direction'], dtype=float)
    lam = np.asarray(q['multipliers'], dtype=float)
    faces = q['working_set']
    if (type(bounded) is not int or not 0 < bounded <= n or x.shape != (n,) or
            g.shape != (n,) or d.shape != (n,) or c.ndim != 1 or a.shape != (len(c), n) or
            not isinstance(faces, list) or any(type(f) is not int or not 0 <= f < 2*bounded for f in faces) or
            len(set(f//2 for f in faces)) != len(faces) or lam.shape != (len(c)+len(faces),)):
        raise ValueError('INVALID_QP_DIMENSIONS')
    if not all(np.isfinite(v).all() for v in (x, g, c, a, d, lam)) or np.any((x[:bounded] < 0) | (x[:bounded] > 1)):
        raise ValueError('INVALID_QP_VALUES')
    pairs = []
    for p in value['history']:
        s, y = np.asarray(p['s'], dtype=float), np.asarray(p['y'], dtype=float)
        if s.shape != (n,) or y.shape != (n,) or not np.isfinite(s).all() or not np.isfinite(y).all():
            raise ValueError('INVALID_HISTORY')
        pairs.append((s, y))
    stationarity = g+hessian_action(d, pairs)+a.T@lam[:len(c)]
    face_error = dual_error = complementarity = 0.
    for j, f in enumerate(faces):
        multiplier = lam[len(c)+j]
        gap = d[f//2]-(f%2-x[f//2])
        stationarity[f//2] += multiplier
        face_error = max(face_error, abs(gap))
        dual_error = max(dual_error, -(1 if f%2 else -1)*multiplier)
        complementarity = max(complementarity, abs(multiplier*gap))
    norm = lambda v: float(np.max(np.abs(v), initial=0.))
    result = dict(stationarity=norm(stationarity), equality_error=norm(a@d+c),
                  box_error=max(0., float(np.max(-x[:bounded]-d[:bounded], initial=0.)),
                                float(np.max(d[:bounded]-(1-x[:bounded]), initial=0.))),
                  face_error=float(face_error), dual_error=float(dual_error),
                  complementarity=float(complementarity))
    if (not all(np.isfinite(v) for v in result.values()) or
            max(result['stationarity'], result['equality_error'], face_error, complementarity) > 1e-9 or
            max(result['box_error'], dual_error) > 1e-10):
        raise ValueError('INDEPENDENT_QP_KKT_FAILED')
    return dict(result, status='FIXED_QP_FIRST_ORDER_VERIFIED', nonlinear_convergence='NOT_CHECKED')
