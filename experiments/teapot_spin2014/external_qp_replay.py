"""ProxQP comparison of a saved SLQP subproblem, never a paper solver profile.

The immutable replay supplies x, g, c, A and BFGS secants, NOT an initial
solution. ProxQP solves min g.d + .5*d.B.d, A.d=-c, -x<=d<=1-x on
the bounded prefix; remaining variables are free. Independent direct-BFGS
KKT checks retain the native 1e-9 residual and 1e-10 box/dual gates.
Dense expansion is intentionally limited: this is a reference experiment,
not a scalable replacement for the matrix-free native implementation.
"""
import argparse
import hashlib
import importlib.metadata
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from direction_replay import hessian_action


def problem(value, bounded):
    t = value['terminal']
    x, g, c, a = (np.asarray(t[k], dtype=float) for k in
                  ('state', 'gradient', 'constraints', 'constraint_jacobian'))
    if (x.ndim != 1 or not len(x) or type(bounded) is not int or not 0 < bounded <= len(x)
            or g.shape != x.shape or c.ndim != 1 or a.shape != (len(c), len(x))
            or not all(np.isfinite(v).all() for v in (x, g, c, a))
            or np.any((x[:bounded] < 0) | (x[:bounded] > 1))):
        raise ValueError('INVALID_QP_INPUT')
    pairs = []
    for pair in value['history']:
        s, y = np.asarray(pair['s'], dtype=float), np.asarray(pair['y'], dtype=float)
        if (s.shape != x.shape or y.shape != x.shape or
                not all(np.isfinite(v).all() for v in (s, y)) or
                not np.isfinite(s @ y) or s @ y <= 1e-12 * np.linalg.norm(s) * np.linalg.norm(y)):
            raise ValueError('INVALID_BFGS_HISTORY')
        pairs.append((s, y))
    return x, g, c, a, pairs


def dense_hessian(n, pairs):
    # Independently assemble the ordinary BFGS recursion; checker uses actions.
    scale = 1. if not pairs else np.clip((pairs[-1][0] @ pairs[-1][1]) /
                                       (pairs[-1][1] @ pairs[-1][1]), 1e-6, 1e6)
    h = np.eye(n) / scale
    for s, y in pairs:
        bs = h @ s
        sbs, sy = s @ bs, s @ y
        if not np.isfinite([sbs, sy]).all() or min(sbs, sy) <= 0:
            raise ValueError('INVALID_BFGS_CURVATURE')
        h += np.outer(y, y) / sy
        h -= np.outer(bs, bs) / sbs
    if not np.isfinite(h).all():
        raise ValueError('NONFINITE_HESSIAN')
    return h


def check_external_solution(value, bounded, direction, equality_dual, box_dual):
    x, g, c, a, pairs = problem(value, bounded)
    d, y, z = (np.asarray(v, dtype=float) for v in (direction, equality_dual, box_dual))
    if (d.shape != x.shape or y.shape != c.shape or z.shape != x.shape or
            not all(np.isfinite(v).all() for v in (d, y, z))):
        raise ValueError('INVALID_EXTERNAL_SOLUTION')
    norm = lambda v: float(np.max(np.abs(v), initial=0.))
    lower_gap = d[:bounded] + x[:bounded]
    upper_gap = 1 - x[:bounded] - d[:bounded]
    result = dict(
        stationarity=norm(g + hessian_action(d, pairs) + a.T @ y + z),
        equality_error=norm(a @ d + c),
        box_error=max(0., float(np.max(-lower_gap, initial=0.)),
                      float(np.max(-upper_gap, initial=0.))),
        free_dual_error=norm(z[bounded:]),
        complementarity=max(norm(np.maximum(z[:bounded], 0) * upper_gap),
                            norm(np.minimum(z[:bounded], 0) * lower_gap)))
    # Signed box dual: negative on a lower face, positive on an upper face.
    # A non-negligible multiplier must belong to its own face, not merely
    # satisfy a small product when the coefficient happens to be tiny.
    wrong_lower = (z[:bounded] < -1e-10) & (np.abs(lower_gap) > 1e-9)
    wrong_upper = (z[:bounded] > 1e-10) & (np.abs(upper_gap) > 1e-9)
    if (not all(np.isfinite(v) for v in result.values()) or
            max(result['stationarity'], result['equality_error'], result['complementarity']) > 1e-9 or
            max(result['box_error'], result['free_dual_error']) > 1e-10 or
            np.any(wrong_lower | wrong_upper)):
        raise ValueError('INDEPENDENT_QP_KKT_FAILED:' + json.dumps(result, allow_nan=False))
    return dict(result, status='FIXED_QP_FIRST_ORDER_VERIFIED')


def solve_replay(value, *, bounded, maximum_dimension=9500):
    x, g, c, a, pairs = problem(value, bounded)
    if type(maximum_dimension) is not int or len(x) > maximum_dimension:
        raise ValueError('DENSE_DIMENSION_LIMIT')
    import proxsuite
    start = time.perf_counter()
    h = dense_hessian(len(x), pairs)
    assembled = time.perf_counter()
    lower, upper = np.full(len(x), -np.inf), np.full(len(x), np.inf)
    lower[:bounded], upper[:bounded] = -x[:bounded], 1 - x[:bounded]
    qp = proxsuite.proxqp.dense.QP(len(x), len(c), 0, True)
    qp.settings.eps_abs = 1e-12
    qp.settings.eps_rel = 0.
    qp.settings.max_iter = 1000
    qp.settings.compute_timings = True
    qp.settings.initial_guess = proxsuite.proxqp.InitialGuess.NO_INITIAL_GUESS
    qp.init(H=h, g=g, A=a, b=-c, l_box=lower, u_box=upper)
    initialized = time.perf_counter()
    qp.solve()
    stopped = time.perf_counter()
    d, y, z = qp.results.x.copy(), qp.results.y.copy(), qp.results.z.copy()
    result = dict(schema='spin2014-external-qp-comparison-v1',
                  scope='EXTERNAL_FIXED_QP_COMPARISON_ONLY',
                  nonlinear_convergence='NOT_CHECKED', paper_reproduction='NOT_CHECKED',
                  solver='ProxQP', version=proxsuite.__version__, dimension=len(x), bounded=bounded,
                  history_pairs=len(pairs), initial_guess='NO_INITIAL_GUESS',
                  eps_abs=1e-12, eps_rel=0., maximum_iterations=1000,
                  solver_status=str(qp.results.info.status), iterations=qp.results.info.iter,
                  direction=d.tolist(), equality_dual=y.tolist(), box_dual=z.tolist(),
                  seconds=dict(hessian=assembled-start, init=initialized-assembled,
                               solve=stopped-initialized, total=stopped-start))
    try:
        result['verification'] = check_external_solution(value, bounded, d, y, z)
    except ValueError as error:
        result['verification'] = dict(status='UNVERIFIED', reason=str(error))
    # Library and independent status are separately retained; no auto-polishing.
    return result


def file_hash(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''): h.update(chunk)
    return h.hexdigest()


def check_snapshot(snapshot):
    for path, expected in snapshot.items():
        if not Path(path).is_file() or file_hash(path) != expected:
            raise ValueError('STALE_RESULT:' + path)


def mapped_shared_objects():
    """Linux runtime inventory includes transitive BLAS/OpenMP, not just modules."""
    paths = set()
    for line in Path('/proc/self/maps').read_text().splitlines():
        fields = line.split(None, 5)
        if len(fields) == 6 and fields[5].startswith('/') and '.so' in Path(fields[5]).name:
            path = Path(fields[5])
            if not path.is_file(): raise ValueError('UNREADABLE_LOADED_LIBRARY:' + str(path))
            paths.add(str(path.resolve()))
    if not paths: raise ValueError('EMPTY_RUNTIME_BINARY_INVENTORY')
    return sorted(paths)


def verify_scientific_result(result, value):
    if (result.get('schema') != 'spin2014-external-qp-comparison-v1' or
            result.get('scope') != 'EXTERNAL_FIXED_QP_COMPARISON_ONLY' or
            result.get('nonlinear_convergence') != 'NOT_CHECKED' or
            result.get('paper_reproduction') != 'NOT_CHECKED'):
        raise ValueError('INVALID_COMPARISON_SCOPE')
    return check_external_solution(value, result['bounded'], result['direction'],
                                   result['equality_dual'], result['box_dual'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--bounded', type=int, required=True)
    args = parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    thread_environment = {k: os.environ.get(k) for k in
                          ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS')}
    if any(v != '1' for v in thread_environment.values()):
        raise ValueError('REFERENCE_RUN_REQUIRES_SINGLE_THREAD_ENVIRONMENT')
    raw = args.input.read_bytes()
    value = json.loads(raw)
    root = Path(__file__).resolve().parents[2]
    sources = [Path(__file__), Path(__file__).with_name('direction_replay.py'),
               root/'python/make_it_stand/verification/source_trace_archive.py']
    snapshot = {str(p.resolve()): file_hash(p) for p in sources}
    result = solve_replay(value, bounded=args.bounded)
    scientific = {k: v for k, v in result.items() if k != 'seconds'}
    result['scientific_sha256'] = hashlib.sha256(json.dumps(scientific, sort_keys=True,
                                            allow_nan=False).encode()).hexdigest()
    loaded = mapped_shared_objects()
    result['provenance'] = dict(input_path=str(args.input.resolve()),
        input_sha256=hashlib.sha256(raw).hexdigest(),
        sources=snapshot,
        loaded_binaries={p: file_hash(p) for p in sorted(loaded)},
        thread_environment=thread_environment,
        python_executable=str(Path(sys.executable).resolve()),
        python_sha256=file_hash(Path(sys.executable).resolve()),
        versions={p: importlib.metadata.version(p) for p in ('numpy', 'scipy', 'proxsuite', 'cmeel')},
        git_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
        tracked_diff_sha256=hashlib.sha256(subprocess.check_output(['git', 'diff', 'HEAD', '--binary'], cwd=root)).hexdigest(),
        maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if file_hash(args.input) != result['provenance']['input_sha256']:
        raise ValueError('INPUT_CHANGED_DURING_RUN')
    check_snapshot(snapshot)
    with args.output.open('x') as stream:
        json.dump(result, stream, sort_keys=True, allow_nan=False, indent=2)
        stream.write('\n')
    print(json.dumps({k: result[k] for k in ('solver_status', 'verification', 'seconds', 'scientific_sha256')}))
    return 0 if result['verification']['status'] == 'FIXED_QP_FIRST_ORDER_VERIFIED' else 2


if __name__ == '__main__':
    raise SystemExit(main())
