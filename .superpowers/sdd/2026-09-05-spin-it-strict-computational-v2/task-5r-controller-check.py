"""Independent point-frame/tetrahedron/inertia oracle; NumPy and stdlib only."""
from collections import defaultdict
import sys
import numpy as np


def indexed(records, dtype=float):
    assert records and [int(r[0]) for r in records] == list(range(len(records)))
    data = np.asarray(records, dtype=dtype)
    assert np.isfinite(data).all()
    return data[:, 1:]


def integral(points, faces, coefficients):
    # Exact signed tetrahedron monomial moments; not a surface-normal kernel.
    a, b, c = (points[faces[:, i]] for i in range(3))
    volumes = np.einsum('ij,ij->i', a, np.cross(b, c)) * coefficients / 6
    total = a + b + c
    first = np.sum(volumes[:, None] * total, axis=0) / 4
    second = sum(np.einsum('ni,nj,n->ij', x, x, volumes)
                 for x in (total, a, b, c)) / 20
    return np.array([volumes.sum(), *first, second[0, 0], second[1, 1],
                     second[2, 2], second[1, 2], second[0, 2], second[0, 1]])


def verify(records):
    fixture, mode = map(int, records['CASE'][0][:2])
    density = float(records['CASE'][0][2])
    nb, nh = map(int, records['CASE'][0][3:])
    assert nb == 8 and nh == 26 and mode in range(5)
    assert list(map(int, records['FLAGS'][0])) == [1, 0, 0, 0, 0, 0, 0]
    vec = lambda key: np.asarray(records[key][0], dtype=float)
    rest = indexed(records['POINT'])
    centers = indexed(records['HANDLE'])
    nodes, cells = indexed(records['NODE']), indexed(records['CELL'], int)
    indices = cells[:, 3:]
    lo, hi = nodes[indices[:, 0], :3], nodes[indices[:, 7], :3]
    weights = []
    for p in rest:
        candidates = np.flatnonzero(np.all(p >= lo, axis=1) & np.all(p <= hi, axis=1))
        cell = min(candidates, key=lambda k: tuple(cells[k, :3]))
        u = (p - lo[cell]) / (hi[cell] - lo[cell])
        basis = [np.prod([u[k] if corner & (1 << k) else 1-u[k]
                          for k in range(3)]) for corner in range(8)]
        weights.append(np.asarray(basis) @ nodes[indices[cell], 3:])
    weights = np.asarray(weights)
    assert weights.shape == (len(rest), nh)
    faces, cavity_faces = indexed(records['FACE'], int), indexed(records['TRI'], int)
    coefficients = np.zeros((len(cavity_faces), nb))
    covered = set()
    for row in records['PATCH']:
        neg, pos, begin, count = map(int, row)
        ids = set(range(begin, begin+count))
        assert count >= 4 and begin >= 0 and begin+count <= len(cavity_faces)
        assert -1 <= neg < nb and -1 <= pos < nb and not covered.intersection(ids)
        covered.update(ids)
        if neg >= 0: coefficients[begin:begin+count, neg] -= 1
        if pos >= 0: coefficients[begin:begin+count, pos] += 1
    assert len(covered) == len(cavity_faces)
    all_faces = np.vstack((faces, cavity_faces))
    assert all_faces.min() >= 0 and all_faces.max() < len(rest)
    edges = np.asarray(records['EDGE'], int)
    assert edges.ndim == 2 and edges.shape[1] == 2
    assert edges.min() >= 0 and edges.max() < nb
    assert len({tuple(sorted(e)) for e in edges}) == len(edges)
    origin, rotation = vec('ORIGIN'), indexed(records['R'])
    assert np.max(np.abs(rotation @ rotation.T-np.eye(3))) < 1e-10
    assert abs(np.linalg.det(rotation)-1) < 1e-10
    gi, gc, gl = vec('WEIGHTS')
    has_angle = mode != 4
    state = np.concatenate((vec('BETA'), vec('PHI') if has_angle else [], vec('PARAM')))
    n = nb + int(has_angle) + 4*nh
    assert len(state) == n and np.isfinite(state).all()

    def calculate(x):
        beta, phi = x[:nb], x[nb] if has_angle else 0.
        parameters = x[nb+int(has_angle):].reshape(nh, 4)
        scales = parameters[:, 0]
        moved = rest * (weights @ scales)[:, None] + weights @ (
            (1-scales)[:, None] * centers + parameters[:, 1:])
        signed = np.concatenate((np.ones(len(faces)), coefficients @ beta))
        world = density * integral(moved, all_faces, signed)
        raw = density * integral((moved-origin) @ rotation.T, all_faces, signed)
        mass, lx, ly, lz, xx, yy, zz, yz, xz, xy = raw
        second = np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
        c, s = np.cos(phi), np.sin(phi)
        transverse = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
        rotated_q = transverse @ second @ transverse.T
        inertia = np.trace(rotated_q) * np.eye(3) - rotated_q
        constraints = [lx, ly]
        spin = 0. * mass
        if has_angle:
            # Paper top's constrained axial height, not recentering at arbitrary CoM.
            axial_shift = 0. if mode == 0 else lz*lz/mass
            ia, ib, ic = inertia[0, 0]-axial_shift, inertia[1, 1]-axial_shift, inertia[2, 2]
            spin = gi * ((ia/ic)**2 + (ib/ic)**2)
            if mode == 0: constraints.append(lz)
            constraints += [xz, yz, rotated_q[0, 1]]
        if mode in (2, 4): spin += gc*(lz/mass)**2
        if mode == 3: spin += gc*lz**2
        fill = .5*gl*np.sum((beta[edges[:, 0]]-beta[edges[:, 1]])**2)
        values = np.concatenate(([spin+fill, spin, fill], raw, constraints))
        return values, world, moved

    expected, world, moved = calculate(state)
    actual = np.concatenate((vec('OBJ'), vec('SPIN'), vec('CONSTRAINT')))
    assert np.isfinite(actual).all() and actual.shape == expected.shape
    length = np.linalg.norm(np.ptp(rest[:8], axis=0))
    raw_scales = density * length**np.array([3, 4, 4, 4, 5, 5, 5, 5, 5, 5])
    constraint_scales = raw_scales[[1, 2] + ([3] if mode == 0 else []) + ([8, 7, 9] if has_angle else [])]
    objective_scale = 1 + gi + gc*(length**2 + density**2*length**8) + gl*nb
    output_scales = np.concatenate((np.full(3, objective_scale), raw_scales, constraint_scales))
    value_error = np.max(np.abs(actual-expected)/output_scales)
    assert value_error <= 1e-10
    assert np.max(np.abs(world-vec('WORLD'))/raw_scales) <= 1e-10
    assert np.max(np.abs(moved-indexed(records['MOVED']))) <= 1e-12
    for left, right in [('OBJ','VOBJ'), ('SPIN','VSPIN'), ('CONSTRAINT','VCONSTRAINT')]:
        assert np.array_equal(vec(left),vec(right))
    assert np.array_equal(indexed(records['MOVED']), indexed(records['VMOVED']))
    expected_j = np.empty((len(expected), n))
    for column in range(n):
        perturbed = state.astype(complex);perturbed[column] += 1j*1e-30
        expected_j[:, column] = calculate(perturbed)[0].imag/1e-30
    actual_j = np.vstack((vec('GRAD'), indexed(records['JM']).T, indexed(records['JC']).T))
    selection = [0] + list(range(3, len(expected)))
    oracle_j, scales = expected_j[selection], output_scales[selection]
    assert actual_j.shape == oracle_j.shape and np.isfinite(actual_j).all()
    difference = np.abs(actual_j-oracle_j)
    assert np.all(difference <= 1e-10*scales[:, None] + 1e-7*np.abs(oracle_j))
    if has_angle: assert np.array_equal(actual_j[1:11, nb], np.zeros(10))
    derivative_error = np.max(difference/scales[:, None])
    nonzero = np.abs(oracle_j) > 1e-7
    relative_error = np.max(difference[nonzero]/np.abs(oracle_j[nonzero]))
    nearzero_error = np.max(difference[~nonzero]) if np.any(~nonzero) else 0.
    assert relative_error <= 1e-6 and nearzero_error <= 1e-11
    print(f'PASS fixture={fixture} mode={mode} columns={n} value_scaled={value_error:.3e} '
          f'joint_derivative_scaled={derivative_error:.3e} nonzero_relative={relative_error:.3e} '
          f'nearzero_absolute={nearzero_error:.3e}')
    return fixture, mode


def main(path):
    case, seen = None, set()
    with open(path, encoding='utf-8') as stream:
        for line in stream:
            pieces = line.split();assert pieces
            if pieces[0] == 'CASE': assert case is None;case = defaultdict(list)
            assert case is not None
            if pieces[0] == 'END':
                key = verify(case);assert key not in seen;seen.add(key);case = None
            else: case[pieces[0]].append(pieces[1:])
    assert case is None and seen == {(f, m) for f in range(2) for m in range(5)}
    print('10 actual source cases verified; synthetic integration, not a joint solver or paper reproduction.')


if __name__ == '__main__':
    main(sys.argv[1])
