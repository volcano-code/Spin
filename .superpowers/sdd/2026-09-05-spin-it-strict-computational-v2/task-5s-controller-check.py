"""Independent per-ring QR fit + point-frame/tetrahedron oracle; NumPy only."""
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
    gi, gc, gl, gd, gv = vec('WEIGHTS')
    has_angle = mode != 4
    state = np.concatenate((vec('BETA'), vec('PHI') if has_angle else [], vec('PARAM')))
    n = nb + int(has_angle) + 4*nh
    assert len(state) == n and np.isfinite(state).all()

    outer_count = int(vec('OUTERCOUNT')[0])
    assert outer_count == 8 and faces.min() >= 0 and faces.max() < outer_count
    neighbors = [set() for _ in range(outer_count)]
    for triangle in faces:
        for vertex in triangle:
            neighbors[vertex].update(int(j) for j in triangle if j != vertex)
    local = []
    # Raw-coordinate QR fits: not the production normalized SVD or sparse LD.
    for vertex, ring in enumerate(neighbors):
        ids = [vertex] + sorted(ring)
        assert len(ring) >= 2
        a = []
        for xx, yy, zz in rest[ids]:
            a.extend([[xx, 0., zz, -yy, 1., 0., 0.],
                      [yy, -zz, 0., xx, 0., 1., 0.],
                      [zz, yy, -xx, 0., 0., 0., 1.]])
        q, r = np.linalg.qr(a, mode='reduced')
        assert np.linalg.matrix_rank(r) == 7
        delta = rest[vertex]-rest[sorted(ring)].mean(axis=0)
        local.append((ids, q, r, delta))

    def shape_residual(moved):
        residuals = []
        for ids, q, r, delta in local:
            fit = np.linalg.solve(r, q.T @ moved[ids].reshape(-1))
            transformed = fit[0]*delta + np.cross(fit[1:4], delta)
            residual = transformed-(moved[ids[0]]-moved[ids[1:]].mean(axis=0))
            residuals.append(residual)
        return np.asarray(residuals)

    def shape_value(moved):
        residuals = shape_residual(moved)
        energy = np.sum(residuals*residuals)
        displacement = np.sum((moved[:outer_count]-rest[:outer_count])**2)
        return gd*energy, gv*displacement

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
        shape, displacement = shape_value(moved)
        values = np.concatenate(([spin+fill+shape+displacement, spin, fill, shape, displacement], raw, constraints))
        return values, world, moved

    expected, world, moved = calculate(state)
    actual = np.concatenate((vec('OBJ'), vec('SHAPE'), vec('SPIN'), vec('CONSTRAINT')))
    assert np.isfinite(actual).all() and actual.shape == expected.shape
    length = np.linalg.norm(np.ptp(rest[:8], axis=0))
    raw_scales = density * length**np.array([3, 4, 4, 4, 5, 5, 5, 5, 5, 5])
    constraint_scales = raw_scales[[1, 2] + ([3] if mode == 0 else []) + ([8, 7, 9] if has_angle else [])]
    objective_scale = 1 + gi + gc*(length**2 + density**2*length**8) + gl*nb + (gd+gv)*outer_count*length**2
    output_scales = np.concatenate((np.full(5, objective_scale), raw_scales, constraint_scales))
    value_error = np.max(np.abs(actual-expected)/output_scales)
    assert value_error <= 1e-10
    assert np.max(np.abs(world-vec('WORLD'))/raw_scales) <= 1e-10
    assert np.max(np.abs(moved-indexed(records['MOVED']))) <= 1e-12
    for left, right in [('OBJ','VOBJ'), ('SHAPE','VSHAPE'), ('SPIN','VSPIN'), ('CONSTRAINT','VCONSTRAINT')]:
        assert np.array_equal(vec(left),vec(right))
    assert np.array_equal(indexed(records['MOVED']), indexed(records['VMOVED']))
    expected_j = np.empty((len(expected), n))
    for column in range(n):
        perturbed = state.astype(complex);perturbed[column] += 1j*1e-30
        expected_j[:, column] = calculate(perturbed)[0].imag/1e-30
    actual_j = np.vstack((vec('GRAD'), indexed(records['JM']).T, indexed(records['JC']).T))
    selection = [0] + list(range(5, len(expected)))
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
    shape_vertex_gradient = np.empty((outer_count, 3))
    for vertex in range(outer_count):
        for axis in range(3):
            moved_step = moved.astype(complex)
            moved_step[vertex, axis] += 1j*1e-30
            shape_vertex_gradient[vertex, axis] = sum(shape_value(moved_step)).imag/1e-30
    actual_vertex_gradient = indexed(records['SHAPEVGRAD'])
    assert actual_vertex_gradient.shape == shape_vertex_gradient.shape
    shape_difference = np.abs(actual_vertex_gradient-shape_vertex_gradient)
    assert np.all(shape_difference <= 1e-10*length + 1e-7*np.abs(shape_vertex_gradient))
    assert list(map(int, records['SHAPEFLAGS'][0])) == [1, 0, 0, 0]
    assert np.max(np.abs(indexed(records['SHAPERES'])-shape_residual(moved))) <= 1e-10*length
    shape_joint = expected_j[3] + expected_j[4]
    handle_begin = nb + int(has_angle)
    assert np.array_equal(shape_joint[:handle_begin], np.zeros(handle_begin))
    assert np.all(np.abs(vec('SHAPEHGRAD')-shape_joint[handle_begin:]) <=
                  1e-10*objective_scale+1e-7*np.abs(shape_joint[handle_begin:]))
    assert np.all(np.abs(vec('GRAD')-vec('BASEGRAD')-shape_joint) <=
                  1e-10*objective_scale+1e-7*np.abs(shape_joint))
    # Additional diagnostics and explicit near-zero gate; no threshold relaxation.
    shape_actual = np.concatenate((actual_vertex_gradient.reshape(-1), vec('SHAPEHGRAD')))
    shape_expected = np.concatenate((shape_vertex_gradient.reshape(-1), shape_joint[handle_begin:]))
    shape_error = np.abs(shape_actual-shape_expected)
    nz_shape = np.abs(shape_expected) > 1e-7
    shape_relative = np.max(shape_error[nz_shape]/np.abs(shape_expected[nz_shape]))
    shape_absolute = np.max(shape_error[~nz_shape]) if np.any(~nz_shape) else 0.
    assert shape_relative <= 1e-6 and shape_absolute <= 1e-10*length
    rest_shape, _ = shape_value(rest)
    translation = np.array([.17, -.13, .21])
    small_axis = np.array([.03, -.04, .02])
    linear_motion = 1.2*rest + np.cross(small_axis, rest) + translation
    null_shape, _ = shape_value(linear_motion)
    assert abs(rest_shape) < 1e-20*length**2
    assert abs(null_shape) < 1e-20*length**2
    print(f'PASS fixture={fixture} mode={mode} columns={n} value_scaled={value_error:.3e} '
          f'joint_derivative_scaled={derivative_error:.3e} nonzero_relative={relative_error:.3e} '
          f'nearzero_absolute={nearzero_error:.3e} shape_relative={shape_relative:.3e} '
          f'shape_nearzero_absolute={shape_absolute:.3e}')
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
    print('10 actual source cases verified; synthetic shape integration, not a joint solver or paper reproduction.')


if __name__ == '__main__':
    main(sys.argv[1])
