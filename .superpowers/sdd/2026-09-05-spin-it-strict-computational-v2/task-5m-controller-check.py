"""Independent signed-tetra / matrix-inertia / complex-step joint check.

Consumes the native export, imports no make-it-balance module. Unit-density
mesh-functional check only: not a geometry seal or paper-case certificate.
"""
import json
import sys
import numpy as np


def moments(vertices, faces):
    total = np.zeros(10, dtype=vertices.dtype)
    for indices in faces:
        a, b, c = vertices[indices]
        volume = np.dot(a, np.cross(b, c)) / 6
        point_sum = a + b + c
        first = volume * point_sum / 4
        second = volume * (np.outer(point_sum, point_sum) +
                           np.outer(a, a) + np.outer(b, b) + np.outer(c, c)) / 20
        total += np.array([volume, *first, second[0, 0], second[1, 1],
                           second[2, 2], second[1, 2], second[0, 2], second[0, 1]])
    return total


def deform(vertices, weights, pivots, rotations, handle_state):
    result = np.zeros(vertices.shape, dtype=handle_state.dtype)
    for index, pivot in enumerate(pivots):
        translation = handle_state[4 * index:4 * index + 3]
        scale = handle_state[4 * index + 3]
        transformed = pivot + translation + scale * ((vertices - pivot) @ rotations[index].T)
        result += weights[:, index, None] * transformed
    return result


def independent(data, state, mode):
    beta_count = len(data['cells'])
    angle = mode != 'STATIC_LOWERING'
    beta = state[:beta_count]
    phi = state[beta_count] if angle else 0
    handles = state[beta_count + int(angle):]
    pivots = np.asarray(data['pivots'])
    rotations = np.asarray(data['rotations'])
    meshes = [data['outer'], *data['cells']]
    deformed = [deform(np.asarray(m['vertices']), np.asarray(m['weights']),
                       pivots, rotations, handles) for m in meshes]
    raw = moments(deformed[0], np.asarray(meshes[0]['faces']))
    for value, vertices, mesh in zip(beta, deformed[1:], meshes[1:]):
        raw -= value * moments(vertices, np.asarray(mesh['faces']))
    mass = raw[0]
    first = raw[1:4]
    second = np.array([[raw[4], raw[9], raw[8]],
                       [raw[9], raw[5], raw[7]],
                       [raw[8], raw[7], raw[6]]])
    inertia = np.trace(second) * np.eye(3) - second
    constraints = [first[0], first[1]]
    spin = 0
    if angle:
        c, s = np.cos(phi), np.sin(phi)
        rotation = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
        rotated = rotation @ inertia @ rotation.T
        if mode != 'YOYO':
            rotated = rotated - np.diag([first[2]**2 / mass, first[2]**2 / mass, 0])
        spin = data['gamma_i'] * ((rotated[0, 0] / inertia[2, 2])**2 +
                                  (rotated[1, 1] / inertia[2, 2])**2)
        if mode == 'YOYO':
            constraints.append(first[2])
        constraints.extend([second[0, 2], second[1, 2], -rotated[0, 1]])
    if mode == 'TOP_MASS_AND_LOWERING':
        spin += data['gamma_c'] * first[2]**2
    elif mode in ('TOP_LOWERING_ONLY', 'STATIC_LOWERING'):
        spin += data['gamma_c'] * (first[2] / mass)**2
    fill = 0.5 * data['gamma_l'] * beta @ np.asarray(data['fill_laplacian']) @ beta
    rest = np.asarray(data['outer']['vertices'])
    graph = np.zeros((len(rest), len(rest)))
    edges = {tuple(sorted((int(face[i]), int(face[(i + 1) % 3]))))
             for face in data['outer']['faces'] for i in range(3)}
    for a, b in edges:
        graph[a, a] += 1
        graph[b, b] += 1
        graph[a, b] -= 1
        graph[b, a] -= 1
    laplacian_positions = graph @ deformed[0]
    shape = data['gamma_dl'] * np.sum(laplacian_positions**2)
    displacement = data['gamma_dv'] * np.sum((deformed[0] - rest)**2)
    return np.array([spin + fill + shape + displacement, *constraints]), raw, deformed


def require_close(actual, expected, label, atol=1e-11, rtol=1e-10):
    actual, expected = np.asarray(actual), np.asarray(expected)
    assert actual.shape == expected.shape, (label, actual.shape, expected.shape)
    assert np.isfinite(actual).all() and np.isfinite(expected).all(), label
    residual = np.abs(actual - expected)
    assert np.all(residual <= atol + rtol * np.abs(expected)), (label, residual.max())
    return float(residual.max(initial=0))


def main(path):
    with open(path, encoding='utf8') as source:
        data = json.load(source)
    assert data['evidence_class'] == 'SYNTHETIC_REGRESSION'
    assert len(data['pivots']) > 16
    assert len(data['cases']) == 5 and len({c['mode'] for c in data['cases']}) == 5
    assert data['bbw_success'] and len(data['bbw_output_hash']) == 64
    summary = {}
    for case in data['cases']:
        mode = case['mode']
        state = np.asarray(case['state'], dtype=float)
        assert len(state) == len(data['cells']) + int(mode != 'STATIC_LOWERING') + 4 * len(data['pivots'])
        values, raw, vertices = independent(data, state, mode)
        value_error = require_close(case['values'], values, mode + ':values')
        moment_error = require_close(case['moments'], raw, mode + ':moments')
        geometry_error = require_close(case['deformed_outer'], vertices[0], mode + ':outer')
        derivative = np.zeros((len(values), len(state)))
        raw_derivative = np.zeros((10, len(state)))
        for j in range(len(state)):
            point = state.astype(complex)
            point[j] += 1e-26j
            complex_values, complex_raw, _ = independent(data, point, mode)
            derivative[:, j] = complex_values.imag / 1e-26
            raw_derivative[:, j] = complex_raw.imag / 1e-26
        derivative_error = require_close(case['derivatives'], derivative, mode + ':derivatives',
                                         atol=2e-10, rtol=1e-6)
        raw_derivative_error = require_close(case['moment_derivatives'], raw_derivative,
                                             mode + ':moment_derivatives', atol=1e-11, rtol=1e-9)
        central_errors = []
        for step in (1e-3, 1e-4, 1e-5):
            central = np.zeros_like(derivative)
            for j in range(len(state)):
                plus, minus = state.copy(), state.copy()
                plus[j] += step
                minus[j] -= step
                central[:, j] = (independent(data, plus, mode)[0] -
                                 independent(data, minus, mode)[0]) / (2 * step)
            central_errors.append(float(np.max(np.abs(central - derivative))))
        assert central_errors[1] < central_errors[0] / 20, (mode, central_errors)
        assert central_errors[2] < central_errors[0] / 100, (mode, central_errors)
        summary[mode] = dict(state_count=len(state), value_error=value_error,
                             moment_error=moment_error, geometry_error=geometry_error,
                             derivative_error=derivative_error,
                             raw_derivative_error=raw_derivative_error, central_errors=central_errors)
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main(sys.argv[1])
