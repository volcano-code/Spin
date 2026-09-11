"""Independent geometry-derived small BBW reference, not native solver code.

Literal Freudenthal brick stencil and exact rational nodal boundaries; SciPy
SLSQP is an independent numerical comparison, not a paper-case/release gate.
Native acceptance uses independently recomputed residuals, never these flags.
"""
import argparse
import itertools
import json
import math
from fractions import Fraction as F

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, minimize


def reference(cells_per_axis, rule):
    n = cells_per_axis
    assert n in (2, 4)
    keys = list(itertools.product(range(n + 1), repeat=3))
    index = {key: i for i, key in enumerate(keys)}
    positions = [[F(k + 2, 8) for k in key] for key in keys]
    size = len(keys)
    k = [[F(0) for _ in keys] for _ in keys]
    mass = [F(0)] * size
    # Independent integrated path-tet stencil: twelve brick-edge terms.
    # First six touch diagonal endpoints000/111, so have twice the others.
    edges = [(0, 1, 2), (0, 2, 2), (0, 4, 2), (3, 7, 2), (5, 7, 2), (6, 7, 2),
             (1, 3, 1), (1, 5, 1), (2, 3, 1), (2, 6, 1), (4, 5, 1), (4, 6, 1)]
    for cell in itertools.product(range(n), repeat=3):
        ids = [index[tuple(cell[a] + ((b >> a) & 1) for a in range(3))] for b in range(8)]
        for a, b, factor in edges:
            i, j = ids[a], ids[b]
            value = F(factor, 48)  # physical cell side1/8 times factor/6
            k[i][i] += value
            k[j][j] += value
            k[i][j] -= value
            k[j][i] -= value
        for b, i in enumerate(ids):
            mass[i] += F(1, 512) * (F(1, 8) if rule == "RESTRICTED_GRID_VORONOI"
                                    else F(1, 4) if b in (0, 7) else F(1, 12))
    assert rule in ("RESTRICTED_GRID_VORONOI", "BARYCENTRIC_TET")
    assert sum(mass) == F(n, 8) ** 3
    assert all(sum(row) == 0 for row in k)
    handles = [key for key in itertools.product(range(0, n + 1, 2), repeat=3)
               if any(x in (0, n) for x in key)]
    handle_index = {key: i for i, key in enumerate(handles)}
    boundary = {}
    for cell in itertools.product(range(0, n, 2), repeat=3):
        for axis in range(3):
            axes = [a for a in range(3) if a != axis]
            for direction in (0, 2):
                plane = cell[axis] + direction
                if plane not in (0, n):
                    continue
                for u, v in itertools.product((F(0), F(1, 2), F(1)), repeat=2):
                    key = list(cell)
                    key[axis] = plane
                    key[axes[0]] += int(2 * u)
                    key[axes[1]] += int(2 * v)
                    row = [F(0)] * len(handles)
                    weights = {(0, 0): 1 - max(u, v), (2, 2): min(u, v),
                               (2, 0): max(u - v, 0), (0, 2): max(v - u, 0)}
                    for uv, value in weights.items():
                        corner = list(cell)
                        corner[axis] = plane
                        for a, bit in zip(axes, uv):
                            corner[a] += bit
                        row[handle_index[tuple(corner)]] += value
                    node = index[tuple(key)]
                    assert sum(row) == 1 and min(row) >= 0
                    assert node not in boundary or boundary[node] == row
                    boundary[node] = row
    free = [i for i in range(size) if i not in boundary]
    assert len(free) == (n - 1) ** 3
    assert len(boundary) == size - (n - 1) ** 3
    stiffness = np.array(k, dtype=float)
    mass_array = np.array(mass, dtype=float)
    q = stiffness.T @ (stiffness / mass_array[:, None])
    scale = float(np.max(np.diag(q)))
    fixed = np.zeros((size, len(handles)))
    for i, row in boundary.items():
        fixed[i, :] = np.array(row, dtype=float)
    return dict(nodes=np.array(positions, dtype=float), keys=keys, handles=[index[x] for x in handles],
                stiffness=stiffness, mass=mass_array, q=q, scale=scale, free=free,
                boundary=sorted(boundary), fixed=fixed, rule=rule, cells_per_axis=n)


def kkt(data, weights, coupled):
    assert weights.shape == data["fixed"].shape and np.isfinite(weights).all()
    x = weights[data["free"]]
    gradient = (data["q"] @ weights)[data["free"]] / data["scale"]
    stationarity = 0.
    complementarity = 0.
    residual = 0.
    multipliers = []
    for row, g in zip(x, gradient):
        lam = 0.
        if coupled:
            interior = (row > 1e-9) & (row < 1 - 1e-9)
            if interior.any():
                lam = -float(np.mean(g[interior]))
            else:
                lower = float(np.max(-g[row <= 1e-9])) if np.any(row <= 1e-9) else -math.inf
                upper = float(np.min(-g[row >= 1 - 1e-9])) if np.any(row >= 1 - 1e-9) else math.inf
                lam = min(max(0., lower), upper) if lower <= upper else (lower + upper) / 2
        shifted = g + lam
        # Unit projected gradient detects missing bound-sign conditions and
        # all-at-bound candidates as well as ordinary interior stationarity.
        projected = row - np.clip(row - shifted, 0., 1.)
        stationarity = max(stationarity, float(np.max(np.abs(projected))))
        lower_dual = np.where(row <= 1e-9, np.maximum(shifted, 0.), 0.)
        upper_dual = np.where(row >= 1 - 1e-9, np.maximum(-shifted, 0.), 0.)
        residual = max(residual, float(np.max(np.abs(shifted - lower_dual + upper_dual))))
        complementarity = max(complementarity, float(np.max(np.abs(lower_dual * row))),
                              float(np.max(np.abs(upper_dual * (1 - row)))))
        multipliers.append(lam)
    return dict(boundary=float(np.max(np.abs(weights[data["boundary"]] - data["fixed"][data["boundary"]]))),
                bounds=max(0., -float(weights.min()), float(weights.max()) - 1.),
                partition=float(np.max(np.abs(weights.sum(axis=1) - 1))),
                projected_stationarity=stationarity,
                stationarity=residual, complementarity=complementarity,
                energy=float(.5 * np.sum(weights * (data["q"] @ weights))),
                multipliers=multipliers)


def solve_reference(data, coupled, polish=False):
    free = data["free"]
    h = data["fixed"].shape[1]
    a = data["q"][np.ix_(free, free)] / data["scale"]
    c = (data["q"] @ data["fixed"])[free] / data["scale"]
    initial = np.full((len(free), h), 1 / h)
    history = []

    def objective(z):
        w = z.reshape(len(free), h)
        return float(.5 * np.sum(w * (a @ w)) + np.sum(w * c))

    def gradient(z):
        return (a @ z.reshape(len(free), h) + c).ravel()

    constraints = ()
    if coupled:
        rows = np.kron(np.eye(len(free)), np.ones((1, h)))
        constraints = (LinearConstraint(rows, np.ones(len(free)), np.ones(len(free))),)
    result = minimize(objective, initial.ravel(), jac=gradient, method="SLSQP",
                      bounds=Bounds(0., 1.), constraints=constraints,
                      callback=lambda z: history.append(objective(z)),
                      options={"maxiter": 1000, "ftol": 1e-13, "disp": False})
    raw = data["fixed"].copy()
    raw[free] = result.x.reshape(len(free), h)
    pre_polish = kkt(data, raw, coupled)
    if polish:
        # Independently solve the KKT equations for the active set inferred
        # from this actual SLSQP iterate. No altered objective or silent clamp:
        # the resulting candidate must pass original bounds/dual conditions.
        z = result.x.copy()
        active_low = z <= 1e-10
        active_high = z >= 1 - 1e-10
        active = active_low | active_high
        inactive = ~active
        fixed_z = np.zeros_like(z)
        fixed_z[active_high] = 1.
        hessian = np.kron(a, np.eye(h))
        b = -(c.ravel() + hessian @ fixed_z)[inactive]
        restricted = hessian[np.ix_(inactive, inactive)]
        if coupled:
            e = rows[:, inactive]
            rhs = np.ones(len(free)) - rows @ fixed_z
            keep = np.any(e != 0, axis=1)
            assert np.max(np.abs(rhs[~keep]), initial=0) < 1e-12
            e, rhs = e[keep], rhs[keep]
            system = np.block([[restricted, e.T], [e, np.zeros((len(rhs), len(rhs)))]])
            correction = np.linalg.solve(system, np.concatenate((b, rhs)))[:np.count_nonzero(inactive)]
        else:
            correction = np.linalg.solve(restricted, b)
        fixed_z[inactive] = correction
        raw[free] = fixed_z.reshape(len(free), h)
        refined = kkt(data, raw, coupled)
        assert refined["bounds"] < 1e-12 and refined["projected_stationarity"] < 1e-12, refined
        if coupled:
            assert refined["partition"] < 1e-12, refined
    final = raw.copy() if coupled else raw / raw.sum(axis=1)[:, None]
    assessment = kkt(data, raw, coupled)
    assert result.success, result.message
    # Independent optimizer accuracy is disclosed separately. Native strict
    # acceptance cannot be lowered because this comparison is less accurate.
    assert assessment["projected_stationarity"] < 3e-7, assessment
    if data["cells_per_axis"] == 2:
        expected = np.full(8, 1 / 16)
        expected[[0, 7]] = 5 / 16
        assert np.max(np.abs(raw[free[0]] - expected)) < 1e-7
    return dict(scope="INDEPENDENT_SYNTHETIC_REFERENCE_NOT_NATIVE_VERIFICATION",
                rule=data["rule"], cells_per_axis=data["cells_per_axis"],
                protocol="COUPLED_PARTITION" if coupled else "INDEPENDENT_NORMALIZED",
                nodes=data["nodes"].tolist(), handle_nodes=data["handles"],
                backend="scipy.optimize.SLSQP", iterations=int(result.nit),
                options={"maxiter": 1000, "ftol": 1e-13},
                active_set_polish=polish, pre_polish_assessment=pre_polish,
                objective_history=history, raw_weights=raw.tolist(), final_weights=final.tolist(),
                raw_assessment=assessment, final_coupled_assessment=kkt(data, final, True),
                normalization_change=float(np.max(np.abs(final - raw))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", type=int, choices=(2, 4), default=2)
    parser.add_argument("--rule", choices=("RESTRICTED_GRID_VORONOI", "BARYCENTRIC_TET"),
                        default="RESTRICTED_GRID_VORONOI")
    parser.add_argument("--protocol", choices=("COUPLED_PARTITION", "INDEPENDENT_NORMALIZED"),
                        default="COUPLED_PARTITION")
    parser.add_argument("--polish", action="store_true")
    args = parser.parse_args()
    data = reference(args.cells, args.rule)
    print(json.dumps(solve_reference(data, args.protocol == "COUPLED_PARTITION", args.polish), sort_keys=True, allow_nan=False))
