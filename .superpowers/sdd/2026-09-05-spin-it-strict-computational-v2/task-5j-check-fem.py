"""Independent rational small-grid FEM oracle; no project solver imports.

Default: derive exact reference for [1/4,1/2]^3 subdivided 2^3.
Optional JSON: compare actual public native API export against that reference.
This is SYNTHETIC_REGRESSION, not BBW minimization or paper-case evidence.
"""
import itertools
import json
import math
import sys
from collections import Counter
from fractions import Fraction as F

errors = []


def inverse(a):
    n = len(a)
    rows = [list(a[i]) + [F(i == j) for j in range(n)] for i in range(n)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if rows[i][j])
        rows[j], rows[pivot] = rows[pivot], rows[j]
        divisor = rows[j][j]
        rows[j] = [x / divisor for x in rows[j]]
        for i in range(n):
            if i != j:
                factor = rows[i][j]
                rows[i] = [x - factor * y for x, y in zip(rows[i], rows[j])]
    return [r[n:] for r in rows]


def determinant3(a):
    return sum(a[0][j] *
               (a[1][(j + 1) % 3] * a[2][(j + 2) % 3] -
                a[1][(j + 2) % 3] * a[2][(j + 1) % 3]) for j in range(3))


def reference(rule="BARYCENTRIC_TET"):
    keys = list(itertools.product(range(3), repeat=3))
    nodes = [[F(1, 4) + F(k, 8) for k in key] for key in keys]
    indices = {k: i for i, k in enumerate(keys)}
    n = len(nodes)
    stiffness = [[F(0) for _ in range(n)] for _ in range(n)]
    mass = [F(0) for _ in range(n)]
    tets = []
    face_count = Counter()
    for cell in itertools.product(range(2), repeat=3):
        for order in itertools.permutations(range(3)):
            path = [cell]
            for axis in order:
                next_key = list(path[-1])
                next_key[axis] += 1
                path.append(tuple(next_key))
            ids = [indices[k] for k in path]
            pts = [nodes[i] for i in ids]
            edge = [[pts[i + 1][a] - pts[0][a] for a in range(3)] for i in range(3)]
            volume = abs(determinant3(edge)) / 6
            assert volume == F(1, 3072)
            # Solve four affine interpolation equations exactly (independent
            # Gauss-Jordan path, not the production element-gradient routine).
            coefficients = inverse([[F(1)] + p for p in pts])
            gradients = [[coefficients[a + 1][j] for a in range(3)] for j in range(4)]
            for a, i in enumerate(ids):
                mass[i] += volume / 4
                for b, j in enumerate(ids):
                    stiffness[i][j] += volume * sum(x * y for x, y in zip(gradients[a], gradients[b]))
            tets.append(tuple(sorted(ids)))
            for face in itertools.combinations(ids, 3):
                face_count[tuple(sorted(face))] += 1
    assert len(tets) == 48 and sum(mass) == F(1, 64)
    assert set(face_count.values()) == {1, 2}
    boundary = sorted({i for face, count in face_count.items() if count == 1 for i in face})
    assert len(boundary) == 26 and set(range(27)) - set(boundary) == {13}
    assert all(sum(row) == 0 for row in stiffness)
    if rule == "RESTRICTED_GRID_VORONOI":
        # Independent integration of Cartesian Voronoi boxes clipped at the
        # exterior boundary, not production's per-brick V/8 accumulation.
        widths = {0: F(1, 16), 1: F(1, 8), 2: F(1, 16)}
        mass = [math.prod(widths[k] for k in key) for key in keys]
        assert sum(mass) == F(1, 64)
    else:
        assert rule == "BARYCENTRIC_TET"
    q = [[sum(stiffness[k][i] * stiffness[k][j] / mass[k] for k in range(n))
          for j in range(n)] for i in range(n)]
    field = [F((i * 7) % 11 - 5, 8) for i in range(n)]
    action = [sum(a * b for a, b in zip(row, field)) for row in q]
    energy = sum(a * b for a, b in zip(field, action)) / 2
    assert energy > 0
    # Piecewise triangular linear boundary rows for full cube corner handles.
    corners = list(itertools.product((0, 2), repeat=3))
    constraints = {}
    for i in boundary:
        key = keys[i]
        candidates = []
        for axis in range(3):
            if key[axis] not in (0, 2):
                continue
            axes = [a for a in range(3) if a != axis]
            u, v = [F(key[a], 2) for a in axes]
            values = {(0, 0): 1 - max(u, v), (2, 2): min(u, v),
                      (2, 0): max(u - v, 0), (0, 2): max(v - u, 0)}
            row = [F(0)] * 8
            for uv, value in values.items():
                corner = [0] * 3
                corner[axis] = key[axis]
                for a, coordinate in zip(axes, uv):
                    corner[a] = coordinate
                row[corners.index(tuple(corner))] = value
            candidates.append(row)
        assert all(row == candidates[0] for row in candidates)
        constraints[i] = candidates[0]
    return nodes, sorted(tets), mass, stiffness, field, action, energy, constraints, corners


def close(got, expected, what):
    assert math.isfinite(got), (what, "nonfinite")
    errors.append(abs(got - float(expected)))
    assert abs(got - float(expected)) <= 1e-11 * max(1., abs(float(expected))), (what, got, str(expected))


def check(path):
    errors.clear()
    with open(path) as stream:
        data = json.load(stream)
    ref = reference(data["mass_rule"])
    nodes, tets, mass, stiffness, field, action, energy, constraints, corners = ref
    assert data["scope"] == "SYNTHETIC_FEM_NOT_BBW_SOLVE"
    for name in ("config_sha256", "output_sha256"):
        assert len(data[name]) == 64 and all(c in "0123456789abcdef" for c in data[name])
    assert data["nodes"] == [[float(x) for x in row] for row in nodes]
    assert sorted(tuple(sorted(t)) for t in data["tetrahedra"]) == tets
    assert len(data["mass"]) == len(mass)
    for i, expected in enumerate(mass):
        close(data["mass"][i], expected, ("mass", i))
    assert len(data["stiffness"]) == len(stiffness)
    for i, row in enumerate(stiffness):
        assert len(data["stiffness"][i]) == len(row)
        for j, expected in enumerate(row):
            close(data["stiffness"][i][j], expected, ("K", i, j))
    assert data["field"] == [float(x) for x in field]
    assert len(data["action"]) == len(action)
    for i, expected in enumerate(action):
        close(data["action"][i], expected, ("Qx", i))
    close(data["energy"], energy, "energy")
    # Map actual handle ordering through coordinates, not expected row numbers.
    handle_keys = [tuple(round((x - .25) * 8) for x in nodes[i]) for i in data["handle_nodes"]]
    assert sorted(handle_keys) == corners
    rows = {row["node"]: row["weights"] for row in data["boundary"]}
    assert len(rows) == len(data["boundary"]) == 26 and set(rows) == set(constraints)
    for i, expected in constraints.items():
        assert len(rows[i]) == 8
        for j, key in enumerate(handle_keys):
            close(rows[i][j], expected[corners.index(key)], ("boundary", i, j))
    print(json.dumps({"file": path, "passed": True, "mass_rule": data["mass_rule"], "nodes": 27, "tetrahedra": 48,
                      "boundary_nodes": 26, "exact_energy": str(energy), "max_absolute_error": max(errors),
                      "scope": "SYNTHETIC_REGRESSION"}, sort_keys=True))


if __name__ == "__main__":
    if len(sys.argv) == 1:
        for rule in ("BARYCENTRIC_TET", "RESTRICTED_GRID_VORONOI"):
            oracle = reference(rule)
            print(json.dumps({"reference_only": True, "mass_rule": rule, "nodes": 27, "tetrahedra": 48,
                              "volume": str(sum(oracle[2])), "energy": str(oracle[6]),
                              "interior_mass": str(oracle[2][13]),
                              "corner_mass": str(oracle[2][0])}, sort_keys=True))
    for filename in sys.argv[1:]:
        check(filename)
    if len(sys.argv) == 5:
        observed = []
        for path in sys.argv[1:]:
            with open(path) as stream:
                observed.append(json.load(stream))
        assert observed[0] == observed[1] and observed[2] == observed[3]
        assert observed[0]["mass_rule"] != observed[2]["mass_rule"]
        for name in ("config_sha256", "output_sha256", "mass", "energy"):
            assert observed[0][name] != observed[2][name], name
        for name in ("stiffness", "boundary", "tetrahedra", "nodes", "handle_nodes"):
            assert observed[0][name] == observed[2][name], name
        print("PASS mass rules remain distinct; geometry, stiffness and boundary constraints agree")
