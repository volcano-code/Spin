"""Independent synthetic binding check; no project/native module imports."""
import ctypes
import copy
import itertools
import sys
from collections import defaultdict
import numpy as np

fma = ctypes.CDLL("libm.so.6").fma
fma.argtypes = [ctypes.c_double] * 3
fma.restype = ctypes.c_double

def box_moments(lo, hi):
    volume = np.prod(hi - lo)
    mean = (hi + lo) / 2
    second = (hi * hi + hi * lo + lo * lo) / 3
    return volume * np.r_[1, mean, second, mean[1]*mean[2], mean[0]*mean[2], mean[0]*mean[1]]

def indexed(rows, count):
    assert len(rows) == count
    assert [int(r[0]) for r in rows] == list(range(count)), "missing/duplicate/out-of-order id"
    a = np.array(rows, dtype=float)
    assert np.isfinite(a).all()
    return a

def verify(c):
    fixture, protocol, depth, nl, np_, nq, nh, nn, nc, no = map(int, c["CASE"][0])
    assert nh == 26 and no == 8
    root = np.array(c["ROOT"][0], dtype=float)
    origin, side = root[:3], root[3]
    points = indexed(c["POINT"], nq)[:, 1:]
    outer = indexed(c["OUTER"], no)[:, 1:]
    assert points[:no].tobytes() == outer.tobytes(), "outer order/bytes changed"
    leaf = indexed(c["LEAF"], nl)
    prim = indexed(c["PRIM"], np_)
    nodes = indexed(c["NODE"], nn)
    cells = indexed(c["GRIDCELL"], nc)
    embeddings = indexed(c["EMBED"], nq)
    weights = indexed(c["WEIGHT"], nq)[:, 1:]
    keys = indexed(c["KEY"], nq-no)[:, 1:].astype(int)
    lo_key = tuple(map(int, c["EXPECTED_LO"][0]))
    hi_key = tuple(map(int, c["EXPECTED_HI"][0]))
    planes = lambda key: np.array([fma(float(key[a])/(1<<depth), side, origin[a]) for a in range(3)])
    sum_by_leaf = np.zeros((nl, 10))
    all_corners, occupied = set(), set()
    max_error = 0.
    owners = {(int(l[1]), *map(int,l[2:5])): int(l[0]) for l in leaf}
    for p in prim:
        owner, level = map(int, p[1:3])
        key = tuple(map(int, p[3:6]))
        ancestors = [owners[(l, *(k >> (level-l) for k in key))]
                     for l in range(level+1) if (l, *(k >> (level-l) for k in key)) in owners]
        assert ancestors == [owner], "material beta owner mismatch"
        low = tuple(k << (depth-level) for k in key)
        high = tuple((k+1) << (depth-level) for k in key)
        assert np.array_equal(p[6:9], planes(low)) and np.array_equal(p[9:12], planes(high))
        for unit in itertools.product(*(range(low[a], high[a]) for a in range(3))):
            assert unit not in occupied, "overlapping source Full primitive"
            occupied.add(unit)
        for b, point_id in enumerate(p[12:].astype(int)):
            corner = tuple(high[a] if b & (1<<a) else low[a] for a in range(3))
            all_corners.add(corner)
            assert point_id >= no and np.array_equal(points[point_id], planes(corner))
        sum_by_leaf[owner] += box_moments(p[6:9], p[9:12])
    assert occupied == set(itertools.product(*(range(lo_key[a],hi_key[a]) for a in range(3))))
    assert list(map(tuple, keys)) == sorted(all_corners)
    assert np.array_equal(points[no:], np.array([planes(k) for k in keys]))
    moment_error = np.max(np.abs(sum_by_leaf-leaf[:,5:]))
    analytic = box_moments(planes(lo_key), planes(hi_key))
    total_error = np.max(np.abs(sum_by_leaf.sum(axis=0)-analytic))
    beta = np.linspace(.13, .77, nl)
    weighted_error = np.max(np.abs(beta@sum_by_leaf-beta@leaf[:,5:]))
    assert max(moment_error,total_error,weighted_error) < 1e-10
    # Search all actual small fixture grid cells independently. Closed incidents
    # choose lexicographically smallest key, rather than trusting embed.cellIndex.
    cell_nodes = cells[:,4:].astype(int)
    cell_lows = nodes[cell_nodes[:,0],1:4]
    cell_highs = nodes[cell_nodes[:,7],1:4]
    for i, e in enumerate(embeddings):
        incident = np.flatnonzero(np.all(points[i]>=cell_lows,axis=1) & np.all(points[i]<=cell_highs,axis=1))
        selected = min(incident, key=lambda j: tuple(cells[j,1:4]))
        assert int(e[1]) == selected and np.array_equal(e[2:5],cells[selected,1:4])
        assert np.array_equal(e[5:13],cell_nodes[selected])
        t = (points[i]-cell_lows[selected])/(cell_highs[selected]-cell_lows[selected])
        alpha = np.array([np.prod([t[a] if b&(1<<a) else 1-t[a] for a in range(3)]) for b in range(8)])
        assert np.max(np.abs(alpha-e[13:])) < 1e-13
        expected = alpha@nodes[cell_nodes[selected],4:]
        max_error = max(max_error,np.max(np.abs(expected-weights[i])))
    assert max_error < 1e-13
    assert np.max(np.abs(weights.sum(axis=1)-1)) < 1e-8
    assert weights.min() >= -1e-8 and weights.max() <= 1+1e-8
    if fixture == 0:
        assert nl == 8 and np_ == 496 and analytic[0] == .030517578125
    print(f"PASS fixture={fixture} protocol={protocol} leaves={nl} primitives={np_} points={nq} handles={nh} moment_error={moment_error:.3e} total_error={total_error:.3e} weighted_error={weighted_error:.3e} weight_error={max_error:.3e}")

count = 0
with open(sys.argv[1], encoding="utf-8") as stream:
    c = None
    for line in stream:
        fields = line.split()
        assert fields, "unexpected blank record"
        if fields[0] == "CASE":
            assert c is None
            c = defaultdict(list)
        assert c is not None
        if fields[0] == "END":
            verify(c)
            if "--negative-checks" in sys.argv[2:]:
                for target, record, column in [("PRIM",0,1),("POINT",0,1),("WEIGHT",0,1),("LEAF",0,5),("EMBED",0,13),("KEY",0,1)]:
                    bad = copy.deepcopy(c)
                    old = float(bad[target][record][column])
                    bad[target][record][column] = str(old + 1)
                    try:
                        verify(bad)
                    except AssertionError:
                        pass
                    else:
                        raise AssertionError(f"checker accepted corrupted {target}")
                print("PASS six independently corrupted export records rejected")
            count += 1; c = None
        else:
            c[fields[0]].append(fields[1:])
assert c is None and count == 4
print("SYNTHETIC_REGRESSION: binding/interpolation only; not deformation geometry or paper reproduction")
