"""Independent rational geometry oracle for exported synthetic source audits.

No native solver/checker imports. Distance via simplex features and explicit
segment/triangle intersection; moments via signed origin tetrahedra. This is
independent verification of the exported fixtures, not arbitrary mesh admission.
"""
from fractions import Fraction as F
from functools import lru_cache
import json
import math
import sys


def vec(x):
    assert all(math.isfinite(float(a)) for a in x)
    return tuple(F(float(a)) for a in x)


def sub(a, b): return tuple(x-y for x, y in zip(a, b))
def add(a, b): return tuple(x+y for x, y in zip(a, b))
def scale(a, s): return tuple(x*s for x in a)
def dot(a, b): return sum(x*y for x, y in zip(a, b))
def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def norm2(a): return dot(a, a)
def clamp(a): return max(F(0), min(F(1), a))


def point_segment(p, a, b):
    e = sub(b, a)
    assert norm2(e) > 0
    return norm2(sub(p, add(a, scale(e, clamp(dot(sub(p, a), e)/norm2(e))))))


def inside_triangle_projection(p, tri):
    a, b, c = tri
    u, v, w = sub(b, a), sub(c, a), sub(p, a)
    uu, uv, vv = dot(u, u), dot(u, v), dot(v, v)
    den = uu*vv-uv*uv
    assert den > 0
    s, t = (dot(w, u)*vv-dot(w, v)*uv)/den, (dot(w, v)*uu-dot(w, u)*uv)/den
    return s >= 0 and t >= 0 and s+t <= 1


def point_triangle(p, tri):
    n = cross(sub(tri[1], tri[0]), sub(tri[2], tri[0]))
    assert norm2(n) > 0
    if inside_triangle_projection(p, tri):
        return dot(sub(p, tri[0]), n)**2/norm2(n)
    return min(point_segment(p, tri[i], tri[(i+1)%3]) for i in range(3))


def segment_segment(a, b, c, d):
    u, v, w = sub(b, a), sub(d, c), sub(a, c)
    aa, bb, cc, dd, ee = dot(u, u), dot(u, v), dot(v, v), dot(u, w), dot(v, w)
    candidates = [point_segment(a, c, d), point_segment(b, c, d),
                  point_segment(c, a, b), point_segment(d, a, b)]
    den = aa*cc-bb*bb
    if den:
        s, t = (bb*ee-cc*dd)/den, (aa*ee-bb*dd)/den
        if 0 <= s <= 1 and 0 <= t <= 1:
            candidates.append(norm2(sub(add(a, scale(u, s)), add(c, scale(v, t)))))
    return min(candidates)


def pierces(a, b, tri):
    n = cross(sub(tri[1], tri[0]), sub(tri[2], tri[0]))
    da, db = dot(sub(a, tri[0]), n), dot(sub(b, tri[0]), n)
    if da == db:
        return da == 0 and (inside_triangle_projection(a, tri) or inside_triangle_projection(b, tri))
    t = da/(da-db)
    return 0 <= t <= 1 and inside_triangle_projection(add(a, scale(sub(b, a), t)), tri)


def triangle_distance(a, b):
    if any(pierces(a[i], a[(i+1)%3], b) or pierces(b[i], b[(i+1)%3], a) for i in range(3)):
        return F(0)
    return min([point_triangle(p, b) for p in a]+[point_triangle(p, a) for p in b]+
               [segment_segment(a[i], a[(i+1)%3], b[j], b[(j+1)%3]) for i in range(3) for j in range(3)])


def boundary_distance(points, fa, fb):
    return _boundary_distance(tuple(points), tuple(map(tuple, fa)), tuple(map(tuple, fb)))


@lru_cache(maxsize=256)
def _boundary_distance(points, fa, fb):
    best = None
    for ta in fa:
        a = [points[i] for i in ta]
        amin = tuple(min(p[k] for p in a) for k in range(3))
        amax = tuple(max(p[k] for p in a) for k in range(3))
        for tb in fb:
            b = [points[i] for i in tb]
            lower = sum(max(F(0), amin[k]-max(p[k] for p in b), min(p[k] for p in b)-amax[k])**2 for k in range(3))
            if best is not None and lower >= best: continue
            d = triangle_distance(a, b)
            best = d if best is None else min(best, d)
            if best == 0: return best
    assert best is not None
    return best


def moments(points, faces):
    edges = {}
    result = [F(0)]*10
    for face in faces:
        assert len(set(face)) == 3 and all(0 <= i < len(points) for i in face)
        a, b, c = [points[i] for i in face]
        assert norm2(cross(sub(b, a), sub(c, a))) > 0
        for i in range(3):
            u, v = face[i], face[(i+1)%3]
            edges.setdefault(tuple(sorted((u, v))), []).append((u, v))
        vol = dot(a, cross(b, c))/6
        total = add(add(a, b), c)
        q = [[vol*(total[i]*total[j]+sum(p[i]*p[j] for p in (a,b,c)))/20 for j in range(3)] for i in range(3)]
        part = [vol]+[vol*x/4 for x in total]+[q[0][0], q[1][1], q[2][2], q[1][2], q[0][2], q[0][1]]
        result = [x+y for x,y in zip(result, part)]
    assert all(len(e)==2 and e[0] == e[1][::-1] for e in edges.values())
    return result


def self_test():
    tri = [vec(x) for x in [(0,0,0),(1,0,0),(0,1,0)]]
    assert point_triangle(vec((1,1,1)), tri) == F(3,2)
    assert triangle_distance(tri, [add(p, vec((0,0,2))) for p in tri]) == 4
    assert triangle_distance(tri, [vec(x) for x in [(.25,.25,-1),(.25,.25,1),(.75,.25,0)]]) == 0
    assert segment_segment(vec((0,0,0)),vec((1,0,0)),vec((0,1,1)),vec((1,1,1))) == 2


def verify(c):
    assert c['scope'] == 'SYNTHETIC_REGRESSION' and c['feasibility'] == 'UNKNOWN'
    assert c['sealed'] is False and c['converged'] is False
    assert c['id'] in ('identity', 'affine', 'beta_phi')
    assert c['status'] == 'CHECKS_PASSED'
    points = [vec(p) for p in c['points']]
    rest = [vec(p) for p in c['rest_points']]
    assert len(rest) == len(points) and len(c['state']) == 113 and len(c['beta']) == 8
    factor = F(2) if c['id'] == 'affine' else F(1)
    offset = vec((-.125,.25,.5)) if c['id'] == 'affine' else vec((0,0,0))
    assert max(abs(float(a-b)) for p,r in zip(points,rest) for a,b in zip(p,add(scale(r,factor),offset))) <= 1e-12
    beta = [F(float(x)) for x in c['beta']]
    outer = moments(points, c['outer_faces'])
    assert outer[0] > 0
    total = outer[:]
    for i, p in enumerate(c['primitives']):
        assert p['id'] == i and 0 <= p['owner'] < len(beta)
        local = moments(points, p['faces'])
        assert local[0] > 0
        distance = boundary_distance(points, c['outer_faces'], p['faces'])
        assert distance == F(p['distance_squared_exact']), (c['id'], i, distance, p['distance_squared_exact'])
        assert p['stage'] == 'PASSED'
        assert distance > 0 and distance >= F(float(c['minimum_clearance']))**2
        total = [a-beta[p['owner']]*b for a,b in zip(total,local)]
    measured = [F(float(x)) for x in c['world_moments']]
    expected = [F(float(c['density']))*x for x in total]
    error = max(abs(float(a-b)) for a,b in zip(measured,expected))/max(1e-30,max(abs(float(x)) for x in expected))
    assert error <= 1e-10, error
    assert len(c['primitives']) == 8
    n = len(c['primitives'])
    assert c['pair_count'] == c['pairs_visited'] == n*(n-1)//2
    pairs = [(e['first'], e['second']) for e in c['events'] if e['stage'] == 6]
    assert pairs == [(i,j) for i in range(n) for j in range(i+1,n)]
    assert len(c['events']) == 3+3*n+n*(n-1)//2
    absent = (1 << 64)-1
    expected_events = [(i, absent, absent) for i in (0,1,2)]
    expected_events += [(stage, i, absent) for i in range(n) for stage in (3,4,5)]
    expected_events += [(6,i,j) for i in range(n) for j in range(i+1,n)]
    assert [(e['stage'],e['first'],e['second']) for e in c['events']] == expected_events
    assert all(e['outcome'] == 1 if e['stage'] != 6 else e['outcome'] in (1,5) for e in c['events'])
    print(f"PASS {c['id']}: {len(c['primitives'])} exact clearances, independent moments error={error:.3e}; no seal")


if __name__ == '__main__':
    self_test()
    if len(sys.argv) == 1:
        print('PASS independent distance primitive self-tests')
    else:
        cases = [json.loads(line) for line in open(sys.argv[1], encoding='utf-8')]
        assert len(cases) == 3 and {c['id'] for c in cases} == {'identity', 'affine', 'beta_phi'}
        for case in cases: verify(case)
