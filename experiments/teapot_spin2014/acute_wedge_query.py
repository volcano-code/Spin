"""Conditional material-ball-cut radius infimum at a convex acute edge.

This is a one-sided continuum certificate, not a complete epsilon-map or a
manufacturing wall metric. A separately authenticated embedded, outward-material
boundary is required. Local topology checks here do not replace that audit.

For midpoint m on an acute edge, centers x(t)=m+t*u on one incident triangle
and projections q(t)=m+t*w on the other give two disjoint boundary disks in a
ball of squared radius t^2*r. Exact inequalities certify this for ALL
0<t<=t_max, not merely for a sample. The convex material chord identifies the
middle ball region as material, giving the local ball-cut Euler change -1.
Thus the material ball-cut radius infimum is zero near this edge. We name the
operator explicitly: the reference's formal +/- map definitions and its later
material-side terminology do not agree. This certificate concerns subtraction
of an OPEN ball, never union with a ball, and never modifies the source mesh.
"""
from fractions import Fraction as Q

from .boundary_cap_query import BoundaryCapQuery
from .near_feature_witness import add, sub, mul, dot, cross, point_triangle


def _strictly_inside(point, triangle, normal):
    return (dot(normal, sub(point, triangle[0])) == 0 and
            all(dot(cross(sub(triangle[(i+1) % 3], triangle[i]),
                          sub(point, triangle[i])), normal) > 0 for i in range(3)))


class AcuteWedgeQuery:
    def __init__(self, caps):
        if not isinstance(caps, BoundaryCapQuery):
            raise TypeError('computed closed-manifold cap query required')
        self.caps = caps

    def _edge(self, edge):
        if (not isinstance(edge, (list, tuple)) or len(edge) != 2 or
                any(type(i) is not int for i in edge) or edge[0] == edge[1]):
            raise ValueError('two distinct integer edge vertices required')
        edge = tuple(sorted(edge))
        if edge not in self.caps.edges:
            raise ValueError('source edge is absent')
        return edge

    def _acute(self, edge):
        b = self.caps.boundary
        a, c = self.caps.edges[edge]
        ta, na = b._triangle(a)
        tc, nc = b._triangle(c)
        va = next(int(v) for v in b.faces[a] if int(v) not in edge)
        vc = next(int(v) for v in b.faces[c] if int(v) not in edge)
        return (dot(na, nc) < 0 and
                dot(na, sub(self.caps.points[vc], ta[0])) < 0 and
                dot(nc, sub(self.caps.points[va], tc[0])) < 0)

    def convex_acute_edges(self):
        # Exact normals and side predicates on EVERY edge, no float prefilter.
        return [list(e) for e in sorted(self.caps.edges) if self._acute(e)]

    def certificate(self, edge, threshold_squared, *, maximum_halvings):
        edge = self._edge(edge)
        if isinstance(threshold_squared, bool):
            raise ValueError('boolean radius')
        limit = Q(threshold_squared)
        if limit <= 0:
            raise ValueError('positive squared radius required')
        if type(maximum_halvings) is not int or maximum_halvings <= 0:
            raise ValueError('positive integer halving limit required')
        row = {'schema': 'conditional-convex-acute-wedge-v1', 'edge': list(edge),
               'operator': 'MATERIAL_MINUS_OPEN_BALL',
               'threshold_squared_exact': str(limit), 'maximum_halvings': maximum_halvings,
               'independent_embedded_material_audit_required': True,
               'all_surface_centers_checked': False, 'wall_metric_adopted': False,
               'complete_wall_thickness': 'UNKNOWN', 'wall_thickness_admission': 'NOT_GRANTED'}
        if not self._acute(edge):
            return {**row, 'state': 'NOT_CONVEX_ACUTE'}
        b = self.caps.boundary
        a, c = self.caps.edges[edge]
        ta, na = b._triangle(a)
        tc, nc = b._triangle(c)
        start, end = (self.caps.points[i] for i in edge)
        direction = sub(end, start)
        midpoint = mul(add(start, end), Q(1, 2))
        opposite = next(self.caps.points[int(v)] for v in b.faces[a] if int(v) not in edge)
        into_source = sub(opposite, start)
        u = sub(into_source, mul(direction, dot(into_source, direction) / dot(direction, direction)))
        w = sub(u, mul(nc, dot(nc, u) / dot(nc, nc)))
        u2 = dot(u, u)
        distance_coefficient = dot(sub(u, w), sub(u, w))
        radius_coefficient = (distance_coefficient + u2) / 2
        if not 0 < distance_coefficient < radius_coefficient < u2:
            raise ValueError('acute wedge metric inequalities failed')
        # Positive distance at a relative-interior edge point from every other
        # closed triangle. This includes all nonshared triangle edges/endpoints.
        clearance, clearance_face = None, None
        for i in range(len(b.faces)):
            if i in (a, c):
                continue
            triangle, normal = b._triangle(i)
            delta = sub(midpoint, point_triangle(midpoint, triangle, normal))
            candidate = dot(delta, delta)
            if clearance is None or candidate < clearance:
                clearance, clearance_face = candidate, i
        if clearance is None or clearance <= 0:
            raise ValueError('edge patch not isolated; embedded material precondition missing')

        t = Q(1)
        for halvings in range(maximum_halvings + 1):
            x, q = add(midpoint, mul(u, t)), add(midpoint, mul(w, t))
            if (_strictly_inside(x, ta, na) and _strictly_inside(q, tc, nc) and
                    dot(na, sub(q, x)) < 0 and dot(nc, sub(x, q)) < 0 and
                    4*t*t*u2 < clearance and t*t*radius_coefficient < limit):
                # Convexity preserves strict face containment for 0<t<=t_max.
                # |x-m|+r < 2t|u| < distance(m, every other triangle), so the
                # whole ball sees only these two planes. r<t|u| keeps it off the
                # shared edge, while r>distance(x,target plane) reaches both.
                return {**row, 'state': 'CONDITIONAL_MATERIAL_BALL_CUT_INFIMUM_ZERO',
                        'incident_faces': [a, c], 'edge_midpoint_exact': list(map(str, midpoint)),
                        'source_direction_exact': list(map(str, u)),
                        'target_direction_exact': list(map(str, w)),
                        'maximum_scale_exact': str(t), 'halvings': halvings,
                        'all_scales_domain': '0 < t <= maximum_scale_exact',
                        'squared_source_to_edge_coefficient_exact': str(u2),
                        'squared_source_to_target_coefficient_exact': str(distance_coefficient),
                        'radius_squared_coefficient_exact': str(radius_coefficient),
                        'other_faces_clearance_squared_exact': str(clearance),
                        'other_faces_clearance_face': clearance_face,
                        'ball_isolation_upper_squared_exact': str(4*t*t*u2),
                        'material_ball_cut_radius_infimum_exact': '0', 'local_cut_euler_change': -1}
            t /= 2
        return {**row, 'state': 'RESOURCE_LIMIT',
                'reason': 'no admissible rational scale within halving budget'}
