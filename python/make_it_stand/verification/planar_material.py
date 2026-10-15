"""Independent integration of saved planar line/arc extrusions.

This checker does not call native polygon/disk integration or optimize again.
It uses Green boundary integrals evaluated at three Gauss--Legendre orders,
then integrates the extrusion coordinate analytically. Agreement of quadrature
orders is numerical evidence, not a rigorous interval error bound. Boundary
closure and parameter consistency are checked; simplicity, nesting, Boolean
identity, inter-material overlap and wall thickness are NOT certified here.
Thus matching moments never creates a SealedSpinArtifact or paper approval.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np

from .spin_common import ORDER, digest, parse, require, scalar

SCOPE = "SAVED_PLANAR_MOMENTS_AND_BALANCE_NOT_FULL_MATERIAL_OR_PAPER_APPROVAL"
ORDERS = (32, 64, 128)
EXPONENTS = ((0, 0), (1, 0), (0, 1), (2, 0), (0, 2), (1, 1))


def _vector(value, size, label):
    require(isinstance(value, (list, tuple)) and len(value) == size,
            label + ": wrong vector size")
    return [scalar(x, label) for x in value]


@lru_cache(maxsize=3)
def _rule(order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return .5*(nodes+1.), .5*weights


def _cycles(entity):
    representation = entity.get("geometry_representation")
    if representation == "exact_polygonal_prism":
        loops = entity.get("loops")
        require(isinstance(loops, list), "polygon loops missing")
        cycles = []
        for loop in loops:
            require(isinstance(loop, list) and len(loop) >= 3, "polygon loop is incomplete")
            points = [_vector(point, 2, "polygon vertex") for point in loop]
            cycles.append([{"type": "LINE_SEGMENT", "start": a, "end": b}
                           for a, b in zip(points, points[1:]+points[:1])])
        return cycles
    require(representation == "exact_line_arc_profile_swept_between_planes",
            "unsupported saved planar entity representation")
    cycles = entity.get("boundary_cycles")
    require(isinstance(cycles, list), "line/arc boundary cycles missing")
    return cycles


def integrate_planar_entity(entity, *, origin, length_scale):
    """Return independent density-weighted 3D moments, preserving saved curves.

    origin and length_scale must come from the input, not the optimized object.
    The native representation names contain 'exact'; these saved parameters are
    still binary64, not authoritative rational coordinates or certified CAD.
    """
    origin = np.asarray(_vector(origin, 2, "input origin"))
    length = scalar(length_scale, "input length", positive=True)
    density = scalar(entity.get("volume_density"), "volume density", positive=True)
    low, high = _vector(entity.get("z_interval"), 2, "extrusion interval")
    require(high > low, "extrusion interval must have positive height")
    cycles = _cycles(entity)
    require(cycles or entity.get("role") == "editable_design", "fixed material cannot be empty")
    prepared, maximum_gap = [], 0.
    for cycle in cycles:
        require(isinstance(cycle, list) and len(cycle) > 0, "empty boundary cycle")
        edges = []
        for row in cycle:
            start = np.asarray(_vector(row.get("start"), 2, "boundary start"))
            end = np.asarray(_vector(row.get("end"), 2, "boundary end"))
            kind = row.get("type")
            if kind == "LINE_SEGMENT":
                require(np.any(end != start), "zero-length boundary segment")
                edges.append((kind, start, end, None))
            else:
                require(kind == "CIRCULAR_ARC", "unsupported boundary primitive")
                center = np.asarray(_vector(row.get("center"), 2, "arc center"))
                radius = scalar(row.get("radius"), "arc radius", positive=True)
                angle = scalar(row.get("start_angle"), "arc start angle")
                sweep = scalar(row.get("sweep"), "arc sweep")
                require(0 < abs(sweep) <= 2*math.pi*(1+16*np.finfo(float).eps)
                        and abs(angle) <= 4*math.pi, "unsupported arc angular range")
                for point, theta in ((start, angle), (end, angle+sweep)):
                    expected = center + radius*np.array([math.cos(theta), math.sin(theta)])
                    gap = float(np.linalg.norm(point-expected))
                    maximum_gap = max(maximum_gap, gap)
                    envelope = 64*np.finfo(float).eps*max(length, radius, float(np.max(np.abs(center))))
                    require(gap <= envelope, "arc endpoint differs from saved parameters")
                edges.append((kind, start, end, (center, radius, angle, sweep)))
        for edge, following in zip(edges, edges[1:]+edges[:1]):
            gap = float(np.linalg.norm(edge[2]-following[1]))
            maximum_gap = max(maximum_gap, gap)
            envelope = 64*np.finfo(float).eps*max(length,
                float(np.max(np.abs(edge[2]))), float(np.max(np.abs(following[1]))))
            require(gap <= envelope, "boundary closure mismatch")
        prepared.extend(edges)

    approximations = []
    for order in ORDERS:
        t, weights = _rule(order)
        contributions = []
        for kind, start, end, arc in prepared:
            if kind == "LINE_SEGMENT":
                points = (start[None, :]-origin + t[:, None]*(end-start))/length
                dy = np.full(order, (end[1]-start[1])/length)
            else:
                center, radius, angle, sweep = arc
                theta = angle+sweep*t
                points = ((center-origin)[None, :] + radius*np.column_stack(
                    (np.cos(theta), np.sin(theta))))/length
                dy = radius*sweep*np.cos(theta)/length
            x, y = points[:, 0], points[:, 1]
            # d/dx [x^(p+1)y^q/(p+1)] = x^p y^q, signed Green integral.
            contributions.append([float(np.dot(weights, x**(p+1)*y**q*dy)/(p+1))
                                  for p, q in EXPONENTS])
        approximations.append([math.fsum(row[i] for row in contributions) for i in range(6)])
    values = np.asarray(approximations)
    require(np.isfinite(values).all(), "independent planar integration is nonfinite")
    agreement = float(np.max(np.abs(values[:-1]-values[-1])/np.maximum(1., np.abs(values[-1]))))
    require(agreement <= 1e-12, "independent planar quadrature did not agree")
    area, mx, my, qxx, qyy, qxy = values[-1]
    require(area > 0. or (not cycles and area == 0.), "nonpositive signed planar area")
    ox, oy = origin
    area *= length**2
    mx0, my0 = mx*length**3, my*length**3
    mx, my = math.fsum((mx0, ox*area)), math.fsum((my0, oy*area))
    qxx = math.fsum((qxx*length**4, 2*ox*mx0, ox*ox*area))
    qyy = math.fsum((qyy*length**4, 2*oy*my0, oy*oy*area))
    qxy = math.fsum((qxy*length**4, ox*my0, oy*mx0, ox*oy*area))
    height, zmean = high-low, .5*(high+low)
    zsecond = (low*low+low*high+high*high)/3.
    weight = density*height
    raw = [float(weight*x) for x in (area, mx, my, area*zmean, qxx, qyy,
                             area*zsecond, my*zmean, mx*zmean, qxy)]
    require(all(math.isfinite(x) for x in raw), "extruded moments overflowed")
    return {"raw_mass_moments": raw, "quadrature_orders": list(ORDERS),
            "maximum_normalized_quadrature_difference": agreement,
            "maximum_boundary_endpoint_gap": maximum_gap,
            "full_topology_verified": False}


def _input_entities(request):
    layers = request.get("fixed_extrusion_layers")
    require(isinstance(layers, list) and bool(layers), "actual fixed extrusion layers required")
    design = {"geometry_representation": "exact_polygonal_prism", "role": "editable_design",
              "loops": request.get("design_region_loops"),
              "z_interval": request.get("extrusion_interval"),
              "volume_density": request.get("extrusion_volume_density")}
    fixed = [{"geometry_representation": "exact_polygonal_prism", "role": "fixed_material",
              "loops": layer.get("loops"), "z_interval": layer.get("z_interval"),
              "volume_density": layer.get("volume_density")} for layer in layers]
    return [design, *fixed]


def _sum(rows):
    return [math.fsum(row[i] for row in rows) for i in range(10)]


def _properties(raw):
    mass = raw[0]
    require(mass > 0., "final material must retain positive mass")
    center = np.asarray(raw[1:4])/mass
    q = np.array([[raw[4], raw[9], raw[8]], [raw[9], raw[5], raw[7]], [raw[8], raw[7], raw[6]]])
    inertia = np.trace(q)*np.eye(3)-q
    central = inertia-mass*(np.dot(center, center)*np.eye(3)-np.outer(center, center))
    return {"mass": mass, "center_of_mass": center.tolist(),
            "inertia_at_origin": inertia.tolist(), "inertia_at_center": central.tolist()}


def audit_planar_case_material(case):
    """Check actual saved layered material; caller's passed fields are ignored."""
    require(case.get("solver_family") == "PLANAR_DOGLEG_V1", "not a planar dogleg application")
    request = case.get("planar_problem")
    require(isinstance(request, dict), "saved planar request is missing")
    inputs = _input_entities(request)
    points = np.asarray([_vector(point, 2, "input vertex")
                         for loop in inputs[0]["loops"] for point in loop])
    require(len(points) > 0, "empty input design")
    length = float(np.linalg.norm(np.ptp(points, axis=0)))
    require(math.isfinite(length) and length > 0., "invalid input length")
    origin = (.5*(np.max(points, axis=0)+np.min(points, axis=0))).tolist()
    integrate = lambda entity: integrate_planar_entity(entity, origin=origin, length_scale=length)
    initial_parts = [integrate(entity)["raw_mass_moments"] for entity in inputs]
    initial = _sum(initial_parts)
    initial_mass = scalar(initial[0], "initial geometry mass", positive=True)
    denominators = np.array([initial_mass*length**degree
                            for degree in (0, 1, 1, 1, 2, 2, 2, 2, 2, 2)])
    require(np.isfinite(denominators).all() and np.all(denominators > 0.), "normalization overflow")
    design_height = inputs[0]["z_interval"][1]-inputs[0]["z_interval"][0]
    areal = scalar(request.get("design_areal_density"), "design areal density", positive=True)
    require(math.isclose(areal, design_height*inputs[0]["volume_density"], rel_tol=1e-12),
            "planar and extrusion density mismatch")
    require(not request.get("fixed_material_loops"),
            "saved checker currently requires unambiguous layer-driven fixed material")

    def compare(claim, actual, label):
        claim = np.asarray(_vector(claim, 10, label))
        error = float(np.max(np.abs(claim-np.asarray(actual))/denominators))
        require(error <= 1e-6, label + " differs from independent saved-geometry moments")
        return error

    attempts = case.get("attempts")
    require(isinstance(attempts, list) and bool(attempts), "saved attempts missing")
    audited = []
    for attempt in attempts:
        result = attempt.get("solver_result")
        identity = {key: attempt.get(key) for key in ("run", "seed", "restart", "status")}
        if result is None:
            require(attempt.get("status") == "EVALUATION_FAILED", "missing returned material")
            audited.append({**identity, "moments_match": False, "balance_constraints_pass": False,
                            "reason": "NO_MATERIAL_RETURNED_BY_FAILED_EVALUATION"})
            continue
        entities = result.get("extruded_material_entities")
        require(isinstance(entities, list) and len(entities) == len(inputs), "missing or extra material entity")
        parts, errors, quadratures = [], [], []
        for index, (entity, expected) in enumerate(zip(entities, inputs)):
            require(entity.get("role") == expected["role"], "material role mismatch")
            for key in ("z_interval", "volume_density"):
                require(entity.get(key) == expected[key], "material input " + key + " mismatch")
            if index:
                require(entity.get("geometry_representation") == expected["geometry_representation"],
                        "fixed material representation mismatch")
                require(entity.get("loops") == expected["loops"], "fixed material geometry mismatch")
            else:
                require(entity.get("geometry_representation") ==
                        "exact_line_arc_profile_swept_between_planes", "missing curved final design")
                center = _vector(result.get("center"), 2, "solved center")
                radius = scalar(result.get("radius"), "solved radius")
                require(math.isclose(radius, math.hypot(*center), rel_tol=1e-14, abs_tol=0.),
                        "fidget circle does not pass through target axis")
                for cycle in entity.get("boundary_cycles", []):
                    for edge in cycle:
                        if edge.get("type") == "CIRCULAR_ARC":
                            require(edge.get("center") == center and edge.get("radius") == radius,
                                    "saved arc differs from solved carving circle")
            measured = integrate(entity)
            parts.append(measured["raw_mass_moments"])
            quadratures.append(measured)
            errors.append(compare(entity.get("raw_mass_moments"), parts[-1], "entity moments"))
        total = _sum(parts)
        errors.append(compare(result.get("extruded_complete_material_raw_moments"), total, "total moments"))
        # Independent physical constraint uses the INPUT's initial mass/length,
        # never a scale determined by the amount of retained material.
        residual = [total[i]/(initial_mass*length) for i in (1, 2)]
        audited.append({**identity, "moments_match": True,
                        "material_sha256": digest(entities),
                        "balance_constraints_pass": attempt.get("status") == "CONVERGED"
                            and result.get("status") == "CONVERGED" and math.hypot(*residual) <= 1e-6,
                        "maximum_normalized_moment_error": max(errors),
                        "independent_normalized_first_moment_residual": residual,
                        "independent_total_raw_mass_moments": total,
                        "independent_properties": _properties(total),
                        "entity_integrations": quadratures})
    return {"id": case.get("id"), "normalization": {"origin": origin, "length": length,
                "initial_mass": initial_mass, "source": "INPUT_GEOMETRY_BEFORE_CARVING"},
            "attempts": audited, "request_sha256": digest(request)}


def audit_saved_planar_application(path):
    path = Path(path).resolve(strict=True)
    contents = path.read_bytes()
    payload = parse(contents.decode("utf-8"))
    cases = payload.get("cases")
    require(isinstance(cases, list) and bool(cases), "saved application cases missing")
    require(len({row.get("id") for row in cases}) == len(cases), "duplicate case identity")
    rows = [audit_planar_case_material(case) for case in cases]
    return {"schema": "saved-planar-material-moment-audit-v1", "scope": SCOPE,
            "input": {"path": str(path), "sha256": hashlib.sha256(contents).hexdigest()},
            "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "mass_moment_order": list(ORDER), "cases": rows,
            "full_topology_verified": False, "material_partition_verified": False,
            "manufacturing_approved": False, "original_paper_reproduction": False}
