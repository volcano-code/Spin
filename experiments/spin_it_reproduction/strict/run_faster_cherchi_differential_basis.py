#!/usr/bin/env python3
"""Check every tangent direction of one actual rational Cherchi moment map.

This is a fixed-input local differential check, not a Newton run, author case,
proof of global differentiability or manufacturing certification. Every sample
uses the production source, provenance and exact-material kernels. Independent
signed-tetrahedron moments are compared with the primary surface integrals.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import time

import numpy as np

import run_faster_cherchi_newton_integration as kernel
from make_it_stand.verification.run_evidence import (
    begin_run, finish_run, EXPLICIT_EFFECTIVE_SOURCE_POLICY)
from make_it_stand.verification.spin_common import digest, parse, require

SCOPE = "FULL_TANGENT_BASIS_AT_ONE_FIXED_CHERCHI_MAP_NOT_NEWTON_OR_PAPER_APPROVAL"
STRATUM_SCOPE = "PARABOLOID_FIXED_STRATUM_10x8_NOT_FULL_NEWTON_OR_PAPER_APPROVAL"
HEIGHT_SCOPE = "PURE_HEIGHT_ACTIVE_SUBSPACE_10x5_NOT_ENDPOINT_NEWTON_OR_PAPER_APPROVAL"
PLANE_SCOPE = "PLANE_FIXED_STRATUM_10x3_NOT_FULL_TOP_NEWTON_OR_PAPER_APPROVAL"
FIXED_NORMAL_MAP = "BOUNDED_FINITE_NORMAL_PARABOLOID_V1"


def _write(path, value):
    with path.open("x", encoding="ascii") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def _settings(extension, problem):
    paper = problem.get("paper_projective_hyperboloid", False)
    require(type(paper) is bool, "paper_projective_hyperboloid must be a boolean")
    algebraic = problem.get("algebraic_finite_map", False)
    require(type(algebraic) is bool, "algebraic_finite_map must be a boolean")
    fixed = problem.get("fixed_normal_paraboloid", False)
    require(type(fixed) is bool, "fixed_normal_paraboloid must be a boolean")
    if fixed:
        require(problem['schema']=='cherchi-paraboloid-basis-method-v1' and not paper and not algebraic,
                'fixed-normal paraboloid requires its explicit reduced profile')
        require(problem.get('maximum_icosphere_subdivisions',-1)==-1,
                'fixed-normal paraboloid cannot adapt normal resolution')
        level=problem.get('icosphere_subdivisions')
        require(type(level) is int and 0<=level<=7, 'invalid fixed-normal resolution')
    radius = problem.get("truncation_radius", 0.0)
    require(type(radius) in (int, float) and np.isfinite(radius) and radius >= 0.,
            "truncation radius must be finite and nonnegative; zero selects automatic radius")
    if problem["schema"] == "cherchi-plane-basis-method-v1":
        require(not algebraic and not paper, "plane stratum does not use full-rank source modes")
        settings = extension.SpinTriangulationSettings()
        settings.truncation_radius = radius
        settings.maximum_radius_doublings = problem["maximum_radius_doublings"]
        return settings
    if problem["schema"] in ("cherchi-paraboloid-basis-method-v1", "cherchi-pure-height-basis-method-v1"):
        require(not algebraic, "algebraic complete rank-two source is not implemented")
        require(not paper, "paper hyperboloid map cannot select a paraboloid stratum")
        settings = extension.SpinTriangulationSettings()
        settings.paraboloid_grid_resolution = problem["paraboloid_grid_resolution"]
        settings.truncation_radius = problem["truncation_radius"]
        settings.maximum_radius_doublings = problem["maximum_radius_doublings"]
        if fixed:
            require(hasattr(settings,'fixed_normal_paraboloid'), 'native fixed-normal source interface unavailable')
            settings.fixed_normal_paraboloid=True
            settings.icosphere_subdivisions=level
        return settings
    settings = kernel._make_triangulation_configuration(extension,
        icosphere_subdivisions=problem["icosphere_subdivisions"],
        maximum_icosphere_subdivisions=problem["icosphere_subdivisions"],
        maximum_radius_doublings=problem["maximum_radius_doublings"],
        paper_projective_hyperboloid=paper, algebraic_finite_map=algebraic)[0]
    settings.truncation_radius = radius
    return settings


def _fixed_normal_features(source):
    require(source.get('fixed_normal_paraboloid') is True and source['surface_type']=='PARABOLOID',
            'actual source did not use the requested fixed-normal paraboloid')
    normals=np.asarray(source['fixed_normals'],dtype=float)
    ids=np.asarray(source['source_normal_vertices'])
    faces=np.asarray(source['faces'])
    count=len(source['vertices'])
    require(normals.ndim==2 and normals.shape[1]==3 and len(normals)>0 and
            np.isfinite(normals).all() and np.max(abs(np.linalg.norm(normals,axis=1)-1))<1e-12,
            'invalid fixed-normal directions')
    require(ids.shape==(count,) and np.issubdtype(ids.dtype,np.integer) and
            np.min(ids)>=0 and np.max(ids)<len(normals), 'invalid source normal identities')
    require(faces.ndim==2 and faces.shape[1]==3 and np.issubdtype(faces.dtype,np.integer) and
            np.min(faces)>=0 and np.max(faces)<count,'invalid source faces')
    require(source['initial_icosphere_subdivisions']==source['effective_icosphere_subdivisions'],
            'fixed-normal resolution changed')
    return {'schema':'finite-normal-paraboloid-features-v1','surface_type':'PARABOLOID',
            'normal_vertices':ids.tolist(),'normals':normals.tolist(),'faces':faces.tolist(),
            'interface_vertex_count':int(source['interface_vertex_count'])}


def _retract_paraboloid(q, direction, step):
    """Rank-two spectral projection and unit normalization, for FD only.

    At a regular rank-two point and tangent direction this has derivative d.
    Coefficients 7,8,9 are matrix yz,xz,xy (polynomial factors are twice these).
    This is not a retraction added to the paper's Newton solver.
    """
    p = q+step*direction
    matrix = np.array([[p[4],p[9],p[8]], [p[9],p[5],p[7]], [p[8],p[7],p[6]]])
    values, vectors = np.linalg.eigh(matrix)
    values[np.argmin(abs(values))] = 0.
    matrix = (vectors*values)@vectors.T
    p[4:] = [matrix[0,0],matrix[1,1],matrix[2,2],matrix[1,2],matrix[0,2],matrix[0,1]]
    return p/np.linalg.norm(p)


def stationary_boundary_signature(identities, rational_vertices, faces, design_vertices, design_faces):
    """Ignore only exact retriangulation inside stationary design planes.

    Preserve every provenance-labelled vertex, all non-design-plane triangles,
    and the directed boundary of each stationary planar patch. Material validity
    is checked separately before this call. No approximate coplanarity tolerance
    or moving-quadric face merging is used. This does not prove a neighborhood
    differentiable; it removes an irrelevant fixed-plane diagonal choice from
    the local comparison identity.
    """
    points = [tuple(Fraction(x) for x in row) for row in rational_vertices]
    design = [tuple(Fraction(x) for x in row) for row in design_vertices]
    require(len(identities) == len(points) and len(set(identities)) == len(points),
            "stationary patch identity needs unique vertex provenance")
    planes = set()
    for face in design_faces:
        a,b,c = (design[int(i)] for i in face)
        u,v = [b[i]-a[i] for i in range(3)], [c[i]-a[i] for i in range(3)]
        n = [u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]]
        require(any(n), "degenerate input plane")
        equation = n+[-sum(n[i]*a[i] for i in range(3))]
        scale = next(x for x in n if x)
        planes.add(tuple(x/scale for x in equation))
    planes = sorted(planes)
    memberships = [{j for j,p in enumerate(planes) if sum(p[i]*v[i] for i in range(3))+p[3] == 0}
                   for v in points]
    groups, free = defaultdict(Counter), []
    for face in faces:
        indices = tuple(int(i) for i in face)
        require(len(indices) == 3 and len(set(indices)) == 3 and min(indices) >= 0
                and max(indices) < len(points), "invalid material triangle")
        common = set.intersection(*(memberships[i] for i in indices))
        require(len(common) <= 1, "degenerate material triangle on multiple planes")
        if not common:
            free.append(tuple(sorted(identities[i] for i in indices)))
            continue
        edges = groups[next(iter(common))]
        for a,b in zip(indices, indices[1:]+indices[:1]):
            edges[(identities[a], identities[b])] += 1
    patches = []
    for plane, edges in sorted(groups.items()):
        require(all(count == 1 for count in edges.values()), "duplicate directed patch edge")
        boundary = sorted(edge for edge in edges if (edge[1],edge[0]) not in edges)
        patches.append({"plane": [str(x) for x in planes[plane]], "boundary": boundary})
    return digest({"schema": "exact-stationary-design-patches-v1", "vertices": sorted(identities),
                   "patches": patches, "free_triangles": sorted(free)})


def _vertex_identities(construction):
    # The production raw signature validates uniqueness and source resolution
    # before these same descriptor fields are used to identify stationary patches.
    identities = []
    for r in construction["records"]:
        operands = [(o["source_vertices"][0]["input_object"], o["source_vertices"][0]["input_vertex"])
                    for o in r["operands"]]
        if r["point_type"] == "LPI":
            # A line-plane intersection is invariant to its line endpoint order.
            # Keep line/plane roles and all plane provenance. Do NOT apply this
            # to LNC: swapping its endpoints would also require replacing t.
            operands[:2] = sorted(operands[:2])
        identities.append(json.dumps((r["point_type"], operands,
            [(p["source_triangles"][0]["input_object"], p["source_triangles"][0]["input_triangle"])
             for p in r.get("surface_planes", [])],
            r.get("parameter_t") if r["point_type"] == "LNC" else None), separators=(",", ":")))
    return identities


def _plane_frame_chart(coefficients):
    """Identity of existing native makeQuadricPlaneFrame, not a new chart.

    Match its x-reference/y-reference switch and its seam rejection. Actual
    source motion is still computed and admitted by the native implementation.
    """
    q = np.asarray(coefficients, dtype=float)
    require(q.shape == (10,) and np.isfinite(q).all() and np.all(q[4:] == 0.),
            'plane frame needs finite plane coefficients')
    length = float(np.linalg.norm(q[1:4]))
    require(np.isfinite(length) and length > 0., 'plane frame needs nonzero normal')
    nx = abs(float(q[1])/length)
    require(abs(nx-.9) > 1e-10, 'plane frame is on the native chart boundary')
    return 0 if nx < .9 else 1


def _validate_plane_protocol(protocol):
    """Require the complete frozen three-direction stencil before replay."""
    q = np.asarray(protocol['normalized_coefficients'], dtype=float)
    chart = _plane_frame_chart(q)
    basis = np.asarray(protocol['parameter_basis'], dtype=float)
    require(abs(np.linalg.norm(q)-1.) <= 1e-12 and basis.shape == (10,3) and
            np.isfinite(basis).all() and np.all(basis[4:] == 0.) and
            np.max(np.abs(basis.T@basis-np.eye(3))) <= 1e-12 and
            np.max(np.abs(q@basis)) <= 1e-12, 'invalid plane tangent basis')
    require(type(protocol.get('plane_frame_chart')) is int and
            protocol['plane_frame_chart'] == chart, 'plane protocol frame identity mismatch')
    steps = protocol['steps']
    require(len(steps) == 3 and len(set(steps)) == 3 and
            all(type(h) is float and np.isfinite(h) and h > 0. for h in steps),
            'plane comparison requires three distinct positive steps')
    directions = protocol['directions']
    require(len(directions) == 3, 'plane comparison requires all three directions')
    expected = [{'label': 'base', 'coefficients': q.tolist()}]
    for i, direction in enumerate(directions):
        require(direction['label'] == f'tangent-{i}' and
                np.array_equal(direction['vector'],basis[:,i]),
                'plane direction does not match ordered parameter basis')
        for j, h in enumerate(steps):
            for sign, label in ((1,'plus'),(-1,'minus')):
                perturbed = q+sign*h*basis[:,i]
                expected.append({'label': f'tangent-{i}-{j}-{label}',
                    'coefficients': (perturbed/np.linalg.norm(perturbed)).tolist()})
    require(protocol['samples'] == expected, 'plane samples do not match frozen stencil')
    require(protocol['relative_tolerance'] == 1e-6 and
            protocol['near_zero_absolute_tolerance'] == 1e-10 and
            protocol['near_zero_sensitivity'] == 1e-8,
            'plane comparison cannot relax scientific tolerances')


def _protocol(problem):
    require(problem["schema"] in ("cherchi-full-basis-method-v1", "cherchi-paraboloid-basis-method-v1",
                                  "cherchi-pure-height-basis-method-v1", "cherchi-plane-basis-method-v1") and
            problem["classification"] == "SYNTHETIC_REGRESSION", "invalid method input")
    q = np.asarray(problem["coefficients"], dtype=float)
    require(q.shape == (10,) and np.isfinite(q).all() and np.linalg.norm(q) > 0.,
            "invalid coefficient vector")
    q = q / np.linalg.norm(q)
    # Stable Householder complement, fixed by input, with no eigenvector choice.
    v = q.copy()
    v[0] += 1. if q[0] >= 0. else -1.
    v /= np.linalg.norm(v)
    basis = (np.eye(10)-2*np.outer(v, v))[:, 1:]
    vertices = np.asarray(problem["design_vertices"], dtype=float)
    faces = np.asarray(problem["design_faces"], dtype=np.int64)
    density = float(problem["density"])
    require(np.isfinite(density) and density > 0., "invalid density")
    domain = kernel.MeshSnapshot((vertices, faces))
    raw = kernel._raw_moments(domain, density)
    length = float(np.max(np.ptp(vertices, axis=0)))
    scales = [float(raw[0])*length**degree for degree in (0,1,1,1,2,2,2,2,2,2)]
    require(all(np.isfinite(x) and x > 0 for x in scales), "invalid input normalization")
    height = problem["schema"] == "cherchi-pure-height-basis-method-v1"
    plane = problem["schema"] == "cherchi-plane-basis-method-v1"
    reduced = plane or height or problem["schema"] == "cherchi-paraboloid-basis-method-v1"
    if plane:
        require(np.all(q[4:] == 0.) and np.linalg.norm(q[1:4]) > 0.,
                "plane stratum requires nonzero linear and zero quadratic coefficients")
    if height:
        # c3, l1, l2, Ic13 and Ic23 depend only on moments 0,1,2,3,7,8.
        # Eq.14 with mu>0 therefore restricts p to these SAME six positions
        # in native and dual coordinates. This is not a general 10D Jacobian.
        require(np.all(q[[4,5,6,9]] == 0.), "pure-height active support requires zero diagonal and xy coefficients")
        active = [0,1,2,3,7,8]
        h = q[active].copy()
        h[0] += 1. if h[0] >= 0. else -1.
        h /= np.linalg.norm(h)
        height_basis = np.zeros((10,5))
        height_basis[active] = (np.eye(6)-2*np.outer(h,h))[:,1:]
    if reduced:
        extension = kernel.spin._native_module()
        source = extension.evaluate_spin_it_faster_stratum_source_motion(q, vertices, _settings(extension, problem))
        require(source["defined"], "reduced source unavailable: "+source["diagnostic"])
        basis = np.asarray(source["coefficient_basis"])
        dimension = 3 if plane else 8
        require(source['surface_type'] == ('PLANE' if plane else 'PARABOLOID') and
                basis.shape == (10,dimension) and np.isfinite(basis).all(),
                "unexpected reduced source type or parameter basis dimension")
        if height:
            require(np.max(np.abs(basis@(basis.T@height_basis)-height_basis)) <= 1e-12,
                    "pure-height active directions leave the available native stratum")
            basis = height_basis
    steps = problem["steps"]
    require(len(steps) == 3 and len(set(steps)) == 3 and
            all(type(h) is float and np.isfinite(h) and h > 0 for h in steps),
            "three distinct positive difference steps required")
    directions = [{"label": f"tangent-{i}", "vector": basis[:, i].tolist()}
                  for i in range(basis.shape[1])]
    if not reduced:
        directions.append({"label": "radial-null", "vector": q.tolist()})
    samples = [{"label": "base", "coefficients": q.tolist()}]
    for direction in directions:
        d = np.asarray(direction["vector"])
        for index, h in enumerate(steps):
            for sign, label in ((1, "plus"), (-1, "minus")):
                samples.append({"label": f"{direction['label']}-{index}-{label}",
                                "coefficients": ((q+sign*h*d)/np.linalg.norm(q+sign*h*d) if height or plane else
                                                 _retract_paraboloid(q, d, sign*h)
                                                 if reduced else q+sign*h*d).tolist()})
    protocol = {"normalized_coefficients": q.tolist(), "parameter_basis": basis.tolist(),
            "directions": directions, "steps": steps, "samples": samples,
            "moment_scales": scales, "normalization": "INPUT_DOMAIN_MASS_AND_MAX_SIDE",
            "derivative_coordinates": "REDUCED_PLANE_BASIS" if plane else "PURE_HEIGHT_ACTIVE_BASIS" if height else
                                      "REDUCED_PARABOLOID_BASIS" if reduced else "FULL_COEFFICIENTS",
            "retraction": "PLANE_LINEAR_PERTURBATION_AND_NORMALIZATION" if plane else
                          "ACTIVE_LINEAR_PERTURBATION_AND_NORMALIZATION" if height else
                          "RANK_TWO_SPECTRAL_PROJECTION_AND_NORMALIZATION" if reduced else "LINEAR_PERTURBATION",
            "branch_identity": "EXACT_STATIONARY_DESIGN_PATCHES_AND_UNORDERED_LPI_LINE_V1",
            "relative_tolerance": 1e-6, "near_zero_absolute_tolerance": 1e-10,
            "near_zero_sensitivity": 1e-8, "input_sha256": digest(problem)}
    if problem.get("paper_projective_hyperboloid", False):
        protocol["source_map"] = "PAPER_PROJECTIVE_HYPERBOLOID_V1"
    if problem.get('fixed_normal_paraboloid',False):
        protocol['source_map']=FIXED_NORMAL_MAP
        protocol['fixed_normal_source_features_sha256']=digest(_fixed_normal_features(source))
    algebraic=problem.get("algebraic_finite_map",False)
    require(type(algebraic) is bool,"algebraic_finite_map must be a boolean")
    if algebraic:
        protocol["finite_point_map"]="ALGEBRAIC_V1"
    if plane:
        protocol['plane_frame_chart'] = _plane_frame_chart(q)
        _validate_plane_protocol(protocol)
    return protocol


def check_records(protocol, records):
    """Recompute comparisons from retained calls; no hand-written pass field."""
    height = protocol.get("derivative_coordinates") == "PURE_HEIGHT_ACTIVE_BASIS"
    plane = protocol.get("derivative_coordinates") == "REDUCED_PLANE_BASIS"
    reduced = plane or height or protocol.get("derivative_coordinates") == "REDUCED_PARABOLOID_BASIS"
    if protocol.get('source_map') == FIXED_NORMAL_MAP:
        require(protocol.get('derivative_coordinates') == 'REDUCED_PARABOLOID_BASIS' and
                [direction['label'] for direction in protocol['directions']] ==
                [f'tangent-{i}' for i in range(8)] and len(protocol['steps']) == 3,
                'fixed-normal paraboloid requires all eight reduced directions and three steps')
    require(len(records) == len(protocol["samples"]) == (19 if plane else 31 if height else 49 if reduced else 61), "missing or extra basis evaluation")
    if plane:
        _validate_plane_protocol(protocol)
    for expected, actual in zip(protocol["samples"], records):
        require(actual["label"] == expected["label"] and
                actual["coefficients"] == expected["coefficients"], "evaluation identity mismatch")
    base = records[0]
    if any(row["status"] != "VALID_EXACT_MATERIAL" for row in records):
        return {"status": "FAILED_EVALUATION", "checks": [],
                "rejections": [row for row in records if row["status"] != "VALID_EXACT_MATERIAL"]}
    jacobian = np.asarray(base["jacobian"], dtype=float)
    require(jacobian.shape == (10, 3 if plane else 5 if height else 8 if reduced else 10) and np.isfinite(jacobian).all(), "invalid base Jacobian")
    scales = np.asarray(protocol["moment_scales"])
    values = []
    for row in records:
        fixed_normal_record = ('fixed_normal_source_features' in row or
                               row['branch'].get('source_map') == FIXED_NORMAL_MAP)
        require(not fixed_normal_record or protocol.get('source_map') == FIXED_NORMAL_MAP,
                'fixed-normal source cannot be checked under a different protocol')
        if plane:
            require(row['branch']['surface_type'] == 'PLANE',
                    'plane differential record changed surface type')
            require(type(row['branch'].get('plane_frame_chart')) is int and
                    row['branch']['plane_frame_chart'] == _plane_frame_chart(row['coefficients']),
                    'plane differential record frame identity mismatch')
        expected_map=protocol.get("finite_point_map","LEGACY_V1")
        require(expected_map in ("LEGACY_V1","ALGEBRAIC_V1") and
                row["branch"].get("finite_point_map","LEGACY_V1")==expected_map,
                "actual finite point map differs from registered protocol")
        if protocol.get("source_map") == "PAPER_PROJECTIVE_HYPERBOLOID_V1":
            features = row.get("paper_source_features")
            require(isinstance(features, dict), "missing paper source provenance")
            require(kernel._paper_source_feature_digest(features) ==
                    row["branch"].get("paper_source_features_sha256"),
                    "paper source provenance disagrees with derivative branch")
            require(features["surface_type"] == row["branch"]["surface_type"],
                    "paper source type disagrees with derivative branch")
        if protocol.get('source_map')==FIXED_NORMAL_MAP:
            features=row.get('fixed_normal_source_features')
            require(isinstance(features,dict) and row['branch'].get('source_map')==FIXED_NORMAL_MAP,
                    'missing or different fixed-normal source map')
            require(features.get('surface_type') == row['branch'].get('surface_type') == 'PARABOLOID',
                    'fixed-normal source type disagrees with derivative branch')
            require(digest(features)==row['branch'].get('fixed_normal_source_features_sha256')==
                    protocol['fixed_normal_source_features_sha256'],
                    'fixed-normal provenance changed or differs from registered protocol')
        raw = np.asarray(row["raw_moments"], dtype=float)
        independent = np.asarray(row["independent_raw_moments"], dtype=float)
        require(raw.shape == independent.shape == (10,) and np.isfinite(raw).all()
                and np.isfinite(independent).all(), "nonfinite or invalid actual moments")
        require(np.max(np.abs(raw-independent)/scales) <= 1e-10,
                "primary moments disagree with independent geometry integration")
        values.append(raw)
    checks = []
    all_pass = True
    for direction_index, direction in enumerate(protocol["directions"]):
        expected = (jacobian[:, direction_index] if reduced else
                    jacobian@np.asarray(direction["vector"]))/scales
        sensitivity = float(np.max(np.abs(expected)))
        rows = []
        for step_index, h in enumerate(protocol["steps"]):
            index = 1+direction_index*6+step_index*2
            same = all(records[j]["branch"] == base["branch"] for j in (index, index+1))
            difference = (values[index]-values[index+1])/(2*h)/scales
            error = float(np.max(np.abs(difference-expected)))
            relative = error/sensitivity if sensitivity > protocol["near_zero_sensitivity"] else None
            numeric_pass = (relative <= protocol["relative_tolerance"] if relative is not None
                            else error <= protocol["near_zero_absolute_tolerance"])
            if direction["label"] == "radial-null":
                numeric_pass = numeric_pass and sensitivity <= protocol["near_zero_absolute_tolerance"]
            rows.append({"step": h, "same_branch": same,
                         "normalized_expected_derivative": expected.tolist(),
                         "normalized_finite_difference": difference.tolist(),
                         "normalized_absolute_error": error, "relative_error": relative,
                         "passed": bool(same and numeric_pass)})
        all_pass = all_pass and all(row["passed"] for row in rows)
        checks.append({"direction": direction["label"], "steps": rows})
    return {"status": "PASSED" if all_pass else "FAILED_DIFFERENTIAL_OR_BRANCH", "checks": checks}


def capture_comparison(protocol, records):
    """A rejected comparison is a scientific failure, not a lost run receipt."""
    try:
        return check_records(protocol, records)
    except Exception as error:
        return {"status": "FAILED_COMPARISON", "checks": [],
                "exception_type": type(error).__name__, "diagnostic": str(error)}


def invoke_boolean(command, output, label, *, timeout):
    def save(stdout, stderr, **details):
        streams = {}
        for name, value in (("stdout", stdout), ("stderr", stderr)):
            raw = value.encode("utf-8") if isinstance(value, str) else (value or b"")
            path = output/f"boolean-{label}.{name}"
            path.write_bytes(raw)
            streams[name] = raw.decode("utf-8", errors="backslashreplace")
            streams[name+"_sha256"] = kernel._sha(path)
        _write(output/f"boolean-{label}.json", {**details, **streams})
        return streams
    try:
        completed = subprocess.run(command, capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired as error:
        save(error.stdout, error.stderr, returncode=None, status="TIMEOUT",
             exception_type=type(error).__name__, diagnostic=str(error))
        raise
    except OSError as error:
        save(b"", b"", returncode=None, status="SPAWN_FAILED",
             exception_type=type(error).__name__, diagnostic=str(error))
        raise
    streams = save(completed.stdout, completed.stderr,
                   returncode=completed.returncode, status="PROCESS_EXITED")
    return subprocess.CompletedProcess(completed.args, completed.returncode,
                                       streams["stdout"], streams["stderr"])


def run(output, *, fixture, probe):
    output, fixture, probe = Path(output), Path(fixture).resolve(strict=True), Path(probe).resolve(strict=True)
    output.mkdir(parents=True, exist_ok=False)
    problem = parse(fixture.read_text())
    protocol = _protocol(problem)
    _write(output/"protocol.json", protocol)  # Fixed before any sampled Boolean.
    extension = kernel.spin._native_module()
    sources = kernel._terminal_audit_checker_paths()
    tick = time.perf_counter()
    begin = begin_run(repository_root=kernel.ROOT, runner_files=[Path(__file__)],
        checker_files=sources, asset_files=[fixture, probe],
        build_configuration_files=[kernel.ROOT/"CMakeLists.txt"],
        parameters=protocol, seed_protocol={"seed": 0, "protocol": "FIXED_INPUT_FULL_BASIS_NO_RANDOMNESS"},
        dependency_roots={"numpy": Path(np.__file__).resolve().parent},
        native_modules=[extension], pure_python=False, source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY)
    _write(output/"begin.json", begin.as_dict())
    design = np.asarray(problem["design_vertices"], dtype=float)
    faces = np.asarray(problem["design_faces"], dtype=np.int64)
    density = float(problem["density"])
    design_path = output/"design.obj"
    kernel._write_obj(design_path, design, faces)
    settings = _settings(extension, problem)
    height = protocol["derivative_coordinates"] == "PURE_HEIGHT_ACTIVE_BASIS"
    plane = protocol["derivative_coordinates"] == "REDUCED_PLANE_BASIS"
    reduced = plane or height or protocol["derivative_coordinates"] == "REDUCED_PARABOLOID_BASIS"
    evaluate_source = (extension.evaluate_spin_it_faster_stratum_source_motion if reduced else
                       extension.evaluate_spin_it_faster_source_motion)
    records = []
    for sample in protocol["samples"]:
        record = dict(sample)
        label = sample["label"]
        try:
            source = evaluate_source(sample["coefficients"], design, settings)
            require(source["defined"], "source differential unavailable: " + source["diagnostic"])
            source_vertices = np.asarray(source["vertices"], dtype=float)
            source_faces = np.asarray(source["faces"], dtype=np.int64)
            operation = kernel._source_side_selection_details(source)["operation"]
            source_path, material_path = output/f"source-{label}.obj", output/f"material-{label}.obj"
            kernel._write_obj(source_path, source_vertices, source_faces)
            completed = invoke_boolean([str(probe), operation, str(design_path), str(source_path),
                str(material_path), "--dump-exact-rational"], output, label, timeout=30)
            require(completed.returncode == 0, "Cherchi Boolean failed: " + completed.stderr)
            raw_vertices, raw_faces = kernel._read_obj(material_path)
            snapshot, derivation, orientation = kernel._prepare_exact_motion_material(
                raw_vertices, raw_faces, Path(str(material_path)+".exact.tsv"))
            vertices, oriented = snapshot.arrays()
            construction = kernel._summarize_vertex_construction_provenance(
                (design_path, source_path), material_path, Path(str(material_path)+".construction.tsv"))
            moments = kernel.integrate_rational_polyhedral_raw_moments(
                vertices, oriented, density=density, rational_vertices=snapshot.rational_vertices)
            independent = np.asarray(kernel._raw_moments(snapshot, density))
            raw_topology = kernel._construction_topology_signature(construction, oriented)
            source_branch, source_features = kernel._source_derivative_branch(
                source, paper_projective_hyperboloid=problem.get("paper_projective_hyperboloid", False),
                algebraic_finite_map=problem.get("algebraic_finite_map", False))
            branch = {"operation": operation,
                "construction_topology": stationary_boundary_signature(
                    _vertex_identities(construction), snapshot.rational_vertices, oriented,
                    design.tolist(), faces),
                **source_branch}
            if source_features is not None:
                record["paper_source_features"] = source_features
            if problem.get('fixed_normal_paraboloid',False):
                features=_fixed_normal_features(source)
                record['fixed_normal_source_features']=features
                branch['source_map']=FIXED_NORMAL_MAP
                branch['fixed_normal_source_features_sha256']=digest(features)
            if plane:
                branch['plane_frame_chart'] = _plane_frame_chart(sample['coefficients'])
            record.update(status="VALID_EXACT_MATERIAL", branch=branch,
                raw_triangle_construction_topology=raw_topology,
                raw_moments=moments.tolist(), independent_raw_moments=independent.tolist(),
                independent_moment_pass=bool(np.max(np.abs(moments-independent)/protocol["moment_scales"]) <= 1e-10),
                geometry_identity=snapshot.identity, sidecar_sha256=derivation["sidecar_sha256"],
                provenance_types=sorted({row["point_type"] for row in construction["records"]}),
                vertices=len(vertices), faces=len(oriented))
            if label == "base":
                if reduced and not height:
                    require(np.array_equal(source["coefficient_basis"], protocol["parameter_basis"]),
                            "base coefficient basis differs from frozen protocol")
                velocity=np.asarray(source["vertex_velocity"]).reshape(-1,3,3 if plane else 8 if reduced else 10)
                if height:
                    source_basis=np.asarray(source['coefficient_basis'])
                    target_basis=np.asarray(protocol['parameter_basis'])
                    change=source_basis.T@target_basis
                    require(source_basis.shape == (10,8) and change.shape == (8,5) and
                            np.max(np.abs(source_basis@change-target_basis)) <= 1e-12,
                            'pure-height motion direction is unsupported by the native chart')
                    velocity=velocity@change
                motion = kernel.recover_construction_vertex_motion(construction["records"],
                    design, source_vertices, velocity, vertices)
                record["jacobian"] = kernel.differentiate_polyhedral_raw_moments(vertices, oriented,
                    motion, density=density, rational_vertices=snapshot.rational_vertices).tolist()
        except Exception as error:
            record.update(status="REJECTED", exception_type=type(error).__name__, diagnostic=str(error))
        records.append(record)
        with (output/"evaluations.jsonl").open("a", encoding="ascii") as stream:
            stream.write(json.dumps(record, sort_keys=True, allow_nan=False)+"\n")
    checked = capture_comparison(protocol, records)
    report = {"schema": "cherchi-full-basis-differential-v1", "scope": PLANE_SCOPE if plane else HEIGHT_SCOPE if height else STRATUM_SCOPE if reduced else SCOPE, **checked,
        "protocol": protocol, "normalized_coefficients": protocol["normalized_coefficients"],
        "parameter_basis": protocol["parameter_basis"], "evaluation_count": len(records),
        "optimizer_executed": False, "scientific_evaluations_sha256": digest(records)}
    _write(output/"diagnostic.json", report)
    events = [{"iteration": i, "status": "ok" if row["status"] == "VALID_EXACT_MATERIAL" else "failed",
               "data": row} for i, row in enumerate(records)]
    events.append({"iteration": len(records), "status": "ok" if checked["status"] == "PASSED" else "failed",
                   "data": {"kind": "COMPLETE_BASIS_COMPARISON", **checked}})
    evidence = finish_run(begin, events=events,
        # Empty raw streams remain on disk and are bound inside Boolean JSON
        # receipts; the shared scientific-output schema rejects empty files.
        output_files=sorted(path for path in output.iterdir() if path.is_file()
                            and path.name != "begin.json" and path.stat().st_size > 0),
        measurements={"elapsed_seconds": time.perf_counter()-tick})
    _write(output/"evidence.json", evidence.as_dict())
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.output, fixture=args.fixture, probe=args.probe)
    print(json.dumps({key: report[key] for key in ("status", "scope", "evaluation_count")}))
    return 0 if report["status"] == "PASSED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
