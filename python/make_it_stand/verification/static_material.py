"""Reintegrate saved static-plane material with independent tetrahedron sums.

No optimizer or native moment integrator is loaded. Oriented manifoldness is
checked by mass_moments; self-intersections, nesting, exact Boolean identity,
complete material partition and thickness are not certified. A successful
moment comparison therefore cannot grant material/manufacturing approval.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from .mass_moments import TriangleMesh, compute_mass_properties
from .spin_common import ORDER, parse, require, scalar

SCOPE = "SAVED_STATIC_MOMENTS_AND_BALANCE_NOT_FULL_MATERIAL_OR_PAPER_APPROVAL"
SPIN2014_STATIC_SCOPE = "INDEPENDENT_MOMENTS_SQUARE_KKT_AND_ACTUAL_SUPPORT_NOT_COMPLETE_PARTITION"


def _array(value, shape, label):
    result = np.asarray(value, dtype=float)
    require(result.shape == shape and np.isfinite(result).all(), label+": invalid shape or nonfinite")
    return result


def _moments(vertices, faces, density):
    properties = compute_mass_properties(TriangleMesh(vertices, faces), density)
    q = properties.second_moment
    return np.array([properties.mass, *properties.first_moment,
                     q[0][0], q[1][1], q[2][2], q[1][2], q[0][2], q[0][1]])


def _hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def audit_static_case_material(case):
    require(case.get("solver_family") == "STATIC_PLANE_FAMILY_V1", "not a static-plane record")
    problem = case.get("static_problem")
    require(isinstance(problem, dict) and _hash(problem) == case.get("static_problem_sha256"),
            "saved static input identity mismatch")
    require(problem.get("fixed_material_definition") == "OUTER_MINUS_DESIGN_SAME_DENSITY",
            "unsupported fixed-material definition")
    origin = _array(problem["target_origin"], (3,), "origin")
    axes = _array(problem["target_axes"], (3,3), "axes")
    require(np.max(np.abs(axes.T@axes-np.eye(3))) <= 1e-12 and np.linalg.det(axes) > 0,
            "invalid target frame")
    density = scalar(problem["material_density"], "density", positive=True)
    measure = problem.get("static_objective_measure","DENSITY_WEIGHTED_MASS_V1")
    require(measure in ("DENSITY_WEIGHTED_MASS_V1","UNIFORM_DENSITY_VOLUME_V1"),"invalid static objective measure")
    mass_scale = 1/density if measure == "UNIFORM_DENSITY_VOLUME_V1" else 1.
    outer = (np.asarray(problem["outer_vertices"], dtype=float)-origin)@axes
    design = (np.asarray(problem["design_vertices"], dtype=float)-origin)@axes
    length = float(np.max(np.ptp(outer, axis=0)))
    recorded_length = scalar(problem["target_length"], "input length", positive=True)
    require(np.isfinite(length) and length > 0 and abs(length-recorded_length) <=
            16*np.finfo(float).eps*max(length, recorded_length),
            "normalization is not fixed by input geometry")
    outer_raw = _moments(outer, problem["outer_faces"], density)
    design_raw = _moments(design, problem["design_faces"], density)
    fixed = outer_raw-design_raw
    require(fixed[0] >= 0., "negative fixed shell mass")
    scales = density*length**np.array([3,4,4,4,5,5,5,5,5,5])
    require(np.isfinite(scales).all() and (scales > 0).all(), "invalid normalization scales")
    attempts = case.get("attempts")
    require(isinstance(attempts, list) and attempts, "missing saved static attempts")
    rows = []
    for attempt in attempts:
        require(attempt.get("static_objective_measure","DENSITY_WEIGHTED_MASS_V1") == measure,
                "attempt objective measure differs from input")
        for key in ("controls", "target_length", "material_density", "static_objective", "static_pareto_weight"):
            require(attempt.get(key) == problem[key], "attempt differs from input: "+key)
        result = attempt["solver_result"]
        require(result.get("status") == attempt.get("status"), "terminal status mismatch")
        identity = {"variable_material_vertices": result["variable_material_vertices"],
                    "variable_material_faces": result["variable_material_faces"],
                    "plane_coefficients": result["plane_coefficients"], "density": density}
        if result.get("boundary_optimality") is not None:
            identity["boundary_optimality"] = result["boundary_optimality"]
            identity["material_region"] = result.get("material_region")
        trace = result.get("trace")
        require(isinstance(trace, list) and trace, "missing static evaluation trace")
        if not any(row.get("evaluation_succeeded") is True for row in trace):
            require(attempt["status"] == "GEOMETRY_EVALUATION_FAILED" and len(trace) == 1 and
                    trace[0].get("evaluation") == 0 and trace[0].get("accepted") is False and
                    not identity["variable_material_vertices"] and not identity["variable_material_faces"],
                    "invalid unavailable terminal geometry")
            rows.append({"run": attempt["run"], "restart": attempt["restart"],
                         "status": attempt["status"], "moments_match": False,
                         "terminal_numerics_qualified": False, "material_sha256": _hash(identity),
                         "geometry_status": "UNAVAILABLE_INITIAL_EVALUATION_FAILED",
                         "diagnostic": trace[0].get("diagnostic")})
            continue
        faces, vertices = result["variable_material_faces"], result["variable_material_vertices"]
        if not faces:
            require(not vertices, "empty material retains orphan vertices")
            variable = np.zeros(10)
        else:
            variable = _moments(vertices, faces, density)
        total = variable+fixed
        errors = []
        for actual, field in ((variable, "variable_raw_moments"), (total, "total_raw_moments")):
            recorded = _array(result[field], (10,), field)
            errors.append(float(np.max(np.abs(actual-recorded)/scales)))
        require(max(errors) <= 1e-6, "saved geometry does not reproduce solver mass moments")
        require(total[0] > 0., "nonpositive total material mass")
        plane = _array(result["plane_coefficients"], (4,), "plane")
        require(abs(np.linalg.norm(plane)-1.) <= 1e-12, "unnormalized terminal plane")
        if vertices:
            values = plane[0]+np.asarray(vertices)@plane[1:]/length
            require(float(np.min(values)) >= -1e-12,
                    "saved retained material crosses its terminal plane")
        normalized = total[1:3]/scales[1]
        tolerance = scalar(problem["controls"]["constraint_tolerance"], "tolerance", positive=True)
        require(tolerance <= 1e-10, "static tolerance was relaxed")
        objective, weight = problem["static_objective"], problem["static_pareto_weight"]
        require(objective in ("MAXIMUM_MASS", "LOWER_CENTER_OF_MASS", "PARETO_TRADEOFF"), "invalid objective")
        require(type(weight) in (int,float) and 0 <= weight <= 1, "invalid Pareto weight")
        maximum = objective == "MAXIMUM_MASS" or (objective == "PARETO_TRADEOFF" and weight == 1.)
        boundary_verified = False
        if maximum:
            stationary = plane[0] > 0 and abs(plane[3]) <= 1e-12
            residual = abs(float(plane[3]))
        else:
            weight = 0. if objective == "LOWER_CENTER_OF_MASS" else weight
            residual = abs(float(-plane[0]*(1-weight)/total[0]-plane[3]/length*(
                weight*mass_scale+(1-weight)*total[3]/total[0]**2)))
            stationary = plane[3] < -1e-12 and residual <= tolerance
        witness = result.get("boundary_optimality")
        if not maximum and result.get("material_region") == "STRICT_FULL" and witness is None:
            require(attempt["status"] != "CONVERGED",
                    "typed full-domain state requires its boundary witness")
            stationary = False
        if witness is not None:
            require(not maximum and isinstance(witness,dict) and
                    result.get("material_region") == "STRICT_FULL" and
                    witness.get("kind") == "STRICT_FULL_DOMAIN_ZERO_EQUALITY_MULTIPLIERS_V1" and
                    witness.get("objective_measure") == measure,
                    "invalid full-domain witness identity")
            recorded_design = _array(attempt.get("design_target_vertices"), design.shape, "design bounds")
            require(np.allclose(recorded_design,design,rtol=0,atol=16*np.finfo(float).eps*
                               max(length,float(np.max(np.abs(design))))), "design bounds mismatch")
            # The saved mesh is independently integrated. This moment agreement
            # plus strict full-halfspace sign is numerical supporting evidence,
            # NOT a full set-equality/topology/partition certificate.
            require(np.max(np.abs(variable-design_raw)/scales) <= 1e-10 and
                    float(np.min(plane[0]+design@plane[1:]/length)) > 1e-12,
                    "full-domain witness disagrees with saved material or actual design")
            a = weight*mass_scale+(1-weight)*total[3]/total[0]**2
            bz = -(1-weight)/total[0]
            values = a+bz*design[:,2]
            guard = 64*np.finfo(float).eps*float(np.max(abs(a)+abs(bz*design[:,2])))
            expected = np.array([a,0.,0.,length*bz])
            require(np.isfinite(expected).all() and np.linalg.norm(expected) > 0,
                    "invalid full-domain negative gradient")
            expected /= np.linalg.norm(expected)
            claimed = _array(witness["plane_coefficients"],(4,),"optimality plane")
            minimum = scalar(witness["minimum_gradient_polynomial"],"gradient minimum",positive=True)
            require(float(np.min(values)) > guard and
                    abs(minimum-float(np.min(values))) <= 1e-12*abs(float(np.min(values))) and
                    np.allclose(claimed,expected,rtol=0,atol=1e-13),
                    "full-domain negative-gradient witness failed independent recomputation")
            boundary_verified = True
            stationary = True
        qualified = bool(np.linalg.norm(normalized) <= tolerance and stationary)
        if attempt["status"] == "CONVERGED":
            require(qualified, "saved converged material fails independent static constraints")
        rows.append({"run": attempt["run"], "restart": attempt["restart"], "status": attempt["status"],
                     "material_sha256": _hash(identity),
                     "moments_match": True, "normalized_moment_error": max(errors),
                     "independent_total_raw_moments": total.tolist(),
                     "independent_variable_raw_moments": variable.tolist(),
                     "normalized_horizontal_residual": float(np.linalg.norm(normalized)),
                     "full_domain_gradient_sign_verified": boundary_verified,
                     "stationarity_residual": residual, "terminal_numerics_qualified": qualified})
    return {"id": case["id"], "input_sha256": _hash(problem), "attempts": rows}


def audit_saved_static_application(path):
    path = Path(path).resolve(strict=True)
    contents = path.read_bytes()
    payload = parse(contents.decode("utf-8"))
    cases = payload.get("cases")
    require(isinstance(cases, list) and cases, "missing static cases")
    require(len({case["id"] for case in cases}) == len(cases), "duplicate static case identity")
    return {"schema": "saved-static-material-moment-audit-v1", "scope": SCOPE,
            "input": {"path": str(path), "sha256": hashlib.sha256(contents).hexdigest()},
            "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "mass_moment_order": list(ORDER), "cases": [audit_static_case_material(case) for case in cases],
            "full_topology_verified": False, "material_partition_verified": False,
            "manufacturing_approved": False, "original_paper_reproduction": False}


def audit_spin2014_static_directory(directory):
    """Reintegrate saved 2014 square-objective material, without native calls.

    This extends the static checker, not the material seal. Topology beyond
    oriented manifoldness and exact complete partition remain separate gates.
    All numerical expectations derive from saved geometry and input units;
    neither the candidate's old pass flag nor the new flag grants acceptance.
    """
    from dataclasses import asdict
    from ..stability_oracle import _convex_hull, evaluate_quasi_static_tilt
    root=Path(directory).resolve(strict=True)
    input_path,run_path=root/"frozen-input.json",root/"native-run.json"
    inputs=parse(input_path.read_text());record=parse(run_path.read_text())
    request,objective=inputs["request"],inputs["objective"]
    require(request["mode"]=="STATIC_BALANCE" and objective["objective"]=="STATIC_LOWERING",
            "saved request is not the 2014 static square objective")
    require(record["objective_profile"]=="SPIN2014_OBJECTIVES_QUADRIC_V1" and
            record["objective_settings"]==objective,"saved objective identity mismatch")
    measured=record["measured_request"]
    for key,value in request.items():
        if key!="fixed_shell_moments":
            require(measured[key]==value,"saved request/geometry identity mismatch: "+key)
    with np.load(root/"input-geometry.npz",allow_pickle=False) as geometry:
        for field,key in (("outer_vertices","outer_vertices"),("outer_faces","outer_faces"),
                          ("design_vertices","design_vertices"),("design_faces","design_faces")):
            dtype=np.int64 if field.endswith("faces") else np.float64
            one,two=np.asarray(geometry[field],dtype=dtype),np.asarray(request[key],dtype=dtype)
            require(one.shape==two.shape and one.tobytes()==two.tobytes(),"saved input geometry mismatch: "+field)
    origin=_array(request["target_origin"],(3,),"origin")
    axes=_array(request["target_axes"],(3,3),"axes")
    require(np.max(abs(axes.T@axes-np.eye(3)))<=1e-12 and np.linalg.det(axes)>0,"invalid target axes")
    fixed_density=scalar(inputs["fixed_material_density"],"fixed density",positive=True)
    variable_density=scalar(request["material_density"],"design density",positive=True)
    require(record["fixed_material_density"]==fixed_density,"saved fixed density mismatch")
    mass_unit=scalar(objective["solver_mass_unit_in_declared_units"],"mass unit",positive=True)
    length_unit=scalar(objective["solver_length_unit_in_declared_units"],"length unit",positive=True)
    weight=scalar(objective["center_weight"],"center weight")
    require(weight>=0,"negative center weight")
    outer=np.asarray(request["outer_vertices"],dtype=float)
    design=np.asarray(request["design_vertices"],dtype=float)
    outer_target=(outer-origin)@axes;design_target=(design-origin)@axes
    length=float(np.max(np.ptp(outer_target,axis=0)))
    require(abs(length-record["target_length"])<=16*np.finfo(float).eps*length,
            "normalization length does not come from the input")
    require(float(np.min(outer_target[:,2]))>=0.,"outer geometry lies below support")
    candidate,checked=record["candidate"],record["squared_height"]
    # Transform raw model moments analytically instead of moving outer vertices
    # then integrating rounded coordinates (which may collapse thin features).
    def target_raw(vertices,faces,density):
        raw=_moments(vertices,faces,density);mass=raw[0]
        first=raw[1:4]
        second=np.array([[raw[4],raw[9],raw[8]],[raw[9],raw[5],raw[7]],
                         [raw[8],raw[7],raw[6]]])
        q=axes.T@(second-np.outer(origin,first)-np.outer(first,origin)+mass*np.outer(origin,origin))@axes
        l=axes.T@(first-mass*origin)
        return np.array([mass,*l,q[0,0],q[1,1],q[2,2],q[1,2],q[0,2],q[0,1]])
    outer_raw=target_raw(outer,request["outer_faces"],fixed_density)
    design_fixed=target_raw(design,request["design_faces"],fixed_density)
    design_raw=target_raw(design,request["design_faces"],variable_density)
    fixed=outer_raw-design_fixed
    faces,vertices=candidate["variable_material_faces"],candidate["variable_material_vertices"]
    require(bool(faces)==bool(vertices),"saved material has orphan vertices or faces")
    variable=_moments(vertices,faces,variable_density) if faces else np.zeros(10)
    total=fixed+variable;full=fixed+design_raw
    scales=max(fixed_density,variable_density)*length**np.array([3,4,4,4,5,5,5,5,5,5])
    require(np.isfinite(scales).all() and (scales>0).all(),"unrepresentable fixed normalization")
    errors=[float(np.max(abs(actual-_array(candidate[field],(10,),field))/scales))
            for actual,field in ((variable,"variable_raw_moments"),(total,"total_raw_moments"))]
    # Fixed moments in the measured request are in model coordinates, not target.
    fixed_model=_moments(outer,request["outer_faces"],fixed_density)-_moments(
        design,request["design_faces"],fixed_density)
    errors.append(float(np.max(abs(fixed_model-_array(measured["fixed_shell_moments"],(10,),"fixed moments"))/scales)))
    require(max(errors)<=1e-6,"saved geometry/density does not reproduce recorded mass moments")
    mass=scalar(float(total[0]),"final mass",positive=True)
    h=float(total[3]/mass)
    require(h>=0,"negative center height")
    factor=weight*length_unit**2
    objective_value=factor*h*h
    gradient=np.zeros(10);gradient[0]=-2*factor*h*h/mass;gradient[3]=2*factor*h/mass
    hessian=np.zeros((10,10));hessian[0,0]=6*factor*h*h/mass**2
    hessian[0,3]=hessian[3,0]=-4*factor*h/mass**2;hessian[3,3]=2*factor/mass**2
    for value,expected,label in ((checked["objective"],np.asarray(objective_value),"objective"),
        (checked["objective_gradient"],gradient,"objective gradient"),
        (checked["objective_hessian"],hessian,"objective Hessian")):
        actual=_array(value,expected.shape,label)
        require(float(np.max(abs(actual-expected)/np.maximum(1.,abs(expected))))<=1e-10,
                "saved squared "+label+" mismatch")
    plane=_array(candidate["plane_coefficients"],(4,),"plane")
    require(abs(np.linalg.norm(plane)-1)<=1e-12,"terminal plane is not normalized")
    physical=plane.copy();physical[1:]/=length
    equality=_array(checked["equality_multipliers"],(2,),"equality multipliers")
    residual=gradient[:4]+float(checked["multiplier"])*physical
    residual[1:3]+=mass_unit*length_unit*equality
    stationarity=float(np.max(abs(residual)))
    constraint_scale=scalar(candidate["constraint_scale"],"constraint scale",positive=True)
    horizontal=float(np.max(abs(total[1:3]))/constraint_scale)
    tolerance=1e-10
    boundary=False
    if candidate["material_region"]=="STRICT_FULL":
        negative=-gradient[:4];values=negative[0]+design_target@negative[1:]
        guard=64*np.finfo(float).eps*float(np.max(abs(negative[0])+abs(design_target)@abs(negative[1:])))
        boundary=bool(float(np.min(values))>guard and horizontal<=tolerance and
                      np.max(abs(variable-design_raw)/scales)<=1e-10)
        if checked["boundary_witness_verified"]:
            witness=negative.copy();witness[1:]*=length
            multiplier=float(np.linalg.norm(witness))
            require(multiplier>0 and boundary,"invalid saved full-domain square witness")
            require(np.allclose(witness/multiplier,checked["boundary_witness_coefficients"],atol=1e-13,rtol=0),
                    "saved boundary witness uses a different objective")
    terminal=bool(weight>0 and h>0 and horizontal<=tolerance and
        (boundary or (plane[3]<0 and checked["multiplier"]>0 and stationarity<=tolerance)))
    require(bool(checked["squared_height_conditions_verified"])==terminal,
            "saved square-objective flag disagrees with independent conditions")
    # Verify that the displayed support is the hull of actual support faces.
    points={tuple(outer[i]) for face in request["outer_faces"]
            if all(outer_target[i,2]==0 for i in face) for i in face}
    require(len(points)>=3,"actual outer geometry has no bottom support faces")
    actual=np.asarray(sorted(points));projected=(actual-origin)@axes
    hull=_convex_hull(projected[:,:2]);hull_world=origin+np.c_[hull,np.zeros(len(hull))]@axes.T
    support=_array(inputs["actual_support_hull"],hull_world.shape,"support hull")
    require(np.allclose(np.asarray(sorted(map(tuple,support))),np.asarray(sorted(map(tuple,hull_world))),
                        rtol=0,atol=16*np.finfo(float).eps*length),"saved support is not the actual bottom-face hull")
    center_model=origin+axes@(total[1:4]/mass)
    support_check=asdict(evaluate_quasi_static_tilt(support,center_model,-axes[:,2],
        tilt_degrees=0.,azimuth_degrees=0.)) if h>0 else {"inside_support":False,"diagnostic":"ZERO_HEIGHT_DEGENERATE"}
    # The oracle's comparison with a numpy tolerance returns numpy.bool_.
    # Persist a real strict-JSON boolean without changing its truth value.
    support_check["inside_support"]=bool(support_check["inside_support"])
    for field in ("projected_center","tilted_gravity"):
        if field in support_check:
            support_check[field]=list(support_check[field])
    material_change=float(np.max(abs(variable-design_raw)/scales))
    initial_objective=factor*float(full[3]/full[0])**2
    nontrivial=bool(candidate["accepted_steps"]>0 and material_change>1e-10 and objective_value<initial_objective)
    return {"schema":"saved-spin2014-static-square-audit-v1",
        "scope":SPIN2014_STATIC_SCOPE,
        "inputs":{"frozen_input_sha256":hashlib.sha256(input_path.read_bytes()).hexdigest(),
                  "native_run_sha256":hashlib.sha256(run_path.read_bytes()).hexdigest()},
        "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "independent_total_raw_moments":total.tolist(),"normalized_moment_error":max(errors),
        "independent_objective":objective_value,"initial_full_material_objective":initial_objective,
        "independent_stationarity_residual":stationarity,"normalized_horizontal_residual":horizontal,
        "terminal_numerics_qualified":terminal,"nontrivial_optimization":nontrivial,
        "normalized_material_change":material_change,"support":support_check,
        "material_partition_verified":False,"manufacturing_approved":False,
        "original_paper_reproduction":False}
