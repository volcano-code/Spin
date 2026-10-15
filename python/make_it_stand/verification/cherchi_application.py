"""Rejoin saved application records to original runs and re-load exact material.

No optimizer is called. Ledger parsing is shared with the Table adapter;
material approval instead runs the existing independent SealedSpinArtifact
geometry/integration checks on persisted payloads. Failures remain auditable
but do not grant successful material or original-paper evidence.
"""
import importlib.util
import json
import math
from pathlib import Path

from .run_evidence import PaperRunEvidenceV2, check_freshness
from .spin_bundle import load_spin_artifact

SCOPE="SAVED_CHERCHI_APPLICATION_AND_EXACT_MATERIAL_NOT_PAPER_OR_MANUFACTURING"

def _table():
    path=Path(__file__).resolve().parents[3]/"experiments/spin_it_reproduction/faster/run_paper_table1.py"
    spec=importlib.util.spec_from_file_location("cherchi_table_saved_checker",path)
    table=importlib.util.module_from_spec(spec);spec.loader.exec_module(table)
    return table

def _require(value,message):
    if not value:raise ValueError(message)

def _recheck_registered_plan(table,case,manifest):
    """Re-read every registered mesh, not only the manifest's JSON bytes."""
    plan=case["initialization_plan"]
    relative=Path(plan["problem_identity"]["path"])
    _require(not relative.is_absolute() and relative.parts and ".." not in relative.parts,
             "registered manifest path must be canonical and relative")
    manifest=Path(manifest).resolve(strict=True)
    asset_root=manifest
    for _ in relative.parts:asset_root=asset_root.parent
    _require((asset_root/relative).resolve()==manifest,"registered manifest path differs from original run")
    rebuilt=table.cherchi_case_plan({"id":case["id"],"mode":case["mode"],
        "cherchi_problem":plan["problem_identity"],"cherchi_settings":plan["settings"],
        "center_of_mass_weight":plan["center_of_mass_weight"],
        "spinning_top_objective":plan.get("spinning_top_objective","INERTIA_AND_LOWERING")},len(case["seed_sequence"]),
        plan["runs"][0]["maximum_restarts"],asset_root=asset_root,initial_seeds=tuple(case["seed_sequence"]))
    if "objective_definition" not in plan and "spinning_top_objective" not in plan:
        # Old records had only the finite objective. Never infer a pure endpoint.
        rebuilt.pop("objective_definition");rebuilt.pop("spinning_top_objective")
    _require(rebuilt==plan,"registered Cherchi plan differs from actual input assets or protocol")

def audit_saved_cherchi_case(case, *, bundle_root):
    table=_table()
    _require(isinstance(case,dict) and case.get("solver_family")=="CHERCHI_RAW_EQ14_V1",
             "saved application is not a Cherchi run")
    _require(case.get("scientific_hash")==table.scientific_hash({k:v for k,v in case.items() if k!="scientific_hash"}),
             "saved Cherchi application digest mismatch")
    root=Path(bundle_root).resolve(strict=True)
    plan=case["initialization_plan"];seeds=case["seed_sequence"]
    _require(plan["id"]==case["id"] and plan["mode"]==case["mode"],"application plan identity mismatch")
    table.registered_seeds(len(seeds),seeds)
    _require(len(case["runs"])==len(seeds)==len(plan["runs"]),"missing or duplicate top-level run")
    recomputed=[];reports=[];qualified=0
    for i,(row,planned,seed) in enumerate(zip(case["runs"],plan["runs"],seeds)):
        _require(row["run"]==planned["run"]==i and row["seed"]==planned["seed"]==seed,
                 "top-level run/seed differs from registered plan")
        relative=Path(row["relative_output"])
        directory=(root/relative).resolve()
        _require(str(relative)==f"{case['id']}/run-{i:06d}" and directory.is_relative_to(root),
                 "saved Cherchi run path escapes or substitutes the application")
        record_path=directory/"result.json"
        record=json.loads(record_path.read_text())
        _require(record["status"]==row["status"] and record["solver_status"]==row["solver_status"],
                 "application status differs from native run")
        if record["status"]=="INFEASIBLE_CERTIFIED":
            from .faster_infeasibility import verify_axis_aligned_yoyo_moment_capacity
            certificate=json.loads((directory/"preflight-feasibility-report.json").read_text())
            # A preflight-only case has no Newton evidence. Never fabricate a
            # solver attempt or call a failure an infeasibility certificate.
            assets=certificate["continuous_certificate"]["asset_sha256"]
            _require(isinstance(assets,list) and len(assets)==4,"certificate asset binding is incomplete")
            manifest=Path(assets[0][0])
            _recheck_registered_plan(table,case,manifest)
            _require(assets[0][1]==plan["problem_identity"]["sha256"] and
                     json.loads(manifest.read_text())==plan["problem"],"certificate belongs to another problem")
            fresh=verify_axis_aligned_yoyo_moment_capacity(manifest)
            _require(fresh.to_dict()==certificate and fresh.status=="INFEASIBLE_CERTIFIED",
                     "infeasibility certificate did not independently reproduce")
            _require(fresh.input_hash==record["input_hash"]==row["preflight_input_hash"] and
                     row["function_evaluations"]==record["function_evaluations"]==0,
                     "certificate input or zero-cost record mismatch")
            reports.append({"run":i,"seed":seed,"status":"INFEASIBLE_CERTIFIED","material_verified":False})
            continue
        evidence=PaperRunEvidenceV2.from_json((directory/"evidence.json").read_text())
        _require(check_freshness(evidence).status=="CURRENT","STALE_RESULT: original Cherchi run evidence")
        bound={Path(x["path"]).resolve() for x in evidence.as_dict()["scientific"]["outputs"]}
        _require(record_path.resolve() in bound and (directory/"inputs.json").resolve() in bound,
                 "run record or request is not bound to evidence")
        inputs=json.loads((directory/"inputs.json").read_text())
        _require(evidence.as_dict()["begin"]["request"]["parameters"]==inputs,
                 "saved request differs from evidence")
        expected_map=plan["settings"].get("paper_projective_hyperboloid",False)
        actual_map=inputs.get("paper_projective_hyperboloid",False)
        _require(type(expected_map) is bool and type(actual_map) is bool and actual_map==expected_map,
                 "saved paper_projective source map differs from registered plan")
        expected_algebraic=plan["settings"].get("algebraic_finite_map",False)
        actual_algebraic=inputs.get("algebraic_finite_map",False)
        _require(type(expected_algebraic) is bool and type(actual_algebraic) is bool and
                 actual_algebraic==expected_algebraic,
                 "saved algebraic finite map differs from registered plan")
        selector=plan.get("spinning_top_objective","INERTIA_AND_LOWERING")
        _require(inputs.get("spinning_top_objective","INERTIA_AND_LOWERING")==selector,
                 "saved objective selector differs from plan")
        if "objective_definition" in plan or "objective_definition" in inputs:
            from .faster_material import rotational_objective_definition
            expected=rotational_objective_definition(case["mode"],plan["center_of_mass_weight"],selector)
            _require(plan.get("objective_definition")==inputs.get("objective_definition")==expected,
                     "saved objective definition differs from registered formula")
        for key,expected in (("mode",case["mode"]),("deterministic_seed",seed),
            ("initialization_mode","deterministic_mode_seed"),("initial_coefficients",None),
            ("maximum_restarts",planned["maximum_restarts"]),
            ("maximum_function_evaluations_per_restart",plan["settings"]["maximum_function_evaluations"]),
            ("center_of_mass_weight",plan["center_of_mass_weight"]),
            ("quadric_icosphere_subdivisions",plan["settings"]["icosphere_subdivisions"]),
            ("maximum_icosphere_subdivisions",-1),
            ("maximum_radius_doublings",plan["settings"]["maximum_radius_doublings"]),
            ("cherchi_boolean_timeout_seconds",plan["settings"]["cherchi_timeout_seconds"]),
            ("equation_profile","DUAL_MONOMIAL_RAW_V1")):
            _require(inputs.get(key)==expected,"saved Cherchi request differs from plan: "+key)
        manifest=inputs["problem_manifest"]
        _recheck_registered_plan(table,case,manifest["path"])
        _require(manifest["sha256"]==plan["problem_identity"]["sha256"] and
                 json.loads(Path(manifest["path"]).read_text())==plan["problem"],
                 "saved locked material differs from application plan")
        _require(row["locked_problem_digest"]==record["locked_problem_digest"]==inputs["locked_problem_digest"],
                 "locked problem identity mismatch")
        runner=table._cherchi_runner()
        _require(record["scientific_event_hash_v2"]==runner._scientific_event_digest_v2(record["events"])==row["scientific_event_hash_v2"],
                 "event identity differs from original ledger")
        parsed=table.cherchi_attempt_rows(record,run=i,seed=seed,
            maximum_restarts=planned["maximum_restarts"],
            maximum_evaluations=plan["settings"]["maximum_function_evaluations"])
        table.validate_attempt_ledger([{**r,"run":0} for r in parsed],expected_runs=1,
                                     maximum_restarts=planned["maximum_restarts"])
        recomputed.extend(parsed)
        _require(row["function_evaluations"]==sum(r["function_evaluations"] for r in parsed) and
                 row["geometry_callback_evaluations"]==len(record["evaluations"]),"saved evaluation cost mismatch")
        geometry=None
        if record["status"]=="CONVERGED":
            claimed=record["sealed_artifact"]
            bundle=directory/"sealed-material-artifact"
            _require(Path(claimed["bundle_path"]).resolve()==bundle.resolve(),"substituted material bundle")
            artifact=load_spin_artifact(bundle,expected_problem_digest=inputs["locked_problem_digest"])
            _require(artifact.scientific_digest==claimed["scientific_digest"]==row["sealed_material_digest"],
                     "saved material identity mismatch")
            _require(parsed[-1]["numerical_acceptance"]["passed"],"material is not a qualified Eq14 solution")
            geometry={"scientific_digest":artifact.scientific_digest,"facts":artifact.facts}
            qualified+=1
        else:
            _require(row["sealed_material_digest"]==record.get("sealed_artifact",{}).get("scientific_digest"),
                     "failed run changed material identity")
        reports.append({"run":i,"seed":seed,"status":record["status"],"material_verified":geometry is not None,
            "function_evaluations":record["function_evaluations"],
            "scientific_event_hash_v2":record["scientific_event_hash_v2"],"geometry":geometry})
    _require(recomputed==case["attempts"],"application attempts differ from saved native ledger")
    summary=table.summarize_attempts(recomputed,len(seeds))
    _require(summary==case["summary"],"application summary differs from complete population")
    _require(case["passed"]==(qualified==len(seeds) and summary["strict_terminal_numerics"]["all_runs_numerically_qualified"]),
             "application passed flag differs from independently checked runs")
    return {"id":case["id"],"scope":SCOPE,"runs":reports,"qualified_material_runs":qualified,
            "expected_runs":len(seeds),"total_function_evaluations":sum(r["function_evaluations"] for r in recomputed),
            "formal_paper_admission":"NOT_GRANTED_BY_METHOD_CHECKER"}


def _top_material_metrics(facts, target, weight):
    """Only call with independently reloaded/integrated material facts.

    The native top objective uses linear target-frame height in input-normalized
    length units, not squared height and not the model's z coordinate.
    """
    import numpy as np
    from .faster_material import _target_axes
    physical=facts["physical"]
    inertia=np.asarray(physical["central_inertia"],dtype=float)
    center=np.asarray(physical["center_of_mass"],dtype=float)
    origin=np.asarray(target["origin"],dtype=float)
    axis=np.asarray(target["axis"],dtype=float)
    length=float(facts["normalization"]["L"])
    mass=float(physical["mass"])
    _require(inertia.shape==(3,3) and center.shape==origin.shape==axis.shape==(3,),"invalid material metrics shape")
    _require(all(np.isfinite(x).all() for x in (inertia,center,origin,axis)) and
             math.isfinite(length) and length>0 and math.isfinite(mass) and mass>0 and
             np.max(np.abs(axis))>0 and type(weight) in (int,float) and math.isfinite(weight) and weight>=0,
             "invalid material metrics or target")
    axes=_target_axes(axis)
    rotated=axes.T@inertia@axes
    _require(rotated[2,2]>0,"top inertia about target axis is nonpositive")
    fyo=float((rotated[0,0]**2+rotated[1,1]**2+2*rotated[0,1]**2)/rotated[2,2]**2)
    height=float(axes[:,2]@(center-origin)/length)
    values={"physical_mass":mass,"f_yo":fyo,"center_height_input_normalized":height,
            "objective":fyo+weight*height}
    _require(all(math.isfinite(x) for x in values.values()),"nonfinite material metrics")
    return values


def _audit_finite_pareto(value, root, reports):
    """Rejoin the frozen sweep, each native request and the full population."""
    table=_table()
    _require(value.get("scientific_hash")==table.scientific_hash({k:v for k,v in value.items() if k!="scientific_hash"}),
             "Pareto digest mismatch")
    weights=value.get("weights",[])
    _require(weights and all(type(w) in (int,float) and math.isfinite(w) and w>=0 for w in weights)
             and len(set(weights))==len(weights),"invalid Pareto weights")
    _require(value.get("mode")=="SPINNING_TOP" and value.get("solver_family")=="CHERCHI_RAW_EQ14_V1",
             "wrong Pareto mode or solver")
    _require(value.get("warm_started") is False and
             value.get("initialization_protocol")=="INDEPENDENT_INPUT_SEEDS_PER_WEIGHT_V1", "invalid Pareto starts")
    endpoint=value.get("endpoint_status",{})
    _require(endpoint.get("pure_lowering") in ("BLOCKED_MISSING_OBJECTIVE",
             "BLOCKED_MISSING_GEOMETRY_ENDPOINT_VALIDATION") and
             endpoint=={"pure_inertia":"REQUESTED" if 0 in weights else "NOT_REQUESTED",
                        "pure_lowering":endpoint.get("pure_lowering")}, "unsupported Pareto endpoint claim")
    _require(value.get("formal_paper_admission")=="NOT_GRANTED_BY_METHOD_RUNNER", "unverified paper admission")
    plans=value["initialization_plans"];cases=value["cases"]
    _require(len(weights)==len(plans)==len(cases)==len(reports),"missing or extra Pareto weight runs")
    frozen=json.loads((root/"pareto-plan.json").read_text())
    _require(frozen=={k:v for k,v in value.items() if k not in ("cases","cost_ledger","scientific_hash")},
             "Pareto request differs from frozen plan")
    reference=None;points=[]
    for i,(weight,plan,case,report) in enumerate(zip(weights,plans,cases,reports)):
        _require(case["id"]==f"{value['case_id']}-pareto-{i:06d}" and case["mode"]=="SPINNING_TOP",
                 "Pareto case identity mismatch")
        _require(plan==case["initialization_plan"] and plan["center_of_mass_weight"]==weight,
                 "Pareto weight differs from actual solver input")
        _require(plan.get("spinning_top_objective","INERTIA_AND_LOWERING")=="INERTIA_AND_LOWERING",
                 "pure endpoint substituted for finite-weight Pareto objective")
        # Every objective definition was already re-derived against its saved
        # native request. It includes the scanned weight, so is not invariant.
        common={k:v for k,v in plan.items() if k not in ("id","center_of_mass_weight","objective_definition")}
        if reference is None:reference=common
        _require(common==reference,"Pareto inputs or starts changed between weights")
        _require(json.loads((root/f"pareto-case-{i:06d}.json").read_text())==case,
                 "Pareto case differs from persisted application")
        for run in report["runs"]:
            point={"weight_index":i,"weight":weight,"run":run["run"],"seed":run["seed"],
                   "status":run["status"],"material_verified":run["material_verified"]}
            if run["material_verified"]:
                point["metrics"]=_top_material_metrics(run["geometry"]["facts"],plan["problem"]["target"],weight)
                point["material_digest"]=run["geometry"]["scientific_digest"]
            points.append(point)
    expected={"reference_attempts":0,"weight_solves":len(cases),
              "top_level_runs":sum(r["expected_runs"] for r in reports),
              "function_evaluations_all_solves":sum(r["total_function_evaluations"] for r in reports)}
    _require(value.get("cost_ledger")==expected,"Pareto cost excludes or changes actual work")
    return {"weights":weights,"endpoint_status":value["endpoint_status"],"cost_ledger":expected,"points":points,
            "original_figure12_complete":False}


def audit_saved_cherchi_application(path, *, bundle_root):
    value=json.loads(Path(path).read_text())
    cases=value.get("cases")
    _require(isinstance(cases,list) and cases,"saved application has no cases")
    identities=[case.get("id") for case in cases]
    _require(all(isinstance(i,str) and i for i in identities) and len(set(identities))==len(identities),
             "saved application has duplicate/invalid cases")
    selected=[case for case in cases if case.get("solver_family")=="CHERCHI_RAW_EQ14_V1"]
    _require(selected,"saved application has no Cherchi cases")
    report={"schema":"saved-cherchi-application-audit-v1","scope":SCOPE,
        "cases":[audit_saved_cherchi_case(case,bundle_root=bundle_root) for case in selected],
        "excluded_case_ids":[case["id"] for case in cases if case not in selected]}
    if value.get("schema")=="cherchi-top-finite-pareto-v1":
        report["pareto"]=_audit_finite_pareto(value,Path(bundle_root).resolve(strict=True),report["cases"])
    return report
