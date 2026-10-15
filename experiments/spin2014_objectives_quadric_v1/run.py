"""Frozen fixed-shape function runs. Never paper-original or manufacture claims.

No history/answer input, parameter overrides, resumed budget or backend fallback.
Each invocation owns a new directory; failed outputs remain available.
"""
from __future__ import annotations
import argparse
from dataclasses import asdict
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = Path(__file__).with_name("protocol.json")
PROTOCOL_SHA256 = "133279c176793316f35503ac33c522670a3f921258e2b485daef515809d2498e"


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _source_module(relative, expected_sha):
    path = ROOT / relative
    if _sha(path) != expected_sha:
        raise ValueError("frozen generator identity mismatch: " + relative)
    spec = importlib.util.spec_from_file_location("frozen_"+path.parent.name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_frozen_rotation(case_id, objective_token, seed, *, protocol_path=PROTOCOL):
    from make_it_stand.spin import SpinItFasterRequest, SpinSolverConfig, Spin2014QuadricObjectiveConfig
    if _sha(protocol_path) != PROTOCOL_SHA256:
        raise ValueError("frozen protocol identity mismatch")
    protocol = json.loads(Path(protocol_path).read_text())
    if type(seed) is not int or seed not in protocol["seeds"]:
        raise ValueError("unregistered initialization seed")
    matches = [x for x in protocol["cases"] if x["id"] == case_id]
    if len(matches) != 1:
        raise ValueError("unregistered case")
    case = matches[0]
    choices = [x for x in case["objectives"] if x["objective"] == objective_token]
    if len(choices) != 1 or objective_token == "STATIC_LOWERING":
        raise ValueError("unregistered rotational strategy; static needs its plane adapter")
    source_path = case["generator"].split(":")[0]
    source = _source_module(source_path, case["generator_sha256"])
    if case_id == "heart-original-axis-three-objectives":
        outer_v, outer_f = source.heart_mesh(case["samples"], case["layers"])
        center = np.asarray(case["design_center"])
        design_v = (outer_v-center)*case["design_scale"]+center
        design_f = outer_f.copy()
    elif case_id == "cuboid-offaxis-pairs-yoyo":
        import trimesh
        angle = math.radians(case["outer_rotation_y_degrees"])
        outer = source.rotated_box([0.,0.,0.], case["outer_extents"], angle)
        rotation = trimesh.transformations.rotation_matrix(angle, [0.,1.,0.])
        parts = [source.rotated_box((rotation @ np.r_[c,1.])[:3],
                                   case["design_box_extents"], angle)
                 for c in case["local_design_centers"]]
        design = trimesh.util.concatenate(parts)
        outer_v, outer_f, design_v, design_f = outer.vertices, outer.faces, design.vertices, design.faces
    else:
        raise ValueError("unsupported registered rotational source")
    config = Spin2014QuadricObjectiveConfig(**choices[0],
        **{name:case[name] for name in ("declared_length_unit", "declared_mass_unit",
            "solver_mass_unit_in_declared_units", "solver_length_unit_in_declared_units")})
    tuples = lambda array, scalar: tuple(tuple(scalar(x) for x in row) for row in array)
    request = SpinItFasterRequest(tuples(outer_v,float), tuples(outer_f,int),
        tuples(design_v,float), tuples(design_f,int),
        target_origin=case["target_origin"], target_axes=case["target_axes"],
        mode="YO_YO" if objective_token == "YOYO" else "SPINNING_TOP",
        evidence="RECONSTRUCTED_ASSET", material_density=case["material_density"],
        settings=SpinSolverConfig(equation_profile=protocol["newton_equation_profile"],
            equation_tolerance=protocol["raw_equation_tolerance"],
            maximum_restarts=protocol["maximum_attempts"],
            maximum_function_evaluations_per_restart=protocol["maximum_evaluations_per_attempt"],
            deterministic_seed=seed, capture_evaluation_payload=True))
    return request, config, case


def load_frozen_static(case_id, objective_token, seed, *, protocol_path=PROTOCOL):
    from make_it_stand.spin import SpinItFasterRequest, Spin2014QuadricObjectiveConfig
    from make_it_stand.stability_oracle import _convex_hull
    if _sha(protocol_path) != PROTOCOL_SHA256:
        raise ValueError("frozen protocol identity mismatch")
    protocol=json.loads(Path(protocol_path).read_text())
    if case_id != "static-faceted-shell-lowering" or objective_token != "STATIC_LOWERING" or \
            type(seed) is not int or seed not in protocol["seeds"]:
        raise ValueError("unregistered static input/strategy/seed")
    case=next(x for x in protocol["cases"] if x["id"]==case_id)

    def mesh(field, hash_field):
        path=ROOT/case[field]
        if _sha(path)!=case[hash_field]:
            raise ValueError("frozen static geometry identity mismatch: "+field)
        vertices,faces=[],[]
        for line in path.read_text().splitlines():
            row=line.split()
            if row and row[0]=="v":
                if len(row)!=4: raise ValueError("frozen OBJ requires xyz vertices")
                vertices.append(tuple(float(x) for x in row[1:]))
            elif row and row[0]=="f":
                if len(row)!=4: raise ValueError("frozen OBJ requires triangles")
                faces.append(tuple(int(x)-1 for x in row[1:]))
        return tuple(vertices),tuple(faces)

    outer,faces=mesh("outer","outer_sha256")
    design,inner_faces=mesh("design_domain","design_sha256")
    config=Spin2014QuadricObjectiveConfig(**case["objectives"][0],
        **{name:case[name] for name in ("declared_length_unit","declared_mass_unit",
            "solver_mass_unit_in_declared_units","solver_length_unit_in_declared_units")})
    request=SpinItFasterRequest(outer,faces,design,inner_faces,mode="STATIC_BALANCE",
        target_origin=case["target_origin"],target_axes=case["target_axes"],
        material_density=case["material_density"],evidence="RECONSTRUCTED_ASSET")
    # Support comes from complete actual coplanar bottom triangles, never AABB.
    bottom=case["support_plane_model_z"]
    points={outer[i] for face in faces if all(outer[i][2]==bottom for i in face) for i in face}
    hull=_convex_hull(np.asarray(sorted(points))[:,:2])
    support=tuple((float(x),float(y),float(bottom)) for x,y in hull)
    return request,config,case,support


def _write(path, value):
    """Publish an exclusive durable attachment before later checks proceed.

    A persistence error propagates; a partially written file remains diagnostic
    evidence and cannot be overwritten by a retry.
    """
    from make_it_stand.verification.faster_material import _json_value
    path = Path(path)
    with path.open("x", encoding="utf-8") as stream:
        json.dump(_json_value(value), stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def durable_rotational_event_sink(path):
    """Persist each actual native event before returning to the solver.

    This is an observation attachment, not a checkpoint or numerical pass.
    A cut may leave an incomplete trailing record; no missing event is inferred.
    Existing journals are never reopened as a new run or overwritten.
    """
    from make_it_stand.verification.faster_material import _json_value
    path = Path(path)
    with path.open("x", encoding="utf-8") as stream:
        stream.flush()
        os.fsync(stream.fileno())
    directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)
    sequence = 0
    def persist(event):
        nonlocal sequence
        row = dict(schema="spin2014-solver-event-journal-v1", sequence=sequence,
                   event=_json_value(event))
        encoded = json.dumps(row, sort_keys=True, allow_nan=False) + "\n"
        with path.open("a", encoding="utf-8") as stream:
            stream.write(encoded)
            stream.flush()
            os.fsync(stream.fileno())
        sequence += 1
    return persist


def frozen_rotation_evaluator(*, exact_source_boundary=False, exact_material_boundary=False,
                              rational_source_coordinates=False,
                              compact_rational_source_coordinates=False):
    """Explicit audit policy, fixed level and unchanged historical default.

    A true boundary option identifies a NEW run; it never upgrades an old
    failed result or changes its model, objective, units or numerical budget.
    """
    from make_it_stand.spin import SpinEvaluatorConfig
    return SpinEvaluatorConfig(derivative_method="DISCRETE_FULL_RANK_QUADRIC",
        icosphere_subdivisions=2, maximum_icosphere_subdivisions=2,
        exact_binary64_self_intersection=True,
        exact_binary64_source_boundary=exact_source_boundary,
        exact_binary64_material_boundary=exact_material_boundary,
        rational_source_coordinates=rational_source_coordinates,
        compact_rational_source_coordinates=compact_rational_source_coordinates)


def frozen_run_inputs(request, objective, case, seed, loaded, evaluator):
    """The exact strict-JSON parameters supplied to the existing evidence API."""
    from make_it_stand.verification.faster_material import _json_value
    return _json_value(dict(schema="spin2014-fixed-shape-function-run-v1", case=case["id"],
        protocol_sha256=PROTOCOL_SHA256, seed=seed, request=asdict(request),
        objective=asdict(objective), fixed_material_density=case["fixed_material_density"],
        evaluator=None if evaluator is None else evaluator.identity_dict(),
        native_path=str(loaded), native_sha256=_sha(loaded),
        initialization=("DETERMINISTIC_FROM_INPUT" if objective.objective=="STATIC_LOWERING"
                        else "REQUEST_AND_SEED_ONLY_NO_EXPLICIT_QUADRICS"),
        physical_validation="NOT_IN_SCOPE", manufacturing_admission="NOT_VERIFIED"))


def _scientific(value, *, canonical_cavity_cycles=False):
    """Retain science and failures; optionally normalize triangle start only.

    Legacy defaults remain unchanged. The explicit oriented-cycle policy never
    welds vertices, reorders faces, changes winding or removes geometry. Saved
    raw outputs remain separate evidence from this derived representation.
    """
    if isinstance(value, dict):
        result = {key:_scientific(x, canonical_cavity_cycles=canonical_cavity_cycles)
                  for key,x in value.items()
                  if key not in {"evaluation_timings", "evaluator_wall_seconds",
                    "residual_evaluation_seconds", "derivative_evaluation_seconds",
                    "event_hash", "output_hash"}}
        if canonical_cavity_cycles and result.get("cavity_faces") is not None:
            faces = result["cavity_faces"]
            if not isinstance(faces, list) or any(
                    not isinstance(face, list) or len(face) != 3 or
                    any(type(index) is not int or index < 0 for index in face) or
                    len(set(face)) != 3 for face in faces):
                raise ValueError("scientific cavity cycles require valid oriented triangles")
            result["cavity_faces"] = [list(min(tuple(face[k:]+face[:k])
                                             for k in range(3))) for face in faces]
        return result
    if isinstance(value, (tuple,list)):
        return [_scientific(x, canonical_cavity_cycles=canonical_cavity_cycles) for x in value]
    return value


def compare_authority_scientific_outputs(first_directory, repeat_directory):
    """Compare authoritative saved science; NEVER grant final repeat admission.

    Frozen inputs, exact coefficient strings, failed attempts and durable events
    are retained. Only the existing named runtime-field projection is allowed.
    Final numerical/material admission, build freshness and distinct process
    receipts must still be checked by the enclosing final audit.
    """
    from make_it_stand.verification.quadric_point_jobs import _read
    from make_it_stand.verification._spin2014_authority import check_authority_event_journal
    roots=[Path(first_directory).resolve(),Path(repeat_directory).resolve()]
    if roots[0]==roots[1]:
        raise ValueError('two distinct run directories required')
    records=[]
    for root in roots:
        paths=[root/name for name in ('frozen-input.json','native-run.json',
                                      'scientific-output.json','solver-events.jsonl')]
        before={str(path):_sha(path) for path in paths}
        frozen,native,scientific=(_read(path) for path in paths[:3])
        if frozen.get('coefficient_policy')!='COEFFICIENT_MPFR1024_V1' or \
           frozen.get('coefficient_precision_bits')!=1024 or \
           native.get('schema')!='spin2014-authority-restarts-v1' or \
           native.get('protocol_scope')!='CONFIGURED_FULL_ATTEMPT_BUDGET' or \
           scientific!=_scientific(native):
            raise ValueError('saved authority science or policy mismatch')
        events=native.get('events')
        if not isinstance(events,list) or not events or not isinstance(events[0],dict):
            raise ValueError('missing returned authority events')
        check_authority_event_journal(paths[3],run=native,
            expected_request_identity=native.get('request_identity'),
            expected_context_identity=native.get('context_identity'),
            expected_objective_identity=events[0].get('objective_identity'))
        if any(_sha(path)!=before[str(path)] for path in paths):
            raise ValueError('saved repeat inputs changed during comparison')
        records.append((frozen,scientific,before))
    canonical=lambda value: json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)
    return dict(schema='spin2014-authority-scientific-comparison-v1',
        scope='SCIENTIFIC_EQUALITY_ONLY_NOT_FINAL_REPEAT_ADMISSION',
        frozen_inputs_identical=records[0][0]==records[1][0],
        scientific_outputs_identical=records[0][1]==records[1][1],
        scientific_sha256=[hashlib.sha256(canonical(row[1]).encode()).hexdigest() for row in records],
        source_files=[row[2] for row in records],
        final_numeric_material_admission='NOT_PERFORMED',
        independent_process_identity='NOT_CHECKED',sealed_artifact=False)


def reinspect_saved_authority_complete_material(output_directory):
    """Read-only replay of full-run numerical and saved-material gates.

    Saved passed flags are not inputs. Request-only starts and the actual initial
    material are reconstructed, then the existing independent checkers run.
    This does not by itself certify repeats, build provenance or all five goals.
    """
    from make_it_stand.verification.quadric_point_jobs import _read
    from make_it_stand.spin import _native_module, _native_request, _native_evaluator_settings
    from make_it_stand.verification._spin2014_authority import check_authority_complete_run_evidence
    from make_it_stand.verification.spin_bundle import _load_spin2014_authority_point_attachment
    from make_it_stand.verification.spin_common import digest
    root=Path(output_directory).resolve()
    paths=[root/name for name in ('frozen-input.json','native-run.json','solver-events.jsonl')]
    if any(path.is_symlink() or not path.is_file() for path in paths[:2]):
        raise ValueError('regular saved full-run inputs required')
    before={str(path):_sha(path) for path in paths[:2]}
    frozen,payload=(_read(path) for path in paths[:2])
    if payload.get('schema')!='spin2014-authority-restarts-v1':
        raise ValueError('diagnostic point is not a complete authority run')
    try:
        request,objective,case=load_frozen_rotation(frozen['case'],frozen['objective']['objective'],frozen['seed'])
    except (KeyError,TypeError) as error:
        raise ValueError('malformed frozen full-run request') from error
    evaluator=frozen_rotation_evaluator(exact_source_boundary=True,exact_material_boundary=True,
        rational_source_coordinates=True,compact_rational_source_coordinates=True)
    native=_native_module();loaded=Path(native.__file__).resolve()
    expected=frozen_run_inputs(request,objective,case,frozen['seed'],loaded,evaluator)
    expected.update(coefficient_policy='COEFFICIENT_MPFR1024_V1',coefficient_precision_bits=1024)
    if frozen!=expected or payload.get('fixed_material_density')!=frozen['fixed_material_density']:
        raise ValueError('saved request differs from frozen inputs or loaded native product')
    if paths[2].is_symlink() or not paths[2].is_file():
        raise ValueError('durable whole-run journal required')
    before[str(paths[2])]=_sha(paths[2])
    agreement=audit_authority_event_attempt_agreement(payload)
    ledger=audit_authority_seed_ledger(request,objective,payload)
    if not ledger['full_run_candidate']:
        raise ValueError('actual optimization updates required')
    native,nr=_native_request(request)
    starts=native.spin2014_authority_initial_sequence(nr.mode,nr.settings)
    start=starts[payload['stationary_attempt_index']]
    settings=native.Spin2014QuadricObjectiveSettings()
    for key,value in asdict(objective).items():
        setattr(settings,key,getattr(native.Spin2014Objective,value) if key=='objective' else value)
    initial=native.evaluate_spin2014_authority_quadric_point(nr,settings,
        frozen['fixed_material_density'],start['geometry_authority'],start['dual_authority'],
        _native_evaluator_settings(native,evaluator))
    # Match the existing persisted initial-point JSON contract (native arrays
    # versus saved lists), without writing or changing any authority strings.
    from make_it_stand.verification.faster_material import _json_value
    initial=json.loads(json.dumps(_json_value(initial),allow_nan=False))
    numerical=check_authority_complete_run_evidence(payload,settings=asdict(objective),
        expected_initial_sequence=starts,journal_path=paths[2],
        expected_request_identity=initial['request_identity'],expected_context_identity=initial['context_identity'],
        expected_objective_identity=initial['objective_identity'],initial_point=initial)
    _,material=_load_spin2014_authority_point_attachment(root/'material-geometry',
        expected_point_digest=digest(payload['stationary_point']),
        expected_request_identity=initial['request_identity'],expected_context_identity=initial['context_identity'])
    if any(_sha(path)!=before[str(path)] for path in paths):
        raise ValueError('saved full-run files changed during reinspection')
    return dict(schema='spin2014-authority-saved-complete-reinspection-v1',
        scope='NUMERICAL_AND_SAVED_MATERIAL_NOT_REPEAT_OR_BUILD_ADMISSION',
        source_sha256=before,event_attempt_agreement=agreement,numerical=numerical,
        material=material,sealed_artifact=False,manufacturing_approved=False)


def audit_authority_repeat_processes(first_job, repeat_job):
    """Check saved whole-workflow receipts, not scientific/material admission.

    Distinct output directories alone are insufficient. This checks the existing
    source-job layout, completed bounded processes, guardian/config bindings and
    the fresh authority invocation. Receipts are same-user observations, not
    cryptographic authentication or compiled-source attestation.
    """
    import sys
    from make_it_stand.verification.quadric_point_jobs import _read
    from make_it_stand.verification.source_checkpoint_store import _fingerprints
    roots=[Path(path).resolve() for path in (first_job,repeat_job)]
    if roots[0]==roots[1]: raise ValueError('two distinct complete job roots required')
    observations=[]
    for root in roots:
        paths={name:root/name for name in ('config.json','result.json','complete.json','kernel.json')}
        paths['frozen-input.json']=root/'output/frozen-input.json'
        if any(path.is_symlink() or not path.is_file() for path in paths.values()):
            raise ValueError('regular complete process receipt files required')
        hashes={name:_sha(path) for name,path in paths.items()}
        try:
            config,result,complete,kernel,frozen=(_read(paths[name]) for name in
                ('config.json','result.json','complete.json','kernel.json','frozen-input.json'))
            memory=24*1024**3;wall=3600_000_000_000
            service=result['service'];unit=config['unit'];boot=config['boot_id']
            elapsed=result['elapsed_ns'];began=config['started_boottime_ns']
            deadline=config['deadline_boottime_ns'];observed=kernel['observed_boottime_ns']
            if config['schema']!='source-job-v1' or result['schema']!='source-job-result-v1' or \
               type(unit) is not str or not unit.strip() or type(boot) is not str or not boot.strip() or \
               result['unit']!=unit or result['status']!='PROCESS_EXITED' or \
               type(result['returncode']) is not int or result['returncode']!=0 or \
               complete['result_sha256']!=hashes['result.json'] or \
               any(type(value) is not int for value in (elapsed,began,deadline,observed)) or \
               not 0<=elapsed<wall or deadline-began!=wall or not began<=observed<=began+elapsed or \
               config['limits']['maximum_wall_ns']!=wall or \
               config['limits']['maximum_memory_bytes']!=memory or \
               service['Result']!='success' or service['SubState']!='exited' or \
               service['ExecMainCode']!='1' or service['ExecMainStatus']!='0' or \
               service['MemoryMax']!=str(memory) or service['MemorySwapMax']!='0' or \
               service['MemorySwapPeak']!='0' or kernel['memory_max']!=memory or \
               type(kernel['swap_max']) is not int or kernel['swap_max']!=0 or kernel['boot_id']!=boot or \
               kernel['absolute_deadline_boottime_ns']!=deadline or \
               type(kernel['guardian_pid']) is not int or kernel['guardian_pid']<=0 or \
               kernel['cgroup']!=service['ControlGroup'] or \
               Path(kernel['cgroup']).name!=unit:
                raise ValueError('incomplete or mismatched whole-workflow resource receipt')
            if _fingerprints(config['identity_files'])!=config['identity']:
                raise ValueError('STALE_RESULT: process effective dependencies changed')
            argv=config['argv'];runner=str(Path(__file__).resolve())
            if not isinstance(argv,list) or any(type(arg) is not str for arg in argv) or \
               argv.count(runner)!=1 or any(flag in argv for flag in
                   ('--saved-static-run','--saved-rotational-run')):
                raise ValueError('fresh authority function invocation required')
            prefix=argv[:argv.index(runner)+1]
            ordinary=[sys.executable,'-B',runner]
            preloaded=['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',*ordinary]
            if prefix not in (ordinary,preloaded):
                raise ValueError('receipt did not invoke the authority Python runner')
            expected={'--case':frozen['case'],'--objective':frozen['objective']['objective'],
                      '--seed':str(frozen['seed']),'--coefficient-policy':'COEFFICIENT_MPFR1024_V1',
                      '--output':str(root/'output')}
            for flag,value in expected.items():
                if argv.count(flag)!=1 or argv.index(flag)+1>=len(argv) or argv[argv.index(flag)+1]!=value:
                    raise ValueError('process invocation differs from frozen '+flag)
            if frozen['coefficient_policy']!='COEFFICIENT_MPFR1024_V1' or \
               type(frozen['coefficient_precision_bits']) is not int or frozen['coefficient_precision_bits']!=1024:
                raise ValueError('wrong authority process policy')
            bound=config['identity']
            for path,digest in ((runner,_sha(runner)),(frozen['native_path'],frozen['native_sha256'])):
                matches=[item for item in bound.values() if item.get('path')==path and item.get('sha256')==digest]
                if len(matches)!=1 or _sha(path)!=digest:
                    raise ValueError('runner or native product not bound to process')
        except (KeyError,TypeError,IndexError) as error:
            raise ValueError('malformed authority process receipt') from error
        if any(_sha(path)!=hashes[name] for name,path in paths.items()):
            raise ValueError('process receipts changed during inspection')
        observations.append(dict(job_directory=str(root),unit=unit,boot_id=boot,
            started_boottime_ns=began,elapsed_ns=elapsed,receipt_sha256=hashes,
            native_path=frozen['native_path'],native_sha256=frozen['native_sha256']))
    if (observations[0]['boot_id'],observations[0]['unit'])== \
       (observations[1]['boot_id'],observations[1]['unit']):
        raise ValueError('cloned process unit is not an independent repeat')
    if observations[0]['native_path']!=observations[1]['native_path'] or \
       observations[0]['native_sha256']!=observations[1]['native_sha256']:
        raise ValueError('repeat uses a different native product')
    return dict(schema='spin2014-authority-repeat-process-audit-v1',
        scope='PROCESS_RECEIPTS_ONLY_NOT_FINAL_REPEAT_ADMISSION',
        independent_bounded_jobs=True,processes=observations,
        scientific_material_admission='NOT_PERFORMED',
        same_environment_admission='NOT_PERFORMED',sealed_artifact=False)


def export_saved_static_material(source_directory, output_directory):
    """Reaudit/export the existing frozen candidate, without invoking a solver."""
    from make_it_stand.spin import _native_module
    from make_it_stand.verification.spin_common import parse, canonical
    from make_it_stand.verification.faster_material import _json_value, prepare_spin2014_static_material_output
    from make_it_stand.verification.spin_bundle import _static_source_run, save_spin2014_static_material_output
    from make_it_stand.verification.static_material import audit_spin2014_static_directory
    from make_it_stand.verification.run_evidence import begin_run, finish_run, native_build_configuration_files
    source=Path(source_directory).absolute()
    inputs=parse((source/"frozen-input.json").read_bytes())
    request,objective,case,support=load_frozen_static(inputs["case"],inputs["objective"]["objective"],inputs["seed"])
    if (inputs["protocol_sha256"]!=PROTOCOL_SHA256 or
        canonical(inputs["request"])!=canonical(_json_value(asdict(request))) or
        canonical(inputs["objective"])!=canonical(_json_value(asdict(objective))) or
        inputs["fixed_material_density"]!=case["fixed_material_density"]):
        raise ValueError("saved static run does not match the frozen function input")
    output=Path(output_directory).absolute()
    output.mkdir(exist_ok=False)
    native=_native_module()
    loaded=Path(native.__file__).resolve()
    expected=ROOT/"build/spin2014-objectives-quadric-20261007/python"
    if loaded.parent!=expected:
        raise RuntimeError("wrong loaded extension: "+str(loaded))
    assets=[source/name for name in ("frozen-input.json","native-run.json","input-geometry.npz")]
    parameters=dict(schema="spin2014-static-material-reaudit-input-v1",
        case=case["id"],source_files={p.name:_sha(p) for p in assets},
        optimizer_invocations=0,solver_freshness="HISTORICAL_RUN_NOT_REEXECUTED",
        protocol_sha256=PROTOCOL_SHA256)
    begin=begin_run(repository_root=ROOT,runner_files=[Path(__file__)],
        checker_files=[ROOT/"python/make_it_stand/verification"/name for name in
            ("static_material.py","mass_moments.py","faster_material.py","spin_problem.py","spin_bundle.py")]+
            [ROOT/"python/make_it_stand/stability_oracle.py"],asset_files=[*assets,PROTOCOL],
        build_configuration_files=native_build_configuration_files(expected.parent,loaded),
        parameters=parameters,seed_protocol={"seed":inputs["seed"],"protocol":"no optimization; saved-state geometry reaudit"},
        dependency_roots={"eigen":ROOT/"libigl/external/eigen","cgal":ROOT/"libigl/external/cgal"},
        native_modules=[native],pure_python=False)
    (output/"begin-evidence.json").write_text(begin.to_json()+"\n")
    audit=audit_spin2014_static_directory(source)
    if not audit["terminal_numerics_qualified"] or not audit["support"]["inside_support"]:
        raise ValueError("saved static numerical/support conditions are not independently qualified")
    source_run=_static_source_run(parse((source/"native-run.json").read_bytes()))
    geometry=prepare_spin2014_static_material_output(source_run)
    persisted=save_spin2014_static_material_output(geometry,source_run,output/"material-geometry")
    _write(output/"independent-numerical-audit.json",audit)
    summary=dict(schema="spin2014-static-material-reaudit-summary-v1",case=case["id"],
        optimizer_invocations=0,source_run_sha256=persisted.identity["source_run_sha256"],
        material_geometry_scientific_sha256=persisted.scientific_digest,
        material_export_complete=True,material_geometry_attachment_verified=True,
        solver_freshness="HISTORICAL_RUN_NOT_REEXECUTED",goal_complete=False,
        sealed_artifact=False,manufacturing_approved=False)
    _write(output/"summary.json",summary)
    outputs=[output/"summary.json",output/"independent-numerical-audit.json",
             *sorted(p for p in (output/"material-geometry").iterdir() if p.is_file())]
    end=finish_run(begin,events=[dict(iteration=0,status="ok",data=dict(summary=summary,audit=audit))],
                   output_files=outputs,measurements={"optimizer_invocations":0})
    (output/"run-evidence.json").write_text(end.to_json()+"\n")
    print(json.dumps(summary,sort_keys=True,allow_nan=False),flush=True)
    return 0


def _save_rotational_material(run, evaluator, output):
    """Shared actual material export for new solves and saved-state reaudit."""
    from make_it_stand.verification.spin_bundle import prepare_and_save_spin2014_quadric_material_output
    persisted=prepare_and_save_spin2014_quadric_material_output(run,output/"material-geometry",
        evaluator_settings=evaluator)
    _write(output/"independent-material-audit.json",persisted.identity)
    return persisted


def export_saved_rotational_material(source_directory, output_directory, *, exact_source_boundary=False,
                                    exact_material_boundary=False, rational_source_coordinates=False,
                                    compact_rational_source_coordinates=False):
    """Export/reload one frozen saved endpoint, never initialize an optimizer."""
    from make_it_stand.spin import SpinEvaluatorConfig, _native_module
    from make_it_stand.verification.spin_common import parse, canonical
    from make_it_stand.verification.faster_material import _json_value
    from make_it_stand.verification.spin_bundle import _rotational_source_run
    from make_it_stand.verification.run_evidence import begin_run, finish_run, native_build_configuration_files
    source=Path(source_directory).absolute()
    inputs=parse((source/"frozen-input.json").read_bytes())
    request,objective,case=load_frozen_rotation(inputs["case"],inputs["objective"]["objective"],inputs["seed"])
    evaluator=frozen_rotation_evaluator(exact_source_boundary=exact_source_boundary,
        exact_material_boundary=exact_material_boundary,
        rational_source_coordinates=rational_source_coordinates,
        compact_rational_source_coordinates=compact_rational_source_coordinates)
    if (inputs["protocol_sha256"]!=PROTOCOL_SHA256 or
        canonical(inputs["request"])!=canonical(_json_value(asdict(request))) or
        canonical(inputs["objective"])!=canonical(_json_value(asdict(objective))) or
        inputs["fixed_material_density"]!=case["fixed_material_density"] or
        canonical(inputs["evaluator"])!=canonical(_json_value(evaluator.identity_dict())) or
        inputs["initialization"]!="REQUEST_AND_SEED_ONLY_NO_EXPLICIT_QUADRICS"):
        raise ValueError("saved rotational run does not match the frozen function input")
    run=_rotational_source_run(parse((source/"native-run.json").read_bytes()))
    if (run.objective_settings!=objective or run.fixed_material_density!=case["fixed_material_density"] or
        any(canonical(_json_value(asdict(run.measured_request)[key]))!=canonical(value)
            for key,value in inputs["request"].items() if key!="fixed_shell_moments")):
        raise ValueError("saved numerical run differs from the frozen request/objective")
    with np.load(source/"input-geometry.npz",allow_pickle=False) as geometry:
        names=("outer_vertices","outer_faces","design_vertices","design_faces")
        if set(geometry.files)!=set(names):
            raise ValueError("saved input geometry fields mismatch")
        for name in names:
            dtype=np.int64 if name.endswith("faces") else np.float64
            one,two=np.asarray(geometry[name],dtype=dtype),np.asarray(getattr(request,name),dtype=dtype)
            if one.shape!=two.shape or one.tobytes()!=two.tobytes():
                raise ValueError("saved input geometry mismatch: "+name)
    output=Path(output_directory).absolute()
    output.mkdir(exist_ok=False)
    native=_native_module();loaded=Path(native.__file__).resolve()
    expected=ROOT/"build/spin2014-objectives-quadric-20261007/python"
    if loaded.parent!=expected:
        raise RuntimeError("wrong loaded extension: "+str(loaded))
    assets=[source/name for name in ("frozen-input.json","native-run.json","input-geometry.npz")]
    parameters=dict(schema="spin2014-rotational-material-reaudit-input-v1",case=case["id"],
        source_files={path.name:_sha(path) for path in assets},optimizer_invocations=0,
        solver_freshness="HISTORICAL_RUN_NOT_REEXECUTED",protocol_sha256=PROTOCOL_SHA256)
    begin=begin_run(repository_root=ROOT,runner_files=[Path(__file__)],
        checker_files=[ROOT/"python/make_it_stand/verification"/name for name in
            ("mass_moments.py","faster_material.py","spin_problem.py","spin_bundle.py")],
        asset_files=[*assets,PROTOCOL],
        build_configuration_files=native_build_configuration_files(expected.parent,loaded),parameters=parameters,
        seed_protocol={"seed":inputs["seed"],"protocol":"saved-state reaudit; no optimization"},
        dependency_roots={"eigen":ROOT/"libigl/external/eigen","cgal":ROOT/"libigl/external/cgal"},
        native_modules=[native],pure_python=False)
    (output/"begin-evidence.json").write_text(begin.to_json()+"\n")
    persisted=_save_rotational_material(run,evaluator,output)
    summary=dict(schema="spin2014-rotational-material-reaudit-summary-v1",case=case["id"],
        objective=objective.objective,seed=inputs["seed"],optimizer_invocations=0,
        source_run_sha256=persisted.identity["source_run_sha256"],
        material_geometry_scientific_sha256=persisted.scientific_digest,
        material_geometry_attachment_verified=True,material_export_complete=True,
        accepted_updates=persisted.audit["accepted_updates"],
        independent_raw_residual=persisted.audit["independent_raw_residual"],
        solver_freshness="HISTORICAL_RUN_NOT_REEXECUTED",goal_complete=False,
        sealed_artifact=False,manufacturing_approved=False)
    _write(output/"summary.json",summary)
    outputs=[output/"summary.json",output/"independent-material-audit.json",
        *sorted(path for path in (output/"material-geometry").iterdir() if path.is_file())]
    end=finish_run(begin,events=[dict(iteration=0,status="ok",data=summary)],
        output_files=outputs,measurements={"optimizer_invocations":0})
    (output/"run-evidence.json").write_text(end.to_json()+"\n")
    print(json.dumps(summary,sort_keys=True,allow_nan=False),flush=True)
    return 0


def solver_run_payload(run):
    """Preserve authority strings; serialization does not grant admission."""
    from make_it_stand.verification.faster_material import _json_value
    if isinstance(run,dict):
        if run.get("schema") != "spin2014-authority-restarts-v1":
            raise ValueError("unsupported solver run mapping")
        return _json_value(run)
    return _json_value(asdict(run))


def authority_audit_outputs(output):
    names=('independent-event-journal-audit.json','independent-seed-ledger-audit.json',
           'independent-event-attempt-agreement.json',
           'material-change-moment-witness.json','authority-initial-point.json',
           'complete-authority-numerical-admission.json')
    return [output/name for name in names if (output/name).is_file()]


def audit_authority_event_attempt_agreement(payload):
    """Reconcile observation records with returned attempts, not convergence.

    Journal byte/identity validation and independent numerical/material admission
    remain separate required checks. This guard does not authenticate a solver.
    """
    events=payload.get('events'); attempts=payload.get('attempts')
    if not isinstance(events,list) or not isinstance(attempts,list) or not attempts:
        raise ValueError('complete events and attempts required')
    cursor=0
    def take(kind,restart,index,**expected):
        nonlocal cursor
        if cursor>=len(events): raise ValueError('event prefix is incomplete')
        event=events[cursor]
        required=dict(type=kind,restart=restart,evaluation_index=index,sequence=cursor,**expected)
        if any(event.get(key)!=value for key,value in required.items()):
            raise ValueError('event/attempt disagreement at '+str(cursor)+' ('+kind+')')
        for key in ('restart','evaluation_index','sequence'):
            if type(event.get(key)) is not int:
                raise ValueError('event index must be an integer')
        cursor+=1
    total=0
    final_accepted_material=None
    for restart,attempt in enumerate(attempts):
        trials=attempt['trials']; accepted=attempt['accepted_state_identities']
        count=attempt['function_evaluations']; updates=attempt['accepted_updates']
        if type(count) is not int or count!=1+len(trials) or \
           type(updates) is not int or updates!=len(accepted):
            raise ValueError('attempt evaluation/update accounting mismatch')
        take('attempt_start',restart,0,
             geometry_authority=attempt['initial_geometry_authority'],
             dual_authority=attempt['initial_dual_authority'])
        take('initial_evaluation_start',restart,1)
        initial_present=bool(attempt['initial_state_identity'])
        initial_status='' if initial_present else 'initial_geometry_or_derivative_failure'
        take('initial_evaluation_terminal',restart,1,
             state_identity=attempt['initial_state_identity'],
             material_identity=attempt['initial_material_identity'],
             raw_residual=attempt['initial_raw_residual'],status=initial_status,
             diagnostic='' if initial_present else attempt['diagnostic'])
        committed=0
        base_state=attempt['initial_state_identity']
        for index,trial in enumerate(trials,2):
            point=dict(geometry_authority=trial['geometry_authority'],
                       dual_authority=trial['dual_authority'],arc_length=trial['arc_length'],
                       state_identity=base_state)
            take('trial_evaluation_start',restart,index,**point)
            if type(trial['accepted']) is not bool:
                raise ValueError('trial acceptance must be a boolean')
            take('trial_evaluation_terminal',restart,index,**point,
                 status=trial['reason'],diagnostic=trial['diagnostic'],
                 raw_residual=trial['raw_equation_residual'],
                 trial_identity=trial['trial_identity'],accepted=trial['accepted'])
            if trial['accepted']:
                if committed>=len(accepted): raise ValueError('accepted trial omitted from attempt')
                take('state_accepted',restart,index,
                     geometry_authority=trial['geometry_authority'],
                     dual_authority=trial['dual_authority'],state_identity=accepted[committed],
                     raw_residual=trial['raw_equation_residual'],accepted=True)
                final_accepted_material=events[cursor-1].get('material_identity')
                base_state=accepted[committed]
                committed+=1
        if committed!=updates: raise ValueError('accepted event/update count mismatch')
        take('attempt_terminal',restart,count,status=attempt['status'],
             diagnostic=attempt['diagnostic'],state_identity=attempt['final_state_identity'],
             raw_residual=attempt['final_raw_residual'])
        total+=count
    if type(payload.get('function_evaluations')) is not int or payload['function_evaluations']!=total:
        raise ValueError('run evaluation total mismatch')
    take('run_terminal',len(attempts),total,status=payload['status'])
    if cursor!=len(events): raise ValueError('unexpected events after run terminal')
    if payload['status']=='stationary_candidate' and attempts[-1]['accepted_updates']:
        terminal=payload.get('stationary_point')
        if not isinstance(terminal,dict) or \
           final_accepted_material!=terminal.get('material_scientific_hash'):
            raise ValueError('last accepted material differs from stationary material')
    return dict(scope='EVENT_ATTEMPT_AGREEMENT_NOT_FINAL_ADMISSION',
                events_match_attempts=True,event_count=cursor,attempt_count=len(attempts),
                function_evaluations=total,sealed_artifact=False)


def record_authority_seed_audit(request, objective, payload, output, summary):
    """Persist structural observations once; never grant final admission here."""
    try:
        ledger_check=audit_authority_seed_ledger(request,objective,payload)
        if ledger_check['full_run_candidate']:
            from make_it_stand.verification._spin2014_authority import check_authority_material_change
            change=check_authority_material_change(
                payload['attempts'][payload['stationary_attempt_index']],payload['stationary_point'])
            _write(output/'material-change-moment-witness.json',change)
            summary['material_change_moment_witness']=change['material_change_witness']
    except Exception as error:
        summary['seed_and_restart_ledger_verified']=False
        _write(output/'independent-seed-ledger-audit.json',dict(
            scope='STRUCTURAL_RESTART_LEDGER_NOT_FINAL_ADMISSION',
            exception_type=type(error).__name__,diagnostic=str(error),sealed_artifact=False))
    else:
        # Publish success only after downstream observations succeed. Keep
        # exclusive creation: an existing audit must never be overwritten.
        _write(output/'independent-seed-ledger-audit.json',ledger_check)
        summary['seed_and_restart_ledger_verified']=True


def export_authority_complete_material(request, objective, evaluator, payload, output):
    """Reconstruct initial material, independently admit, then save/reload terminal.

    This is an algorithm material attachment, never manufacturing approval.
    Only a full returned candidate can invoke this expensive path.
    """
    if payload.get('status')!='stationary_candidate' or payload.get('protocol_scope')!='CONFIGURED_FULL_ATTEMPT_BUDGET':
        raise ValueError('stationary full run required')
    if json.loads((output/'native-run.json').read_text())!=payload:
        raise ValueError('saved native run differs from candidate')
    source_hash=_sha(output/'native-run.json')
    frozen_path=output/'frozen-input.json'
    frozen_hash=_sha(frozen_path)
    frozen=json.loads(frozen_path.read_text())
    returned_density=payload.get('fixed_material_density')
    frozen_density=frozen.get('fixed_material_density')
    try:
        density=Fraction(frozen_density)
        if isinstance(frozen_density,bool) or isinstance(returned_density,bool) or density<=0 or Fraction(returned_density)!=density:
            raise ValueError('density mismatch')
    except (TypeError,ValueError,OverflowError,ZeroDivisionError) as error:
        raise ValueError('fixed density differs from frozen request') from error
    audit_authority_event_attempt_agreement(payload)
    ledger=audit_authority_seed_ledger(request,objective,payload)
    if not ledger['full_run_candidate']:
        raise ValueError('actual optimization updates required')
    from make_it_stand.spin import _native_request, _native_evaluator_settings
    from make_it_stand.verification._spin2014_authority import check_authority_complete_run_evidence
    from make_it_stand.verification.spin_bundle import _save_spin2014_authority_point_attachment
    native,nr=_native_request(request)
    starts=native.spin2014_authority_initial_sequence(nr.mode,nr.settings)
    start=starts[payload['stationary_attempt_index']]
    settings=native.Spin2014QuadricObjectiveSettings()
    for key,value in asdict(objective).items():
        setattr(settings,key,getattr(native.Spin2014Objective,value) if key=='objective' else value)
    initial=native.evaluate_spin2014_authority_quadric_point(nr,settings,
        frozen_density,start['geometry_authority'],start['dual_authority'],
        _native_evaluator_settings(native,evaluator))
    _write(output/'authority-initial-point.json',initial)
    initial=json.loads((output/'authority-initial-point.json').read_text())
    admission=check_authority_complete_run_evidence(payload,settings=asdict(objective),
        expected_initial_sequence=starts,journal_path=output/'solver-events.jsonl',
        expected_request_identity=initial['request_identity'],expected_context_identity=initial['context_identity'],
        expected_objective_identity=initial['objective_identity'],initial_point=initial)
    _write(output/'complete-authority-numerical-admission.json',admission)
    regions,audit=_save_spin2014_authority_point_attachment(payload['stationary_point'],asdict(objective),
        output/'material-geometry',expected_request_identity=initial['request_identity'],
        expected_context_identity=initial['context_identity'],arithmetic_backend='GMPY2_2_2_1_EXACT_V1')
    if _sha(output/'native-run.json')!=source_hash or _sha(frozen_path)!=frozen_hash:
        raise ValueError('saved run or frozen request changed during material admission')
    return audit


def audit_authority_seed_ledger(request, objective, payload):
    """Reconstruct request-only starts; structural evidence never grants a seal."""
    from make_it_stand.spin import _native_request
    from make_it_stand.verification._spin2014_authority import check_authority_restart_ledger
    native,native_request=_native_request(request)
    starts=native.spin2014_authority_initial_sequence(native_request.mode,native_request.settings)
    return check_authority_restart_ledger(payload,objective_settings=asdict(objective),
        expected_initial_sequence=starts,evaluation_budget_per_attempt=60,
        expected_request_identity=payload['request_identity'],
        expected_context_identity=payload['context_identity'])


def authority_run_accounting(payload):
    """Report returned cost, never infer zero work from an interrupted return."""
    events=payload.get("events",[])
    complete=payload.get("schema")=="spin2014-authority-restarts-v1" and \
        bool(events) and events[-1].get("type")=="run_terminal"
    if not complete:
        return dict(native_return_complete=False,function_evaluations=None,
                    restart_function_evaluations=None,accepted_steps=None)
    attempts=payload["attempts"]
    return dict(native_return_complete=True,function_evaluations=payload["function_evaluations"],
                restart_function_evaluations=[a["function_evaluations"] for a in attempts],
                accepted_steps=sum(a["accepted_updates"] for a in attempts))


def authority_event_outcome(event):
    """Keep rejected trials and terminal failures in the evidence ledger."""
    if event["type"] in ("attempt_terminal","run_terminal"):
        return "ok" if event["status"]=="stationary_candidate" else "failed"
    if event["type"]=="trial_evaluation_terminal":
        return "ok" if event["status"] in ("residual_decreased","topology_transition_relinearized") else "failed"
    if event["status"]=="initial_geometry_or_derivative_failure":
        return "failed"
    return "ok"


def main():
    from make_it_stand.spin import (SpinEvaluatorConfig, _native_module,
        solve_spin2014_quadric_hollowing, solve_spin2014_static_quadric_hollowing)
    from make_it_stand.verification.run_evidence import begin_run, finish_run, native_build_configuration_files
    from make_it_stand.verification.faster_material import _json_value
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case")
    parser.add_argument("--objective")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--coefficient-policy", choices=("COEFFICIENT_MPFR1024_V1",),
        help="Explicit authority solver; old invocation retains historical behavior")
    parser.add_argument("--saved-static-run", type=Path)
    parser.add_argument("--saved-rotational-run", type=Path)
    parser.add_argument("--exact-source-boundary", action="store_true",
        help="Explicit source-only complete exact audit; records a distinct evaluator identity")
    parser.add_argument("--exact-material-boundary", action="store_true",
        help="Explicit complete exact signed-shell solver-material proof; not manufacturing approval")
    parser.add_argument("--rational-source-coordinates", action="store_true",
        help="Explicit enclosed rational production source profile; changes evaluator/run identity")
    parser.add_argument("--compact-rational-source-coordinates", action="store_true",
        help="Enclosed dyadic representative; requires --rational-source-coordinates")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    authority = args.coefficient_policy is not None
    if authority and (args.saved_static_run is not None or args.saved_rotational_run is not None or
                      args.objective == "STATIC_LOWERING"):
        parser.error("coefficient policy requires a fresh rotational run")
    if args.compact_rational_source_coordinates and not args.rational_source_coordinates:
        parser.error("compact rational source requires --rational-source-coordinates")
    rotational_policy = (args.exact_source_boundary or args.exact_material_boundary or
                         args.rational_source_coordinates or args.compact_rational_source_coordinates)
    if args.saved_static_run is not None and args.saved_rotational_run is not None:
        parser.error("saved-state export requires exactly one source")
    if args.saved_rotational_run is not None:
        if any(x is not None for x in (args.case,args.objective,args.seed)):
            parser.error("saved-state export cannot override case, objective or seed")
        return export_saved_rotational_material(args.saved_rotational_run,args.output,
            exact_source_boundary=args.exact_source_boundary,
            exact_material_boundary=args.exact_material_boundary,
            rational_source_coordinates=args.rational_source_coordinates,
            compact_rational_source_coordinates=args.compact_rational_source_coordinates)
    if args.saved_static_run is not None:
        if rotational_policy:
            parser.error("exact source/material audit policy applies only to rotational runs")
        if any(x is not None for x in (args.case,args.objective,args.seed)):
            parser.error("saved-state export cannot override case, objective or seed")
        return export_saved_static_material(args.saved_static_run,args.output)
    if any(x is None for x in (args.case,args.objective,args.seed)):
        parser.error("optimization requires --case, --objective and --seed")
    static=args.objective=="STATIC_LOWERING"
    if static and rotational_policy:
        parser.error("exact source/material audit policy applies only to rotational runs")
    support=None
    if static:
        request,objective,case,support=load_frozen_static(args.case,args.objective,args.seed)
    else:
        request, objective, case = load_frozen_rotation(args.case, args.objective, args.seed)
    output = args.output.absolute()
    output.mkdir(exist_ok=False)
    native = _native_module()
    loaded = Path(native.__file__).resolve()
    expected = ROOT/"build/spin2014-objectives-quadric-20261007/python"
    if loaded.parent != expected:
        raise RuntimeError("wrong loaded extension: " + str(loaded))
    evaluator = None if static else frozen_rotation_evaluator(
        exact_source_boundary=args.exact_source_boundary or authority,
        exact_material_boundary=args.exact_material_boundary or authority,
        rational_source_coordinates=args.rational_source_coordinates or authority,
        compact_rational_source_coordinates=args.compact_rational_source_coordinates or authority)
    inputs = frozen_run_inputs(request, objective, case, args.seed, loaded, evaluator)
    if authority:
        inputs.update(coefficient_policy=args.coefficient_policy,coefficient_precision_bits=1024)
    if support is not None:
        inputs["actual_support_hull"]=[list(x) for x in support]
    _write(output/"frozen-input.json", inputs)
    np.savez(output/"input-geometry.npz", outer_vertices=request.outer_vertices,
        outer_faces=request.outer_faces, design_vertices=request.design_vertices, design_faces=request.design_faces)
    begin = begin_run(repository_root=ROOT, runner_files=[Path(__file__)],
        checker_files=[ROOT/"python/make_it_stand/verification"/name for name in
            ("mass_moments.py", "faster_material.py", "spin_problem.py", "spin_bundle.py",
             "_spin2014_authority.py", "rational_text.py")],
        asset_files=[PROTOCOL, output/"frozen-input.json", output/"input-geometry.npz"],
        build_configuration_files=native_build_configuration_files(expected.parent, loaded),
        parameters=inputs, seed_protocol={"seed":args.seed,"protocol":"frozen request-only deterministic quadrics"},
        dependency_roots={"eigen":ROOT/"libigl/external/eigen", "cgal":ROOT/"libigl/external/cgal"},
        native_modules=[native], pure_python=False)
    (output/"begin-evidence.json").write_text(begin.to_json()+"\n")
    print(json.dumps(dict(stage="START_SOLVE", case=args.case, objective=args.objective,
                         seed=args.seed, explicit_initial_guesses=0)), flush=True)
    journal = None if static else output/"solver-events.jsonl"
    event_sink = None if static else durable_rotational_event_sink(journal)
    started = time.monotonic()
    try:
        if static:
            run=solve_spin2014_static_quadric_hollowing(request,objective,
                fixed_material_density=case["fixed_material_density"])
        else:
            run = solve_spin2014_quadric_hollowing(request, objective,
                fixed_material_density=case["fixed_material_density"], evaluator_settings=evaluator,
                event_sink=event_sink, coefficient_policy=args.coefficient_policy)
        payload = solver_run_payload(run)
    except Exception as error:
        # Preserve admission/construction failures without calling them infeasible.
        payload = dict(schema="spin2014-quadric-admission-failure-v1",
            exception_type=type(error).__name__, diagnostic=str(error),
            feasibility="UNKNOWN", material_admission="NOT_VERIFIED")
    elapsed = time.monotonic()-started
    _write(output/"native-run.json", payload)
    scientific = _scientific(payload)
    _write(output/"scientific-output.json", scientific)
    authority_return = authority and payload.get("schema") == "spin2014-authority-restarts-v1"
    solve = payload if authority_return else payload.get("candidate" if static else "solve", {})
    summary = dict(case=args.case, objective=args.objective, seed=args.seed,
        status=solve.get("status", "INPUT_OR_CONSTRUCTION_FAILED"),
        diagnostic=solve.get("diagnostic",payload.get("diagnostic")),
        function_evaluations=solve.get("function_evaluations",0),
        restart_function_evaluations=solve.get("restart_function_evaluations",[]),
        final_geometry_verified=solve.get("final_geometry_verified",False),
        independent_raw_residual=solve.get("independent_raw_residual"),
        material_export_complete=payload.get("material_export_complete",False),
        material_seal="PENDING_INDEPENDENT_PARTITION_AND_RELOAD",
        goal_complete=False, elapsed_seconds=elapsed,
        scientific_sha256=_sha(output/"scientific-output.json"))
    if authority:
        summary.update(coefficient_policy=args.coefficient_policy,
            solver_family="MPFR1024_ACTUAL_MATERIAL_NEWTON_V1",
            material_geometry_attachment_verified=False,manufacturing_approved=False,
            material_seal="PENDING_COMPLETE_AUTHORITY_RUN_ADMISSION")
        summary.update(authority_run_accounting(payload))
        if authority_return:
            record_authority_seed_audit(request,objective,payload,output,summary)
            from make_it_stand.verification._spin2014_authority import check_authority_event_journal
            try:
                agreement=audit_authority_event_attempt_agreement(payload)
                _write(output/'independent-event-attempt-agreement.json',agreement)
                journal_check=check_authority_event_journal(journal,run=payload,
                    expected_request_identity=payload["request_identity"],
                    expected_context_identity=payload["context_identity"],
                    expected_objective_identity=payload["events"][0]["objective_identity"])
                _write(output/"independent-event-journal-audit.json",journal_check)
                summary["durable_event_stream_verified"]=True
            except Exception as error:
                summary["durable_event_stream_verified"]=False
                _write(output/"independent-event-journal-audit.json",dict(
                    scope="DURABLE_OBSERVATION_NOT_FINAL_ADMISSION",exception_type=type(error).__name__,
                    diagnostic=str(error),completed_stream_matches=False,sealed_artifact=False))
            if summary.get('durable_event_stream_verified') and summary.get('seed_and_restart_ledger_verified') and solve.get('status')=='stationary_candidate':
                try:
                    audit=export_authority_complete_material(request,objective,evaluator,payload,output)
                    summary.update(material_geometry_attachment_verified=True,
                        material_seal='COMPLETE_RUN_EXACT_ATTACHMENT_NOT_MANUFACTURING_SEAL',
                        material_geometry_scientific_sha256=audit['source_point_digest'])
                except Exception as error:
                    failure=dict(schema='spin2014-authority-complete-material-failure-v1',
                        exception_type=type(error).__name__,diagnostic=str(error),
                        material_verified=False,manufacturing_approved=False)
                    _write(output/'material-export-diagnostic.json',failure)
                    summary['material_export_diagnostic']=failure
    elif static:
        trace=solve.get("trace",[])
        summary.update(solver_family="PLANE_FAMILY_FRESH_2014_SQUARE_V1",
            function_evaluations=max((x["evaluation"] for x in trace),default=-1)+1,
            accepted_steps=solve.get("accepted_steps",0),
            fresh_squared_height_verified=payload.get("squared_height",{}).get(
                "squared_height_conditions_verified",False),
            fresh_squared_height_diagnostic=payload.get("squared_height",{}).get("diagnostic"))
        summary.update(material_geometry_attachment_verified=False, manufacturing_approved=False)
        if solve.get("status")=="CONVERGED" and summary["fresh_squared_height_verified"]:
            try:
                from make_it_stand.verification.faster_material import prepare_spin2014_static_material_output
                from make_it_stand.verification.spin_bundle import save_spin2014_static_material_output
                geometry=prepare_spin2014_static_material_output(run)
                persisted=save_spin2014_static_material_output(geometry,run,output/"material-geometry")
                summary.update(material_export_complete=True, material_geometry_attachment_verified=True,
                    material_geometry_scientific_sha256=persisted.scientific_digest,
                    material_seal="GEOMETRY_ATTACHMENT_NOT_MANUFACTURING_SEAL")
            except Exception as error:
                # Never suppress the numerical outcome or replace a failed
                # geometry export with a repaired mesh. Preserve any partial
                # exclusive bundle and the actual exception as separate evidence.
                failure=dict(schema="spin2014-static-material-export-failure-v1",
                    exception_type=type(error).__name__,diagnostic=str(error),
                    manufacturing_approved=False,material_verified=False)
                _write(output/"material-export-diagnostic.json",failure)
                summary["material_export_diagnostic"]=failure
    else:
        summary.update(material_geometry_attachment_verified=False,manufacturing_approved=False)
        if solve.get("status")=="CONVERGED" and solve.get("final_geometry_verified"):
            try:
                persisted=_save_rotational_material(run,evaluator,output)
                summary.update(material_geometry_attachment_verified=True,
                    material_geometry_scientific_sha256=persisted.scientific_digest,
                    material_seal="GEOMETRY_ATTACHMENT_NOT_MANUFACTURING_SEAL")
            except Exception as error:
                failure=dict(schema="spin2014-rotational-material-export-failure-v1",
                    exception_type=type(error).__name__,diagnostic=str(error),
                    manufacturing_approved=False,material_verified=False)
                _write(output/"material-export-diagnostic.json",failure)
                summary["material_export_diagnostic"]=failure
    _write(output/"summary.json", summary)
    if static:
        rows=scientific.get("candidate",{}).get("trace",[])
        events=[dict(iteration=i,status="ok" if row["evaluation_succeeded"] else "failed",data=row)
                for i,row in enumerate(rows)]
    elif authority:
        rows = scientific.get("events", [])
        events = [dict(iteration=i,status=authority_event_outcome(row),data=row)
                  for i,row in enumerate(rows)]
    else:
        rows = scientific.get("solve", {}).get("events", [])
        events = [dict(iteration=i, status="failed" if row["outcome"] in ("failed","rejected") else "ok",
                       data=row) for i,row in enumerate(rows)]
    if not events:
        events = [dict(iteration=0,status="failed",data=scientific)]
    outputs=[output/"native-run.json",output/"scientific-output.json",output/"summary.json"]
    if journal is not None:
        outputs.append(journal)
    outputs.extend(authority_audit_outputs(output))
    if (output/"material-geometry").is_dir():
        outputs.extend(sorted(p for p in (output/"material-geometry").iterdir() if p.is_file()))
    if (output/"material-export-diagnostic.json").is_file():
        outputs.append(output/"material-export-diagnostic.json")
    if (output/"independent-material-audit.json").is_file():
        outputs.append(output/"independent-material-audit.json")
    end = finish_run(begin, events=events,
        output_files=outputs,
        measurements={"elapsed_seconds":elapsed})
    (output/"run-evidence.json").write_text(end.to_json()+"\n")
    print(json.dumps(summary, sort_keys=True, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
