"""Run the frozen case through the existing seeded formal Cherchi entry."""
import hashlib
import json
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent


def registered_solver(settings,seed):
    expected=dict(maximum_iterations=60,maximum_restarts=8,
        maximum_function_evaluations_per_restart=60,equation_tolerance=1e-12,
        equation_profile='DUAL_MONOMIAL_RAW_V1',capture_evaluation_payload=True)
    allowed=set(expected)|{'physical_relative_tolerance','deterministic_seed'}
    if seed not in (0,1) or set(settings)-allowed or any(settings.get(k)!=v for k,v in expected.items()):
        raise ValueError('FROZEN_CASE_SOLVER_PROTOCOL_CHANGED')
    return dict(settings,deterministic_seed=seed)


def main(output,seed):
    import numpy as np
    import trimesh
    from make_it_stand import spin
    from make_it_stand.verification.spin_problem import SpinArtifactProblem
    from make_it_stand.verification.spin_bundle import save_spin_artifact
    from make_it_stand.verification.cherchi_inprocess import solve_cherchi_inprocess_and_seal_faster
    from noodles_formal_case import SEAL_WALL_M
    load=lambda p:json.loads(p.read_text())
    prepared=load(HERE/'noodles-formal-prepare-job-01/output/report.json')
    if not prepared['passed']:raise ValueError('FORMAL_CASE_PREPARATION_NOT_PASSED')
    request_path=HERE/'noodles-balanced-request-freeze-job-01/output/request.json'
    frozen=load(request_path)
    for name,digest in frozen['input_sources'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest()!=digest:raise ValueError('FROZEN_INPUT_CHANGED')
    old=frozen['requests'][seed];solver=spin.SpinSolverConfig(**registered_solver(old['settings'],seed))
    paths={Path(p).name:Path(p) for p in frozen['input_sources']}
    meshes={r:trimesh.load(paths[r+'-si-candidate.obj'],force='mesh',process=False) for r in ('outer','design','fixed')}
    spec=load(next(Path(p) for p in frozen['input_sources'] if Path(p).parent.name=='noodles-dense-v5'))['spec']
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic_ns()
    problem=SpinArtifactProblem(original=(meshes['outer'].vertices,meshes['outer'].faces),
        fixed_regions=[('fixed-tubes-handle-tip',(meshes['fixed'].vertices,meshes['fixed'].faces),spec['density_kg_m3'])],
        design_regions=[('design',(meshes['design'].vertices,meshes['design'].faces),spec['density_kg_m3'])],
        length_unit='m',mass_unit='kg',minimum_outer_wall=SEAL_WALL_M,
        target=dict(mode='SPINNING_TOP',axis=spec['target_axis'],
                    origin=(np.asarray(spec['target_origin_model'])*spec['model_unit_m']).tolist(),support_points=[]))
    if problem.input_digest!=prepared['input_digest']:raise ValueError('FORMAL_PREPARED_INPUT_MISMATCH')
    evaluator=spin.SpinEvaluatorConfig(**dict(frozen['evaluator_settings'],algebraic_finite_map=True))
    (output/'registration.json').write_text(json.dumps(dict(seed=seed,input_digest=problem.input_digest,
        request_sha256=hashlib.sha256(request_path.read_bytes()).hexdigest(),
        classification='RECONSTRUCTED_METHOD_VALIDATION',backend='CHERCHI',
        initialization='REQUEST_SEED_ONLY_NO_EXPLICIT_QUADRICS',seal_wall_m=SEAL_WALL_M,
        weight_contract='FASTER_NORMALIZED_REQUEST',center_of_mass_weight=old['center_of_mass_weight'])))
    result=solve_cherchi_inprocess_and_seal_faster(problem,evaluator_settings=evaluator,solver_settings=solver,
        output_directory=output/'run',center_of_mass_weight=old['center_of_mass_weight'],
        spinning_top_objective=old['spinning_top_objective'],rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1',
        solver_boundary_profile='EXACT_RATIONAL_FIXED_SUBSET_V1',durability_profile='EVALUATION_CHECKPOINT_V1',
        validation_profile='LAYERED_SAME_CALL_V1')
    returned=time.monotonic_ns()
    events=[json.loads(line) for line in (output/'run/events.jsonl').read_text().splitlines()]
    begin=next(e for e in events if e.get('kind')=='STAGE_STARTED' and e.get('stage')=='native_newton_solver')
    end=next(e for e in events if e.get('kind')=='STAGE_RETURNED' and e.get('stage')=='native_newton_solver')
    native=result.identity['native_run']
    report=dict(seed=seed,status=result.status,seal_error=result.seal_error,
        normal_solver_seconds=end['elapsed_ns']/1e9,input_preparation_seconds=(begin['monotonic_ns']-started)/1e9,
        terminal_sealing_seconds=(returned-end['monotonic_ns'])/1e9,
        raw_equation_residual=native.get('raw_equation_residual'),scientific_digest=result.scientific_digest,
        material_closed_loop_passed=False,case_complete=False)
    if result.status=='CONVERGED' and result.artifact is not None:
        verification=load(output/'run/record-verification.json')
        if native['raw_equation_residual']>1e-12 or not verification['record_integrity_valid'] or verification['accepted_updates']<1:
            raise ValueError('CONVERGED_RUN_EVIDENCE_INVALID')
        save_start=time.monotonic();reloaded=save_spin_artifact(result.artifact,output/'material')
        report.update(material_closed_loop_passed=True,material_digest=reloaded.scientific_digest,
            terminal_saved_reload_seconds=time.monotonic()-save_start,
            accepted_updates=verification['accepted_updates'],scientific_event_digest=verification['scientific_event_digest'])
    report['total_seconds']=(time.monotonic_ns()-started)/1e9
    (output/'report.json').write_text(json.dumps(report));print(json.dumps(report),flush=True)
    if not report['material_closed_loop_passed']:raise SystemExit(1)


if __name__=='__main__':main(Path(sys.argv[1]),int(sys.argv[2]))
