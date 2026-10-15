"""One frozen scale and seed through the formal Cherchi material chain."""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
import time

from make_it_stand import spin
from make_it_stand.verification.spin_problem import SpinArtifactProblem
from make_it_stand.verification.spin_bundle import save_spin_artifact
from make_it_stand.verification.cherchi_inprocess import solve_cherchi_inprocess_and_seal_faster

HERE = Path(__file__).resolve().parent


def main(output, scale, seed, arithmetic_backend='PYTHON_FRACTION_V1',
         boundary_profile='EXACT_NEF_REGION_V1',durability_profile='PER_EVENT_FSYNC_V1',
         validation_profile='INDEPENDENT_EACH_EVALUATION_V1', input_validation_profile='RECHECK_INPUT_V1',
         terminal_region_profile='INDEXED_BOUNDARY_REBUILD_V1'):
    if scale not in ('small', 'faces-12288', 'faces-98304') or seed not in range(10):
        raise ValueError('UNREGISTERED_SCALE_OR_SEED')
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic_ns()
    inputs = HERE / 'core-benchmark-inputs-job-01/output'
    manifest = json.loads((inputs / 'manifest.json').read_text())
    row = next(x for x in manifest['rows'] if x['name'] == scale)
    raw = (inputs / row['configuration_file']).read_bytes()
    if hashlib.sha256(raw).hexdigest() != row['sha256']:
        raise ValueError('INPUT_CONFIGURATION_CHANGED')
    config = json.loads(raw)
    r = config['request']
    r['settings']['deterministic_seed'] = seed
    fixed = config['fixed_material']
    problem = SpinArtifactProblem(original=(r['outer_vertices'], r['outer_faces']),
        fixed_regions=[('shell', (fixed['vertices'], fixed['faces']), fixed['density'])],
        design_regions=[('design', (r['design_vertices'], r['design_faces']), r['material_density'])],
        length_unit='m', mass_unit='kg', minimum_outer_wall=0.05,
        target=dict(mode='SPINNING_TOP', axis=[0., 0., 1.], origin=r['target_origin'], support_points=[]))
    evaluator = spin.SpinEvaluatorConfig(**config['evaluator'])
    solver = spin.SpinSolverConfig(**r['settings'])
    if (not evaluator.paper_projective_hyperboloid or solver.maximum_restarts != 8
            or solver.maximum_function_evaluations_per_restart != 60 or solver.equation_tolerance != 1e-12):
        raise ValueError('BENCHMARK_PROTOCOL_MISMATCH')
    (output / 'registration.json').write_text(json.dumps(dict(
        scale=scale, seed=seed, input_sha256=row['sha256'], input_digest=problem.input_digest,
        evaluator=asdict(evaluator), solver=asdict(solver), initialization='REQUEST_SEED_ONLY',
        rational_arithmetic_backend=arithmetic_backend,
        solver_boundary_profile=boundary_profile,
        durability_profile=durability_profile,
        validation_profile=validation_profile,
        input_validation_profile=input_validation_profile,
        terminal_region_profile=terminal_region_profile,
        native_module_sha256=hashlib.sha256(Path(spin._native_module().__file__).read_bytes()).hexdigest(),
        core_admission=False), sort_keys=True))
    result = solve_cherchi_inprocess_and_seal_faster(problem, evaluator_settings=evaluator,
        solver_settings=solver, output_directory=output / 'run', center_of_mass_weight=1.,
        rational_arithmetic_backend=arithmetic_backend,solver_boundary_profile=boundary_profile,
        durability_profile=durability_profile,validation_profile=validation_profile,
        input_validation_profile=input_validation_profile, terminal_region_profile=terminal_region_profile)
    solved_at = time.monotonic_ns()
    events = [json.loads(line) for line in (output / 'run/events.jsonl').read_text().splitlines()]
    native_start = next(e for e in events if e.get('kind') == 'STAGE_STARTED' and e.get('stage') == 'native_newton_solver')
    native_finish = next(e for e in events if e.get('kind') == 'STAGE_RETURNED' and e.get('stage') == 'native_newton_solver')
    native = result.identity['native_run']
    report = dict(scale=scale, seed=seed, status=result.status, seal_error=result.seal_error,
        rational_arithmetic_backend=arithmetic_backend,
        solver_boundary_profile=boundary_profile,
        durability_profile=durability_profile,
        validation_profile=validation_profile,
        input_validation_profile=input_validation_profile,
        terminal_region_profile=terminal_region_profile,
        normal_solver_seconds=native_finish['elapsed_ns'] / 1e9,
        input_preparation_seconds=(native_start['monotonic_ns'] - started) / 1e9,
        terminal_sealing_seconds=(solved_at - native_finish['monotonic_ns']) / 1e9,
        raw_equation_residual=native.get('raw_equation_residual'),
        scientific_digest=result.scientific_digest, core_admission=False)
    if result.status == 'CONVERGED' and result.artifact is not None:
        if native['raw_equation_residual'] > 1e-12:
            raise ValueError('RAW_RESIDUAL_GATE_FAILED')
        verification = json.loads((output / 'run/record-verification.json').read_text())
        if not verification['record_integrity_valid'] or verification['accepted_updates'] < 1:
            raise ValueError('NO_VALID_ACCEPTED_UPDATES')
        save_start = time.monotonic()
        reloaded = save_spin_artifact(result.artifact, output / 'material')
        report.update(terminal_saved_reload_seconds=time.monotonic() - save_start,
                      material_digest=reloaded.scientific_digest,
                      accepted_updates=verification['accepted_updates'],
                      scientific_event_digest=verification['scientific_event_digest'],
                      material_closed_loop_passed=True)
    else:
        report['material_closed_loop_passed'] = False
    report['total_seconds'] = (time.monotonic_ns() - started) / 1e9
    (output / 'report.json').write_text(json.dumps(report, sort_keys=True, allow_nan=False))
    if not report['material_closed_loop_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main(Path(sys.argv[1]), sys.argv[2], int(sys.argv[3]),
         sys.argv[4] if len(sys.argv)>4 else 'PYTHON_FRACTION_V1',
         sys.argv[5] if len(sys.argv)>5 else 'EXACT_NEF_REGION_V1',
         sys.argv[6] if len(sys.argv)>6 else 'PER_EVENT_FSYNC_V1',
         sys.argv[7] if len(sys.argv)>7 else 'INDEPENDENT_EACH_EVALUATION_V1',
         sys.argv[8] if len(sys.argv)>8 else 'RECHECK_INPUT_V1',
         sys.argv[9] if len(sys.argv)>9 else 'INDEXED_BOUNDARY_REBUILD_V1')
