"""Frozen Heart seed through the common Cherchi/Newton/material chain."""
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent


def main(output, seed, *, certified=False):
    from make_it_stand import spin
    from make_it_stand.verification.spin_problem import SpinArtifactProblem
    from make_it_stand.verification.spin_bundle import save_spin_artifact
    from make_it_stand.verification.cherchi_inprocess import solve_cherchi_inprocess_and_seal_faster
    from case_coordinate_contract import normalized_faster_weight
    from run_noodles_formal_case import registered_solver
    load = lambda path: json.loads(path.read_text())
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    admission_path = HERE/('heart-certified-solve-admission-01/report.json' if certified else
                          'heart-guarded-solve-admission-01/report.json')
    admission = load(admission_path)
    if (not admission['guarded_solve_admission'] or not admission['runtime_branch_regression_passed'] or
            admission['old_five_point_gate_passed'] or admission['neighborhood_certification']):
        raise ValueError('GUARDED_HEART_ADMISSION_NOT_VALID')
    for record in admission['bindings'].values():
        if sha(Path(record['path'])) != record['sha256']:
            raise ValueError('GUARDED_HEART_ADMISSION_DEPENDENCY_CHANGED:'+record['path'])
    if sha(Path(spin._native_module().__file__)) != admission['bindings']['native']['sha256']:
        raise ValueError('LOADED_HEART_NATIVE_MISMATCH')
    frozen_path = HERE.parent/'faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json'
    frozen = load(frozen_path)
    if sha(frozen_path) != admission['request_sha256']:
        raise ValueError('FROZEN_HEART_REQUEST_CHANGED')
    if frozen['units'] != dict(model_length_m=.025, model_mass_kg=.01875,
            density_kg_m3=1200., objective_coordinate_contract='INPUT_MODEL_UNITS_TARGET_ORIGIN_ZERO_V1'):
        raise ValueError('FROZEN_HEART_UNITS_CHANGED')
    old = frozen['requests'][seed]
    solver = spin.SpinSolverConfig(**registered_solver(old['settings'], seed))
    paths = {Path(path).name: Path(path) for path in frozen['input_sources']}
    meshes = {role: load(paths[role+'.binary64.json']) for role in ('outer', 'fixed', 'design')}
    pair = lambda role: (meshes[role]['vertices'], meshes[role]['faces'])
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic_ns()
    problem = SpinArtifactProblem(original=pair('outer'),
        fixed_regions=[('fixed-shell-handle-pin', pair('fixed'), 1.)],
        design_regions=[('design', pair('design'), 1.)],
        length_unit='MODEL_UNIT_0.025_M', mass_unit='MODEL_UNIT_0.01875_KG', minimum_outer_wall=.015,
        target=dict(mode='SPINNING_TOP', axis=[0., 0., 1.], origin=[0., 0., 0.], support_points=[]))
    if problem.input_digest != admission['input_digest']:
        raise ValueError('FROZEN_HEART_INPUT_DIGEST_CHANGED')
    evaluator = spin.SpinEvaluatorConfig(**load(HERE/'heart-cherchi-coordinate-map-job-01/output/registration.json')['settings'])
    weight = normalized_faster_weight(old['center_of_mass_weight'], problem.normalization['L'])
    (output/'registration.json').write_text(json.dumps(dict(seed=seed, input_digest=problem.input_digest,
        request_sha256=sha(frozen_path), guarded_admission_sha256=sha(admission_path),
        classification='RECONSTRUCTED_METHOD_VALIDATION', backend='CHERCHI',
        initialization='REQUEST_SEED_ONLY_NO_EXPLICIT_QUADRICS', original_five_point_gate_passed=False,
        neighborhood_certification=False, center_of_mass_weight=weight,
        si_conversion=dict(length_m_per_model_unit=.025, mass_kg_per_model_unit=.01875,
                           inertia_kg_m2_per_model_unit=.01875*.025**2))))
    result = solve_cherchi_inprocess_and_seal_faster(problem, evaluator_settings=evaluator, solver_settings=solver,
        output_directory=output/'run', center_of_mass_weight=weight,
        spinning_top_objective=old['spinning_top_objective'], rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1',
        solver_boundary_profile=('EXACT_RATIONAL_FIXED_CLASSIFIED_V1' if certified else 'EXACT_RATIONAL_FIXED_SUBSET_V1'),
        durability_profile='EVALUATION_CHECKPOINT_V1', validation_profile='LAYERED_SAME_CALL_V1',
        input_validation_profile=('CONSTRUCTOR_PROOF_REUSE_V1' if certified else 'RECHECK_INPUT_V1'),
        terminal_region_profile=('OWNED_INTERIOR_DIFFERENCE_V1' if certified else 'INDEXED_BOUNDARY_REBUILD_V1'),
        integral_profile=('CERTIFIED_MPFR256_GMP_FALLBACK_V1' if certified else 'EXACT_GMP_V1'))
    returned = time.monotonic_ns()
    events = [json.loads(line) for line in (output/'run/events.jsonl').read_text().splitlines()]
    start = next(row for row in events if row.get('kind') == 'STAGE_STARTED' and row.get('stage') == 'native_newton_solver')
    end = next(row for row in events if row.get('kind') == 'STAGE_RETURNED' and row.get('stage') == 'native_newton_solver')
    native = result.identity['native_run']
    report = dict(seed=seed, status=result.status, seal_error=result.seal_error,
        normal_solver_seconds=end['elapsed_ns']/1e9,
        input_preparation_seconds=(start['monotonic_ns']-started)/1e9,
        terminal_sealing_seconds=(returned-end['monotonic_ns'])/1e9,
        raw_equation_residual=native.get('raw_equation_residual'), multiplier_mu=native.get('multiplier_mu'),
        material_closed_loop_passed=False, case_complete=False,
        old_five_point_gate_passed=False, neighborhood_certification=False)
    if result.status == 'CONVERGED' and result.artifact is not None:
        verified = load(output/'run/record-verification.json')
        if native['raw_equation_residual'] > 1e-12 or native['multiplier_mu'] <= 0 or not verified['record_integrity_valid']:
            raise ValueError('HEART_TERMINAL_NOT_QUALIFIED')
        save_start = time.monotonic()
        reloaded = save_spin_artifact(result.artifact, output/'material')
        report.update(material_closed_loop_passed=True, material_digest=reloaded.scientific_digest,
            terminal_saved_reload_seconds=time.monotonic()-save_start, accepted_updates=verified['accepted_updates'],
            scientific_event_digest=verified['scientific_event_digest'])
    report['total_seconds'] = (time.monotonic_ns()-started)/1e9
    (output/'report.json').write_text(json.dumps(report, sort_keys=True, allow_nan=False))
    print(json.dumps(report), flush=True)
    if not report['material_closed_loop_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main(Path(sys.argv[1]), int(sys.argv[2]),certified='--certified-integral' in sys.argv)
