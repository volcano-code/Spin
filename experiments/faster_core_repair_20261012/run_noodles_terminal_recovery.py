"""Independent postprocessing of the interrupted, genuinely accepted terminal.

Never calls Newton or constructs a new source quadric. The original resource
job remains failed; this output cannot complete its missing run record.
"""
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
from make_it_stand import spin
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator, _json_value
from make_it_stand.verification.cherchi_run_evidence import load_cherchi_terminal_prefix
from make_it_stand.verification.cherchi_run_journal import RunJournal
from make_it_stand.verification import faster_material as fm
from make_it_stand.verification.spin_problem import SpinArtifactProblem
from make_it_stand.verification.spin_bundle import save_spin_artifact


def load_problem(frozen, expected_digest):
    import trimesh
    from noodles_formal_case import SEAL_WALL_M
    for name, sha in frozen['input_sources'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != sha:
            raise ValueError('FROZEN_INPUT_CHANGED')
    paths = {Path(path).name: Path(path) for path in frozen['input_sources']}
    meshes = {role: trimesh.load(paths[role+'-si-candidate.obj'], force='mesh', process=False)
              for role in ('outer', 'fixed', 'design')}
    spec = json.loads(next(Path(p) for p in frozen['input_sources']
                          if Path(p).parent.name == 'noodles-dense-v5').read_text())['spec']
    problem = SpinArtifactProblem(original=(meshes['outer'].vertices, meshes['outer'].faces),
        fixed_regions=[('fixed-tubes-handle-tip', (meshes['fixed'].vertices, meshes['fixed'].faces), spec['density_kg_m3'])],
        design_regions=[('design', (meshes['design'].vertices, meshes['design'].faces), spec['density_kg_m3'])],
        length_unit='m', mass_unit='kg', minimum_outer_wall=SEAL_WALL_M,
        target=dict(mode='SPINNING_TOP', axis=spec['target_axis'], support_points=[],
                    origin=(np.asarray(spec['target_origin_model'])*spec['model_unit_m']).tolist()))
    if problem.input_digest != expected_digest:
        raise ValueError('SAVED_INPUT_IDENTITY_MISMATCH')
    return problem


def preflight_saved_backend(config, frozen):
    """Reject changed operands/backends before constructing expensive regions."""
    import _cherchi_core_boolean as bridge
    backend = 'cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:' + hashlib.sha256(Path(bridge.__file__).read_bytes()).hexdigest()
    if backend != config['boolean_backend_identity']:
        raise ValueError('RECOVERY_NATIVE_BACKEND_CHANGED')
    if config.get('rational_arithmetic_backend') == 'NATIVE_GMP_SURFACE_V1':
        import _cherchi_rational_moments as arithmetic
        identity = dict(backend='NATIVE_GMP_SURFACE_V1', algorithm=arithmetic.arithmetic_identity,
                        module_sha256=hashlib.sha256(Path(arithmetic.__file__).read_bytes()).hexdigest())
        if identity != config.get('rational_arithmetic_identity'):
            raise ValueError('RECOVERY_NATIVE_BACKEND_CHANGED')
    for name, sha in frozen['input_sources'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != sha:
            raise ValueError('FROZEN_INPUT_CHANGED')


def independent_terminal_numerics(material, terminal):
    from independent_eq14 import formulas, dual
    f, df, _, g, dg, _ = formulas()
    r = np.asarray(material['raw_moments'], dtype=float)
    p = dual(terminal['quadric_coefficients'])
    grad = np.asarray(df(*r)).ravel(); G = np.asarray(dg(*r))
    constraints = np.asarray(g(*r)).ravel()
    multipliers = np.linalg.lstsq(np.column_stack((-G.T, p)), -grad, rcond=None)[0]
    mu = float(multipliers[-1])
    raw = float(np.linalg.norm(np.r_[grad+mu*p-G.T@multipliers[:-1], constraints]))
    report = dict(objective=float(f(*r)), multiplier_mu=mu, raw_equation_residual=raw,
                  maximum_constraint=float(np.max(np.abs(constraints))),
                  objective_error=abs(float(f(*r))-terminal['objective']),
                  mu_error=abs(mu-terminal['multiplier_mu']),
                  source='Independent symbolic formulas and NumPy least squares')
    if (not np.isfinite(list(report[k] for k in report if k != 'source')).all() or
            mu <= 0 or raw > 1e-12 or report['maximum_constraint'] > 1e-6 or
            report['objective_error'] > 1e-10 or report['mu_error'] > 1e-10):
        raise ValueError('INDEPENDENT_TERMINAL_NUMERICS_FAILED:'+json.dumps(report))
    return report


def main(output, source, *, classified_boundary=False, reuse_input=False, owned_terminal=False):
    if classified_boundary:
        raise ValueError('OFFLINE_RECOVERY_REQUIRES_INDEPENDENT_BOUNDARY: '
                         'saved terminals have no live same-call production validator')
    output.mkdir(parents=True, exist_ok=False)
    with RunJournal(output, durability_profile='PER_EVENT_FSYNC_V1', json_value=_json_value) as journal:
        def stage(name, operation):
            start = time.monotonic_ns()
            journal.emit(dict(kind='STAGE_STARTED', stage=name, monotonic_ns=start))
            print('STAGE_STARTED', name, flush=True)
            try:
                value = operation()
            except BaseException as error:
                journal.emit(dict(kind='STAGE_FAILED', stage=name, error=f'{type(error).__name__}: {error}',
                                  elapsed_ns=time.monotonic_ns()-start))
                raise
            journal.emit(dict(kind='STAGE_RETURNED', stage=name, elapsed_ns=time.monotonic_ns()-start))
            print('STAGE_RETURNED', name, (time.monotonic_ns()-start)/1e9, flush=True)
            return value
        prefix = stage('load_verified_accepted_prefix', lambda: load_cherchi_terminal_prefix(source/'output/run'))
        registration = json.loads((source/'output/registration.json').read_text())
        config = prefix['configuration']; old_request = config['request']
        frozen_path = HERE/'noodles-balanced-request-freeze-job-01/output/request.json'
        frozen = json.loads(frozen_path.read_text())
        if hashlib.sha256(frozen_path.read_bytes()).hexdigest() != registration['request_sha256']:
            raise ValueError('FROZEN_REQUEST_IDENTITY_MISMATCH')
        if old_request['center_of_mass_weight'] != 1. or old_request['spinning_top_objective'] != 'INERTIA_AND_LOWERING':
            raise ValueError('INDEPENDENT_FORMULA_SCOPE_MISMATCH')
        stage('preflight_saved_backend_and_inputs', lambda: preflight_saved_backend(config, frozen))
        journal.publish('registration.json', dict(schema='cherchi-independent-terminal-recovery-v1',
            source_job=str(source), source_job_config_sha256=hashlib.sha256((source/'config.json').read_bytes()).hexdigest(),
            accepted_material=prefix['terminal_reference'], input_digest=registration['input_digest'],
            event_log_sha256=prefix['event_log_sha256'], configuration_sha256=prefix['configuration_sha256'],
            optimization_invocations=0, original_run_complete=False, case_complete=False,
            classification='RECONSTRUCTED_METHOD_VALIDATION',
            validation_profile='INDEPENDENT_RATIONAL_BOUNDARY_AND_EXACT_SEAL_V1'))
        journal.publish('accepted-terminal-event.json', prefix['terminal_event'])
        problem = stage('locked_problem_from_frozen_input', lambda: load_problem(frozen, registration['input_digest']))
        evaluator_settings = spin.SpinEvaluatorConfig(**config['evaluator'])
        solver_settings = spin.SpinSolverConfig(**old_request['settings'])
        prepared = stage('prepare_locked_request_once', lambda:
                         fm._prepare_faster_problem(problem, evaluator_settings, solver_settings,
                             input_validation_profile=('CONSTRUCTOR_PROOF_REUSE_V1' if reuse_input else 'RECHECK_INPUT_V1')))
        problem = prepared.problem
        request = replace(prepared.request, center_of_mass_weight=old_request['center_of_mass_weight'],
                          spinning_top_objective=old_request['spinning_top_objective'])
        actual = json.loads(json.dumps(asdict(request)))
        expected = dict(old_request)
        fixed_error = float(np.max(np.abs(np.asarray(actual.pop('fixed_shell_moments'))-
                                         np.asarray(expected.pop('fixed_shell_moments')))))
        if actual != expected or fixed_error > 1e-10:
            raise ValueError('RECOVERY_REQUEST_CHANGED')
        journal.publish('prepared-request.json', asdict(request))
        producer = InProcessCherchiEvaluator(request, evaluator_settings,
            normalization=problem.normalization, design_vertices_physical=problem.design_regions[0].mesh.arrays()[0],
            rational_arithmetic_backend=config.get('rational_arithmetic_backend', 'PYTHON_FRACTION_V1'),
            stage_observer=stage, fixed_boundary_reuse=classified_boundary)
        if producer.backend_identity != config['boolean_backend_identity'] or producer.arithmetic_identity != config.get('rational_arithmetic_identity'):
            raise ValueError('RECOVERY_NATIVE_BACKEND_CHANGED')
        fixed = stage('prepare_independent_fixed_material', lambda:
                      fm._prepare_external_fixed_material(problem, request, stage_observer=stage))
        retained = prefix['material']; coefficients = retained['quadric_coefficients']
        retained_map = stage('accepted_material_independent_audit', lambda:
            fm._external_geometry_moment_map(problem, request, retained, expected_coefficients=coefficients,
                solver_boundary_profile=('EXACT_RATIONAL_FIXED_CLASSIFIED_V1' if classified_boundary else
                                         'EXACT_RATIONAL_BOUNDARY_V1'), prepared_fixed_material=fixed,
                stage_observer=stage))
        numeric = independent_terminal_numerics(dict(retained, raw_moments=retained_map['raw_moments']),
                                               prefix['terminal_event'])
        journal.publish('independent-terminal-numerics.json', numeric)
        carved = stage('recover_complement_from_saved_source', lambda:
            producer.complement_of_saved_source(retained,
                packet_sink=lambda packet: journal.publish('terminal-complement-raw.json', packet)))
        journal.publish('terminal-complement.json', carved)
        carved_map = stage('complement_independent_audit', lambda:
            fm._external_geometry_moment_map(problem, request, carved, require_derivative=False,
                expected_coefficients=coefficients,
                solver_boundary_profile=('EXACT_RATIONAL_FIXED_CLASSIFIED_V1' if classified_boundary else
                                         'EXACT_RATIONAL_BOUNDARY_V1'),
                prepared_fixed_material=fixed, stage_observer=stage))
        journal.publish('independent-material-moments.json', dict(
            retained=retained_map['raw_moments'], complement=carved_map['variable_raw_moments'],
            fixed_preparation_identity=fixed.identity, fixed_roundoff=fixed_error))
        def model_mesh(raw):
            exact = fm._target_rationals_to_normalized_model(raw['variable_rational_vertices'], request)
            vertices = np.asarray([[float(fm.fraction_from_rational_text(x)) for x in row] for row in exact])
            return vertices, np.asarray(raw['variable_faces'], dtype=np.int32), exact
        artifact = stage('seal_actual_terminal_partition', lambda:
            fm._seal(problem, request, model_mesh(retained), model_mesh(carved),
                (request.outer_vertices, request.outer_faces), retained_map['raw_moments'], stage_observer=stage,
                reuse_problem_validation=reuse_input,
                terminal_region_profile=('OWNED_INTERIOR_DIFFERENCE_V1' if owned_terminal else
                                          'INDEXED_BOUNDARY_REBUILD_V1')))
        journal.publish('terminal-seal-facts.json', artifact.facts)
        reloaded = stage('save_and_independently_reload_material', lambda:
                         save_spin_artifact(artifact, output/'material'))
        if producer.call_counts['source'] != 0:
            raise ValueError('RECOVERY_RECONSTRUCTED_SOURCE')
        journal.publish('report.json', dict(status='RECOVERED_MATERIAL_VALIDATED',
            material_closed_loop_passed=True, material_digest=reloaded.scientific_digest,
            independent_terminal_numerics=numeric, accepted_updates=prefix['accepted_updates'],
            optimization_invocations=0, source_constructions=0, original_run_complete=False,
            case_complete=False, reuse_input_validation=reuse_input,
            classified_boundary=classified_boundary, owned_terminal=owned_terminal,
            remaining=['FRESH_FULL_RUN_AND_REPEAT', 'QUALIFIED_OBJ', 'DYNAMICS']))
        journal.emit(dict(kind='POSTPROCESSING_TERMINATED', status='RECOVERED_MATERIAL_VALIDATED'))


if __name__ == '__main__':
    main(Path(sys.argv[1]), Path(sys.argv[2]),
         classified_boundary='--classified-boundary' in sys.argv[3:],
         reuse_input='--reuse-input-proof' in sys.argv[3:],
         owned_terminal='--owned-terminal' in sys.argv[3:])
