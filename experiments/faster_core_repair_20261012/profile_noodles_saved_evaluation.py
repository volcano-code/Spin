"""Actual accepted noodles map, measured without rebuilding an input Nef."""
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import trimesh
from make_it_stand import spin
from make_it_stand.verification import cherchi_inprocess as adapter
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator, _json_value
from make_it_stand.verification.cherchi_run_evidence import load_cherchi_terminal_prefix
from make_it_stand.verification.cherchi_run_journal import RunJournal
from make_it_stand.verification.spin_problem import rational_solver_candidate_snapshot
from independent_tetra_moments import integrate
from fractions import Fraction

HERE = Path(__file__).resolve().parent


def main(output):
    run = HERE/'noodles-formal-cherchi-seed-0-job-01/output/run'
    frozen_path = HERE/'noodles-balanced-request-freeze-job-01/output/request.json'
    frozen = json.loads(frozen_path.read_text())
    for name, sha in frozen['input_sources'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != sha:
            raise ValueError('FROZEN_INPUT_CHANGED')
    paths = {Path(name).name: Path(name) for name in frozen['input_sources']}
    output.mkdir(parents=True, exist_ok=False)
    with RunJournal(output, durability_profile='EVALUATION_CHECKPOINT_V1', json_value=_json_value) as journal:
        timings = []
        def stage(name, operation):
            start = time.monotonic_ns()
            journal.emit(dict(kind='STAGE_STARTED', stage=name, monotonic_ns=start))
            try:
                result = operation()
            except Exception as error:
                journal.emit(dict(kind='STAGE_FAILED', stage=name,
                    elapsed_ns=time.monotonic_ns()-start, error=f'{type(error).__name__}: {error}'))
                raise
            elapsed = time.monotonic_ns()-start
            timings.append(dict(stage=name, elapsed_ns=elapsed))
            journal.emit(dict(kind='STAGE_RETURNED', **timings[-1]))
            print(name, elapsed/1e9, flush=True)
            return result
        prefix = stage('verify_original_accepted_prefix', lambda: load_cherchi_terminal_prefix(run))
        config, old = prefix['configuration'], prefix['material']
        data = dict(config['request']); data['settings'] = spin.SpinSolverConfig(**data['settings'])
        request = spin.SpinItFasterRequest(**data)
        outer = trimesh.load(paths['outer-si-candidate.obj'], force='mesh', process=False)
        design = trimesh.load(paths['design-si-candidate.obj'], force='mesh', process=False)
        low, high = outer.vertices.min(0), outer.vertices.max(0)
        normalization = dict(origin=(low+(high-low)*.5).tolist(), L=float(max(high-low)))
        producer = stage('prepare_same_input_producer', lambda: InProcessCherchiEvaluator(request,
            spin.SpinEvaluatorConfig(**config['evaluator']), normalization=normalization,
            design_vertices_physical=design.vertices, rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1',
            fixed_boundary_reuse=True, stage_observer=stage))
        if producer.fixed_boundary_proof.reference_identity != config['fixed_boundary_proof_identity']:
            raise ValueError('ORIGINAL_FIXED_BOUNDARY_CHANGED')
        q = old['quadric_coefficients']
        original_orientation = adapter.normalize_closed_component_orientation
        original_target_view = producer._target_view
        def timed_orientation(*args, **kwargs):
            return stage('production_orientation', lambda: original_orientation(*args, **kwargs))
        def timed_target_view(*args, **kwargs):
            return stage('exact_coordinate_conversion', lambda: original_target_view(*args, **kwargs))
        adapter.normalize_closed_component_orientation = timed_orientation
        producer._target_view = timed_target_view
        try:
            mapped = stage('production_geometry_evaluation', lambda: producer(q))
        finally:
            adapter.normalize_closed_component_orientation = original_orientation
            producer._target_view = original_target_view
        actual = producer.material(q)
        fields = ('variable_rational_vertices', 'variable_faces', 'moment_jacobian',
                  'quadric_source_identity', 'geometry_hash', 'topology_hash', 'material_proof')
        checks = {key: actual[key] == old[key] for key in fields}
        if not all(checks.values()):
            raise ValueError('ORIGINAL_PRODUCTION_MAPPING_CHANGED:'+json.dumps(checks))
        stage('evaluation_checkpoint', lambda: journal.publish('accepted-material.json', actual))
        saved = json.loads((output/'accepted-material.json').read_text())
        stage('independent_exact_boundary', lambda: rational_solver_candidate_snapshot(
            saved['variable_rational_vertices'], saved['variable_faces'],
            fixed_boundary_proof=producer.fixed_boundary_proof, allow_signed_classification=True))
        independent = stage('independent_saved_integral', lambda:
            integrate([[Fraction(x) for x in row] for row in saved['variable_rational_vertices']],
                saved['variable_faces'], request.material_density))
        error = float(np.max(np.abs(independent-np.asarray(mapped['variable_raw_moments']))))
        if error > 1e-10:
            raise ValueError('INDEPENDENT_MATERIAL_INTEGRAL_CHANGED')
        journal.publish('report.json', dict(schema='noodles-accepted-map-stage-profile-v1',
            passed=True, checks=checks, independent_error=error, timings=timings,
            arithmetic_identity=producer.arithmetic_identity,
            boolean_backend_identity=producer.backend_identity,
            evaluation_map_identity=producer.context_identity,
            optimization_invocations=0, case_complete=False, performance_gate_passed=False,
            scope='One actual accepted state; excludes cold problem construction and terminal seal',
            original_material_reference=prefix['terminal_reference']))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
