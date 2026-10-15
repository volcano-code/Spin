"""Recompute four actual Heart states using the common production evaluator.

No Newton or answer initialization. Compare source, Boolean material and
Jacobian with original checkpoints, and independently integrate the material.
"""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

from make_it_stand import spin
from make_it_stand.verification.cherchi_inprocess import InProcessCherchiEvaluator
from make_it_stand.verification.spin_problem import rational_solver_candidate_snapshot
from make_it_stand.verification.faster_material import _raw_moments

HERE = Path(__file__).resolve().parent


def main(output):
    output.mkdir(parents=True, exist_ok=False)
    load = lambda p: json.loads(p.read_text())
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    run = HERE/'heart-guarded-formal-seed-0-job-01/output/run'
    diagnostic = load(HERE/'heart-stalled-prefix-independent-algebra-02.json')
    if sha(run/'configuration.json') != diagnostic['configuration_sha256']:
        raise ValueError('CONFIGURATION_CHANGED')
    config = load(run/'configuration.json')
    data = dict(config['request'])
    data['settings'] = spin.SpinSolverConfig(**data['settings'])
    request = spin.SpinItFasterRequest(**data)
    settings = spin.SpinEvaluatorConfig(**config['evaluator'])
    frozen_path = HERE.parent/'faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json'
    frozen = load(frozen_path)
    paths = {Path(path).name: Path(path) for path in frozen['input_sources']}
    outer, design = load(paths['outer.binary64.json']), load(paths['design.binary64.json'])
    vertices = np.asarray(outer['vertices'])
    low, high = vertices.min(0), vertices.max(0)
    spans = high-low
    normalization = dict(origin=(low+spans*.5).tolist(), L=float(max(spans)))
    producer = InProcessCherchiEvaluator(request, settings, normalization=normalization,
        design_vertices_physical=design['vertices'], fixed_boundary_reuse=True,
        rational_arithmetic_backend='NATIVE_GMP_SURFACE_V1')
    if producer.fixed_boundary_proof.reference_identity != config['fixed_boundary_proof_identity']:
        raise ValueError('FROZEN_FRAME_OR_INPUT_CHANGED')
    fields = ('quadric_coefficients', 'variable_rational_vertices', 'variable_faces',
              'variable_vertices', 'variable_vertices_model', 'variable_raw_moments',
              'moment_jacobian', 'quadric_source_vertices', 'quadric_source_faces',
              'quadric_source_identity', 'geometry_hash', 'topology_hash')
    results = []
    for state in diagnostic['diagnostic_states']:
        reference = state['checkpoint']
        path = run/reference['evaluation_file']
        if sha(path) != reference['sha256']:
            raise ValueError('ORIGINAL_CHECKPOINT_CHANGED')
        old = load(path)
        started = time.monotonic()
        producer(state['coefficients'])
        actual = producer.material(state['coefficients'])
        evaluation_seconds = time.monotonic()-started
        checks = {key: json.dumps(actual[key], sort_keys=True) == json.dumps(old[key], sort_keys=True)
                  for key in fields}
        if not all(checks.values()):
            raise ValueError('PRODUCTION_MAPPING_CHANGED:'+json.dumps(checks))
        started = time.monotonic()
        snapshot, mode = rational_solver_candidate_snapshot(actual['variable_rational_vertices'],
            actual['variable_faces'], fixed_boundary_proof=producer.fixed_boundary_proof,
            allow_separated_roots=True)
        snapshot_seconds = time.monotonic()-started
        started = time.monotonic()
        independent = _raw_moments(snapshot, request.material_density)
        independent_seconds = time.monotonic()-started
        error = float(np.max(np.abs(np.asarray(independent)-actual['variable_raw_moments'])))
        if error > 1e-10:
            raise ValueError('INDEPENDENT_MOMENT_MISMATCH')
        name = f"evaluation-{state['evaluation']:06d}.json"
        payload = json.dumps(actual, sort_keys=True, allow_nan=False).encode()
        (output/name).write_bytes(payload)
        row = dict(evaluation=state['evaluation'], original_sha256=reference['sha256'],
            replay_sha256=hashlib.sha256(payload).hexdigest(), checks=checks,
            evaluation_seconds=evaluation_seconds, snapshot_seconds=snapshot_seconds,
            independent_seconds=independent_seconds, independent_moment_error=error,
            snapshot_mode=mode, separation_proof=getattr(snapshot, 'separation_proof', None))
        results.append(row)
        (output/f"comparison-{state['evaluation']:06d}.json").write_text(json.dumps(row, allow_nan=False))
        print(json.dumps(row), flush=True)
    report = dict(schema='heart-four-production-state-replay-v1', rows=results,
        source_module_sha256=sha(Path(spin._native_module().__file__)),
        registered_evaluator=asdict(settings), optimization_invocations=0,
        derivative_finite_difference_gate=False, case_complete=False, passed=True,
        scope='Actual source/Boolean/Jacobian equality at four registered points; not neighborhood derivative certification')
    (output/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
