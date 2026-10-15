"""Replay saved *boundary checks*; never Boolean, Newton or case admission.

Each child process measures one frozen state cold then three times warm.
Reported native phases are exclusive; Python envelope is separate.
"""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from boundary_profile_inputs import select_saved_states
from make_it_stand import spin
from make_it_stand.verification.spin_problem import _canonical_rational_vertices

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reference_coordinates(output):
    config = json.loads((output/'run/configuration.json').read_text())
    scale = json.loads((output/'report.json').read_text())['scale']
    original = json.loads((HERE/'core-benchmark-inputs-job-01/output'/f'{scale}.json').read_text())
    normalization = json.loads((output/'material/problem.json').read_text())['identity']['normalization']
    origin = list(map(Fraction, normalization['origin']))
    length = Fraction(normalization['L'])
    target = list(map(Fraction, config['request']['target_origin']))
    axes = [[Fraction(x) for x in row] for row in config['request']['target_axes']]
    rows = []
    for vertex in original['request']['design_vertices']:
        model = [(Fraction(vertex[i])-origin[i])/length-target[i] for i in range(3)]
        point = [sum(model[i]*axes[i][j] for i in range(3)) for j in range(3)]
        rows.append([f'{x.numerator}/{x.denominator}' for x in point])
    return rows, np.asarray(original['request']['design_faces'], dtype=np.int32), config['fixed_boundary_proof_identity']


def state_run(output, selected, destination, measured=True):
    n = spin._native_module()
    reference, faces, expected = reference_coordinates(output)
    started = time.monotonic()
    proof = n.verify_fixed_rational_boundary(reference, faces)
    preparation = time.monotonic()-started
    if proof.reference_identity != expected:
        raise ValueError('FIXED_PROOF_FRAME_MISMATCH')
    path = Path(selected['path'])
    if sha(path) != selected['sha256']:
        raise ValueError('SAVED_STATE_CHANGED')
    packet = json.loads(path.read_text())
    rows = []
    facts = None
    for repetition in range(4 if measured else 1):
        start = time.monotonic()
        exact = _canonical_rational_vertices(packet['variable_rational_vertices'])
        f = np.asarray(packet['variable_faces'], dtype=np.int32)
        canonical_seconds = time.monotonic()-start
        begin = time.monotonic()
        if measured:
            result = n.audit_rational_boundary_with_fixed_proof(exact, f, 'AnyNonzero', proof, profile=True)
            detail = result.pop('performance')
        else:
            result = n.audit_rational_boundary_with_fixed_proof(exact, f, 'AnyNonzero', proof)
            detail = None
        native_envelope = time.monotonic()-begin
        if facts is not None and result != facts:
            raise ValueError('BOUNDARY_FACTS_NOT_DETERMINISTIC')
        facts = result
        begin = time.monotonic()
        display = np.asarray([[float(Fraction(x)) for x in row] for row in exact], dtype='<f8')
        float_seconds = time.monotonic()-begin
        if display.tobytes() != np.asarray(packet['variable_vertices'], dtype='<f8').tobytes():
            raise ValueError('FLOAT_VIEW_CHANGED')
        begin = time.monotonic()
        serialized = json.dumps(exact, sort_keys=True, separators=(',', ':')).encode('ascii')
        identity = hashlib.sha256(serialized).hexdigest()
        serialization_seconds = time.monotonic()-begin
        rows.append(dict(repetition=repetition, temperature='cold' if repetition == 0 else 'warm',
            canonicalization_seconds=canonical_seconds, native_envelope_seconds=native_envelope,
            float_view_seconds=float_seconds, serialization_seconds=serialization_seconds,
            rational_identity=identity, native=detail))
    destination.write_text(json.dumps(dict(selected=selected, facts=facts, rows=rows,
        fixed_preparation_seconds=preparation, native_path=n.__file__, native_sha256=sha(Path(n.__file__)),
        core_admission=False, scope='saved boundary only; no source/Boolean/derivative/Newton evaluation'), sort_keys=True))


def main(destination, baseline_module, diagonal_candidate=False):
    destination.mkdir(parents=True, exist_ok=False)
    inventory = []
    reports = []
    for scale in ('small', 'faces-12288', 'faces-98304'):
        output = HERE/f'layered-formal-benchmark-{scale}-seed-0-job-05/output'
        for state in select_saved_states(output/'run'):
            inventory.append(dict(scale=scale, **state))
            for mode in ('baseline', 'profile'):
                target = destination/f'{scale}-{state["role"]}-{mode}.json'
                command = [sys.executable, '-B', __file__, '--state', str(output),
                           json.dumps(state), str(target), mode]
                import os
                env = dict(os.environ)
                if mode == 'baseline':
                    env['PYTHONPATH'] = str(Path(baseline_module).parent)+':'+env.get('PYTHONPATH','')
                subprocess.run(command, check=True, env=env)
            baseline = json.loads((destination/f'{scale}-{state["role"]}-baseline.json').read_text())
            profiled = json.loads((destination/f'{scale}-{state["role"]}-profile.json').read_text())
            comparison = dict(baseline_facts_identical=True)
            if diagonal_candidate:
                from boundary_candidate_comparison import compare_boundary_facts
                comparison = compare_boundary_facts(baseline['facts'],profiled['facts'])
            elif baseline['facts'] != profiled['facts']:
                raise ValueError('PROFILE_CHANGED_BASELINE_FACTS')
            reports.append(dict(scale=scale, role=state['role'], **comparison,
                baseline_seconds=baseline['rows'][0]['native_envelope_seconds'],
                measured=profiled['rows']))
            print(scale, state['role'], 'same facts; four measurements saved', flush=True)
    (destination/'registration.json').write_text(json.dumps(inventory, sort_keys=True))
    (destination/'report.json').write_text(json.dumps(dict(rows=reports, passed=True,
        core_admission=False, scope='boundary profiling, not full fixed-state evaluation'), sort_keys=True))


if __name__ == '__main__':
    if sys.argv[1] == '--state':
        state_run(Path(sys.argv[2]), json.loads(sys.argv[3]), Path(sys.argv[4]), sys.argv[5] == 'profile')
    else:
        main(Path(sys.argv[1]), sys.argv[2], len(sys.argv)>3 and sys.argv[3]=='diagonal-candidate')
