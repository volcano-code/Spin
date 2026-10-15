"""Time the actual frozen input constructor and fixed integral; no Newton."""
from dataclasses import asdict, replace
import json
from pathlib import Path
import sys
import time

import numpy as np
from make_it_stand import spin
from make_it_stand.verification import faster_material as fm
from make_it_stand.verification import spin_problem as sp
from make_it_stand.verification.cherchi_inprocess import _json_value
from make_it_stand.verification.cherchi_run_evidence import load_cherchi_terminal_prefix
from make_it_stand.verification.cherchi_run_journal import RunJournal
from run_noodles_terminal_recovery import load_problem

HERE = Path(__file__).resolve().parent


def main(output, source):
    output.mkdir(parents=True, exist_ok=False)
    with RunJournal(output, durability_profile='PER_EVENT_FSYNC_V1', json_value=_json_value) as journal:
        timings = []
        def stage(name, operation, **details):
            start = time.monotonic_ns()
            journal.emit(dict(kind='STAGE_STARTED', stage=name, monotonic_ns=start, **details))
            print('STAGE_STARTED', name, flush=True)
            try:
                value = operation()
            except BaseException as error:
                journal.emit(dict(kind='STAGE_FAILED', stage=name,
                    elapsed_ns=time.monotonic_ns()-start, error=f'{type(error).__name__}: {error}'))
                raise
            elapsed = time.monotonic_ns()-start
            timings.append(dict(stage=name, elapsed_ns=elapsed, **details))
            journal.emit(dict(kind='STAGE_RETURNED', **timings[-1]))
            print('STAGE_RETURNED', name, elapsed/1e9, flush=True)
            return value

        prefix = stage('load_verified_accepted_prefix', lambda:
            load_cherchi_terminal_prefix(source/'output/run'))
        registration = json.loads((source/'output/registration.json').read_text())
        frozen_path = HERE/'noodles-balanced-request-freeze-job-01/output/request.json'
        frozen = json.loads(frozen_path.read_text())
        journal.publish('registration.json', dict(schema='frozen-noodles-input-cost-v1',
            input_digest=registration['input_digest'], optimization_invocations=0,
            case_complete=False, performance_gate_passed=False))

        original_init, original_partition = sp.MeshSnapshot.__init__, sp.partition
        construction_count = 0
        def timed_init(self, mesh, **kwargs):
            nonlocal construction_count
            index = construction_count
            construction_count += 1
            return stage(f'input_mesh_snapshot_{index}', lambda:
                original_init(self, mesh, **kwargs))
        def timed_partition(whole, parts):
            return stage('input_exact_partition', lambda: original_partition(whole, parts))
        sp.MeshSnapshot.__init__, sp.partition = timed_init, timed_partition
        try:
            problem = stage('locked_problem_from_frozen_input', lambda:
                load_problem(frozen, registration['input_digest']))
        finally:
            sp.MeshSnapshot.__init__, sp.partition = original_init, original_partition

        config = prefix['configuration']
        old = config['request']
        prepared = stage('prepare_locked_request_once', lambda:
            fm._prepare_faster_problem(problem, spin.SpinEvaluatorConfig(**config['evaluator']),
                spin.SpinSolverConfig(**old['settings']),
                input_validation_profile='CONSTRUCTOR_PROOF_REUSE_V1'))
        request = replace(prepared.request, center_of_mass_weight=old['center_of_mass_weight'],
            spinning_top_objective=old['spinning_top_objective'])
        actual, expected = json.loads(json.dumps(asdict(request))), dict(old)
        error = float(np.max(np.abs(np.asarray(actual.pop('fixed_shell_moments'))-
            np.asarray(expected.pop('fixed_shell_moments')))))
        if actual != expected or error > 1e-10:
            raise ValueError('FROZEN_REQUEST_CHANGED')
        fixed = stage('prepare_independent_fixed_material', lambda:
            fm._prepare_external_fixed_material(problem, request, stage_observer=stage))
        accepted = prefix['terminal_event']['evaluation_payload']
        previous = np.asarray(accepted['solver_raw_moments'])-np.asarray(accepted['variable_raw_moments'])
        moment_error = float(np.max(np.abs(previous-np.asarray(fixed.moments))))
        if moment_error > 1e-10:
            raise ValueError('ACCEPTED_FIXED_MOMENTS_CHANGED')
        journal.publish('report.json', dict(schema='frozen-noodles-input-cost-v1',
            passed=True, input_digest=problem.input_digest, fixed_moments=fixed.moments,
            request_fixed_error=error, accepted_fixed_error=moment_error,
            input_constructor_calls=construction_count, timings=timings,
            optimization_invocations=0, case_complete=False, performance_gate_passed=False))


if __name__ == '__main__':
    main(Path(sys.argv[1]), Path(sys.argv[2]))
