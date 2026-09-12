"""Mutation audit of actual exported solver observations, not hand-written pass flags."""
import contextlib
import copy
import importlib.util
import io
import json
from pathlib import Path
import sys

spec = importlib.util.spec_from_file_location('trace_checker', Path(__file__).with_name('task-5t-controller-check.py'))
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)
cases = [json.loads(line) for line in Path(sys.argv[1]).read_text().splitlines()]
assert len(cases) == 10
mutations = [
    lambda c: c.update(geometry_status='VERIFIED'),
    lambda c: c.update(feasibility_status='FEASIBLE_WITNESS'),
    lambda c: c.update(algorithm='LEGACY_V3'),
    lambda c: c.update(handle_count=25),
    lambda c: c.update(final_phi=c['final_phi']+.1),
    lambda c: c['final'].__setitem__(0, c['final'][0]+.1),
    lambda c: c['trace'][-1].update(final_state=False),
    lambda c: c['trace'][0]['gradient'].__setitem__(0, c['trace'][0]['gradient'][0]+100),
    lambda c: c['trace'][0]['constraints'].__setitem__(0, 100),
    lambda c: c['trace'][0]['kkt'].update(full=100),
    lambda c: c['trace'][0]['kkt']['alpha'].__setitem__(8, 1),
    lambda c: c['trace'][0].update(phi=c['trace'][0]['phi']+.1),
    lambda c: c['trace'][0].update(membership_hash='0'*64),
    lambda c: c['trace'][0].update(history_before=7),
    lambda c: c['trace'][0].update(merit=float('nan')),
    lambda c: c['trace'][0].update(line_search_evaluations=-1),
    lambda c: c.update(trace=[]),
    lambda c: c['controls'].update(equality_tolerance=1),
    lambda c: c.update(numerical_only=False),
    lambda c: c['calls'][0].update(state_hash='0'*64),
    lambda c: c.update(source_evaluation_count=0),
    lambda c: c['calls'][-1].update(purpose=0),
    lambda c: c['trace'][0]['objective_trials'].append(
        {'step': .5, 'valid': False, 'accepted': True, 'reason': 'forged', 'constraint_l1': 0}),
]
rejected = 0
for case in (cases[0], cases[-1]):
    for index, mutation in enumerate(mutations):
        broken = copy.deepcopy(case)
        mutation(broken)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                checker.verify(broken)
        except (AssertionError, ValueError):
            rejected += 1
        else:
            raise AssertionError(f'corrupted observation admitted: {index}')
print(f'PASS: {rejected}/46 actual trace/call/status/physics-promotion corruptions rejected.')
