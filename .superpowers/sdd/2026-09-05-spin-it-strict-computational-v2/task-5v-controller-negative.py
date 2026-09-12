"""Mutate actual exported evidence; never manufacture a positive result."""
import copy
import importlib.util
import json
from pathlib import Path
import sys

spec = importlib.util.spec_from_file_location('checker', Path(__file__).with_name('task-5v-controller-check.py'))
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


def corrupt(case, label, change):
    altered = copy.deepcopy(case)
    change(altered)
    try:
        checker.verify(altered, quiet=True)
    except (AssertionError, KeyError, ValueError, IndexError, TypeError):
        print('REJECT', label)
    else:
        raise AssertionError('corruption accepted: '+label)


def main(path):
    cases = [json.loads(line) for line in Path(path).read_text().splitlines()]
    checker.verify_all(cases, quiet=True)
    original = cases[0]
    changes = [
        ('scope', lambda c: c.update(scope='AUTHOR_ORIGINAL')),
        ('false seal', lambda c: c.update(sealed=True)),
        ('false feasibility', lambda c: c.update(feasibility=1)),
        ('legacy promoted', lambda c: c.update(legacy_geometry=1)),
        ('mode identity', lambda c: c.update(mode=4)),
        ('missing call', lambda c: c['calls'].pop(1)),
        ('missing guard', lambda c: c['guards'].pop(1)),
        ('missing stage', lambda c: c['stages'].pop(0)),
        ('missing final recheck', lambda c: c.update(final_rechecks=0)),
        ('checker undercount', lambda c: c.update(checker_calls=c['checker_calls']-1)),
        ('value undercount', lambda c: c.update(geometry_values=0)),
        ('gradient undercount', lambda c: c.update(source_evaluations=0)),
        ('incomplete value count', lambda c: c.update(geometry_values_complete=False)),
        ('call budget forged', lambda c: c.update(maximum_guard_calls=1)),
        ('callback identity', lambda c: c['guards'][0].update(callback=1)),
        ('stage offset', lambda c: c['guards'][0].update(first_stage=1)),
        ('stage order', lambda c: c['stages'][0].update(stage=1)),
        ('first pass rejected', lambda c: c['guards'][0].update(outcome=1)),
        ('final pass incomplete', lambda c: c['guards'][-1].update(completed=False)),
        ('final check not invoked', lambda c: c['guards'][-1].update(invoked=False)),
        ('raw state tamper', lambda c: c['calls'][0]['state'].__setitem__(0, .9)),
        ('wrong state digest', lambda c: c['calls'][0].update(state_hash='0'*64)),
        ('wrong value digest', lambda c: c['guards'][0].update(value_hash='0'*64)),
        ('missing trial', lambda c: c['trace'][0]['objective_trials'].pop(0)),
        ('accept invalid trial', lambda c: c['trace'][0]['objective_trials'][0].update(accepted=True)),
        ('mismatched next state', lambda c: c['trace'][1]['state'].__setitem__(0,.9)),
        ('mismatched final', lambda c: c['final'].__setitem__(0,.9)),
        ('wrong KKT', lambda c: c['trace'][0]['kkt'].update(full=0.)),
        ('false convergence', lambda c: c.update(numerical_status='CONVERGED')),
        ('wrong material moment', lambda c: c['material']['moments'].__setitem__(0, 2.)),
        ('wrong density', lambda c: c['material'].update(density=2.)),
        ('primitive omitted', lambda c: c['material']['primitives'].pop()),
        ('owner changed', lambda c: c['material']['primitives'][0].update(owner=1)),
        ('clearance changed', lambda c: c.update(minimum_clearance=1.)),
        ('nonfinite state', lambda c: c['initial'].__setitem__(0,float('nan'))),
        ('nonfinite moment', lambda c: c['material']['moments'].__setitem__(0,float('inf'))),
    ]
    for label, change in changes:
        corrupt(original, label, change)
    resource = next(c for c in cases if c['id'] == 'guard_cap')
    for label, change in [
        ('resource promoted', lambda c: c.update(final_status=1)),
        ('work after resource latch', lambda c: c['guards'][-1].update(invoked=True)),
        ('gradient after latch', lambda c: c['calls'][-1].update(source_invoked=True)),
        ('resource outcome hidden', lambda c: c.update(source_outcome=0)),
    ]:
        corrupt(resource, label, change)
    restoration = next(c for c in cases if c['id'] == 'restoration')
    for label, change in [
        ('restoration direction changed', lambda c: c['trace'][0]['restoration']['direction'].__setitem__(9, 9.)),
        ('restoration accepted hidden', lambda c: c['trace'][0]['restoration'].update(accepted=False)),
        ('restoration step changed', lambda c: c['trace'][0]['restoration']['trials'][-1].update(step=1.)),
        ('restoration residual changed', lambda c: c['trace'][0]['restoration']['trials'][-1].update(constraint_l1=1.)),
    ]:
        corrupt(restoration, label, change)
    exhaustion = next(c for c in cases if c['id'] == 'exhaustion')
    for label, change in [
        ('exhaustion history changed', lambda c: c['trace'][0].update(history_after=1)),
        ('exhaustion converged', lambda c: c.update(numerical_status='CONVERGED')),
        ('exact distance changed', lambda c: c['stages'][c['guards'][-1]['first_stage']+5].update(distance='1/2')),
        ('outer winding changed', lambda c: c['material']['outer_faces'][0].reverse()),
    ]:
        corrupt(exhaustion, label, change)
    for label, altered in [('missing case', cases[:-1]), ('duplicate identity', cases[:-1]+[cases[0]])]:
        try:
            checker.verify_all(altered, quiet=True)
        except AssertionError:
            print('REJECT', label)
        else:
            raise AssertionError('corruption accepted: '+label)
    print(f'PASS {len(changes)+14} actual-export corruptions rejected; no synthetic positive pass injected')


if __name__ == '__main__':
    main(sys.argv[1])
