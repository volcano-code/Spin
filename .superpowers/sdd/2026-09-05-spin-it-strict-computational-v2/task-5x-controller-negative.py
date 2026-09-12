"""Corrupt actual native exports, never manufacture successful evidence."""
import copy
import importlib.util
import json
from pathlib import Path
import sys

spec = importlib.util.spec_from_file_location('checker', Path(__file__).with_name('task-5x-controller-check.py'))
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)
rows = [json.loads(line) for line in Path(sys.argv[1]).read_text().splitlines()]
checker.verify_rows(rows)
cases = {c['id']: c for c in rows}
rejected = 0


def must_reject(label, value, check):
    global rejected
    try:
        check(value)
    except (AssertionError, ValueError, KeyError, TypeError, IndexError):
        rejected += 1
        print('REJECT:', label)
    else:
        raise AssertionError('corrupted evidence accepted: '+label)


def corruption(case, label, edit):
    value = copy.deepcopy(cases[case])
    edit(value)
    must_reject(label, value, checker.verify)


must_reject('missing case', rows[:-1], checker.verify_rows)
must_reject('duplicate ID', rows[:-1]+[rows[0]], checker.verify_rows)
corruption('stationary_multi', 'paper grant', lambda c: c.update(paper_reproduction=True))
corruption('stationary_multi', 'hidden middle pass', lambda c: c['passes'].pop(1))
corruption('stationary_multi', 'reused history', lambda c: c['passes'][1]['guarded']['trace'][0].update(history_before=3))
corruption('stationary_multi', 'false KKT stationarity', lambda c: c['passes'][0]['guarded']['trace'][0]['kkt'].update(stationarity=9.))
corruption('stationary_multi', 'wrong callback hash', lambda c: c['passes'][0]['guarded']['calls'][0].update(state_hash='0'*64))
corruption('stationary_multi', 'wrong intermediate functional', lambda c: c['passes'][1].update(old_functional='1'*64))
corruption('stationary_multi', 'wrong rebound beta', lambda c: c['passes'][0]['rebind']['source']['state'].__setitem__(0,.4))
corruption('stationary_multi', 'wrong raw mass', lambda c: c['passes'][0]['source']['moments'].__setitem__(0,.5))
corruption('stationary_multi', 'forged terminal material', lambda c: c['last_solved']['moments'].__setitem__(0,.5))
corruption('stationary_multi', 'nonfinite trace', lambda c: c['passes'][0]['guarded']['trace'][0].update(objective=float('nan')))
corruption('stationary_pass_limit', 'pending relabeled solved', lambda c: c.update(last_solved=c['pending']))
corruption('nontrivial_stop', 'missing failed callbacks', lambda c: c['passes'][0]['guarded']['calls'].pop())
corruption('rebind_geometry_reject', 'deleted failed rebind', lambda c: c['passes'][0].pop('rebind'))
corruption('aggregate_limit', 'pre-admission failure forged as solved', lambda c: c.update(last_solved=c['passes'][0]['source']))
print(f'PASS: {rejected} actual-export corruptions rejected')
