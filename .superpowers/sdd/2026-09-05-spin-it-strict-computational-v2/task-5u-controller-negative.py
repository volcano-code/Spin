"""Real exported geometry/facts mutations, not user-entered approval booleans."""
import contextlib
import copy
import importlib.util
import io
import json
from pathlib import Path
import sys

spec = importlib.util.spec_from_file_location('geometry_check', Path(__file__).with_name('task-5u-controller-check.py'))
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)
cases = [json.loads(x) for x in open(sys.argv[1], encoding='utf-8')]
assert len(cases) == 3
for case in cases: checker.verify(case)
count = 0
for case in (cases[0], cases[-1]):
    changes = [lambda c: c.update(scope='AUTHOR_ORIGINAL'),
               lambda c: c.update(feasibility='FEASIBLE_WITNESS'),
               lambda c: c.update(sealed=True),
               lambda c: c.update(converged=True),
               lambda c: c.update(status='CONVERGED'),
               lambda c: c.update(id='paper_elephant'),
               lambda c: c['world_moments'].__setitem__(0, c['world_moments'][0]+1),
               lambda c: c['primitives'][0].update(distance_squared_exact='1'),
               lambda c: c['primitives'][0].update(stage='UNKNOWN'),
               lambda c: c.update(pairs_visited=0),
               lambda c: c['primitives'][0].update(owner=len(c['beta'])),
               lambda c: c['primitives'][0].update(id=999),
               lambda c: c['points'][0].__setitem__(0, float('nan')),
               lambda c: c.update(primitives=[]),
               lambda c: c['outer_faces'][0].__setitem__(1, c['outer_faces'][0][0]),
               lambda c: c['events'][3].update(stage=0),
               lambda c: c['events'][4].update(outcome=5),
               lambda c: c['events'].pop(),
               lambda c: c.update(pair_count=0)]
    for change in changes:
        broken = copy.deepcopy(case)
        change(broken)
        try:
            with contextlib.redirect_stdout(io.StringIO()): checker.verify(broken)
        except (AssertionError, ValueError, IndexError):
            count += 1
        else:
            raise AssertionError('corrupted geometry evidence was admitted')
print(f'PASS {count}/38 actual geometry/material/status corruptions rejected')
