"""Corrupt actual exported observations; the independent checker must reject."""
import contextlib
import copy
import importlib.util
import io
from pathlib import Path
import sys
from collections import defaultdict

spec = importlib.util.spec_from_file_location('paper_checker', Path(__file__).with_name('task-5r-controller-check.py'))
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)
cases, current = [], None
with open(sys.argv[1], encoding='utf-8') as stream:
    for line in stream:
        pieces = line.split()
        if pieces[0] == 'CASE': current = defaultdict(list)
        if pieces[0] == 'END': cases.append(current);current = None
        else: current[pieces[0]].append(pieces[1:])
assert len(cases) == 10 and current is None
rejected = 0
for case in (cases[0], cases[-1]):
    for key, row, column in [('OBJ',0,0), ('SPIN',0,4), ('CONSTRAINT',0,0),
                             ('GRAD',0,0), ('JM',0,1), ('JC',0,1),
                             ('VMOVED',0,1), ('VOBJ',0,0), ('FLAGS',0,2),
                             ('R',0,1)]:
        broken = copy.deepcopy(case)
        broken[key][row][column] = str(float(broken[key][row][column])+.1)
        try:
            with contextlib.redirect_stdout(io.StringIO()): checker.verify(broken)
        except (AssertionError, ValueError):
            rejected += 1
        else:
            raise AssertionError(f'corruption incorrectly admitted: {key}')
print(f'PASS: {rejected}/20 objective/constraint/Jacobian/value-only/frame corruptions rejected.')
