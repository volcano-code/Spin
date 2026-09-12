"""Mutate actual exports; these tests assess the independent checker, not native admission."""
import copy
import importlib.util
import json
from pathlib import Path
import sys

spec = importlib.util.spec_from_file_location("check_rebind", Path(__file__).with_name("task-5w-controller-check.py"))
check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)


def set_at(obj, path, value):
    for component in path[:-1]:
        obj = obj[component]
    obj[path[-1]] = value


rows = [json.loads(line) for line in Path(sys.argv[1]).read_text().splitlines()]
assert rows
count = 0
for row in rows:
    check.verify(row)
    mutations = [
        ('stale_mass', ['after', 'moments', 0], row['after']['moments'][0]*1.01),
        ('nonfinite_material', ['after', 'moments', 4], float('nan')),
        ('replaced_bbw', ['after', 'bbw_hash'], '0'*64),
        ('wrong_handle_layout', ['after', 'handle_offset'], row['after']['handle_offset']+1),
        ('omitted_handle', ['after', 'handles'], 25),
        ('wrong_density', ['after', 'density'], 2.),
        ('wrong_point', ['after', 'points', -1, 0], row['after']['points'][-1][0]+.01),
        ('invalid_owner', ['after', 'primitives', 0, 'owner'], row['after']['betas']),
        ('stale_beta_derivative', ['after', 'beta_jacobian', 0, 0], 1.),
        ('wrong_rest_moment', ['after', 'leaves', 0, 'moments', 0], .2),
    ]
    for name, path, value in mutations:
        damaged = copy.deepcopy(row)
        set_at(damaged, path, value)
        try:
            check.verify(damaged)
        except AssertionError:
            count += 1
            print(f'REJECTED {row["id"]}:{name}')
        else:
            raise AssertionError(f'corrupt export admitted: {row["id"]}:{name}')
assert count == 60
print(f'PASS: {count} actual-export corruptions rejected by the independent checker')
