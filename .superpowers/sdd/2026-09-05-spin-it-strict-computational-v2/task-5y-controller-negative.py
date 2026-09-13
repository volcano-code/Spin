"""Corrupt actual native exports, require independent checker rejection."""
import copy
import importlib.util
import json
from pathlib import Path
import sys

spec = importlib.util.spec_from_file_location('binary_oracle', Path(__file__).with_name('task-5y-controller-check.py'))
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)
rows, _ = checker.check_file(sys.argv[1])
by_id = {r['id']: r for r in rows}
mutations = []


def add(name, key, change):
    row = copy.deepcopy(by_id[key])
    change(row)
    mutations.append((name, row))


def change_point(c): c['material']['vertices'][0][0] += .003
def change_face(c): c['material']['faces'][0][0], c['material']['faces'][0][1] = c['material']['faces'][0][1], c['material']['faces'][0][0]
def change_handle(c): c['source']['state'][c['source']['handle_offset']+1] += .01
def change_phi(c): c['source']['state'][c['source']['betas']] += .11
def change_moment(c): c['world_integrated'][7] += .001
def change_cavity(c): c['cavity']['faces'].pop()
def change_rounded(c): c['source']['state'][0] = 0
def change_delta(c): c['delta'][0] = 0
def change_inertia(c): c['central_inertia'][0][1] += .001
def change_constraints(c): c['constraints'][0] = 0
def change_origin(c): c['origin'][0] += .01


add('changed_final_point','full_affine',change_point)
add('reversed_final_triangle','full_affine',change_face)
add('changed_handle','full_affine',change_handle)
add('changed_givens_angle','full_affine',change_phi)
add('wrong_cross_moment','full_affine',change_moment)
add('deleted_cavity_triangle','full_affine',change_cavity)
add('wrong_half_tie','tie_affine',change_rounded)
add('fabricated_delta','tie_affine',change_delta)
add('wrong_nondiagonal_inertia','full_affine',change_inertia)
add('hidden_rounding_constraint','mixed_rounding_failure',change_constraints)
add('changed_axis_origin','full_affine',change_origin)
add('handwritten_validated','mixed_rounding_failure',lambda c:c.update(validated=True))
add('false_empty_cavity','full_affine',lambda c:c.update(empty_cavity=True))
add('nonfinite_world','full_affine',lambda c:c['world_integrated'].__setitem__(0,float('nan')))
add('relaxed_tolerance','mixed_rounding_failure',lambda c:c.update(constraint_tolerance=.1))
add('false_seal','full_affine',lambda c:c.update(sealed=True))

for name, row in mutations:
    try:
        checker.verify(row)
    except (AssertionError, ValueError, IndexError, KeyError):
        print('REJECTED:',name)
    else:
        raise AssertionError('invalid export accepted: '+name)
print(f'PASS: {len(mutations)} actual-export corruptions rejected; not a general release gate')
