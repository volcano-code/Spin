"""Corrupt actual exports; no fabricated positive material or passing flags."""
import copy
import importlib.util
from pathlib import Path
import sys

spec = importlib.util.spec_from_file_location('regions_checker',Path(__file__).with_name('task-5z-controller-check.py'))
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)
rows, _ = checker.check_file(sys.argv[1])


def set_path(row, path, value):
    for item in path[:-1]: row = row[item]
    row[path[-1]] = value


mutations = [
    (1,0,['published'],False),
    (1,0,['same_binary'],False),
    (1,0,['full_thickness'],'PASSED'),
    (1,0,['binary_input_hash'],'0'*64),
    (1,0,['evaluation_hash'],'0'*64),
    (1,0,['material_partition','equal_union'],False),
    (1,0,['material_partition','no_overlap'],False),
    (1,0,['material_partition','missing'],True),
    (1,0,['material_partition','parts'],[]),
    (1,0,['selected_partition','complete_partition'],False),
    (1,0,['selected'],[]),
    (1,0,['selected',0,'owner'],999),
    (1,0,['selected',0,'key'],[0,0,0]),
    (1,0,['selected',0,'mesh_hash'],'0'*64),
    (1,0,['selected',0,'subset'],False),
    (1,0,['material','void_components'],0),
    (1,0,['material','input'],'0'*64),
    (1,0,['clearance','squared_distance'],'1/100'),
    (1,0,['clearance','minimum_exact'],'0'),
    (1,0,['clearance','contact'],True),
    (1,0,['clearance','witness'],[999,999]),
    (1,1,['published'],True),
    (1,2,['published'],True),
    (0,0,['clearance','squared_distance'],'0'),
    (0,0,['selected_partition','parts'],[]),
    (3,0,['published'],True),
]
for number,(row_id,region_id,path,value) in enumerate(mutations,1):
    row = copy.deepcopy(rows[row_id])
    set_path(row['regions'][region_id],path,value)
    try: checker.check_regions(row,row['regions'][region_id],region_id)
    except (AssertionError,IndexError,KeyError,ValueError,TypeError): pass
    else: raise AssertionError(f'corruption {number} accepted: {path}')
print(f'PASS: {len(mutations)} actual-export corruptions rejected; no native or paper admission from JSON')
