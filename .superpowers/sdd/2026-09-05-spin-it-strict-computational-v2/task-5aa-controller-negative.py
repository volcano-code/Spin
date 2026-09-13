"""Mutation checks on actual exports, never a native admission route."""
import copy
import importlib.util
from pathlib import Path
import sys
spec=importlib.util.spec_from_file_location('component_checker',Path(__file__).with_name('task-5aa-controller-check.py'))
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
rows,_=c.check_file(sys.argv[1]);by_id={r['id']:r for r in rows}
changes=[
('two_cavities',['pairs'],[]),('two_cavities',['pair_count'],1),
('two_cavities',['face_components'],[]),('two_cavities',['components',1,'faces'],[]),
('two_cavities',['components',1,'vertices'],[]),('two_cavities',['components',1,'sign'],1),
('two_cavities',['components',1,'reversed'],False),('two_cavities',['components',1,'parent'],-1),
('two_cavities',['components',1,'original_input'],'0'*64),
('two_cavities',['components',1,'checker_input'],'0'*64),
('two_cavities',['pairs',0,'distance'],'1'),('two_cavities',['pairs',0,'parent_witness'],[999,999]),
('two_cavities',['pairs',0,'local_witness'],[999,999]),
('two_cavities',['pairs',0,'components'],[0,0]),('two_cavities',['pairs',0,'completed'],False),
('two_cavities',['pairs',0,'minimum_exact'],'0'),('two_cavities',['distance'],'1'),
('two_cavities',['minimum_components'],[0,0]),('two_cavities',['minimum_parent_witness'],[999,999]),
('two_cavities',['material_input'],'0'*64),('two_cavities',['same_binary'],False),
('two_cavities',['full_thickness'],'PASSED'),('two_cavities',['input_hash'],'0'*64),
('empty',['distance'],'0'),('empty',['applicable'],True),
('two_cavities_rejected',['published'],True),('threshold_rejected',['published'],True),
('resource_limited',['published'],True),('not_admitted',['published'],True),
]
for case,path,value in changes:
    row=copy.deepcopy(by_id[case]);target=row
    for key in path[:-1]:target=target[key]
    target[path[-1]]=value
    try:c.verify(row)
    except (AssertionError,ValueError,IndexError,KeyError,TypeError):pass
    else:raise AssertionError(f'accepted corrupt actual {case}:{path}')
print(f'PASS:{len(changes)} actual-output corruptions rejected; no JSON publication authority')
