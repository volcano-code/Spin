"""Same fixed parts, retain the proved step71/72 rational construction thereafter."""
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import gmpy2
import numpy as np
import trimesh
import _make_it_stand as native

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
demo=ROOT/'experiments/quadric_physics_demo_20261011';case=demo/'noodles-dense-v5'
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
source=HERE/'noodles-rational71-72-diagnosis-job-01/output/original-072-rational.json'
manifest=json.loads((case/'manifest.json').read_text())
packet=json.loads(source.read_text())
if not packet['succeeded']:raise ValueError('FAILED_STARTING_BOOLEAN')
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
with (output/'events.jsonl').open('x') as events:
    def event(status,**values):
        events.write(json.dumps(dict(status=status,seconds=time.monotonic()-started,**values))+'\n')
        events.flush();os.fsync(events.fileno())
    event('RATIONAL_SUFFIX_STARTED',excluded=72,input_sha256=sha(source),
        profile='CGAL_CORE_RATIONAL_FROM71_V1',geometry_admission=False)
    for index,name in enumerate(manifest['construction']['fixed_parts'][72:],73):
        entry=manifest['parts'][name];path=case/entry['model_obj']
        if sha(path)!=entry['sha256']['model']:raise ValueError('FIXED_PART_CHANGED:'+name)
        part=trimesh.load(path,force='mesh',process=False)
        event('BOOLEAN_STARTED',step=index,part=name)
        cutter=[[str(gmpy2.mpq(float(x))) for x in row] for row in part.vertices]
        packet=native.libigl_exact_boolean_construction_from_rationals(packet['vertices'],
            np.asarray(packet['faces'],dtype=np.int32),cutter,part.faces,'difference')
        with (output/f'domain-{index:03d}-rational.json').open('x') as stream:
            json.dump(dict(packet),stream,default=lambda x:x.tolist(),allow_nan=False)
        if not packet['succeeded']:raise ValueError('ACTUAL_BOOLEAN_FAILED')
        event('BOOLEAN_RETURNED_PENDING_BOUNDARY_CHECK',step=index,part=name,
            vertices=len(packet['vertices']),faces=len(packet['faces']))
    with (output/'preparation-report.json').open('x') as stream:
        json.dump(dict(profile='CGAL_CORE_RATIONAL_FROM71_V1',excluded_parts=manifest['construction']['fixed_parts'],
            input_sha256=sha(source),manifest_sha256=sha(case/'manifest.json'),
            seconds=time.monotonic()-started,geometry_admission=False,optimization_admission=False),stream)
    event('ALL81_PARTS_RETURNED_PENDING_BOUNDARY_CHECK')
