"""One registered all-original-input construction, no rounded intermediate operands."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time
import gmpy2
import numpy as np
import trimesh
import _make_it_stand as native

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'experiments/quadric_physics_demo_20261011'))
from core_boolean_mesh import binary64_rational_token,manifest_design_seed
case=ROOT/'experiments/quadric_physics_demo_20261011/noodles-dense-v5'
manifest_path=case/'manifest.json';manifest=json.loads(manifest_path.read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def load(name):
    entry=manifest['parts'][name];path=case/entry['model_obj']
    if sha(path)!=entry['sha256']['model']:raise ValueError('FROZEN_INPUT_CHANGED:'+name)
    mesh=trimesh.load(path,force='mesh',process=False)
    return [[binary64_rational_token(x) for x in row] for row in mesh.vertices],np.asarray(mesh.faces,dtype=np.int32)
if len(manifest['construction']['fixed_parts'])!=81:raise ValueError('EXPECTED_ALL81_FIXED_PARTS')
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
seed_name=manifest_design_seed(manifest)
vertices,faces=load(seed_name);started=time.monotonic()
with (output/'events.jsonl').open('x') as events:
    def event(status,**values):
        events.write(json.dumps(dict(status=status,seconds=time.monotonic()-started,**values))+'\n')
        events.flush();os.fsync(events.fileno())
    event('ALL_ORIGINAL_RATIONAL_CONSTRUCTION_STARTED',profile='CGAL_CORE_RATIONAL_ALL_INPUTS_V2_MANIFEST_DOMAIN',
        domain_seed=seed_name,domain_seed_sha256=manifest['parts'][seed_name]['sha256']['model'],
        manifest_sha256=sha(manifest_path),geometry_admission=False,optimization_admission=False)
    for step,name in enumerate(manifest['construction']['fixed_parts'],1):
        if shutil.disk_usage(output).free<20*1024**3:
            event('RESOURCE_LIMIT_DISK',step=step);raise RuntimeError('DISK_BELOW20_GIB')
        cutter,cut_faces=load(name);event('BOOLEAN_STARTED',step=step,part=name)
        packet=native.libigl_exact_boolean_construction_from_rationals(vertices,faces,cutter,cut_faces,'difference')
        with (output/f'domain-{step:03d}-rational.json').open('x') as stream:
            json.dump(dict(packet),stream,default=lambda x:x.tolist(),allow_nan=False)
        if not packet['succeeded']:raise ValueError('ACTUAL_BOOLEAN_FAILED')
        vertices=packet['vertices'];faces=np.asarray(packet['faces'],dtype=np.int32)
        event('BOOLEAN_RETURNED_PENDING_BOUNDARY_CHECK',step=step,part=name,
            vertices=len(vertices),faces=len(faces))
    with (output/'preparation-report.json').open('x') as stream:
        json.dump(dict(profile='CGAL_CORE_RATIONAL_ALL_INPUTS_V2_MANIFEST_DOMAIN',domain_seed=seed_name,
            excluded_parts=manifest['construction']['fixed_parts'],manifest_sha256=sha(manifest_path),
            native_sha256=sha(Path(native.__file__)),seconds=time.monotonic()-started,
            geometry_admission=False,optimization_admission=False),stream)
    event('ALL81_ORIGINAL_INPUT_PARTS_RETURNED_PENDING_BOUNDARY_CHECK')
