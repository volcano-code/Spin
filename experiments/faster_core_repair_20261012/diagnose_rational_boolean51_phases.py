"""One fixed diagnostic plus native correctness comparison; never continues the input sequence."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
import trimesh
import _make_it_stand as native
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
case=ROOT/'experiments/quadric_physics_demo_20261011/noodles-dense-v5'
original=HERE/'noodles-all-original-rational-job-02'
prefix_path=original/'output/domain-050-rational.json'
manifest_path=case/'manifest.json';manifest=json.loads(manifest_path.read_text())
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
first=json.loads((original/'output/events.jsonl').read_text().splitlines()[0])
if first['manifest_sha256']!=sha(manifest_path):raise ValueError('FROZEN_MANIFEST_CHANGED')
name=manifest['construction']['fixed_parts'][50];entry=manifest['parts'][name];part_path=case/entry['model_obj']
if sha(part_path)!=entry['sha256']['model']:raise ValueError('FIXED_PART_CHANGED')
prefix=json.loads(prefix_path.read_text());part=trimesh.load(part_path,force='mesh',process=False)
cut=[[str(n)+'/'+str(d) for n,d in (float(x).as_integer_ratio() for x in row)] for row in part.vertices]
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
source=HERE/'profile_rational_boolean.cpp';obj=output/'profile.o';exe=output/'profile-rational-boolean'
precompiled=Path(sys.argv[2]) if len(sys.argv)==3 else None
if precompiled is not None:
    origin=json.loads((precompiled.parent.parent/'config.json').read_text())
    if origin['identity']['probe']['sha256']!=sha(source):raise ValueError('PRECOMPILED_SOURCE_CHANGED')
    obj=precompiled
compile=list(json.loads((ROOT/'build/heart-fixed-facet-20261011-v1/commands.json').read_text())[0])
for flag,value in (('-c',str(source)),('-o',str(obj)),('-MF',str(obj)+'.d'),('-MT',str(obj))):
    compile[compile.index(flag)+1]=value
library=ROOT/'build/spin2014-objectives-quadric-20261007/libmake_it_stand_core.a'
receipt=HERE/'exact-core-material-build-job-01/output/receipt.json'
boost=list(dict.fromkeys(x for x in json.loads(receipt.read_text())['commands'][-1]['argv'] if 'libboost_' in x and x.endswith('.a')))
if not boost:raise ValueError('PINNED_BOOST_LINK_INPUTS_MISSING')
link=[compile[0],str(obj),str(library),'-Wl,--start-group',*boost,'-Wl,--end-group',
    '-L/usr/lib/x86_64-linux-gnu','-lgmp','-lmpfr','-lpthread','-o',str(exe)]
commands=[]
with (output/'build.log').open('x') as stream:
    if precompiled is None:
        subprocess.run(compile,cwd=ROOT,stdout=stream,stderr=stream,check=True,timeout=1800);commands.append(compile)
    subprocess.run(link,cwd=ROOT,stdout=stream,stderr=stream,check=True,timeout=300)
    commands.append(link)
input_file=output/'operands.txt';started=time.monotonic()
with input_file.open('x') as stream:
    for vertices,faces in ((prefix['vertices'],prefix['faces']),(cut,part.faces)):
        stream.write(f'{len(vertices)} {len(faces)}\n')
        for row in vertices:stream.write(' '.join(row)+'\n')
        for row in faces:stream.write(' '.join(map(str,row))+'\n')
preparation_seconds=time.monotonic()-started
with (output/'stages.jsonl').open('x') as stream,(output/'native-stderr.log').open('x') as errors:
    subprocess.run([str(exe),str(input_file),str(output/'actual-boolean.json')],stdout=stream,stderr=errors,check=True)
diagnostic=json.loads((output/'actual-boolean.json').read_text());start=time.monotonic()
baseline=native.libigl_exact_boolean_construction_from_rationals(prefix['vertices'],
    np.asarray(prefix['faces'],dtype=np.int32),cut,part.faces,'difference')
baseline_seconds=time.monotonic()-start
with (output/'official-native-boolean.json').open('x') as stream:
    json.dump(dict(baseline),stream,default=lambda x:x.tolist(),allow_nan=False)
def signature(packet):
    points=[tuple(Fraction(x) for x in row) for row in packet['vertices']];result=[]
    for face in packet['faces']:
        tri=tuple(points[int(v)] for v in face);result.append(min(tri,tri[1:]+tri[:1],tri[2:]+tri[:2]))
    return sorted(result)
equal=signature(diagnostic)==signature(baseline)
stages=[json.loads(line) for line in (output/'stages.jsonl').read_text().splitlines()]
timings={};pending={}
for row in stages:
    key=row['stage']
    if key.endswith('-begin'):pending[key[:-6]]=row['seconds']
    elif key.endswith('-end'):timings[key[:-4]]=row['seconds']-pending.pop(key[:-4])
report=dict(role='FIXED_BOOLEAN51_PHASE_DIAGNOSIS_ONLY',part=name,operation='difference',
    source_sha256=sha(source),prefix_sha256=sha(prefix_path),part_sha256=sha(part_path),
    library_sha256=sha(library),native_sha256=sha(Path(native.__file__)),executable_sha256=sha(exe),
    commands=commands,precompiled_object_sha256=sha(obj),
    boost_input_sha256={path:sha(Path(path)) for path in boost},operand_text_preparation_seconds=preparation_seconds,
    phase_seconds=timings,official_native_seconds=baseline_seconds,
    identical_exact_oriented_triangle_boundary=equal,geometry_admission=False,optimization_admission=False)
with (output/'report.json').open('x') as stream:json.dump(report,stream,indent=2)
print(json.dumps(dict(phase_seconds=timings,official_native_seconds=baseline_seconds,
    identical_exact_boundary=equal,optimization_admission=False)),flush=True)
if not equal:raise ValueError('DIAGNOSTIC_BOOLEAN_DIFFERS_FROM_NATIVE')
