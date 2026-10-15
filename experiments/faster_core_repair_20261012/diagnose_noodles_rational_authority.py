"""Distinguish an inherited authority defect from binary64 export defects."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
source=HERE/'inspect_rational_boundary.cpp';obj=output/'inspect.o';exe=output/'inspect-rational-boundary'
template=json.loads((ROOT/'build/heart-fixed-facet-20261011-v1/commands.json').read_text())[0]
command=list(template)
for flag,value in (('-c',str(source)),('-o',str(obj)),('-MF',str(obj)+'.d'),('-MT',str(obj))):
    command[command.index(flag)+1]=value
library=ROOT/'build/spin2014-objectives-quadric-20261007/libmake_it_stand_core.a'
link=[command[0],str(obj),str(library),'-L/usr/lib/x86_64-linux-gnu','-lgmp','-lmpfr','-lpthread','-o',str(exe)]
with (output/'build.log').open('x') as stream:
    subprocess.run(command,cwd=ROOT,stdout=stream,stderr=stream,check=True)
    subprocess.run(link,cwd=ROOT,stdout=stream,stderr=stream,check=True)
packet_path=ROOT/'experiments/quadric_physics_demo_20261011/noodles-dense-v5-core-rational-suffix-v1/domain-job-01/output/domain-081-rational.json'
packet=json.loads(packet_path.read_text())
data=output/'boundary.txt'
with data.open('x') as stream:
    stream.write(f"{len(packet['vertices'])} {len(packet['faces'])}\n")
    for point in packet['vertices']:stream.write(' '.join(point)+'\n')
    for face in packet['faces']:stream.write(' '.join(map(str,face))+'\n')
started=time.monotonic()
with (output/'report.json').open('x') as stream, (output/'stderr.log').open('x') as errors:
    subprocess.run([str(exe),str(data)],stdout=stream,stderr=errors,check=True)
report=json.loads((output/'report.json').read_text())
with (output/'receipt.json').open('x') as stream:
    json.dump(dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
        input_sha256=hashlib.sha256(packet_path.read_bytes()).hexdigest(),
        executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
        commands=[command,link],audit_seconds=time.monotonic()-started),stream)
print(json.dumps(dict(boundary_valid=report['boundary_valid'],vertex_manifold=report['vertex_manifold'],
    intersections=len(report['intersecting_pairs']))),flush=True)
