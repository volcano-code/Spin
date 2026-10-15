"""Re-use the pinned exact audit executable on the new actual rational output."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
HERE=Path(__file__).resolve().parent
source=HERE/'noodles-repaired-rational-suffix-job-01/output/domain-081-rational.json'
if len(sys.argv)==3:source=Path(sys.argv[2])
exe=HERE/'noodles-rational-authority-diagnosis-job-01/output/inspect-rational-boundary'
packet=json.loads(source.read_text())
if not packet['succeeded']:raise ValueError('FAILED_RATIONAL_INPUT')
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
data=output/'boundary.txt'
with data.open('x') as stream:
    stream.write(f"{len(packet['vertices'])} {len(packet['faces'])}\n")
    for point in packet['vertices']:stream.write(' '.join(point)+'\n')
    for face in packet['faces']:stream.write(' '.join(map(str,face))+'\n')
started=time.monotonic()
with (output/'report.json').open('x') as stream,(output/'stderr.log').open('x') as errors:
    subprocess.run([str(exe),str(data)],stdout=stream,stderr=errors,check=True)
report=json.loads((output/'report.json').read_text())
with (output/'receipt.json').open('x') as stream:
    json.dump(dict(input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
        seconds=time.monotonic()-started,optimization_admission=False),stream)
print(json.dumps(dict(boundary_valid=report['boundary_valid'],
    non_manifold_vertices=report['non_manifold_vertices'],
    intersections=len(report['intersecting_pairs']),optimization_admission=False)),flush=True)
