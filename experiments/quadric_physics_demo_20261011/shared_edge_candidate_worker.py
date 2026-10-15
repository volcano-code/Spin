"""Build/test isolated exact-predicate candidates; benchmark the same boundary."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def build(output):
    output.mkdir(parents=True,exist_ok=False)
    cgal=ROOT/'libigl/external/cgal'
    include=[str(ROOT/'include'),str(ROOT/'libigl/external/eigen'),'/home/lin/anaconda3/include']
    include.extend(str(p) for p in sorted(cgal.glob('*/include')))
    flags=['-std=c++17','-O3','-DNDEBUG','-frounding-math','-fno-omit-frame-pointer','-g1']
    flags+=['-DCGAL_HEADER_ONLY=1','-DCGAL_USE_CORE=1']
    flags+=[arg for path in include for arg in ['-I',path]]
    log=(output/'build-and-tests.log').open('x')
    commands=[]
    def run(args):
        commands.append(args)
        log.write(json.dumps(args)+'\n');log.flush()
        subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,check=True)
    try:
        run(['/usr/bin/c++',*flags,str(HERE/'shared_edge_predicate_test.cpp'),'-lgmp','-lmpfr','-o',str(output/'predicate-test')])
        run([str(output/'predicate-test')])
        for name,defines in [('baseline',['-DSPIN_DISABLE_EXACT_SHARED_EDGE_SHORTCUT']),('candidate',[])]:
            obj=output/(name+'-audit.o')
            run(['/usr/bin/c++',*flags,*defines,'-c',str(ROOT/'src/verification/ExactBoundaryAudit.cpp'),'-o',str(obj)])
            run(['/usr/bin/c++',*flags,str(ROOT/'tests/exact_boundary_audit_regression.cpp'),str(obj),'-lgmp','-lmpfr','-o',str(output/(name+'-regression'))])
            run([str(output/(name+'-regression'))])
            run(['/usr/bin/c++',*flags,str(HERE/'exact_boundary_benchmark.cpp'),str(obj),'-lgmp','-lmpfr','-o',str(output/(name+'-benchmark'))])
    finally:log.close()
    binaries={name:hashlib.sha256((output/name).read_bytes()).hexdigest()
              for name in ('baseline-benchmark','candidate-benchmark')}
    (output/'build.json').write_text(json.dumps(dict(commands=commands,binaries=binaries,
        frozen_production_module_changed=False,material_admission=False),indent=2))


def benchmark(binary,source,output):
    output.mkdir(parents=True,exist_ok=False)
    raw=Path(source).read_bytes();packet=json.loads(raw)
    input_file=output/'boundary.txt'
    with input_file.open('x') as stream:
        stream.write(f"{len(packet['vertices'])} {len(packet['faces'])}\n")
        for row in packet['vertices']:stream.write(' '.join(row)+'\n')
        for row in packet['faces']:stream.write(' '.join(str(x) for x in row)+'\n')
    started=time.monotonic()
    with (output/'benchmark.log').open('x') as log:
        subprocess.run([str(binary),str(input_file),str(output/'audit.json')],stdout=log,stderr=subprocess.STDOUT,check=True)
    (output/'identity.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(raw).hexdigest(),
        binary_sha256=hashlib.sha256(Path(binary).read_bytes()).hexdigest(),
        elapsed_seconds=time.monotonic()-started,geometry_admission=False,optimization_admission=False),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('operation',choices=('build','benchmark'))
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--source',type=Path)
    parser.add_argument('--binary',type=Path);args=parser.parse_args()
    if args.operation=='build':build(args.output)
    else:benchmark(args.binary,args.source,args.output)
