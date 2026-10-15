"""Instrumented native call graph on an actual saved boundary, not solver time."""
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys
import numpy as np
from make_it_stand import spin

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def main(output,mode='callgraph'):
    output.mkdir(parents=True,exist_ok=False)
    build=HERE/'formal-native-core-build-job-04/output/build'
    rows=json.loads((build/'compile_commands.json').read_text())
    template=shlex.split(next(r['command'] for r in rows if Path(r['file']).name=='ExactBoundaryAudit.cpp'))
    commands=[]
    for name,source in (('audit',ROOT/'src/verification/ExactBoundaryAudit.cpp'),
                        ('driver',HERE/'profile_exact_boundary.cpp')):
        command=list(template)
        command[command.index('-c')+1]=str(source)
        command[command.index('-o')+1]=str(output/(name+'.o'))
        if mode=='callgraph':command.extend(['-pg','-fno-omit-frame-pointer'])
        elif mode=='phases':command.append('-DMAKE_IT_STAND_PROFILE_EXACT_BOUNDARY=1')
        else:raise ValueError('UNREGISTERED_PROFILE_MODE')
        commands.append(command)
    executable=output/'profile-boundary'
    commands.append(['/usr/bin/c++',*(['-pg'] if mode=='callgraph' else []),'-o',str(executable),str(output/'driver.o'),
                     str(output/'audit.o'),'-lgmp','-lmpfr'])
    with (output/'build.log').open('x') as log:
        for command in commands:
            subprocess.run(command,cwd=ROOT,stdout=log,stderr=log,check=True)
    run=HERE/'formal-core-benchmark-faces-12288-seed-0-job-01/output/run'
    terminal=json.loads((run/'terminal-material-reference.json').read_text())['evaluation_file']
    path=run/terminal
    material=json.loads(path.read_text())
    v,f=material['variable_rational_vertices'],material['variable_faces']
    facts=spin._native_module().audit_exact_rational_triangle_boundary(v,np.asarray(f,dtype=np.int64),'AnyNonzero')
    if not facts['boundary_valid']:raise ValueError('ACTUAL_NATIVE_BASELINE_BOUNDARY_INVALID')
    with (output/'mesh.txt').open('x') as tape:
        tape.write(f'{len(v)} {len(f)}\n')
        for point in v:tape.write(' '.join(point)+'\n')
        for face in f:tape.write(' '.join(map(str,face))+'\n')
    with (output/'driver.log').open('x') as log,(output/'phases.log').open('x') as phases:
        subprocess.run([str(executable),str(output/'mesh.txt'),'10'],cwd=output,
                       stdout=log,stderr=phases,check=True)
    lines=(output/'driver.log').read_text().splitlines()
    samples=[]
    for line in lines:
        fields=line.split()
        if (fields[1]!=facts['input_digest'] or fields[2]!=facts['backend'] or
            int(fields[3])!=facts['exact_pair_tests'] or
            fields[5:]!=facts['component_signed_volumes_exact']):
            raise ValueError('INSTRUMENTED_AUDIT_CHANGED_REGISTERED_EXACT_FACTS')
        samples.append(int(fields[0])/1e9)
    phase_rows=[];pair_rows=[]
    if mode=='callgraph':
        with (output/'gprof.txt').open('x') as profile:
            subprocess.run(['/usr/bin/gprof','-b',str(executable),str(output/'gmon.out')],
                           cwd=output,stdout=profile,stderr=profile,check=True)
    else:
        for line in (output/'phases.log').read_text().splitlines():
            fields=line.split()
            if fields[0]=='AUDIT_PHASE' and len(fields)==3:
                phase_rows.append(dict(stage=fields[1],seconds=int(fields[2])/1e9))
            elif fields[0]=='AUDIT_PAIR' and len(fields)==4:
                pair_rows.append(dict(shared_vertices=int(fields[1]),count=int(fields[2]),seconds=int(fields[3])/1e9))
            else:raise ValueError('UNKNOWN_NATIVE_PROFILE_RECORD')
        if len(phase_rows)!=50:raise ValueError('INCOMPLETE_NATIVE_PHASE_CAPTURE')
    (output/'report.json').write_text(json.dumps(dict(actual_fixed_boundary_agrees=True,faces=len(f),
        mode=mode,phase_rows=phase_rows,pair_rows=pair_rows,
        instrumented_samples_seconds=samples,commands=commands,core_admission=False,
        material_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256((ROOT/'src/verification/ExactBoundaryAudit.cpp').read_bytes()).hexdigest(),
        scope='Instrumented native call graph; shared-library GMP time and instrumentation bias limit attribution; not a production timing gate')))


if __name__=='__main__':main(Path(sys.argv[1]),sys.argv[2] if len(sys.argv)>2 else 'callgraph')
