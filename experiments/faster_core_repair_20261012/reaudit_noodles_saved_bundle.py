"""Reaudit the persisted terminal, never reconstruct its optimization run."""
import hashlib
import json
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def worker(out):
    from make_it_stand.verification import spin_bundle as bundle
    source=HERE/'noodles-terminal-recovery-job-03/output/material'
    manifest=json.loads((source/'manifest.json').read_text())
    problem=json.loads((source/'problem.json').read_text())
    registration=json.loads((source.parent/'registration.json').read_text())
    assert problem['input_digest']==registration['input_digest']
    checked={}
    def check(x):
        if isinstance(x,dict):
            if 'name' in x and 'sha256' in x:
                f=source/x['name']
                assert f.parent==source and not f.is_symlink()
                data=f.read_bytes();assert hashlib.sha256(data).hexdigest()==x['sha256']
                checked[f.name]=len(data)
            for v in x.values():check(v)
        elif isinstance(x,list):
            for v in x:check(v)
    check(manifest);check(problem)
    out=Path(out);out.mkdir(exist_ok=False)
    (out/'integrity.json').write_text(json.dumps(checked,indent=2))
    original=bundle._Reader.mesh
    def timed(reader,value,**kw):
        label=value.get('vertices',{}).get('name','empty');start=time.monotonic()
        print('MESH_START',label,flush=True)
        result=original(reader,value,**kw)
        print('MESH_RETURN',label,time.monotonic()-start,flush=True)
        return result
    bundle._Reader.mesh=timed
    start=time.monotonic()
    artifact,receipt=bundle.reaudit_spin_artifact(source,
        expected_problem_digest=registration['input_digest'],
        expected_scientific_digest=manifest['scientific_digest'])
    (out/'report.json').write_text(json.dumps(dict(receipt=receipt,
        scientific_digest=artifact.scientific_digest,seconds=time.monotonic()-start,
        original_run_complete=False,case_complete=False),indent=2))
    print('REAUDIT_PASSED',flush=True)

def launch():
    sys.path.insert(0,str(ROOT/'python'))
    from make_it_stand.verification.serial_source_job import run_serial_source_job
    native=HERE/'formal-native-core-build-job-12/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    source=HERE/'noodles-terminal-recovery-job-03/output/material'
    files={'worker':Path(__file__),'native':native,'registration':source.parent/'registration.json'}
    files.update({'bundle_'+p.name:p for p in source.iterdir() if p.is_file()})
    files.update({'checker_'+p.name:p for p in (ROOT/'python/make_it_stand/verification').glob('*.py')})
    files['spin']=ROOT/'python/make_it_stand/spin.py'
    result=run_serial_source_job(ROOT,HERE/'noodles-terminal-bundle-reaudit-job-01',
        ['/usr/bin/env','LD_PRELOAD='+str(native.parent.parent/'compiler-runtime/libstdc++.so.6'),
         '/home/lin/anaconda3/bin/python','-B',str(Path(__file__)),'{job}/output'],files,
        environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
    print(result.status,result.returncode,flush=True)

if __name__=='__main__':
    if len(sys.argv)>1:worker(sys.argv[1])
    else:launch()
