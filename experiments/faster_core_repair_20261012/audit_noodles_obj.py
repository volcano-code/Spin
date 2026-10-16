"""Independent binary64 geometry audit of the literal terminal OBJ export."""
from pathlib import Path
import json, sys, time
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def worker(out):
    import numpy as np
    import _make_it_stand as native
    out=Path(out);out.mkdir(exist_ok=False)
    source=HERE/'noodles-terminal-obj-01'
    rows={}
    for name in ('complete','retained','cavity','fixed','original','design'):
        v=[];f=[]
        for line in (source/(name+'.obj')).read_text().splitlines():
            if line.startswith('v '):v.append([float(x) for x in line.split()[1:]])
            elif line.startswith('f '):f.append([int(x)-1 for x in line.split()[1:]])
        start=time.monotonic()
        result=dict(native.audit_exact_triangle_boundary(np.asarray(v),np.asarray(f,dtype=np.int32),'AnyNonzero'))
        rows[name]=dict(result=result,seconds=time.monotonic()-start)
        (out/(name+'.json')).write_text(json.dumps(rows[name],default=lambda x:x.tolist()))
        print(name, result, flush=True)
        if not result['boundary_valid']:
            (out/'report.json').write_text(json.dumps(dict(status='DERIVED_BOUNDARY_REJECTED',part=name,qualified_delivery=False)))
            return
    (out/'report.json').write_text(json.dumps(dict(status='BOUNDARIES_VALID',qualified_delivery=False,pending=['bounded geometry and partition certification'])))

def launch():
    sys.path.insert(0,str(ROOT/'python'))
    from make_it_stand.verification.serial_source_job import run_serial_source_job
    build=HERE/'formal-native-core-build-job-12/output/build'
    files={'worker':Path(__file__),'native':build/'python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'}
    files.update({p.name:p for p in (HERE/'noodles-terminal-obj-01').iterdir()})
    r=run_serial_source_job(ROOT,HERE/'noodles-terminal-obj-audit-job-01',
        ['/usr/bin/env','LD_PRELOAD='+str(build/'compiler-runtime/libstdc++.so.6'),sys.executable,'-B',str(Path(__file__)),'{job}/output'],files,
        environment={'PYTHONPATH':str(build/'python')})
    print(r.status,r.returncode,flush=True)

if __name__=='__main__':
    if len(sys.argv)>1:worker(sys.argv[1])
    else:launch()
