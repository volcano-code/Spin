"""Apply explicitly identified diagnostic candidate to verified saved operands."""
import json,hashlib,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
packet=HERE/'inertia-boolean65-capture-job-07/output/capture/boolean-evaluation65.json'
module=HERE/'cherchi-traversal-diagnostic-build-job-05/output/candidate/build/python'
if len(sys.argv)>1:
    import _cherchi_core_boolean as native
    p=json.loads(packet.read_text());sha=p.pop('operand_packet_sha256')
    assert hashlib.sha256(json.dumps(p,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()==sha
    start=time.monotonic();result=None;error=None
    try:result=native.boolean(*(p[k] for k in ['design_vertices','design_faces','source_vertices','source_faces','operation']))
    except Exception as ex:error=repr(ex)
    report=dict(operand_sha256=sha,candidate_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),seconds=time.monotonic()-start,error=error,result=result,production_admission=False)
    Path(sys.argv[1]).parent.mkdir(parents=True,exist_ok=True)
    Path(sys.argv[1]).write_text(json.dumps(report));print({k:v for k,v in report.items() if k!='result'},flush=True)
else:
    sys.path.insert(0,str(ROOT/'python'))
    from make_it_stand.verification.serial_source_job import run_serial_source_job
    runtime=HERE/'formal-native-core-build-job-13/output/build/compiler-runtime/libstdc++.so.6'
    r=run_serial_source_job(ROOT,HERE/'cherchi-traversal-diagnostic-run-job-05',['/usr/bin/env','LD_PRELOAD='+str(runtime),sys.executable,'-B',str(Path(__file__)),'{job}/output/report.json'],dict(worker=Path(__file__),packet=packet,module=next(module.glob('*.so'))),environment={'PYTHONPATH':str(module)})
    print(r.status,r.returncode,flush=True)
