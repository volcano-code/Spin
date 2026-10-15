from pathlib import Path
import sys
import argparse
from verify_timing_candidate import verify
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
p=argparse.ArgumentParser()
p.add_argument('--native',type=Path,default=HERE/'timing-build-job-03/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so')
p.add_argument('--comparison',default='same-input-timed-job-03')
p.add_argument('--job',default='heart-state-timing-job-02')
p.add_argument('--repetitions',type=int,choices=(1,2,3,4),default=1)
p.add_argument('--perturb-constant',action='store_true')
p.add_argument('--replay-record',type=Path)
args=p.parse_args()
verify(HERE/args.comparison/'output',HERE/'same-input-job-03/output')
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=args.native.resolve(strict=True)
worker=HERE/'time_saved_heart_state.py'
heart=ROOT/'experiments/faster_heart_contact_20261011'
request=heart/'rounded-clearance-v3-faster-frozen-v2/request.json'
import json
frozen=json.loads(request.read_text())
files=dict(worker=worker,launcher=Path(__file__),timing=HERE/'phase_timing.py',native=native,
    core_policy=heart/'core_faster_heart.py',journal=heart/'rounded_solver_journal.py',request=request,
    events=heart/'rounded-clearance-v3-core-faster-v3/seed-0-job-01/output/run/newton-events.jsonl',
    spin=ROOT/'python/make_it_stand/spin.py')
files.update({'input-'+str(i):Path(path) for i,path in enumerate(frozen['input_sources'])})
if args.replay_record is not None:files['replay_record']=args.replay_record.resolve(strict=True)
result=run_source_job(HERE/args.job,
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(worker),'--output','{job}/output','--repetitions',str(args.repetitions),
     *(['--perturb-constant'] if args.perturb_constant else []),
     *(['--replay-record',str(args.replay_record.resolve())] if args.replay_record is not None else [])],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,result.returncode,flush=True)
