"""Compare the two actual boundary audits on a frozen admitted request."""
import argparse
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'python'))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--candidate',type=Path,required=True)
    parser.add_argument('--job',type=Path);parser.add_argument('--output',type=Path)
    args=parser.parse_args();case=args.candidate.resolve()
    if args.output:
        import _make_it_stand as native
        import numpy as np
        args.output.mkdir(parents=True,exist_ok=False)
        for name in ('outer','design'):
            d=json.loads((case/(name+'.binary64.json')).read_text())
            v=np.asarray(d['vertices']);f=np.asarray(d['faces'],dtype=np.int32)
            legacy=native.audit_triangle_boundary(v,f)
            (args.output/(name+'-legacy.json')).write_text(json.dumps(legacy,indent=2)+'\n')
            exact=native.audit_exact_triangle_boundary(v,f,'OutwardPositive')
            (args.output/(name+'-exact.json')).write_text(json.dumps(exact,indent=2)+'\n')
            print(json.dumps(dict(name=name,legacy=legacy.get('strict_boundary_valid'),exact=exact['boundary_valid'])),flush=True)
        return
    from make_it_stand.verification.source_job import run_source_job
    native=ROOT/'build/spin2014-objectives-quadric-20261007/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    files=dict(worker=Path(__file__),native=native,outer=case/'outer.binary64.json',design=case/'design.binary64.json')
    result=run_source_job(args.job.resolve(),['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
        '/home/lin/anaconda3/bin/python3','-B',str(Path(__file__).resolve()),'--candidate',str(case),
        '--output','{job}/output'],files,environment={'PYTHONPATH':str(native.parent)})
    print(result.status,flush=True)


if __name__=='__main__':main()
