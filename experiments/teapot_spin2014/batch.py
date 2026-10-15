"""Sequential 24-GiB/60-minute attempts; failure never disappears from counts."""
import argparse
from pathlib import Path
import os
import signal
import subprocess
import sys
import time
from teapot import STRATEGIES, write


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(exist_ok=False)
    rows=[]
    for name in STRATEGIES:
        print('BATCH START',name,flush=True);start=time.monotonic()
        cmd=['prlimit','--as=25769803776',sys.executable,str(Path(__file__).with_name('run.py')),
            '--input',str(a.input),'--output',str(a.output/name),'--strategy',name]
        with (a.output/(name+'.log')).open('x') as log:
            child=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            try:code=child.wait(timeout=3600);status='EXITED'
            except subprocess.TimeoutExpired:
                os.killpg(child.pid,signal.SIGTERM)
                try:child.wait(timeout=10)
                except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
                code=child.returncode;status='RESOURCE_LIMIT'
        rows.append(dict(strategy=name,command=cmd,returncode=code,status=status,elapsed_seconds=time.monotonic()-start))
        write(a.output/'batch.json',dict(attempts=rows,expected_strategies=list(STRATEGIES),memory_limit_bytes=25769803776))
        print('BATCH END',name,code,status,flush=True)
