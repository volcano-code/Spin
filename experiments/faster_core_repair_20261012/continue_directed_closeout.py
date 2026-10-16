"""Continue explicit approved serial diagnostics after bundle reaudit exits."""
import json
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def main():
    pending=HERE/'noodles-terminal-bundle-reaudit-job-01/result.json'
    while not pending.exists():time.sleep(5)
    result=json.loads(pending.read_text())
    print('BUNDLE_TERMINAL',result['status'],flush=True)
    # Each launcher independently verifies lock, unit and cgroup termination.
    for script in ('launch_inertia_boolean65_capture.py','launch_boolean65_stack.py'):
        code=subprocess.call([sys.executable,'-B',str(HERE/script)],cwd=ROOT)
        if code:raise SystemExit(code)
    print('DIRECTED_DIAGNOSTICS_TERMINATED',flush=True)

if __name__=='__main__':main()
