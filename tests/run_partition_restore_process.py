"""Fresh-process domain reconstruction + next SLQP, not durable run resume."""
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

def main():
    if len(sys.argv)!=2: raise ValueError('native test executable required')
    executable=Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix='spin-partition-process-') as directory:
        state=Path(directory)/'state.txt'
        def run(mode,path):
            return subprocess.run([str(executable),mode,str(path)],capture_output=True,timeout=60)
        baseline=run('--produce',state)
        if baseline.returncode: raise AssertionError(baseline.stderr.decode())
        outcomes=[]
        for _ in range(2):
            restored=run('--consume',state)
            if restored.returncode: raise AssertionError(restored.stderr.decode())
            if baseline.stdout!=restored.stdout: raise AssertionError('cross-process scientific result/trace mismatch')
            outcomes.append(hashlib.sha256(restored.stdout).hexdigest())
        # This is a test fixture reader, not evidence of a production snapshot
        # schema. The production restore factory independently validates keys.
        original=state.read_bytes()
        truncated=Path(directory)/'truncated.txt';truncated.write_bytes(original[:-12])
        if run('--consume',truncated).returncode==0: raise AssertionError('truncated test state accepted')
        tampered=Path(directory)/'tampered.txt'
        altered=bytearray(original)
        lines=original.splitlines(keepends=True)
        for row_index, row in enumerate(lines[3:-1], start=3):
            beta=row.rstrip(b'\n').split()[-1]
            if beta.startswith(b'0.') and len(beta)>2 and 48<=beta[2]<=57:
                position=sum(map(len,lines[:row_index]))+row.rfind(beta)+2
                altered[position]=49 if altered[position]!=49 else 50
                break
        else:
            raise AssertionError('fixture has no fractional beta to tamper')
        tampered.write_bytes(altered)
        rejected=run('--consume',tampered)
        if rejected.returncode==0 or b'corrupt test-state framing' not in rejected.stderr:
            raise AssertionError('tampered beta was not rejected by test-state framing')
        print(json.dumps(dict(scope='SYNTHETIC_PARTITION_PLUS_ONE_REAL_SLQP_STAGE',
            fresh_process_repeats=2,scientific_sha256=outcomes,full_checkpoint_resume='NOT_IMPLEMENTED')))

if __name__=='__main__': main()
