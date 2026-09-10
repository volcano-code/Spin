"""Preserve one observed parent CTest regular-file log unlinked by nested CTest.

Read-only with respect to parent/filesystem inputs; writes an exclusive evidence
copy. Holds a regular-file descriptor, never reads a parent's stdout/stderr pipe.
"""
import hashlib
import json
import os
from pathlib import Path
import stat
import sys
import time

pid = int(sys.argv[1])
descriptor = Path(f'/proc/{pid}/fd/3')
output = Path(sys.argv[2])
assert output.parent == Path('build/task-5g-evidence-vL6KVn')
assert not output.exists()
command = Path(f'/proc/{pid}/cmdline').read_bytes().replace(b'\0', b' ').decode()
assert 'ctest --test-dir build/strict-v2-spin --output-on-failure -j1' in command
target = os.readlink(descriptor)
assert target.endswith('/build/strict-v2-spin/Testing/Temporary/LastTest.log.tmp (deleted)')
with descriptor.open('rb') as held:
    assert stat.S_ISREG(os.fstat(held.fileno()).st_mode)
    print('Holding actual parent regular-file log:',target,flush=True)
    start = time.monotonic()
    while Path(f'/proc/{pid}').exists() and time.monotonic()-start < 3600:
        time.sleep(.5)
    parent_still_exists = Path(f'/proc/{pid}').exists()
    payload = held.read()
with output.open('xb') as sink:
    sink.write(payload)
metadata = {'observed_parent_command':command,'observed_fd_target':target,
            'parent_still_exists_at_capture':parent_still_exists,
            'bytes':len(payload),'sha256':hashlib.sha256(payload).hexdigest()}
with output.with_suffix('.json').open('x') as sink:
    json.dump(metadata,sink,sort_keys=True,indent=2)
print(json.dumps(metadata,sort_keys=True),flush=True)
if parent_still_exists:
    raise SystemExit('capture timed out; retained partial log, not completion evidence')
