"""Opaque durable carrier + actual next-SLQP fixture; NOT full run resume."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from make_it_stand.verification.source_checkpoint_store import StageCheckpointStore


def open_store(directory, identity):
    return StageCheckpointStore(directory, identity, maximum_wall_ns=10**12,
        maximum_trace_bytes=10**8, maximum_payload_bytes=16*1024*1024)


def consume(directory):
    identity=json.loads((directory/'identity.json').read_text())
    state=open_store(directory/'store',identity).load_latest()
    if state.stage!=0:raise AssertionError('unexpected stage')
    # The state decoder below is test-only, not a production state schema.
    with tempfile.TemporaryDirectory(prefix='spin-loaded-state-') as temporary:
        path=Path(temporary)/'state.txt';path.write_bytes(state.payload)
        result=subprocess.run([identity['native'],'--consume',str(path)],
            capture_output=True,timeout=60,check=True)
    sys.stdout.buffer.write(result.stdout)


def verify(executable,directory):
    directory.mkdir(exist_ok=False)
    root=Path(__file__).resolve().parents[1]
    identity={
        'native':str(executable),
        'checker':str(Path(__file__).resolve()),
        'store':str(root/'python/make_it_stand/verification/source_checkpoint_store.py'),
        'native_fixture':str(root/'tests/spin_it_2014_partition_restore_process.cpp'),
        'domain_source':str(root/'src/spin2014/SourceAdaptiveDomain.cpp'),
        'domain_header':str(root/'include/spin2014/SourceAdaptiveDomain.h'),
    }
    (directory/'identity.json').write_text(json.dumps(identity,sort_keys=True))
    baseline=subprocess.run([str(executable),'--produce',str(directory/'state.txt')],
        capture_output=True,timeout=60,check=True)
    (directory/'baseline.stdout').write_bytes(baseline.stdout)
    payload=(directory/'state.txt').read_bytes()
    published=open_store(directory/'store',identity).publish(0,payload=payload,
        cumulative_wall_ns=0,cumulative_trace_bytes=len(baseline.stdout),
        artifacts={'baseline':directory/'baseline.stdout','input_state':directory/'state.txt',
                   'identity_map':directory/'identity.json'})
    hashes=[]
    for repeat in range(2):
        result=subprocess.run([sys.executable,str(Path(__file__).resolve()),
            '--consume-store',str(directory)],capture_output=True,timeout=65,check=True)
        (directory/f'restored-{repeat}.stdout').write_bytes(result.stdout)
        if result.stdout!=baseline.stdout:raise AssertionError('next SLQP/trace mismatch')
        hashes.append(hashlib.sha256(result.stdout).hexdigest())
    report=dict(scope='SYNTHETIC_DURABLE_CARRIER_PLUS_ONE_REAL_SLQP_STAGE',
        fresh_process_repeats=2,scientific_sha256=hashes,
        stored_frame_sha256=published.sha256,
        native_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
        full_checkpoint_resume='NOT_IMPLEMENTED',
        runtime_accounting='NOT_TESTED_COUNTERS_ARE_CALLER_REPORTED',
        native_binary_provenance='FINGERPRINTED_NOT_REBUILT_BY_THIS_CHECKER')
    (directory/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,sort_keys=True))


if __name__=='__main__':
    if len(sys.argv)==3 and sys.argv[1]=='--consume-store':
        consume(Path(sys.argv[2]).resolve())
    elif len(sys.argv) in (2,3):
        executable=Path(sys.argv[1]).resolve(strict=True)
        if len(sys.argv)==3:verify(executable,Path(sys.argv[2]).absolute())
        else:
            with tempfile.TemporaryDirectory(prefix='spin-carrier-check-') as temporary:
                verify(executable,Path(temporary)/'evidence')
    else:raise ValueError('native executable [new evidence directory] required')
