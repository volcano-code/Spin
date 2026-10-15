"""Fresh native processes + durable carrier, SYNTHETIC only, no run supervisor."""
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import time

from make_it_stand.verification.source_checkpoint_store import StageCheckpointStore


def check(executable,root):
    root.mkdir(exist_ok=False)
    repo=Path(__file__).resolve().parents[1]
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    outcomes=[];calls=0;elapsed=0
    def run(mode,fixture,directory,snapshot=None):
        nonlocal calls,elapsed
        command=[str(executable),mode,str(fixture),str(directory)]
        if snapshot is not None:command.append(str(snapshot))
        started=time.monotonic_ns()
        result=subprocess.run(command,capture_output=True,timeout=90,env=env)
        elapsed+=time.monotonic_ns()-started
        (root/f'call-{calls:03d}.stderr').write_bytes(result.stderr)
        calls+=1
        return result
    for fixture in range(4):
        baseline=run('baseline',fixture,root/f'baseline-{fixture}')
        assert baseline.returncode==0,baseline.stderr
        assert baseline.stdout[:1]==b'F'
        (root/f'baseline-{fixture}.science').write_bytes(baseline.stdout)
        expected=hashlib.sha256(baseline.stdout).hexdigest()
        for repeat in range(2):
            work=root/f'fixture-{fixture}-repeat-{repeat}';work.mkdir()
            config=work/'config.json';config.write_text(json.dumps({'fixture':fixture,'protocol':'SYNTHETIC_COMPLETE_NATIVE_STAGE_CHAIN'}))
            identities={'binary':executable,'checker':Path(__file__).resolve(),'config':config,
                'codec':repo/'src/spin2014/SourceCheckpointCodec.cpp',
                'solver':repo/'src/spin2014/SourceAdaptiveSpinSolver.cpp',
                'fixture':repo/'tests/spin_it_2014_checkpoint_codec.cpp',
                'store':repo/'python/make_it_stand/verification/source_checkpoint_store.py'}
            store=StageCheckpointStore(work/'store',identities,maximum_wall_ns=10**12,
                maximum_trace_bytes=2**30,maximum_payload_bytes=64*1024*1024)
            result=run('start',fixture,work/'attempt-0');completed=0
            while True:
                assert result.returncode==0,result.stderr
                if result.stdout[:1]==b'F':break
                assert result.stdout[:1]==b'P'
                traces={str(p.relative_to(work)):p for p in sorted(work.glob('attempt-*/*.slqptrace'))}
                assert traces
                store.publish(completed,payload=result.stdout[1:],cumulative_wall_ns=elapsed,
                    cumulative_trace_bytes=sum(p.stat().st_size for p in traces.values()),artifacts=traces)
                loaded=store.load_latest();assert loaded.stage==completed
                snapshot=work/f'snapshot-{completed}.bin';snapshot.write_bytes(loaded.payload)
                roundtrip=run('roundtrip',fixture,work/'unused-roundtrip',snapshot)
                assert roundtrip.returncode==0,roundtrip.stderr
                assert not (work/'unused-roundtrip').exists()
                completed+=1
                assert completed<8,'global stage counter reset'
                result=run('resume',fixture,work/f'attempt-{completed}',snapshot)
            assert result.stdout==baseline.stdout,'complete resumed output/iteration bytes differ'
            assert completed>0,'fixture did not exercise process restart'
            (work/'final.science').write_bytes(result.stdout)
            outcomes.append({'fixture':fixture,'repeat':repeat,'saved_boundaries':completed,'scientific_sha256':expected})
    # Corrupt the metadata, not the original files; rehash cases must reach the
    # semantic checks, not merely trip the outer checksum.
    original=root/'fixture-0-repeat-0/snapshot-0.bin';data=original.read_bytes()
    position=8
    def field():
        nonlocal position
        n=struct.unpack_from('<Q',data,position)[0];start=position+8;position=start+n;return start,n
    field();budget_offset=position;position+=8;count_offset=position;position+=8
    path=field();field();field();field();solution=field();before=field();after=field()
    negatives=[]
    for kind in ['truncate','trailing','checksum','version','count','length','counter','budget','objective','phi-nan','dual-witness','partition','wrong-request']:
        modified=bytearray(data);rehash=True;fixture=0
        if kind=='truncate':modified=modified[:-1];rehash=False
        elif kind=='trailing':modified+=b'x';rehash=False
        elif kind=='checksum':modified[-1]^=1;rehash=False
        elif kind=='version':modified[0]^=1
        elif kind=='count':struct.pack_into('<Q',modified,count_offset,0)
        elif kind=='length':struct.pack_into('<Q',modified,8,2**64-1)
        elif kind=='counter':struct.pack_into('<Q',modified,len(modified)-72,0)
        elif kind=='budget':struct.pack_into('<Q',modified,budget_offset,2**30+1)
        elif kind=='objective':
            value=struct.unpack_from('<d',modified,solution[0]+8)[0]
            struct.pack_into('<d',modified,solution[0]+8,value+1.)
        elif kind=='phi-nan':struct.pack_into('<d',modified,solution[0]+81,float('nan'))
        elif kind=='dual-witness':
            # Solution record: scalar prefix 89 bytes, kkt.valid byte, then
            # status string and lambda vector (rows,cols,coefficients).
            k=solution[0]+90;status_n=struct.unpack_from('<Q',modified,k)[0];k+=8+status_n
            rows,cols=struct.unpack_from('<QQ',modified,k);k+=16
            assert rows>0 and cols==1 and abs(struct.unpack_from('<d',modified,k)[0])>1e-8
            struct.pack_into('<d',modified,k,0.)
        elif kind=='partition':modified[after[0]]^=1
        elif kind=='wrong-request':fixture=2
        if rehash:modified[-64:]=hashlib.sha256(modified[:-64]).hexdigest().encode()
        snapshot=root/f'bad-{kind}.bin';snapshot.write_bytes(modified)
        attempt=root/f'bad-{kind}-attempt';result=run('resume',fixture,attempt,snapshot)
        assert result.returncode!=0,(kind,'accepted')
        assert not attempt.exists(),(kind,'created output before validation')
        if kind=='objective':assert b'evaluation disagrees' in result.stderr,result.stderr
        if kind=='dual-witness':assert b'KKT witness' in result.stderr,result.stderr
        negatives.append({'case':kind,'returncode':result.returncode,'diagnostic':result.stderr.decode().strip()})
    history=Path(data[path[0]:path[0]+path[1]].decode());original_history=history.read_bytes()
    reopened=run('reopen',0,history);assert reopened.returncode==0,reopened.stderr
    assert history.read_bytes()==original_history,'reopen modified historical bytes'
    visited=run('visit',0,root/'single-pass-visit')
    assert visited.returncode==0,visited.stderr
    capped=run('reopen',0,history,str(history.stat().st_size-1))
    assert capped.returncode!=0 and b'byte budget' in capped.stderr
    alias=root/'history-symlink';alias.symlink_to(history)
    assert run('reopen',0,alias).returncode!=0
    truncated=root/'history-truncated';truncated.write_bytes(original_history[:-1])
    assert run('reopen',0,truncated).returncode!=0
    assert run('reopen',0,root/'missing-history').returncode!=0
    # Modify a trace objective and recompute its per-record checksum. Point a
    # copied metadata payload to that file; original evidence is not modified.
    damaged=bytearray(original_history);record_size=struct.unpack_from('<Q',damaged,24)[0]
    objective_offset=24+8+64+8
    struct.pack_into('<d',damaged,objective_offset,struct.unpack_from('<d',damaged,objective_offset)[0]+1.)
    damaged[32:96]=hashlib.sha256(damaged[96:96+record_size]).hexdigest().encode()
    altered=root/'history-rehashed';altered.write_bytes(damaged)
    newpath=str(altered).encode()
    body=data[:path[0]-8]+struct.pack('<Q',len(newpath))+newpath+data[path[0]+path[1]:-64]
    altered_metadata=root/'bad-rehashed-history.bin';altered_metadata.write_bytes(body+hashlib.sha256(body).hexdigest().encode())
    attempt=root/'bad-rehashed-history-attempt';result=run('resume',0,attempt,altered_metadata)
    assert result.returncode!=0 and b'evaluation disagrees' in result.stderr,result.stderr
    assert not attempt.exists()
    negatives.append({'case':'rehashed-history-objective','returncode':result.returncode,'diagnostic':result.stderr.decode().strip()})
    for kind in ['forge-summary','forge-terminal','forge-count']:
        result=run(kind,0,root/kind)
        assert result.returncode==0,(kind,result.stderr)
        negatives.append({'case':kind,'returncode':result.returncode,
            'diagnostic':'native checker verified rejection after recomputing all hashes'})
    report={'scope':'SYNTHETIC_FRESH_PROCESS_COMPLETE_STAGE_CHAIN','runs':outcomes,
        'negative_cases':negatives,'native_processes':calls,'observed_native_call_wall_ns':elapsed,
        'resource_supervisor':'NOT_IMPLEMENTED','original_paper_case':'NOT_CLAIMED',
        'read_only_reload_checks':5,
        'binary_sha256':hashlib.sha256(executable.read_bytes()).hexdigest()}
    (root/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,sort_keys=True))


if __name__=='__main__':
    if len(sys.argv) not in (2,3):raise ValueError('native executable [new output directory] required')
    binary=Path(sys.argv[1]).resolve(strict=True)
    if len(sys.argv)==3:check(binary,Path(sys.argv[2]).absolute())
    else:
        with tempfile.TemporaryDirectory(prefix='spin-native-codec-') as directory:
            check(binary,Path(directory)/'evidence')
