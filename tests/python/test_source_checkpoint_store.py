"""Storage safety only: no optimizer success or full resume is inferred."""
import hashlib
import importlib
import json
import os
import struct
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.fixture
def api():
    try:
        return importlib.import_module('make_it_stand.verification.source_checkpoint_store')
    except ModuleNotFoundError as exc:
        if exc.name=='make_it_stand.verification.source_checkpoint_store':
            pytest.fail('durable stage checkpoint store is not implemented')
        raise


@pytest.fixture
def identity(tmp_path):
    files={}
    for name in ('source','module','input','config','checker'):
        path=tmp_path/(name+'.dat');path.write_bytes(name.encode());files[name]=path
    return files


def store(api,tmp_path,identity,**options):
    return api.StageCheckpointStore(tmp_path/'checkpoints',identity,
        maximum_wall_ns=100,maximum_trace_bytes=1000,maximum_payload_bytes=4096,**options)


def publish(s,stage=0,**kw):
    args=dict(payload=b'\x00beta\x80phi\x00',cumulative_wall_ns=10+stage,
              cumulative_trace_bytes=20+stage,artifacts={})
    args.update(kw)
    return s.publish(stage,**args)


def test_binary_payload_and_new_process_roundtrip(api,tmp_path,identity):
    s=store(api,tmp_path,identity);assert s.load_latest() is None
    first=publish(s);second=publish(s,1,payload=b'new\x00snapshot')
    assert first.path.read_bytes() and second.parent_sha256==first.sha256
    current=store(api,tmp_path,identity).load_latest()
    assert (current.stage,current.next_stage,current.payload)==(1,2,b'new\x00snapshot')
    program='''import json,sys
from make_it_stand.verification.source_checkpoint_store import StageCheckpointStore
s=StageCheckpointStore(sys.argv[1],json.loads(sys.argv[2]),maximum_wall_ns=100,maximum_trace_bytes=1000,maximum_payload_bytes=4096)
r=s.load_latest();print(r.stage,r.next_stage,r.sha256,r.payload.hex())
'''
    out=subprocess.check_output([sys.executable,'-c',program,str(s.directory),
        json.dumps({k:str(v) for k,v in identity.items()})],text=True)
    assert out.strip()==f'1 2 {second.sha256} {second.payload.hex()}'


@pytest.mark.parametrize('stage',[-1,True,1,2,100000000])
def test_initial_stage_must_be_zero(api,tmp_path,identity,stage):
    s=store(api,tmp_path,identity)
    with pytest.raises(ValueError):publish(s,stage)
    assert s.load_latest() is None


def test_no_duplicate_skipped_or_deleted_stages(api,tmp_path,identity):
    s=store(api,tmp_path,identity);a=publish(s);original=a.path.read_bytes()
    for stage in (0,2):
        with pytest.raises(ValueError):publish(s,stage)
    assert a.path.read_bytes()==original
    b=publish(s,1);a.path.unlink()
    with pytest.raises(ValueError):s.load_latest()
    assert b.path.exists()


@pytest.mark.parametrize('field,value',[
    ('cumulative_wall_ns',9),('cumulative_wall_ns',101),('cumulative_wall_ns',True),
    ('cumulative_wall_ns',float('nan')),('cumulative_trace_bytes',19),
    ('cumulative_trace_bytes',1001),('cumulative_trace_bytes',-1)])
def test_reported_counter_reset_or_limit_rejected(api,tmp_path,identity,field,value):
    s=store(api,tmp_path,identity);a=publish(s)
    with pytest.raises(ValueError):publish(s,1,**{field:value})
    assert s.load_latest().sha256==a.sha256


@pytest.mark.parametrize('kind',['truncate','trailing','checksum','length','duplicate-key'])
def test_corrupt_latest_is_not_silently_rolled_back(api,tmp_path,identity,kind):
    s=store(api,tmp_path,identity);old=publish(s);latest=publish(s,1)
    data=bytearray(latest.path.read_bytes())
    if kind=='truncate':data=data[:-7]
    elif kind=='trailing':data+=b'junk'
    elif kind=='checksum':data[-1]^=1
    elif kind=='length':data[8:16]=b'\xff'*8
    else:
        magic,hn,pn=struct.unpack('<8sQQ',data[:24]);h=bytes(data[24:24+hn])
        h=h[:-1]+b',"stage":1}'
        data=bytearray(struct.pack('<8sQQ',magic,len(h),pn)+h+bytes(data[24+hn:-64]))
        data+=hashlib.sha256(data).hexdigest().encode()
    latest.path.write_bytes(data)
    with pytest.raises(ValueError):s.load_latest()
    assert old.path.exists()


@pytest.mark.parametrize('name',['source','module','input','config','checker'])
def test_effective_file_change_fails_for_same_and_fresh_store(api,tmp_path,identity,name):
    s=store(api,tmp_path,identity);publish(s)
    identity[name].write_bytes(b'changed')
    with pytest.raises(ValueError):s.load_latest()
    with pytest.raises(ValueError):store(api,tmp_path,identity).load_latest()


def test_external_artifacts_bound_and_checked(api,tmp_path,identity):
    artifact=tmp_path/'complete-trace';artifact.write_bytes(b'recorded bytes')
    s=store(api,tmp_path,identity);publish(s,artifacts={'pass-0':artifact})
    artifact.write_bytes(b'corrupted')
    with pytest.raises(ValueError):s.load_latest()


@pytest.mark.parametrize('payload',[b'',b'x'*4097,'text-not-bytes'])
def test_payload_bound(api,tmp_path,identity,payload):
    s=store(api,tmp_path,identity)
    with pytest.raises(ValueError):publish(s,payload=payload)
    assert s.load_latest() is None


@pytest.mark.parametrize('failure',['write','file-sync','directory-sync','link'])
def test_publication_failure_preserves_previous(api,tmp_path,identity,monkeypatch,failure):
    s=store(api,tmp_path,identity);first=publish(s);old=first.path.read_bytes()
    real_fsync=api.os.fsync
    def fail(*args,**kwargs):raise OSError('injected storage fault')
    with monkeypatch.context() as m:
        if failure=='write':m.setattr(api.os,'write',fail)
        elif failure=='link':m.setattr(api.os,'link',fail)
        else:
            import stat
            def sync(fd):
                directory=stat.S_ISDIR(os.fstat(fd).st_mode)
                if directory==(failure=='directory-sync'):fail()
                return real_fsync(fd)
            m.setattr(api.os,'fsync',sync)
        with pytest.raises(OSError):publish(s,1)
    assert first.path.read_bytes()==old
    assert s.load_latest().sha256==first.sha256
    assert list(s.directory.glob('.pending-*'))  # interrupted attempt evidence retained


def test_undeclared_limits_change_rejected(api,tmp_path,identity):
    s=store(api,tmp_path,identity);publish(s)
    changed=api.StageCheckpointStore(s.directory,identity,maximum_wall_ns=101,
        maximum_trace_bytes=1000,maximum_payload_bytes=4096)
    with pytest.raises(ValueError):changed.load_latest()


def test_abandoned_partial_file_not_a_published_snapshot(api,tmp_path,identity):
    s=store(api,tmp_path,identity);a=publish(s)
    (s.directory/'.pending-interrupted').write_bytes(b'partial')
    assert s.load_latest().sha256==a.sha256


def test_symlink_stage_not_followed(api,tmp_path,identity):
    s=store(api,tmp_path,identity);a=publish(s)
    data=a.path.read_bytes();outside=tmp_path/'outside';outside.write_bytes(data)
    a.path.unlink();a.path.symlink_to(outside)
    with pytest.raises((OSError,ValueError)):s.load_latest()


@pytest.mark.parametrize('which',['payload','receipt'])
def test_single_latest_file_deletion_cannot_rewind(api,tmp_path,identity,which):
    s=store(api,tmp_path,identity);first=publish(s);second=publish(s,1)
    target=second.path if which=='payload' else s.directory/'commit-00000001.sha256'
    assert target.exists()
    target.unlink()
    with pytest.raises(ValueError):s.load_latest()
    assert first.path.exists()


def test_identity_symlink_retarget_is_detected(api,tmp_path,identity):
    original=identity['module'];alias=tmp_path/'module-link';alias.symlink_to(original)
    identity['module']=alias;s=store(api,tmp_path,identity);publish(s)
    new=tmp_path/'other-module';new.write_bytes(original.read_bytes())
    alias.unlink();alias.symlink_to(new)
    with pytest.raises(ValueError):s.load_latest()


_CHILD='''import json,os,sys
from pathlib import Path
from make_it_stand.verification.source_checkpoint_store import StageCheckpointStore
s=StageCheckpointStore(sys.argv[1],json.loads(sys.argv[2]),maximum_wall_ns=100,maximum_trace_bytes=1000,maximum_payload_bytes=4096)
mode=sys.argv[3]
link=os.link
def interrupted(src,dst):
    name=Path(dst).name
    if mode=='before-stage' and name.startswith('stage-'):os._exit(73)
    link(src,dst)
    if mode=='after-stage' and name.startswith('stage-'):os._exit(73)
    if mode=='after-receipt' and name.startswith('commit-'):os._exit(73)
os.link=interrupted
if mode=='concurrent':
    print('ready',flush=True);sys.stdin.readline()
try:
    s.publish(1,payload=b'child',cumulative_wall_ns=11,cumulative_trace_bytes=21,artifacts={})
except ValueError as exc:
    if str(exc)!='CHECKPOINT_STAGE_SEQUENCE':raise
    sys.exit(42)
'''


def child_args(s,identity,mode):
    return [sys.executable,'-c',_CHILD,str(s.directory),
            json.dumps({k:str(v) for k,v in identity.items()}),mode]


@pytest.mark.parametrize('mode',['before-stage','after-stage','after-receipt'])
def test_abrupt_exit_is_not_exception_cleanup(api,tmp_path,identity,mode):
    s=store(api,tmp_path,identity);first=publish(s);original=first.path.read_bytes()
    result=subprocess.run(child_args(s,identity,mode),capture_output=True,timeout=10)
    assert result.returncode==73,result.stderr
    assert first.path.read_bytes()==original
    if mode=='after-stage':
        with pytest.raises(ValueError,match='MISSING_STAGE_OR_COMMIT_RECEIPT'):s.load_latest()
    elif mode=='before-stage':assert s.load_latest().sha256==first.sha256
    else:
        # A complete pair may be visible without the writer acknowledging it.
        # This process-crash test does not simulate power-loss durability.
        assert s.load_latest().payload==b'child'


def test_concurrent_writers_cannot_replace_same_stage(api,tmp_path,identity):
    s=store(api,tmp_path,identity);publish(s)
    children=[subprocess.Popen(child_args(s,identity,'concurrent'),stdin=subprocess.PIPE,
              stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True) for _ in range(2)]
    try:
        for p in children:assert p.stdout.readline().strip()=='ready'
        for p in children:p.stdin.write('go\n');p.stdin.flush()
        assert sorted(p.wait(timeout=10) for p in children)==[0,42]
        assert s.load_latest().stage==1
        assert len(list(s.directory.glob('stage-*')))==2
    finally:
        for p in children:
            if p.poll() is None:p.kill()
            p.communicate(timeout=10)


@pytest.mark.parametrize('field,value',[
    ('parent_sha256','0'*64),('stage',0),('next_stage',3),
    ('cumulative_wall_ns',9),('cumulative_trace_bytes',True)])
def test_rehashed_invalid_metadata_is_rejected(api,tmp_path,identity,field,value):
    s=store(api,tmp_path,identity);publish(s);latest=publish(s,1)
    data=latest.path.read_bytes();magic,hn,pn=struct.unpack('<8sQQ',data[:24])
    h=json.loads(data[24:24+hn]);h[field]=value
    header=json.dumps(h,sort_keys=True,separators=(',',':')).encode()
    body=struct.pack('<8sQQ',magic,len(header),pn)+header+data[24+hn:-64]
    digest=hashlib.sha256(body).hexdigest()
    latest.path.write_bytes(body+digest.encode())
    (s.directory/'commit-00000001.sha256').write_text(f'00000001 {digest}\n')
    with pytest.raises(ValueError):s.load_latest()


def test_short_writes_are_completed(api,tmp_path,identity,monkeypatch):
    s=store(api,tmp_path,identity);write=api.os.write
    monkeypatch.setattr(api.os,'write',lambda fd,data:write(fd,data[:7]))
    current=publish(s)
    assert s.load_latest().sha256==current.sha256


def test_zero_write_fails_without_publication(api,tmp_path,identity,monkeypatch):
    s=store(api,tmp_path,identity)
    monkeypatch.setattr(api.os,'write',lambda fd,data:0)
    with pytest.raises(OSError,match='short write'):publish(s)
    assert s.load_latest() is None
