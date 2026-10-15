"""POSIX durable, append-only storage for source-adaptive stage payloads.

This is an integrity carrier, NOT a solver-state decoder, optimizer acceptance
gate, authenticated log or full resume scheduler. Payload bytes remain opaque.
Reported cumulative counters are checked for monotonicity/limits, not measured
here. An experiment supervisor must account for interrupted attempts, enforce
deadlines, and bind actual native stage state before enabling experiment resume.
Use a private directory on a local filesystem with working flock/link/fsync.
External identity/artifact files must be immutable while an operation runs.
"""
from contextlib import contextmanager
from dataclasses import dataclass
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import struct
import uuid

_FRAME=struct.Struct('<8sQQ')
_MAGIC=b'SPCKPT01'
_HEADER_LIMIT=1024*1024
_NAME=re.compile(r'stage-([0-9]{8})\.ckpt\Z')
_RECEIPT=re.compile(r'commit-([0-9]{8})\.sha256\Z')
_FIELDS={'schema','stage','next_stage','parent_sha256','identity','limits',
         'cumulative_wall_ns','cumulative_trace_bytes','artifacts'}


def _integer(value,minimum,maximum):
    if type(value) is not int or not minimum<=value<=maximum:
        raise ValueError('INVALID_CHECKPOINT_INTEGER')
    return value


def _canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode('utf-8')


def _unique(pairs):
    result={}
    for key,value in pairs:
        if key in result:raise ValueError('DUPLICATE_CHECKPOINT_KEY')
        result[key]=value
    return result


@contextmanager
def _regular(path):
    fd=os.open(path,os.O_RDONLY|os.O_NOFOLLOW|os.O_NONBLOCK|os.O_CLOEXEC)
    try:
        if not stat.S_ISREG(os.fstat(fd).st_mode):raise ValueError('REGULAR_FILE_REQUIRED')
        with os.fdopen(fd,'rb',closefd=False) as stream:yield stream
    finally:os.close(fd)


def _fingerprints(files):
    if not isinstance(files,dict):raise ValueError('CHECKPOINT_FILE_MAP_REQUIRED')
    result={}
    for label,path in files.items():
        if type(label) is not str or not label:raise ValueError('CHECKPOINT_FILE_LABEL_REQUIRED')
        anchored=Path(path).resolve(strict=True)
        with _regular(anchored) as stream:
            before=os.fstat(stream.fileno());digest=hashlib.file_digest(stream,'sha256').hexdigest()
            after=os.fstat(stream.fileno())
        if (before.st_size,before.st_mtime_ns,before.st_ctime_ns)!=(after.st_size,after.st_mtime_ns,after.st_ctime_ns):
            raise ValueError('FILE_CHANGED_DURING_CHECKPOINT_HASH')
        result[label]=dict(path=str(anchored),size=before.st_size,sha256=digest)
    return result


def _verify_artifacts(recorded):
    if not isinstance(recorded,dict):raise ValueError('INVALID_ARTIFACT_MAP')
    files={}
    for label,value in recorded.items():
        if not isinstance(value,dict) or set(value)!={'path','size','sha256'}:
            raise ValueError('INVALID_ARTIFACT_RECORD')
        if type(value['path']) is not str or not Path(value['path']).is_absolute():
            raise ValueError('INVALID_ARTIFACT_PATH')
        files[label]=value['path']
    if _canonical(_fingerprints(files))!=_canonical(recorded):
        raise ValueError('STALE_CHECKPOINT_ARTIFACT')


@dataclass(frozen=True)
class StoredStage:
    stage:int
    next_stage:int
    payload:bytes
    sha256:str
    parent_sha256:str|None
    cumulative_wall_ns:int
    cumulative_trace_bytes:int
    path:Path


class StageCheckpointStore:
    def __init__(self,directory,identity_files,*,maximum_wall_ns,
                 maximum_trace_bytes,maximum_payload_bytes=16*1024*1024):
        self.directory=Path(directory).resolve()
        self._limits=dict(maximum_wall_ns=_integer(maximum_wall_ns,1,2**63-1),
            maximum_trace_bytes=_integer(maximum_trace_bytes,1,2**63-1),
            maximum_payload_bytes=_integer(maximum_payload_bytes,1,64*1024*1024))
        self._identity=_fingerprints(identity_files)
        if not self._identity:raise ValueError('CHECKPOINT_IDENTITY_REQUIRED')
        # Anchor relative names but re-resolve symlinks on every observation;
        # otherwise retargeting the running module/input alias goes unnoticed.
        self._identity_files={k:str(Path(v).absolute()) for k,v in identity_files.items()}
        self.directory.mkdir(mode=0o700,exist_ok=True)
        parent=os.open(self.directory.parent,os.O_RDONLY|os.O_DIRECTORY|os.O_CLOEXEC)
        try:os.fsync(parent)
        finally:os.close(parent)

    @contextmanager
    def _lock(self,exclusive):
        fd=os.open(self.directory/'.checkpoint.lock',os.O_RDWR|os.O_CREAT|os.O_NOFOLLOW|os.O_CLOEXEC,0o600)
        try:
            if not stat.S_ISREG(os.fstat(fd).st_mode):raise ValueError('INVALID_CHECKPOINT_LOCK')
            fcntl.flock(fd,fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH)
            yield
        finally:os.close(fd)

    def _identity_current(self):
        if _canonical(_fingerprints(self._identity_files))!=_canonical(self._identity):
            raise ValueError('STALE_CHECKPOINT_IDENTITY')

    def _validate_header(self,h,index,previous):
        if not isinstance(h,dict) or set(h)!=_FIELDS or h['schema']!='source-stage-storage-v1':
            raise ValueError('INVALID_CHECKPOINT_HEADER')
        _integer(h['stage'],0,99999999);_integer(h['next_stage'],1,100000000)
        if h['stage']!=index or h['next_stage']!=index+1:
            raise ValueError('CHECKPOINT_STAGE_SEQUENCE')
        if h['parent_sha256']!=(previous.sha256 if previous else None):
            raise ValueError('CHECKPOINT_PARENT_MISMATCH')
        if _canonical(h['identity'])!=_canonical(self._identity):
            raise ValueError('STALE_CHECKPOINT_IDENTITY')
        if _canonical(h['limits'])!=_canonical(self._limits):
            raise ValueError('CHECKPOINT_LIMITS_CHANGED')
        for field,limit,prior in (
            ('cumulative_wall_ns','maximum_wall_ns',previous.cumulative_wall_ns if previous else 0),
            ('cumulative_trace_bytes','maximum_trace_bytes',previous.cumulative_trace_bytes if previous else 0)):
            _integer(h[field],prior,self._limits[limit])
        _verify_artifacts(h['artifacts'])

    def _read(self,path,index,previous):
        with _regular(path) as stream:
            size=os.fstat(stream.fileno()).st_size
            if size<_FRAME.size+65 or size>_FRAME.size+_HEADER_LIMIT+self._limits['maximum_payload_bytes']+64:
                raise ValueError('CHECKPOINT_FILE_SIZE')
            prefix=stream.read(_FRAME.size)
            if len(prefix)!=_FRAME.size:raise ValueError('TRUNCATED_CHECKPOINT')
            magic,hn,pn=_FRAME.unpack(prefix)
            if magic!=_MAGIC or not 0<hn<=_HEADER_LIMIT or not 0<pn<=self._limits['maximum_payload_bytes']:
                raise ValueError('CHECKPOINT_LENGTH_OR_VERSION')
            if size!=_FRAME.size+hn+pn+64:raise ValueError('CHECKPOINT_LENGTH_MISMATCH')
            header=stream.read(hn);payload=stream.read(pn);digest=stream.read(64)
            computed=hashlib.sha256(prefix+header+payload).hexdigest()
            if len(header)!=hn or len(payload)!=pn or digest!=computed.encode() or stream.read(1):
                raise ValueError('CHECKPOINT_CHECKSUM_OR_EOF')
        h=json.loads(header,object_pairs_hook=_unique)
        if _canonical(h)!=header:raise ValueError('NONCANONICAL_CHECKPOINT_HEADER')
        self._validate_header(h,index,previous)
        return StoredStage(index,index+1,payload,computed,h['parent_sha256'],
            h['cumulative_wall_ns'],h['cumulative_trace_bytes'],path)

    def _latest_locked(self):
        self._identity_current();paths=[];receipts={}
        for path in self.directory.iterdir():
            if path.name.startswith('stage-'):
                match=_NAME.fullmatch(path.name)
                if not match:raise ValueError('INVALID_PUBLISHED_STAGE_NAME')
                paths.append((int(match[1]),path))
            elif path.name.startswith('commit-'):
                match=_RECEIPT.fullmatch(path.name)
                if not match:raise ValueError('INVALID_COMMIT_RECEIPT_NAME')
                receipts[int(match[1])]=path
        if {i for i,_ in paths}!=set(receipts):
            raise ValueError('MISSING_STAGE_OR_COMMIT_RECEIPT')
        previous=None
        for expected,(index,path) in enumerate(sorted(paths)):
            if index!=expected:raise ValueError('MISSING_CHECKPOINT_STAGE')
            previous=self._read(path,index,previous)
            expected_receipt=f'{index:08d} {previous.sha256}\n'.encode()
            with _regular(receipts[index]) as stream:
                if stream.read(len(expected_receipt)+1)!=expected_receipt:
                    raise ValueError('CHECKPOINT_COMMIT_RECEIPT_MISMATCH')
        return previous

    def _write_pending(self,data):
        temporary=self.directory/('.pending-'+uuid.uuid4().hex)
        fd=os.open(temporary,os.O_WRONLY|os.O_CREAT|os.O_EXCL|os.O_CLOEXEC,0o600)
        try:
            offset=0
            while offset<len(data):
                written=os.write(fd,memoryview(data)[offset:])
                if written<=0:raise OSError('checkpoint short write')
                offset+=written
            os.fsync(fd)
        finally:os.close(fd)
        return temporary

    def load_latest(self):
        """Validate the entire published chain. Never fall back past corruption."""
        with self._lock(False):return self._latest_locked()

    def publish(self,stage,*,payload,cumulative_wall_ns,cumulative_trace_bytes,artifacts):
        """Persist opaque bytes. Caller must supply/verify true native stage state.

        Old versions never change. Failed attempt files are retained as .pending-*
        and are never treated as published. Cumulative counters are caller-reported;
        this function cannot measure or charge a killed solver process.
        """
        _integer(stage,0,99999999)
        if type(payload) is not bytes or not 0<len(payload)<=self._limits['maximum_payload_bytes']:
            raise ValueError('INVALID_CHECKPOINT_PAYLOAD')
        with self._lock(True):
            previous=self._latest_locked()
            expected=previous.next_stage if previous else 0
            if stage!=expected:raise ValueError('CHECKPOINT_STAGE_SEQUENCE')
            h=dict(schema='source-stage-storage-v1',stage=stage,next_stage=stage+1,
                parent_sha256=previous.sha256 if previous else None,identity=self._identity,
                limits=self._limits,cumulative_wall_ns=cumulative_wall_ns,
                cumulative_trace_bytes=cumulative_trace_bytes,artifacts=_fingerprints(artifacts))
            self._validate_header(h,stage,previous)
            header=_canonical(h)
            if len(header)>_HEADER_LIMIT:raise ValueError('CHECKPOINT_HEADER_LIMIT')
            body=_FRAME.pack(_MAGIC,len(header),len(payload))+header+payload
            digest=hashlib.sha256(body).hexdigest();data=body+digest.encode()
            target=self.directory/f'stage-{stage:08d}.ckpt'
            receipt=self.directory/f'commit-{stage:08d}.sha256'
            directory_fd=os.open(self.directory,os.O_RDONLY|os.O_DIRECTORY|os.O_CLOEXEC)
            linked_stage=False;linked_receipt=False
            try:
                temporary=self._write_pending(data)
                self._identity_current();_verify_artifacts(h['artifacts'])
                # The durable stage precedes its no-replace commit receipt.
                # Readers hold the same lock and require both. An abrupt death
                # between links leaves an explicit orphan, never a silent rewind.
                os.link(temporary,target);linked_stage=True
                os.fsync(directory_fd)
                receipt_temporary=self._write_pending(f'{stage:08d} {digest}\n'.encode())
                os.link(receipt_temporary,receipt);linked_receipt=True
                os.fsync(directory_fd)
            except BaseException:
                if linked_stage or linked_receipt:
                    # Revoke only this unacknowledged publication, not old stages;
                    # keep its pending file for inspection. Failure here is fatal.
                    if linked_receipt:os.unlink(receipt)
                    if linked_stage:os.unlink(target)
                    os.fsync(directory_fd)
                raise
            finally:os.close(directory_fd)
            # Pending hardlink deliberately retained: no evidence deletion and no
            # post-commit cleanup failure that could disguise successful publication.
            return StoredStage(stage,stage+1,payload,digest,h['parent_sha256'],
                cumulative_wall_ns,cumulative_trace_bytes,target)
