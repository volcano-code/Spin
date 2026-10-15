"""Measured POSIX attempt budgets, NOT a scientific acceptance gate.

Private local filesystem, immutable supplied identities, cooperative processes
that stay in their process group. Disk is sampled (not a hard quota); RLIMIT_AS
is per process (not aggregate RSS). An interrupted supervisor fails closed:
unfinished reservations cannot be automatically refunded or resumed. No shell,
automatic evidence deletion, manual passed flag, or paper-completion inference.
"""
from contextlib import contextmanager
from dataclasses import dataclass
import fcntl
import hashlib
import math
import os
from pathlib import Path
import shutil
import signal
import stat
import subprocess
import sys
import time
import uuid

from .source_checkpoint_store import _canonical, _fingerprints, _integer, _regular, _unique, _verify_artifacts
import json


def _sync(directory):
    fd=os.open(directory,os.O_RDONLY|os.O_DIRECTORY|os.O_CLOEXEC)
    try:os.fsync(fd)
    finally:os.close(fd)


def _write(path,data):
    """Durable no-replace publication; pending hardlinks retained as evidence."""
    body=_canonical(data)
    encoded=_canonical({'data':data,'sha256':hashlib.sha256(body).hexdigest()})
    temporary=path.parent/('.pending-'+uuid.uuid4().hex)
    with temporary.open('xb') as stream:
        stream.write(encoded);stream.flush();os.fsync(stream.fileno())
    os.link(temporary,path);_sync(path.parent)


def _read(path):
    with _regular(path) as stream:
        raw=stream.read(4*1024*1024+1)
    if len(raw)>4*1024*1024:raise ValueError('SUPERVISOR_RECORD_SIZE')
    record=json.loads(raw,object_pairs_hook=_unique)
    if (not isinstance(record,dict) or set(record)!={'data','sha256'} or
        _canonical(record)!=raw or hashlib.sha256(_canonical(record['data'])).hexdigest()!=record['sha256']):
        raise ValueError('SUPERVISOR_RECORD_CORRUPT')
    return record['data'],record['sha256']


def _disk_bytes(root):
    """Account inodes once; permit only links to regular files within root.

    Directory links are never traversed. Library aliases cost their own inode
    size; their regular target is accounted by the ordinary directory walk.
    """
    resolved_root=root.resolve(strict=True)
    info=root.stat();seen={(info.st_dev,info.st_ino)};total=max(info.st_size,info.st_blocks*512)
    for base,dirs,files in os.walk(root,followlinks=False):
        for name in dirs+files:
            path=Path(base)/name
            try:
                info=path.lstat()
            except FileNotFoundError:
                # A live build may rename/unlink an enumerated temporary file.
                # This is sampled accounting, not a filesystem snapshot. Never
                # suppress permission, I/O, or unsupported-file errors.
                continue
            if stat.S_ISLNK(info.st_mode):
                try:
                    target=path.resolve(strict=True)
                    valid=target.is_relative_to(resolved_root) and stat.S_ISREG(target.stat().st_mode)
                except (FileNotFoundError,RuntimeError):
                    valid=False
                except OSError as exc:
                    import errno
                    if exc.errno!=errno.ELOOP:raise
                    valid=False
                if not valid:raise ValueError('UNSUPPORTED_RUN_ARTIFACT')
            elif not (stat.S_ISREG(info.st_mode) or stat.S_ISDIR(info.st_mode)):
                raise ValueError('UNSUPPORTED_RUN_ARTIFACT')
            inode=(info.st_dev,info.st_ino)
            if inode not in seen:
                total+=max(info.st_size,info.st_blocks*512);seen.add(inode)
    return total


@dataclass(frozen=True)
class AttemptResult:
    index:int
    status:str
    returncode:int
    charged_wall_ns:int
    cumulative_wall_ns:int
    peak_disk_bytes:int
    stdout:Path
    stderr:Path


class StageRunSupervisor:
    def __init__(self,directory,identity_files,*,maximum_wall_ns,maximum_disk_bytes,
                 maximum_memory_bytes=24*1024**3,minimum_free_bytes=20*1024**3,
                 poll_seconds=.05):
        self.directory=Path(directory).resolve()
        self.limits=dict(maximum_wall_ns=_integer(maximum_wall_ns,1,2**63-1),
            maximum_disk_bytes=_integer(maximum_disk_bytes,1,2**63-1),
            maximum_memory_bytes=_integer(maximum_memory_bytes,1,2**63-1),
            minimum_free_bytes=_integer(minimum_free_bytes,0,2**63-1))
        if type(poll_seconds) not in (float,int) or not math.isfinite(poll_seconds) or not 0<poll_seconds<=1:
            raise ValueError('INVALID_POLL_INTERVAL')
        self.poll_seconds=float(poll_seconds)
        if not identity_files:raise ValueError('SUPERVISOR_IDENTITY_REQUIRED')
        self.identity_files=dict(identity_files)
        self.launcher=Path(__file__).with_name('_source_supervised_exec.py')
        for key,path in [('supervisor',Path(__file__)),('launcher',self.launcher),('python',Path(sys.executable))]:
            if '__'+key in self.identity_files:raise ValueError('RESERVED_IDENTITY_LABEL')
            self.identity_files['__'+key]=path
        self.identity=_fingerprints(self.identity_files)
        self.directory.mkdir(mode=0o700,exist_ok=True);_sync(self.directory.parent)

    @contextmanager
    def _lock(self):
        fd=os.open(self.directory/'.supervisor.lock',os.O_RDWR|os.O_CREAT|os.O_NOFOLLOW|os.O_CLOEXEC,0o600)
        try:
            if not stat.S_ISREG(os.fstat(fd).st_mode):raise ValueError('INVALID_SUPERVISOR_LOCK')
            try:fcntl.flock(fd,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError:raise ValueError('SUPERVISOR_BUSY') from None
            yield
        finally:os.close(fd)

    def _load(self):
        if _fingerprints(self.identity_files)!=self.identity:raise ValueError('STALE_SUPERVISOR_IDENTITY')
        config=dict(schema='source-run-supervisor-v1',limits=self.limits,
                    poll_seconds=self.poll_seconds,identity=self.identity)
        path=self.directory/'config.json'
        if not path.exists():
            if any(p.name!='.supervisor.lock' for p in self.directory.iterdir()):
                raise ValueError('MISSING_SUPERVISOR_CONFIG')
            _write(path,config)
        saved,parent=_read(path)
        if saved!=config:raise ValueError('SUPERVISOR_IDENTITY_OR_LIMITS_CHANGED')
        reservations=sorted(self.directory.glob('reservation-*'))
        wall=0;peak=0
        for index,path in enumerate(reservations):
            if path.name!=f'reservation-{index:08d}.json':raise ValueError('ATTEMPT_SEQUENCE')
            start,start_hash=_read(path)
            if start.get('index')!=index or start.get('parent')!=parent or start.get('prior_wall_ns')!=wall:
                raise ValueError('ATTEMPT_CHAIN')
            _integer(start.get('reserved_wall_ns'),1,self.limits['maximum_wall_ns']-wall)
            attempt=self.directory/f'attempt-{index:08d}'
            if not (attempt/'result.json').exists() or not (attempt/'complete.json').exists():
                raise ValueError('UNFINISHED_ATTEMPT_RESERVATION_NO_AUTOMATIC_REFUND')
            result,parent=_read(attempt/'result.json')
            receipt,_=_read(attempt/'complete.json')
            if receipt!={'result_sha256':parent}:raise ValueError('ATTEMPT_COMPLETION_RECEIPT')
            if result.get('reservation')!=start_hash or result.get('index')!=index:raise ValueError('ATTEMPT_RESULT_CHAIN')
            _verify_artifacts(result.get('artifacts'))
            recorded=set(result['artifacts'])|{'result.json','complete.json'}
            present={str(p.relative_to(attempt)) for p in attempt.rglob('*')
                     if p.is_file() and not p.name.startswith('.pending-')}
            if recorded!=present:raise ValueError('ATTEMPT_ARTIFACT_SET_CHANGED')
            wall+=_integer(result.get('charged_wall_ns'),1,2**63-1)
            if result.get('cumulative_wall_ns')!=wall:raise ValueError('ATTEMPT_WALL_COUNTER')
            peak=max(peak,_integer(result.get('peak_disk_bytes'),0,2**63-1))
        expected={f'attempt-{i:08d}' for i in range(len(reservations))}
        if {p.name for p in self.directory.glob('attempt-*')}!=expected:raise ValueError('ORPHAN_ATTEMPT')
        return dict(attempts=len(reservations),cumulative_wall_ns=wall,
                    peak_disk_bytes=max(peak,_disk_bytes(self.directory)),parent=parent)

    def inspect(self):
        with self._lock():return self._load()

    def _resource(self,used,peak):
        if used>=self.limits['maximum_wall_ns']:return 'RESOURCE_LIMIT_WALL'
        if peak>=self.limits['maximum_disk_bytes']:return 'RESOURCE_LIMIT_DISK'
        if shutil.disk_usage(self.directory).free<self.limits['minimum_free_bytes']:return 'RESOURCE_LIMIT_FREE_SPACE'
        return None

    def run(self,argv,*,attempt_wall_ns=None):
        """Run argv without a shell. Replace literal {attempt} in arguments.

        Capture binary stdout/stderr to the exclusive attempt directory. A zero
        exit status means only leader exit, never solver convergence or successful
        completion of its descendants (any remaining process group is killed).
        All outputs/checkpoints governed by this budget must stay under directory.
        """
        if (not isinstance(argv,(list,tuple)) or not argv or
            any(type(a) is not str or '\0' in a for a in argv) or not Path(argv[0]).is_absolute()):
            raise ValueError('ABSOLUTE_ARGV_REQUIRED')
        cap=self.limits['maximum_wall_ns'] if attempt_wall_ns is None else _integer(attempt_wall_ns,1,2**63-1)
        with self._lock():
            began=time.monotonic_ns();prior=self._load()
            peak=max(prior['peak_disk_bytes'],_disk_bytes(self.directory))
            reason=self._resource(prior['cumulative_wall_ns']+time.monotonic_ns()-began,peak)
            if reason:raise ValueError(reason)
            index=prior['attempts'];attempt=self.directory/f'attempt-{index:08d}'
            command=[a.replace('{attempt}',str(attempt)) for a in argv]
            reserved=min(cap,self.limits['maximum_wall_ns']-prior['cumulative_wall_ns'])
            start=dict(index=index,parent=prior['parent'],prior_wall_ns=prior['cumulative_wall_ns'],
                       reserved_wall_ns=reserved,argv=command)
            reservation=self.directory/f'reservation-{index:08d}.json'
            _write(reservation,start);_,start_hash=_read(reservation)
            attempt.mkdir(mode=0o700);_sync(self.directory)
            stdout=attempt/'stdout.bin';stderr=attempt/'stderr.bin'
            process=None;status='PROCESS_FAILED';code=-1;error=None
            env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
            try:
                with stdout.open('xb') as out,stderr.open('xb') as err:
                    process=subprocess.Popen([sys.executable,str(self.launcher),str(self.limits['maximum_memory_bytes']),*command],
                        cwd=attempt,env=env,stdout=out,stderr=err,stdin=subprocess.DEVNULL,
                        start_new_session=True,close_fds=True)
                    _write(attempt/'process.json',dict(pid=process.pid))
                    while True:
                        peak=max(peak,_disk_bytes(self.directory))
                        spent=time.monotonic_ns()-began
                        reason=self._resource(prior['cumulative_wall_ns']+spent,peak)
                        if spent>=reserved:reason='RESOURCE_LIMIT_WALL'
                        if reason:status=reason;break
                        # Observe without reaping: retain the child's PID until
                        # the owned process group has been signalled in finally.
                        exited=os.waitid(os.P_PID,process.pid,os.WEXITED|os.WNOHANG|os.WNOWAIT)
                        if exited is not None:
                            status='PROCESS_EXITED' if exited.si_code==os.CLD_EXITED and exited.si_status==0 else 'PROCESS_FAILED';break
                        time.sleep(min(self.poll_seconds,max(0,(reserved-spent)/1e9)))
            except BaseException as exc:
                error=exc;status='SUPERVISOR_ERROR'
            finally:
                if process is not None:
                    try:os.killpg(process.pid,signal.SIGKILL)
                    except ProcessLookupError:pass
                    code=process.wait()
            peak=max(peak,_disk_bytes(self.directory))
            artifacts=_fingerprints({str(p.relative_to(attempt)):p for p in attempt.rglob('*')
                if p.is_file() and not p.name.startswith('.pending-')})
            if _fingerprints(self.identity_files)!=self.identity:
                error=ValueError('STALE_SUPERVISOR_IDENTITY');status='SUPERVISOR_ERROR'
            charged=max(1,time.monotonic_ns()-began)
            reason=self._resource(prior['cumulative_wall_ns']+charged,peak)
            if reason:status=reason
            if charged>=reserved:status='RESOURCE_LIMIT_WALL'
            result=dict(index=index,reservation=start_hash,status=status,returncode=code,
                charged_wall_ns=charged,cumulative_wall_ns=prior['cumulative_wall_ns']+charged,
                peak_disk_bytes=peak,artifacts=artifacts)
            # Conservatively charge this result's own data and two directory
            # entries before publishing it. Also reserve one data block and two
            # directory entries for the final completion receipt. The 256 bytes
            # cover changed scalar widths/checksum/envelope.
            block=max(512,os.statvfs(attempt).f_frsize)
            publication_bytes=((len(_canonical(result))+256+block-1)//block+5)*block
            peak+=publication_bytes;result['peak_disk_bytes']=peak
            if peak>=self.limits['maximum_disk_bytes']:
                result['status']='RESOURCE_LIMIT_DISK'
            _write(attempt/'result.json',result)
            actual=_disk_bytes(self.directory)
            if actual+3*block>peak:
                raise ValueError('PUBLICATION_ALLOCATION_EXCEEDED_ACCOUNTING_BOUND')
            # This is the ONLY completion marker, written after required checks.
            # Death/error after result publication but before here stays blocked.
            _,digest=_read(attempt/'result.json')
            _write(attempt/'complete.json',dict(result_sha256=digest))
            if error is not None:raise error
            return AttemptResult(index,result['status'],code,charged,result['cumulative_wall_ns'],peak,stdout,stderr)
