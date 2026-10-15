"""One-shot whole-workflow user-systemd envelope, Linux cgroup v2 only.

Enforces aggregate job memory, swap=0 and independent service timeout. No sudo,
persistent units, automatic weaker fallback, solver acceptance or restart with
a fresh budget. Disk monitoring is sampled, NOT a filesystem quota. The same
user can deliberately escape/change its cgroup: this is not hostile containment.
"""
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import uuid

from .source_checkpoint_store import _canonical, _fingerprints, _integer
from .source_run_supervisor import _disk_bytes


def _publish(path,data):
    temporary=path.parent/('.pending-'+uuid.uuid4().hex)
    with temporary.open('xb') as out:
        out.write(_canonical(data));out.flush();os.fsync(out.fileno())
    os.link(temporary,path)
    fd=os.open(path.parent,os.O_RDONLY|os.O_DIRECTORY|os.O_CLOEXEC)
    try:os.fsync(fd)
    finally:os.close(fd)


def _show(unit):
    command=['/usr/bin/systemctl','--user','show',unit,
        '--property=LoadState,ActiveState,SubState,Result,ExecMainCode,ExecMainStatus,MainPID,ControlGroup,MemoryPeak,MemorySwapPeak,MemoryMax,MemorySwapMax,RuntimeMaxUSec']
    result=subprocess.run(command,capture_output=True,text=True,timeout=3)
    if result.returncode:raise RuntimeError('JOB_UNIT_OBSERVATION_FAILED: '+result.stderr)
    return dict(line.split('=',1) for line in result.stdout.splitlines() if '=' in line)


@dataclass(frozen=True)
class SourceJobResult:
    directory:Path
    unit:str
    status:str
    returncode:int
    elapsed_ns:int
    service:dict


def run_source_job(directory,argv,identity_files,*,maximum_wall_ns=3600_000_000_000,
                   maximum_memory_bytes=24*1024**3,maximum_disk_bytes=80*1024**3,
                   minimum_start_free_bytes=100*1024**3,minimum_free_bytes=20*1024**3,
                   environment=None):
    """Run a complete workflow once. {job} arguments expand to exclusive root.

    Caller identities must include scripts, assets, settings and loaded native
    modules. This layer additionally pins its machinery and argv executable.
    Retrying within the payload shares the same job cap; reusing root is refused.
    Controller-only preparation/final audit is outside cgroup memory, but its
    elapsed time is included in deadline admission and final status checks.
    """
    limits={name:_integer(value,0 if name.startswith('minimum_') else 1,2**63-1)
        for name,value in dict(maximum_wall_ns=maximum_wall_ns,maximum_memory_bytes=maximum_memory_bytes,
            maximum_disk_bytes=maximum_disk_bytes,minimum_start_free_bytes=minimum_start_free_bytes,
            minimum_free_bytes=minimum_free_bytes).items()}
    if (not isinstance(argv,(list,tuple)) or not argv or not Path(argv[0]).is_absolute() or
        any(type(a) is not str or '\0' in a for a in argv)):
        raise ValueError('ABSOLUTE_JOB_ARGV_REQUIRED')
    if not identity_files or any(k.startswith('__job_') for k in identity_files):raise ValueError('JOB_IDENTITY_REQUIRED')
    env={'PATH':'/usr/bin:/bin','LANG':'C.UTF-8','LC_ALL':'C.UTF-8',
         'PYTHONPATH':str(Path(__file__).resolve().parents[2]),
         'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
    if environment:
        if set(environment)-{'PYTHONPATH','LD_LIBRARY_PATH'}:raise ValueError('UNSUPPORTED_JOB_ENVIRONMENT')
        if any(type(v) is not str or '\0' in v for v in environment.values()):raise ValueError('INVALID_JOB_ENVIRONMENT')
        env.update(environment)
    root=Path(directory).absolute()
    if root.exists():raise FileExistsError(root)
    if shutil.disk_usage(root.parent).free<minimum_start_free_bytes:raise ValueError('RESOURCE_LIMIT_START_FREE_SPACE')
    began=time.clock_gettime_ns(time.CLOCK_BOOTTIME);deadline=began+maximum_wall_ns
    files=dict(identity_files)
    helper=Path(__file__).with_name('_source_job_exec.py')
    for name,path in [('supervisor',Path(__file__)),('worker',helper),('python',Path(sys.executable)),
                      ('store',Path(__file__).with_name('source_checkpoint_store.py')),
                      ('accounting',Path(__file__).with_name('source_run_supervisor.py')),
                      ('systemd_run',Path('/usr/bin/systemd-run')),('systemctl',Path('/usr/bin/systemctl')),
                      ('executable',Path(argv[0]))]:files['__job_'+name]=path
    identity=_fingerprints(files)
    root.mkdir(mode=0o700)
    parent=os.open(root.parent,os.O_RDONLY|os.O_DIRECTORY|os.O_CLOEXEC)
    try:os.fsync(parent)
    finally:os.close(parent)
    unit='spin2014-job-'+uuid.uuid4().hex+'.service'
    config=dict(schema='source-job-v1',unit=unit,identity=identity,
        identity_files={k:str(Path(v).absolute()) for k,v in files.items()},limits=limits,environment=env,
        argv=[a.replace('{job}',str(root)) for a in argv],started_boottime_ns=began,
        deadline_boottime_ns=deadline,boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip())
    _publish(root/'config.json',config)
    for name in ('stdout.log','stderr.log'):(root/name).touch(exist_ok=False)
    remaining=deadline-time.clock_gettime_ns(time.CLOCK_BOOTTIME)
    if remaining<=0:raise ValueError('RESOURCE_LIMIT_WALL_BEFORE_START')
    command=['/usr/bin/systemd-run','--user','--no-ask-password','--quiet',
        '--expand-environment=no','--unit='+unit,'--service-type=exec',
        '--property=RemainAfterExit=yes','--property=KillMode=control-group',
        '--property=TimeoutStopSec=100ms','--property=OOMPolicy=kill',
        '--property=MemoryMax='+str(maximum_memory_bytes),'--property=MemorySwapMax=0',
        '--property=RuntimeMaxSec='+str(max(1,remaining//1_000_000))+'ms',
        '--property=StandardOutput=append:'+str(root/'stdout.log'),
        '--property=StandardError=append:'+str(root/'stderr.log'),
        '--working-directory='+str(root)]
    command+=['--setenv='+key+'='+value for key,value in env.items()]
    command += [sys.executable,'-I',str(helper),str(root/'config.json')]
    service={};status='CONTROLLER_FAILED';code=-1;peak=_disk_bytes(root);failure=None
    try:
        launch=subprocess.run(command,capture_output=True,text=True,timeout=5)
        _publish(root/'launch.json',dict(argv=command,returncode=launch.returncode,stderr=launch.stderr))
        if launch.returncode:raise RuntimeError('JOB_BACKEND_UNAVAILABLE: '+launch.stderr)
        while True:
            service=_show(unit)
            peak=max(peak,_disk_bytes(root))
            now=time.clock_gettime_ns(time.CLOCK_BOOTTIME)
            if now>=deadline:status='RESOURCE_LIMIT_WALL';break
            if peak>=maximum_disk_bytes:status='RESOURCE_LIMIT_DISK';break
            if shutil.disk_usage(root).free<minimum_free_bytes:status='RESOURCE_LIMIT_FREE_SPACE';break
            if service.get('Result')=='oom-kill':status='RESOURCE_LIMIT_MEMORY';break
            if service.get('Result')=='timeout':status='RESOURCE_LIMIT_WALL';break
            if service.get('SubState')=='exited' or service.get('ActiveState') in ('failed','inactive'):
                code=int(service.get('ExecMainStatus','-1'))
                if service.get('ExecMainCode')!='1':code=-code if code else -1
                status='PROCESS_EXITED' if service.get('Result')=='success' and code==0 else 'PROCESS_FAILED'
                break
            time.sleep(.02)
    except BaseException as exc:failure=exc
    finally:
        # Target only this newly generated unit; never stop unrelated user jobs.
        cleanup=subprocess.run(['/usr/bin/systemctl','--user','stop',unit],capture_output=True,text=True,timeout=5)
        _publish(root/'service.json',dict(observed=service,stop_returncode=cleanup.returncode,stop_stderr=cleanup.stderr))
        subprocess.run(['/usr/bin/systemctl','--user','reset-failed',unit],capture_output=True,timeout=3)
    if failure is not None:raise failure
    if cleanup.returncode:raise RuntimeError('JOB_CLEANUP_FAILED')
    if _fingerprints(files)!=identity:raise ValueError('STALE_JOB_IDENTITY')
    if status=='PROCESS_EXITED' and not (root/'kernel.json').is_file():raise RuntimeError('MISSING_KERNEL_CONTAINMENT_EVIDENCE')
    elapsed=time.clock_gettime_ns(time.CLOCK_BOOTTIME)-began
    if elapsed>=maximum_wall_ns and status!='RESOURCE_LIMIT_MEMORY':status='RESOURCE_LIMIT_WALL'
    peak=max(peak,_disk_bytes(root))
    result=dict(schema='source-job-result-v1',unit=unit,status=status,returncode=code,elapsed_ns=elapsed,
        service=service,scientific_validation='NOT_GRANTED',disk_accounting='SAMPLED_WITH_METADATA_RESERVE')
    block=max(512,os.statvfs(root).f_frsize)
    peak+=((len(_canonical(result))+512+block-1)//block+5)*block
    if peak>=maximum_disk_bytes:status='RESOURCE_LIMIT_DISK';result['status']=status
    result['accounted_disk_bytes']=peak
    _publish(root/'result.json',result)
    if _disk_bytes(root)+3*block>peak:raise RuntimeError('JOB_PUBLICATION_ACCOUNTING_BOUND')
    _publish(root/'complete.json',dict(result_sha256=hashlib.sha256((root/'result.json').read_bytes()).hexdigest()))
    return SourceJobResult(root,unit,status,code,elapsed,service)
