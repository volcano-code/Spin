"""Verify real kernel containment before executing a whole scientific workflow."""
import json
import os
from pathlib import Path
import sys
import time
import subprocess

# Invoked with Python -I: payload PYTHONPATH must not select the verifier.
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from make_it_stand.verification.source_checkpoint_store import _fingerprints
from make_it_stand.verification.source_job import _publish


if __name__=='__main__':
    config_path=Path(sys.argv[1]);config=json.loads(config_path.read_text())
    if _fingerprints(config['identity_files'])!=config['identity']:
        raise ValueError('STALE_JOB_IDENTITY')
    relative=Path('/proc/self/cgroup').read_text().strip()
    if not relative.startswith('0::/'):raise RuntimeError('UNIFIED_CGROUP_REQUIRED')
    relative=relative[3:]
    root=Path('/sys/fs/cgroup')/relative.lstrip('/')
    if root.name!=config['unit']:raise RuntimeError('JOB_CGROUP_MISMATCH')
    memory=int((root/'memory.max').read_text());swap=int((root/'memory.swap.max').read_text())
    if memory!=config['limits']['maximum_memory_bytes'] or swap!=0:
        raise RuntimeError('JOB_KERNEL_LIMIT_MISMATCH')
    boot=Path('/proc/sys/kernel/random/boot_id').read_text().strip()
    if boot!=config['boot_id']:raise RuntimeError('JOB_BOOT_CHANGED')
    if time.clock_gettime_ns(time.CLOCK_BOOTTIME)>=config['deadline_boottime_ns']:
        raise RuntimeError('JOB_DEADLINE_BEFORE_EXEC')
    # An absolute BOOTTIME timer, armed before spawning any payload, survives
    # controller death and service activation delay. The guardian remains in
    # this cgroup even if the payload creates another session or exits early.
    # Python 3.13 timerfd + writable delegated cgroup.kill are mandatory.
    timer=os.timerfd_create(time.CLOCK_BOOTTIME,flags=os.TFD_CLOEXEC)
    killer=os.open(root/'cgroup.kill',os.O_WRONLY|os.O_CLOEXEC)
    os.timerfd_settime_ns(timer,flags=os.TFD_TIMER_ABSTIME,initial=config['deadline_boottime_ns'])
    guardian=os.fork()
    if guardian==0:
        try:
            os.read(timer,8)
            os.write(killer,b'1')
        finally:os._exit(124)
    os.close(timer);os.close(killer)
    _publish(config_path.parent/'kernel.json',dict(cgroup=relative,memory_max=memory,
        swap_max=swap,boot_id=boot,observed_boottime_ns=time.clock_gettime_ns(time.CLOCK_BOOTTIME),
        absolute_deadline_boottime_ns=config['deadline_boottime_ns'],guardian_pid=guardian))
    code=subprocess.run(config['argv'],env=config['environment'],close_fds=True).returncode
    raise SystemExit(code if code>=0 else 128-code)
