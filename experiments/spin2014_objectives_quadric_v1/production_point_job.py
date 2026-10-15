"""One frozen production map point, not an optimization or derivative verdict.

Launch sequentially through verification.source_job.run_source_job, which
enforces the absolute 60-minute cgroup deadline. No resume/stitch of optimizers.
Each terminal packet is independently integrated after durable actual capture;
the unified aggregator must re-read all 81 packets and verify the FD protocol.
"""
import argparse
from fractions import Fraction
import json
import os
from pathlib import Path
import shutil
import time
import threading
import signal

from make_it_stand.verification.quadric_point_jobs import (
    _read,validate_point_manifest,verify_point_sources,file_sha256)


def save(path,value):
    with path.open('x') as stream:
        json.dump(value,stream,sort_keys=True,allow_nan=False)
        stream.write('\n');stream.flush();os.fsync(stream.fileno())
    descriptor=os.open(path.parent,os.O_RDONLY)
    try: os.fsync(descriptor)
    finally: os.close(descriptor)


def require_point_resource_envelope(envelope):
    relative=Path('/proc/self/cgroup').read_text().strip()
    root=Path('/sys/fs/cgroup')/relative[3:].lstrip('/')
    if not relative.startswith('0::/') or (root/'memory.max').read_text().strip()!=str(24*1024**3) or \
       (root/'memory.swap.max').read_text().strip()!='0':
        raise RuntimeError('production point requires actual 24 GiB / swap 0 containment')
    # An existing source-job guardian must already enforce the absolute deadline.
    # Its kernel receipt is read, and checked against this process's real cgroup.
    if not envelope: raise RuntimeError('source-job absolute-deadline receipt required')
    try:
        kernel=_read(envelope)
    except ValueError as error:
        raise RuntimeError('SOURCE_JOB_RECEIPT_UNAVAILABLE: kernel/deadline admission failed') from error
    now=time.clock_gettime_ns(time.CLOCK_BOOTTIME)
    if kernel.get('cgroup')!=relative[3:] or kernel.get('memory_max')!=24*1024**3 or \
       kernel.get('swap_max')!=0 or kernel.get('boot_id')!=Path('/proc/sys/kernel/random/boot_id').read_text().strip() or \
       not now < kernel.get('absolute_deadline_boottime_ns',0) <= now+3600_000_000_000:
        raise RuntimeError('source-job kernel/deadline identity mismatch')


def start_disk_watch(directory):
    if shutil.disk_usage(directory).free<20*1024**3:
        raise RuntimeError('RESOURCE_LIMIT: less than 20 GiB disk free')
    done=threading.Event()
    def watch():
        while not done.wait(30):
            if shutil.disk_usage(directory).free<20*1024**3:
                print('RESOURCE_LIMIT: disk below 20 GiB; incomplete evidence retained',flush=True)
                os.kill(os.getpid(),signal.SIGTERM)
    threading.Thread(target=watch,daemon=True).start()
    return done


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--manifest-identity',required=True)
    parser.add_argument('--index',type=int,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--job-envelope',type=Path,
                        help='Existing source_job {job}/kernel.json; no environment whitelist change')
    args=parser.parse_args()
    require_point_resource_envelope(args.job_envelope)
    manifest=_read(args.manifest);validate_point_manifest(manifest)
    if manifest['identity']!=args.manifest_identity or not 0<=args.index<81:
        raise ValueError('wrong frozen manifest or point index')
    settings=manifest['mapping_settings']
    paths=settings['source_paths'];verify_point_sources(manifest,paths)
    from make_it_stand import spin
    from make_it_stand.verification import _spin2014_authority as checker
    from make_it_stand.verification import mass_moments,rational_text,quadric_point_jobs,surface_derivatives
    native=spin._native_module()
    loaded={'native':native.__file__,'runner':__file__,'checker':checker.__file__,
            'mapping_interface':spin.__file__,'mass_moments':mass_moments.__file__,
            'rational_text':rational_text.__file__,'manifest_checker':quadric_point_jobs.__file__,
            'surface_derivatives':surface_derivatives.__file__}
    if any(name not in paths or Path(paths[name]).resolve()!=Path(path).resolve()
           for name,path in loaded.items()):
        raise ValueError('loaded production dependencies differ from frozen paths')
    request=dict(settings['request']);request['settings']=spin.SpinSolverConfig(**request['settings'])
    _,native_request=spin._native_request(spin.SpinItFasterRequest(**request))
    evaluator_config=spin.SpinEvaluatorConfig(**settings['evaluator'])
    if evaluator_config.icosphere_subdivisions!=2 or evaluator_config.maximum_icosphere_subdivisions!=2:
        raise ValueError('actual producer differs from frozen level 2 normal discretization')
    evaluator=spin._native_evaluator_settings(native,evaluator_config)
    density=settings['fixed_material_density']
    raw=manifest['points'][args.index]['raw_coefficients']
    authority=native.normalize_spin2014_authority_coefficients(
        [rational_text.fraction_to_rational_text(Fraction(x)) for x in raw])
    if shutil.disk_usage(args.output.parent).free<20*1024**3:
        raise RuntimeError('RESOURCE_LIMIT: less than 20 GiB disk free')
    args.output.mkdir(exist_ok=False)
    save(args.output/'started.json',dict(schema='quadric-point-job-start-v1',
        manifest_sha256=manifest['identity'],index=args.index,status='STARTED_NOT_COMPLETED',
        raw_coefficients=raw,scope='ONE_PRODUCTION_POINT_NOT_DERIVATIVE_GATE'))
    disk_done=start_disk_watch(args.output)
    begun=time.monotonic();actual_calls=0;packet=dict(index=args.index,raw_coefficients=raw)
    status='FAILED'
    try:
        actual_calls=1
        packet['evaluation']=native.evaluate_spin2014_authority_material_mapping(
            native_request,density,authority['geometry_authority'],authority['dual_authority'],evaluator)
        # Persist actual geometry before independent work; interruption preserves evidence.
        save(args.output/'point.json',packet)
        saved_packet=_read(args.output/'point.json')
        independent=checker.check_independent_authority_material_mapping(saved_packet['evaluation'],
            expected_request_identity=settings['expected_request_identity'],
            expected_context_identity=settings['expected_context_identity'],
            arithmetic_backend='GMPY2_2_2_1_EXACT_V1')
        save(args.output/'independent.json',independent)
        verify_point_sources(manifest,paths)
        if _read(args.manifest)!=manifest:
            raise ValueError('frozen manifest changed during actual point job')
        status='COMPLETE'
    except Exception as error:
        if not (args.output/'point.json').exists():
            packet.update(exception_type=type(error).__name__,exception=str(error))
            save(args.output/'point.json',packet)
        save(args.output/'failure.json',dict(exception_type=type(error).__name__,exception=str(error)))
    save(args.output/'terminal.json',dict(schema='quadric-point-job-terminal-v1',
        manifest_sha256=manifest['identity'],index=args.index,status=status,
        actual_material_evaluations=actual_calls,packet_sha256=file_sha256(args.output/'point.json'),
        independent_sha256=file_sha256(args.output/'independent.json') if status=='COMPLETE' else None,
        elapsed_seconds=time.monotonic()-begun,optimizer_invocations=0,
        scientific_gate_granted=False))
    disk_done.set()
    return 0 if status=='COMPLETE' else 1


if __name__=='__main__': raise SystemExit(main())
