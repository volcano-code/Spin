"""Sample a diagnostic child under GDB; never attach to a production worker.

Samples locate native hotspots, not unbiased wall-time estimates. The original
worker and its material checks run unchanged. No scientific admission is granted.
"""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import re
import selectors
import subprocess
import time


def classify_stack(record):
    # A GMP leaf alone cannot identify its caller's geometry operation.
    if 'PairScan' in record:
        return 'BOUNDARY_INTERSECTION_SCAN'
    if 'auditPoints' in record or 'auditExactTriangleBoundary' in record:
        return 'BOUNDARY_AUDIT'
    if any(name in record for name in ('Nef_polyhedron', 'SNC_', 'Nef::', 'finiteRegularSolid')):
        return 'NEF_CONSTRUCTION'
    if 'ComponentBuilder' in record or 'Polyhedron_incremental_builder' in record:
        return 'COMPONENT_SURFACE_BUILD'
    if '__gmp' in record:
        return 'EXACT_ARITHMETIC_UNRESOLVED_CALLER'
    return 'UNRESOLVED'


def profile_command(command, output, *, interval=60., maximum_samples=24, timeout=3500.):
    if interval <= 0 or maximum_samples < 1 or timeout <= 0:
        raise ValueError('INVALID_PROFILE_BUDGET')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    samples = []
    started = time.monotonic()
    debugger = subprocess.Popen(
        ['/usr/bin/gdb', '--nx', '--quiet', '--interpreter=mi2', '--args', *command],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    selector = selectors.DefaultSelector()
    selector.register(debugger.stdout, selectors.EVENT_READ)
    pending = bytearray()
    records = []
    exit_code = None
    running = False
    raw = (output / 'gdb-mi.log').open('xb')
    journal = (output / 'stack-samples.jsonl').open('x')

    def receive(wait=.1):
        nonlocal exit_code, running
        if time.monotonic() - started > timeout:
            raise TimeoutError('NATIVE_PROFILE_TIMEOUT')
        if selector.select(wait):
            chunk = os.read(debugger.stdout.fileno(), 65536)
            if not chunk:
                selector.unregister(debugger.stdout)
                if exit_code is None:
                    raise RuntimeError('DEBUGGER_EXITED_WITHOUT_INFERIOR_RESULT')
                return []
            raw.write(chunk); raw.flush()
            pending.extend(chunk)
        lines = []
        while b'\n' in pending:
            line, _, rest = pending.partition(b'\n')
            pending[:] = rest
            value = line.decode('utf-8', errors='replace').strip()
            lines.append(value); records.append(value)
            if value.startswith('*running'):
                running = True
            elif value.startswith('*stopped'):
                running = False
                if 'reason="exited-normally"' in value:
                    exit_code = 0
                elif 'reason="exited"' in value:
                    found = re.search(r'exit-code="([0-7]+)"', value)
                    if not found:
                        raise RuntimeError('MISSING_INFERIOR_EXIT_CODE')
                    exit_code = int(found.group(1), 8)
                elif 'reason="exited-signalled"' in value:
                    exit_code = 128
        return lines

    token = 0
    def issue(operation):
        nonlocal token
        token += 1
        prefix = str(token)
        debugger.stdin.write((prefix + operation + '\n').encode()); debugger.stdin.flush()
        while True:
            for record in receive():
                if record.startswith(prefix + '^error'):
                    raise RuntimeError(record)
                if record.startswith(prefix + '^'):
                    return record

    try:
        issue('-gdb-set pagination off')
        issue('-gdb-set mi-async on')
        issue('-gdb-set disable-randomization off')
        issue('-exec-run')
        next_sample = time.monotonic() + interval
        while exit_code is None:
            receive(.05)
            if exit_code is not None:
                break
            if not running:
                raise RuntimeError('UNEXPECTED_INFERIOR_STOP: ' + records[-1])
            if len(samples) < maximum_samples and time.monotonic() >= next_sample:
                issue('-exec-interrupt --all')
                while running and exit_code is None:
                    receive()
                if exit_code is not None:
                    break
                stack = issue('-stack-list-frames 0 30')
                row = dict(seconds=time.monotonic()-started, category=classify_stack(stack), stack=stack)
                samples.append(row)
                journal.write(json.dumps(row)+'\n'); journal.flush(); os.fsync(journal.fileno())
                issue('-exec-continue')
                next_sample = time.monotonic() + interval
        report = dict(schema='native-geometry-stack-profile-v1',
            scope='DIAGNOSTIC_SAMPLES_NOT_PHASE_TIMING_OR_SCIENTIFIC_ACCEPTANCE',
            command=command, inferior_exit_code=exit_code, samples=samples,
            sample_counts=dict(Counter(row['category'] for row in samples)),
            elapsed_seconds=time.monotonic()-started,
            geometry_admission=False, predecessor_admission=False, optimization_admission=False)
        with (output/'profile.json').open('x') as stream:
            json.dump(report, stream, indent=2)
        return report
    finally:
        try:
            debugger.stdin.write(b'-gdb-exit\n'); debugger.stdin.flush()
            debugger.wait(timeout=5)
        except (OSError, subprocess.TimeoutExpired):
            debugger.kill(); debugger.wait()
        debugger.stdin.close(); debugger.stdout.close()
        selector.close(); raw.close(); journal.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--interval', type=float, default=60.)
    parser.add_argument('--maximum-samples', type=int, default=24)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        parser.error('diagnostic command required')
    report = profile_command(command, args.output, interval=args.interval, maximum_samples=args.maximum_samples)
    raise SystemExit(0 if report['inferior_exit_code'] == 0 else 1)
