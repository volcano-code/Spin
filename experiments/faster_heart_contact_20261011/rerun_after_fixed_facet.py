"""One bounded fresh seed-0 run after the registered local diagnostic passes."""
import hashlib
import json
from pathlib import Path
import sys

from make_it_stand.verification.source_job import run_source_job
from make_it_stand.verification.spin_common import digest

ROOT = Path(__file__).resolve().parents[2]
CASE = Path(__file__).parent / 'asymmetric-candidate-v1'
BUILD = ROOT / 'build/heart-fixed-facet-20261011-v1'
NATIVE = BUILD / '_make_it_stand.cpython-313-x86_64-linux-gnu.so'


def main():
    local = CASE / 'fixed-facet-difference-v1'
    receipt = json.loads((local / 'result.json').read_text())
    assert receipt['status'] == 'PROCESS_EXITED' and receipt['returncode'] == 0
    report = json.loads((local / 'difference-report.json').read_text())
    assert report['all_errors_at_most_1e_minus_6']
    for point in local.glob('difference-*-*.json'):
        assert json.loads(point.read_text())['derivative_defined'], point
    assert len(list(local.glob('difference-*-*.json'))) == 6
    # Local direction check is a repair prerequisite, NOT a claim that the
    # complete 81-point protocol or final material verification has passed.
    frozen = json.loads((CASE / 'frozen-requests-v3.json').read_text())
    parent = frozen.pop('frozen_payload_digest')
    assert digest(frozen) == parent
    old_native = ROOT / 'build/heart-exact-input-20261011-v2' / NATIVE.name
    del frozen['input_sources'][str(old_native)]
    frozen['input_sources'][str(NATIVE)] = hashlib.sha256(NATIVE.read_bytes()).hexdigest()
    frozen['schema'] = 'asymmetric-heart-faster-frozen-requests-v4'
    frozen['parent_frozen_payload_digest'] = parent
    frozen['change'] = 'EXACT_STATIONARY_DESIGN_FACET_TRANSPORT_REPAIR'
    frozen['local_diagnostic_scope'] = 'ONE_DIRECTION_THREE_SCALES_NOT_FULL_DERIVATIVE_ACCEPTANCE'
    frozen['frozen_payload_digest'] = digest(frozen)
    target = CASE / 'frozen-requests-v4.json'
    with target.open('x') as out:
        json.dump(frozen, out, indent=2, allow_nan=False)
    previous = json.loads((CASE / 'solve-seed-0-v3/config.json').read_text())
    argv = previous['argv']
    code = argv[argv.index('-c')+1]
    files = {name:path for name,path in previous['identity_files'].items()
             if not name.startswith('__job_')}
    files.update(native=str(NATIVE), frozen=str(target), controller=__file__,
        local_derivative_report=str(local / 'difference-report.json'))
    run = run_source_job(CASE / 'solve-seed-0-v4',
        ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
         'OMP_NUM_THREADS=1', 'OPENBLAS_NUM_THREADS=1', 'MKL_NUM_THREADS=1',
         sys.executable, '-B', '-c', code, str(target), '0', '{job}'],
        files, environment={'PYTHONPATH':str(ROOT / 'python')+':'+str(BUILD)})
    print(run, flush=True)


if __name__ == '__main__':
    main()
