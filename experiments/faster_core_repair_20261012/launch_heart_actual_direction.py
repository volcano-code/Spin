"""Serial bounded launch of the preregistered Heart directional probe."""
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job


def main():
    suffix = '--backtrack-suffix' in sys.argv
    local = '--local-suffix' in sys.argv
    benchmark_args = [x for x in sys.argv if x.startswith('--benchmark-state=')]
    prepared = '--prepared-arithmetic' in sys.argv
    certified = '--certified-integral' in sys.argv
    if benchmark_args and int(benchmark_args[0].split('=')[1]) not in (2,38,45,46):
        raise ValueError('UNREGISTERED_BENCHMARK_STATE')
    build = HERE/'formal-native-core-build-job-12/output/build'
    bridge = HERE/'formal-cherchi-build-job-02/output/build/python'
    arithmetic = HERE/('formal-rational-arithmetic-build-job-05/output/build/python' if certified else
                      'formal-rational-arithmetic-build-job-03/output/build/python' if prepared else
                      'formal-rational-arithmetic-build-job-02/output/build/python')
    run = HERE/'heart-guarded-formal-seed-0-job-01/output/run'
    frozen = HERE.parent/'faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json'
    files = dict(launcher=Path(__file__), worker=HERE/'probe_heart_actual_direction.py',
        registration=HERE/('heart-actual-direction-probe-registration-05.json' if local else
                            'heart-actual-direction-probe-registration-04.json' if suffix else
                            'heart-actual-direction-probe-registration-03.json'), frozen=frozen,
        native=build/'python/_make_it_stand.cpython-313-x86_64-linux-gnu.so',
        runtime=build/'compiler-runtime/libstdc++.so.6', receipt=build.parent/'receipt.json',
        algebra=HERE/'audit_heart_stalled_prefix.py', residual=HERE/'audit_heart_trial_prediction.py',
        integral=HERE/'independent_tetra_moments.py', spin=ROOT/'python/make_it_stand/spin.py')
    for i, path in enumerate(json.loads(frozen.read_text())['input_sources']):
        files['input_'+str(i)] = Path(path)
    for directory in (ROOT/'python/make_it_stand/verification', run, bridge, arithmetic):
        for p in directory.iterdir():
            if p.is_file(): files[str(p.relative_to(ROOT))] = p
    benchmark_job = ('heart-arithmetic-'+('05-certified-integral' if certified else '03' if prepared else '02')+'-state-'+
                     benchmark_args[0].split('=')[1]+'-job-01') if benchmark_args else None
    result = run_serial_source_job(ROOT, HERE/(benchmark_job if benchmark_job else 'heart-actual-direction-local-job-01' if local else
                                              'heart-actual-direction-backtrack-job-01' if suffix else
                                              'heart-actual-direction-probe-job-03'),
        ['/usr/bin/env', 'LD_PRELOAD='+str(files['runtime']), '/home/lin/anaconda3/bin/python', '-B',
         str(files['worker']), '{job}/output']+benchmark_args+(['--certified-integral'] if certified else [])+(['--local-suffix'] if local else ['--backtrack-suffix'] if suffix else []), files,
        environment={'PYTHONPATH': ':'.join(map(str, [ROOT/'python', build/'python', bridge, arithmetic, HERE]))})
    print(result.status, result.returncode, flush=True)
    if result.status != 'PROCESS_EXITED' or result.returncode != 0:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
