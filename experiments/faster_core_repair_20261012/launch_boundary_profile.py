"""Serial, resource-limited gen05 boundary baseline versus instrumented native07."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
from boundary_profile_inputs import select_saved_states

mode = sys.argv[1]
generation = sys.argv[2] if len(sys.argv)>2 else '01'
if generation not in ('01','02'): raise ValueError('UNREGISTERED_PROFILE_GENERATION')
build = '08' if generation=='02' else '07'
basebuild = '07' if generation=='02' else '06'
native = HERE/f'formal-native-core-build-job-{build}/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
baseline = HERE/f'formal-native-core-build-job-{basebuild}/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
files = dict(launcher=Path(__file__), native=native, baseline=baseline,
    receipt=HERE/f'formal-native-core-build-job-{build}/output/receipt.json',
    selector=HERE/'boundary_profile_inputs.py')
for path in (ROOT/'python/make_it_stand').rglob('*.py'):
    files['python_'+str(path.relative_to(ROOT))] = path
env = {'PYTHONPATH': ':'.join(map(str, [ROOT/'python', native.parent, HERE]))}
if mode == 'validate':
    names = ('test_boundary_phase_profile', 'test_verified_fixed_boundary_reuse',
             'test_exact_boundary_audit', 'test_rational_boundary_audit', 'test_rational_solver_snapshot')
    tests = [ROOT/'tests/python'/f'{name}.py' for name in names]
    tests.append(HERE/'test_boundary_profile_inputs.py')
    if generation=='02':
        tests.extend([ROOT/'tests/python/test_diagonal_boundary_candidate.py', HERE/'test_boundary_candidate_comparison.py'])
        files['comparison'] = HERE/'boundary_candidate_comparison.py'
    files.update({str(path): path for path in tests})
    command = ['/home/lin/anaconda3/bin/python3', '-B', '-m', 'pytest', *map(str,tests), '-q']
    job = HERE/f'boundary-phase-validation-job-{generation}'
elif mode == 'measure':
    files['worker'] = HERE/'profile_saved_boundaries.py'
    files['admission'] = HERE/f'boundary-phase-validation-job-{generation}/result.json'
    if generation=='02': files['comparison'] = HERE/'boundary_candidate_comparison.py'
    import json
    if json.loads(files['admission'].read_text())['returncode'] != 0:
        raise ValueError('PROFILE_TESTS_NOT_PASSED')
    for scale in ('small','faces-12288','faces-98304'):
        output = HERE/f'layered-formal-benchmark-{scale}-seed-0-job-05/output'
        for relative in ('report.json','material/problem.json','run/configuration.json',
                         'run/events.jsonl','run/terminal-material-reference.json'):
            files[scale+relative] = output/relative
        files[scale+'input'] = HERE/'core-benchmark-inputs-job-01/output'/f'{scale}.json'
        for selected in select_saved_states(output/'run'):
            files[scale+selected['role']] = Path(selected['path'])
    command = ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output',str(baseline)]
    if generation=='02': command.append('diagonal-candidate')
    job = HERE/f'boundary-phase-measurement-job-{generation}'
else:
    raise ValueError('UNKNOWN_PROFILE_JOB')
command = ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6', *command]
result = run_source_job(job, command, files, environment=env)
print(result.status, result.returncode, flush=True)
