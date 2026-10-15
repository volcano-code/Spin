"""Paired exact audits of frozen same-call materials, in isolated processes."""
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
NATIVE = {
    'before': HERE/'rational-boundary-binding-build-job-01/output',
    'after': HERE/'filtered-boundary-build-job-01/output',
}

def native_directories(profile):
    if profile=='filtered':return NATIVE
    if profile=='ranked':return {
        'before':HERE/'formal-native-core-build-job-04/output/build/python',
        'after':HERE/'ranked-boundary-build-job-01/output'}
    if profile=='sharedvertex':return {
        'before':HERE/'formal-native-core-build-job-04/output/build/python',
        'after':HERE/'shared-vertex-boundary-build-job-01/output'}
    raise ValueError('UNREGISTERED_EXACT_BOUNDARY_PROFILE')


def saved_cases():
    fixed = HERE/'projective-formal-pipeline-job-01/output/fixed'
    cases = [(s, fixed/(s+'-material.json')) for s in ('ellipsoid','one_sheet','two_sheet')]
    run = HERE/'core-benchmark-faces-12288-seed-0-job-01/output/run'
    terminal = json.loads((run/'terminal-material-reference.json').read_text())['evaluation_file']
    return cases + [('medium-initial', run/'evaluation-000001.json'), ('medium-terminal', run/terminal)]


def variant(role, output,profile='filtered'):
    import numpy as np
    import _make_it_stand as native
    loaded = Path(native.__file__).resolve()
    if loaded.parent != native_directories(profile)[role].resolve():
        raise ValueError('WRONG_ACTUALLY_LOADED_NATIVE')
    rows = []
    for name, path in saved_cases():
        material = json.loads(path.read_text())
        vertices = material['variable_rational_vertices']
        faces = np.asarray(material['variable_faces'], dtype=np.int64)
        facts, samples = None, []
        for trial in range(4):
            start = time.perf_counter_ns()
            current = native.audit_exact_rational_triangle_boundary(vertices, faces, 'AnyNonzero')
            samples.append((time.perf_counter_ns()-start)/1e9)
            if facts is not None and current != facts:
                raise ValueError('NONDETERMINISTIC_BOUNDARY_FACTS')
            facts = current
        rows.append(dict(name=name, faces=len(faces), material_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                         facts=facts, cold_seconds=samples[0], warm_seconds=samples[1:],
                         warm_median_seconds=statistics.median(samples[1:])))
        (output/(role+'-rows.json')).write_text(json.dumps(rows, allow_nan=False))
    (output/(role+'.json')).write_text(json.dumps(dict(rows=rows, module=str(loaded),
        module_sha256=hashlib.sha256(loaded.read_bytes()).hexdigest()), allow_nan=False))


def main(output,profile='filtered'):
    output.mkdir(parents=True, exist_ok=False)
    for role in ('before', 'after'):
        env = dict(os.environ, LD_PRELOAD='/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
                   PYTHONPATH=str(ROOT/'python')+':'+str(native_directories(profile)[role]))
        subprocess.run([sys.executable, '-B', __file__, '--variant', role, str(output),profile],
                       env=env, cwd=ROOT, check=True)
    before, after = [json.loads((output/(s+'.json')).read_text()) for s in ('before','after')]
    rows = []
    for old, new in zip(before['rows'], after['rows'], strict=True):
        a, b = dict(old['facts']), dict(new['facts'])
        old_backend, new_backend = a.pop('backend'), b.pop('backend')
        if old['name'] != new['name'] or old['material_sha256'] != new['material_sha256']:
            raise ValueError('PAIRED_INPUT_IDENTITY_MISMATCH')
        rows.append(dict(name=old['name'], faces=old['faces'], material_sha256=old['material_sha256'],
            facts_equal=a==b, boundary_valid=a['boundary_valid'],
            before_backend=old_backend, after_backend=new_backend,
            before_cold_seconds=old['cold_seconds'], after_cold_seconds=new['cold_seconds'],
            before_warm_seconds=old['warm_seconds'], after_warm_seconds=new['warm_seconds'],
            before_warm_median_seconds=old['warm_median_seconds'],
            after_warm_median_seconds=new['warm_median_seconds'],
            speedup=old['warm_median_seconds']/new['warm_median_seconds']))
    passed = all(r['facts_equal'] and r['boundary_valid'] and
                 'Filtered_kernel<Simple_cartesian<Gmpq>>' in r['after_backend'] for r in rows)
    if profile=='ranked':passed=passed and all('exact-coordinate-rank-boxes-v1' in r['after_backend'] for r in rows)
    if profile=='sharedvertex':passed=passed and all('exact-shared-vertex-cone-v1' in r['after_backend'] for r in rows)
    (output/'report.json').write_text(json.dumps(dict(rows=rows, passed=passed, boolean_calls=0,
        core_admission=False,profile=profile, scope='Paired frozen exact boundary audits; not a full solve or performance gate'),
        allow_nan=False))
    if not passed:
        raise ValueError('FILTERED_BOUNDARY_PAIRED_FACTS_MISMATCH')


if __name__ == '__main__':
    if sys.argv[1] == '--variant':
        variant(sys.argv[2], Path(sys.argv[3]),sys.argv[4] if len(sys.argv)>4 else 'filtered')
    else:
        main(Path(sys.argv[1]),sys.argv[2] if len(sys.argv)>2 else 'filtered')
