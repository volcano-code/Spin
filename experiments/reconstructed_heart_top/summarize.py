"""Derive a bounded experiment summary; never grants paper release status."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
from run import write


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(); p.add_argument('run1',type=Path); p.add_argument('run2',type=Path)
    a=p.parse_args(); root=Path(__file__).resolve().parents[2]
    for run in [a.run1,a.run2]:
        summary=json.loads((run/'summary.json').read_text())
        if not summary['sealed'] or not summary['source_unchanged']: raise ValueError('invalid run')
        for name,h in json.loads((run/'source.json').read_text())['files'].items():
            if digest(Path(name))!=h: raise ValueError('STALE_RESULT: '+name)
    compared=['native-run.json','outer.obj','material.obj','cavity.obj','config.json','input.json']
    equality={name:digest(a.run1/name)==digest(a.run2/name) for name in compared}
    if not all(equality.values()): raise ValueError('optimization repeat mismatch')
    base=json.loads((a.run1/'validation-complete/report.json').read_text())
    extra=json.loads((a.run1/'supplemental-rotation/report.json').read_text())
    for report,speeds in [(base,[2,6,15]),(extra,[35,45])]:
        expected=set(itertools.product(['original','carved'],speeds,range(8)))
        ids=[(r['model'],r['hz'],r['direction']) for r in report['scenarios']]
        if len(ids)!=len(set(ids)) or set(ids)!=expected: raise ValueError('incomplete or duplicate matrix')
    executed=a.run1/'validation-complete/executed-evaluator.py'
    if digest(executed)!=base['evaluator_sha256'] or digest(Path(__file__).with_name('evaluate.py'))!=digest(executed):
        raise ValueError('STALE_RESULT: evaluator')
    tests=[]
    for name in ['test_reconstructed_heart_top.py','test_heart_linearization.py']:
        result=subprocess.run([sys.executable,str(root/'tests/python'/name)],capture_output=True,text=True)
        tests.append(dict(test=name,exit_code=result.returncode,stdout=result.stdout,stderr=result.stderr))
        if result.returncode: raise ValueError('test failed: '+name)
    write(a.run1/'tests.json',tests)
    forced=json.loads((a.run1/'forced-step-refinement/report.json').read_text())
    table=[]
    for name,hz in itertools.product(['original','carved'],[2,6,15,35,45]):
        rows=[r for r in base['scenarios']+extra['scenarios'] if r['model']==name and r['hz']==hz]
        table.append(dict(model=name,hz=hz,count=len(rows),
            worst_tilt_degrees=max(r['max_tilt_degrees'] for r in rows),
            directions_with_max_tilt_at_most_15deg=sum(r['max_tilt_degrees']<=15 for r in rows)))
    write(a.run1/'verification-summary.json',dict(status='RECONSTRUCTED_PARTIAL_VALIDATION',
        original_paper_cases_added=0,faster_solver_executed=False,physical_contact=False,
        repeat_comparisons=equality,models=base['models'],rotation_matrix=table,
        original_matrix_refinement_passed=base['all_step_refinements_passed'],
        forced_high_speed_refinement=forced,registered_local_test_cases=7))
    files={}
    for directory in [a.run1,a.run2,Path(__file__).parent]:
        for f in directory.rglob('*'):
            if f.is_file() and '__pycache__' not in f.parts and f.name!='artifact-manifest.json':
                files[str(f.resolve())]=digest(f)
    write(a.run1/'artifact-manifest.json',dict(algorithm='SHA-256',files=files,
        scope='Experiment files and declared source only; not full dependency build attestation'))
    print('Recorded actual test results, repeat comparison, all80 scenarios, and file manifest')


if __name__=='__main__': main()
