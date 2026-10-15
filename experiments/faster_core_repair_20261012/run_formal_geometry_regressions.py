"""Saved failure regressions only; no Heart or Noodles optimization starts."""
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def main(output):
    output.mkdir(parents=True,exist_ok=False)
    commands=[
        [sys.executable,'-m','pytest',str(HERE/'test_saved_heart_projective_source.py'),
         str(HERE/'test_registered_source_side.py'),'-q'],
        [sys.executable,'-B',str(HERE/'check_repaired_heart_states.py'),'--output',str(output/'heart')],
        [sys.executable,'-B',str(HERE/'verify_small_face_candidate.py'),str(output/'noodles')],
        [sys.executable,'-m','pytest',str(HERE/'test_registered_material_normalization.py'),'-q'],
    ]
    for i,command in enumerate(commands):
        with (output/('stage-'+str(i)+'.log')).open('x') as stream:
            subprocess.run(command,cwd=ROOT,stdout=stream,stderr=stream,check=True)
    heart=json.loads((output/'heart/report.json').read_text())['records']
    noodles=json.loads((output/'noodles/report.json').read_text())
    passed=(len(heart)==2 and all(r.get('actual_material_returned') and r.get('derivative_defined') and
        r.get('source_boundary_valid') for r in heart) and noodles.get('evaluation_completed') and
        noodles.get('derivative_defined') and noodles.get('source_boundary_valid'))
    (output/'report.json').write_text(json.dumps(dict(passed=bool(passed),heart=heart,noodles=noodles,
        commands=commands,optimization_invocations=0,core_admission=False,
        scope='Registered failed source/material states on ordinary formal native build'),allow_nan=False))
    if not passed:
        raise ValueError('FORMAL_REGISTERED_GEOMETRY_REGRESSION_FAILED')


if __name__=='__main__':main(Path(sys.argv[1]))
