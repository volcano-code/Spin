"""Full seed regression after actual physical-source admission fix."""
import json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job
c=json.loads((HERE/'spin2014-shared-small-top-inertia-only-seed-1-run-01/config.json').read_text())
argv=c['argv'];argv[5]='{job}/output'
files={str(p.relative_to(ROOT)):p for p in (ROOT/'python/make_it_stand/verification').glob('*.py')}
files.update(worker=HERE/'run_core_benchmark.py',launcher=Path(__file__),spin=ROOT/'python/make_it_stand/spin.py')
for directory in c['environment']['PYTHONPATH'].split(':'):
    for p in Path(directory).glob('*.so'):files[str(p)]=p
r=run_serial_source_job(ROOT,HERE/'spin2014-inertia-actual-source-guard-seed1-job-04',argv,files,environment={'PYTHONPATH':c['environment']['PYTHONPATH']})
print(r.status,r.returncode,flush=True)
