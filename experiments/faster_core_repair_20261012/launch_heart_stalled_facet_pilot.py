"""Run the finite facet witness pilot after the production probe terminates."""
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job


def main():
    generation = sys.argv[1] if len(sys.argv)>1 else '01'
    if generation not in ('01','02','03'):
        raise ValueError('UNREGISTERED_STALLED_FACET_PILOT')
    probe = HERE/'heart-stalled-production-81-job-01/output'
    if not (probe/'report.json').is_file():
        raise ValueError('WAIT_FOR_PRODUCTION_PROBE_TERMINAL')
    frozen = HERE.parent/'faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json'
    files = dict(launcher=Path(__file__),worker=HERE/'classify_heart_stalled_facets.py',frozen=frozen,
        source_manifest=probe/'manifest.json',source_events=probe/'events.jsonl',source_report=probe/'report.json',
        point0=probe/'point-000-material.json',point16=probe/'point-016-material.json',
        spin=ROOT/'python/make_it_stand/spin.py')
    for i,path in enumerate(json.loads(frozen.read_text())['input_sources']):
        files['input_'+str(i)] = Path(path)
    for path in (ROOT/'python/make_it_stand/verification').glob('*.py'):
        files['module_'+path.name] = path
    if generation == '03':
        for path in probe.glob('point-*-material.json'):
            files['sample_'+path.name] = path
    result = run_serial_source_job(ROOT,HERE/('heart-stalled-facet-pilot-job-'+generation),
        ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output']+
        (['--all-samples'] if generation == '03' else []),files,
        environment={'PYTHONPATH':str(ROOT/'python')})
    print(result.status,result.returncode,flush=True)
    if result.status != 'PROCESS_EXITED' or result.returncode != 0:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
