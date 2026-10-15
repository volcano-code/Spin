"""Compare formal seed runs; discard only explicitly named timing fields."""
import json
from pathlib import Path
from compare_pilot_repeats import scientific,digest

HERE=Path(__file__).resolve().parent
paths=[HERE/name/'output' for name in ('formal-cherchi-small-closure-job-04','formal-cherchi-small-repeat-job-01')]
def load(p,name):return json.loads((p/name).read_text())
configs=[load(p,'run/configuration.json') for p in paths]
reports=[load(p,'report.json') for p in paths]
events=[]
for p in paths:
    rows=[json.loads(line) for line in (p/'run/events.jsonl').read_text().splitlines()]
    events.append([scientific(row['event']) for row in rows if row['kind']=='NEWTON_EVENT'])
materials=[{f.name:load(p,'run/'+f.name) for f in sorted((p/'run').glob('evaluation-*.json'))} for p in paths]
terminal=[load(p,'run/terminal-material-reference.json') for p in paths]
results=[load(p,'run/result.json')['identity'] for p in paths]
for identity in results:
    native=identity['native_run']
    for name in ('event_hash','events','evaluation_timings'):native.pop(name,None)
checks=dict(configuration_equal=configs[0]==configs[1],scientific_events_equal=events[0]==events[1],
    all_same_call_materials_equal=materials[0]==materials[1],terminal_references_equal=terminal[0]==terminal[1],
    result_equal=scientific(results[0])==scientific(results[1]),
    sealed_material_equal=reports[0]['material_digest']==reports[1]['material_digest'],
    both_closures_passed=all(r['closure_passed'] for r in reports))
report=dict(checks=checks,passed=all(checks.values()),scientific_event_hashes=[digest(e) for e in events],
    scope='Formal normalized small seed1 repeat only; original derivative gate and core admission remain incomplete')
with (HERE/'formal-cherchi-repeat-comparison-01.json').open('x') as stream:json.dump(report,stream,indent=2)
print(json.dumps(report))
if not report['passed']:raise SystemExit(1)
