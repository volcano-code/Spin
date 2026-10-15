"""Full saved trajectories across the boundary-only optimization.

Build-bound evaluation context is expected to change; every other material
field is compared, including same-call construction proof and Jacobian.
"""
import hashlib
import json
from pathlib import Path
from make_it_stand.verification.cherchi_run_evidence import scientific_projection,verify_cherchi_run_evidence
HERE=Path(__file__).resolve().parent
rows=[]
for scale,seed in (('small',0),('small',1),('faces-12288',0)):
    roots=[HERE/f'layered-formal-benchmark-{scale}-seed-{seed}-job-{g}/output/run' for g in ('05','06')]
    verified=[verify_cherchi_run_evidence(p) for p in roots]
    events=[[scientific_projection(json.loads(line)['event']) for line in (p/'events.jsonl').read_text().splitlines()
             if json.loads(line)['kind']=='NEWTON_EVENT'] for p in roots]
    contexts=[json.loads((p/'configuration.json').read_text())['evaluation_map_identity'] for p in roots]
    strict_event_equality=events[0]==events[1]
    for entries,context in zip(events,contexts):
        for event in entries:
            payload=event.get('evaluation_payload')
            if isinstance(payload,dict) and 'evaluation_map_identity' in payload:
                if payload['evaluation_map_identity']!=context:
                    raise ValueError('EVENT_MAP_DIFFERS_FROM_REGISTERED_CONTEXT')
                payload['evaluation_map_identity']='<REGISTERED_BUILD_BOUND_CONTEXT>'
            if isinstance(event.get('detail'),str):
                event['detail']=event['detail'].replace('evaluation_map='+context+';',
                    'evaluation_map=<REGISTERED_BUILD_BOUND_CONTEXT>;')
    names=[sorted(p.name for p in root.glob('evaluation-*.json')) for root in roots]
    if names[0]!=names[1]:raise ValueError('MATERIAL_EVENT_COUNT_CHANGED')
    changed=[];materials=[]
    for name in names[0]:
        raw=[(p/name).read_bytes() for p in roots]
        a,b=map(json.loads,raw)
        ignored={key:[a.pop(key,None),b.pop(key,None)] for key in ('evaluation_map_identity',)}
        differences=[key for key in sorted(set(a)|set(b)) if a.get(key)!=b.get(key)]
        if differences:changed.append(dict(file=name,fields=differences))
        materials.append(dict(file=name,sha256=[hashlib.sha256(data).hexdigest() for data in raw],
                              registered_context_identities=ignored))
    rows.append(dict(scale=scale,seed=seed,material_count=len(names[0]),materials=materials,
        strict_event_equality=strict_event_equality,registered_context_identities=contexts,
        changed_material_fields=changed,scientific_events_equal=events[0]==events[1],
        both_record_integrity=all(r['record_integrity_valid'] for r in verified),
        passed=not changed and events[0]==events[1] and all(r['record_integrity_valid'] for r in verified)))
report=dict(rows=rows,passed=all(r['passed'] for r in rows),core_admission=False,
    scope='Boundary-only gen05/gen06 comparison for three completed seeded runs; not 30-seed gate')
with (HERE/'gen05-gen06-scientific-comparison-02.json').open('x') as stream:json.dump(report,stream)
print(json.dumps(dict(passed=report['passed'],rows=[{k:v for k,v in r.items() if k!='materials'} for r in rows])))
if not report['passed']:raise SystemExit(1)
