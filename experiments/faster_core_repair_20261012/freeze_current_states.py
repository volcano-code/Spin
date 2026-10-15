"""Freeze saved failure/accepted states as diagnostic replays, never starts."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
heart=ROOT/'experiments/faster_heart_contact_20261011'
run=heart/'rounded-clearance-v3-core-faster-v3/seed-0-job-01/output/run'
events_path=run/'newton-events.jsonl'
mapping_path=run/'mapping-events.jsonl'
events=[json.loads(row) for row in events_path.read_text().splitlines()]
mapping=[json.loads(row) for row in mapping_path.read_text().splitlines()]
accepted=[row for row in events if row.get('outcome')=='accepted' and (row.get('step') or 0)>0]
failures=[row for row in mapping if row.get('status')=='FAILED']
selected={}
if accepted:selected['heart-last-accepted-update']=accepted[-1]
for label,needle in (('heart-first-cap-failure','truncation'),('heart-first-small-face-rejection','degenerate_faces=')):
    row=next((item for item in failures if needle in item.get('exception_message','')),None)
    if row is not None:selected[label]=row
packet=dict(schema='faster-core-saved-state-set-v1',role='DIAGNOSTIC_ONLY_NOT_OPTIMIZER_INITIALIZATION',
    records=selected,source_files={str(path):hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (events_path,mapping_path,heart/'rounded-clearance-v3-faster-frozen-v2/request.json')},
    heart_job_status='RESOURCE_LIMIT_WALL',optimization_admission=False)
with (HERE/'frozen-heart-states.json').open('x') as stream:
    json.dump(packet,stream,allow_nan=False)
print({name:row.get('evaluation_index',row.get('index')) for name,row in selected.items()})
