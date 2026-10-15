"""Read-only extraction of two original failures; no geometry/solver invocation.

Writes a new diagnostic registration, never alters the running job or its
dependencies. The input log prefix is hashed so later appends remain valid.
"""
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
run = HERE / 'noodles-registered-core-seed-0-job-01/output/run'
prefix = (run / 'mapping-events.jsonl').read_bytes()
if not prefix.endswith(b'\n'):
    prefix = prefix[:prefix.rfind(b'\n') + 1]
rows = [json.loads(line) for line in prefix.splitlines()]
config_path = run / 'mapping-config.json'
config = json.loads(config_path.read_text())
failures = {}
markers = {
    'source_side': 'outward quadric interface faces disagree on the interior sign',
    'material_boundary': 'quadric material Boolean is not a positive closed manifold',
}
for label, marker in markers.items():
    record = next(r for r in rows if r['status'] == 'FAILED' and marker in r['exception_message'])
    if record['mapping_identity'] != config['identity']:
        raise ValueError('FAILED_STATE_MAPPING_IDENTITY_MISMATCH')
    stages = [r['stage'] for r in rows if r['index'] == record['index'] and r['status'] == 'STAGE']
    failures[label] = dict(record=record, stages=stages)
if failures['material_boundary']['stages'][-1] != 'material-authority-audit-begin':
    raise ValueError('MATERIAL_FAILURE_STAGE_NOT_CONFIRMED')
if len(failures['material_boundary']['record']['coefficients']) != 10:
    raise ValueError('INCOMPLETE_FAILED_COEFFICIENTS')
report = dict(
    schema='noodles-original-failure-diagnostic-registration-v1',
    mapping_config=str(config_path), mapping_config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest(),
    mapping_log_prefix_bytes=len(prefix), mapping_log_prefix_sha256=hashlib.sha256(prefix).hexdigest(),
    mapping_identity=config['identity'], source_files=config['source_files'], failures=failures,
    geometry_invocations=0, optimization_invocations=0,
    hypothesis='Binary64 conversion/canonicalization may change material incidence; unconfirmed until exact output is captured.',
    competing_explanation='The original exact Boolean may itself be nonmanifold, or the authority audit may reject valid incidence.',
    next_diagnostic='After current job terminates, capture the exact output, birth faces and binary64 view for evaluation 7 with the SAME production source, side and cap branch. Compare original and canonicalized vertex links; do not relax the checker.',
    source_diagnostic='Independently preserve evaluation 3 side-selection failure; do not change its signs or pair roots by proximity.',
    frozen_scope='No new optimizer, geometry change, step scan, threshold relaxation or second concurrent heavy job.',
    admission=False)
output = Path(sys.argv[1])
output.mkdir(parents=True, exist_ok=False)
with (output / 'registration.json').open('x') as stream:
    json.dump(report, stream, indent=2, allow_nan=False)
print(json.dumps(dict(registered={k: v['record']['index'] for k, v in failures.items()}, geometry_invocations=0)))
