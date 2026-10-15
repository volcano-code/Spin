"""Compare every saved event, material and numerical result, excluding known costs."""
import hashlib
import json
from pathlib import Path

TIMING_FIELDS = frozenset(('evaluator_wall_seconds', 'residual_evaluation_seconds',
                          'derivative_evaluation_seconds'))


def scientific(value):
    if isinstance(value, dict):
        return {k: scientific(v) for k, v in value.items() if k not in TIMING_FIELDS}
    if isinstance(value, list): return [scientific(v) for v in value]
    return value


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def main():
    root = Path(__file__).parent
    checks = []
    for seed in (0, 1):
        paths = [root / name / f'seed-{seed}-job-01/output' for name in
                 ('small-top-inprocess-solver-sequence-01', 'small-top-inprocess-solver-repeat-01')]
        records = [[json.loads(line) for line in (p/'events.jsonl').read_text().splitlines()] for p in paths]
        events = [scientific(rows) for rows in records]
        materials = [json.loads((p/'last-material.json').read_text()) for p in paths]
        configs = [json.loads((p/'configuration.json').read_text()) for p in paths]
        results = [json.loads((p/'numerical-result.json').read_text())['result'] for p in paths]
        # Raw event_hash covers timing; exclude it explicitly here, not from scientific events.
        result_views = [{k: scientific(v) for k,v in r.items()
                         if k not in ('event_hash', 'events', 'evaluation_timings')} for r in results]
        check = dict(seed=seed, configuration_equal=configs[0]==configs[1],
                     events_equal=events[0]==events[1], event_count=len(events[0]),
                     scientific_event_hashes=[digest(e) for e in events],
                     material_equal=materials[0]==materials[1],
                     numerical_result_equal=result_views[0]==result_views[1],
                     statuses=[r['status'] for r in results])
        check['passed'] = all(check[k] for k in ('configuration_equal','events_equal',
                                                 'material_equal','numerical_result_equal'))
        checks.append(check)
    report = dict(checks=checks, passed=all(c['passed'] for c in checks),
                  scope='Small regular top pilot only; no geometry or physics admission inferred',
                  excluded_event_fields=sorted(TIMING_FIELDS))
    with (root/'small-top-repeat-comparison-01.json').open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps(report))
    if not report['passed']: raise SystemExit(1)


if __name__ == '__main__': main()
