"""Exact same-version comparison of real saved runs, retaining scientific identity."""
import json
from pathlib import Path
import sys
from make_it_stand.verification.cherchi_run_evidence import verify_cherchi_run_evidence, scientific_projection


def main(first, second, output):
    roots = [first/'output', second/'output']
    verified = [verify_cherchi_run_evidence(p/'run') for p in roots]
    def load(p, name): return json.loads((p/name).read_text())
    configs = [load(p, 'run/configuration.json') for p in roots]
    events = [[scientific_projection(json.loads(line)['event'])
        for line in (p/'run/events.jsonl').read_text().splitlines()
        if json.loads(line)['kind']=='NEWTON_EVENT'] for p in roots]
    materials = [{f.name:json.loads(f.read_text())
        for f in sorted((p/'run').glob('evaluation-*.json'))} for p in roots]
    terminal = [load(p, 'run/terminal-material-reference.json') for p in roots]
    identities = [load(p, 'run/result.json')['identity'] for p in roots]
    for identity in identities:
        for field in ('event_hash', 'events', 'evaluation_timings'):
            identity['native_run'].pop(field, None)
    reports = [load(p, 'report.json') for p in roots]
    checks = dict(configuration_equal=configs[0]==configs[1],
        scientific_events_equal=events[0]==events[1],
        all_same_call_materials_equal=materials[0]==materials[1],
        terminal_references_equal=terminal[0]==terminal[1],
        result_identity_equal=scientific_projection(identities[0])==scientific_projection(identities[1]),
        sealed_material_equal=reports[0]['material_digest']==reports[1]['material_digest'],
        both_material_closures=all(r['material_closed_loop_passed'] for r in reports),
        both_record_checks=all(r['record_integrity_valid'] for r in verified))
    report=dict(checks=checks, passed=all(checks.values()), same_version_repeat=True,
        core_admission=False, scope='Actual registered small seed0 independent repeat; not 30-seed gate')
    with output.open('x') as stream: json.dump(report, stream, sort_keys=True)
    print(json.dumps(report), flush=True)
    if not report['passed']: raise SystemExit(1)


if __name__=='__main__': main(*map(Path, sys.argv[1:]))
