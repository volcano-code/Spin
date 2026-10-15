"""Freeze measurement states from actual events, not arbitrary saved trials."""
import hashlib
import json
from pathlib import Path


def select_saved_states(run):
    run = Path(run)
    events = [json.loads(line) for line in (run/'events.jsonl').read_text().splitlines()]
    maps = {row['evaluation']: row for row in events if row.get('kind') == 'MAP_SAVED'}
    accepted = [row['event'] for row in events if row.get('kind') == 'NEWTON_EVENT'
                and row['event']['type'] == 'evaluation' and row['event']['outcome'] == 'accepted']
    initial = next(row for row in accepted if row['reason'] == 'initial_evaluation_succeeded')
    update = next(row for row in accepted if row['reason'] != 'initial_evaluation_succeeded')
    terminal = json.loads((run/'terminal-material-reference.json').read_text())
    final = next(row for row in maps.values() if row['evaluation_file'] == terminal['evaluation_file'])
    if final['sha256'] != terminal['sha256']:
        raise ValueError('TERMINAL_REFERENCE_CHANGED')
    result = []
    for role, row in (('initial', maps[initial['evaluation_index']]),
                      ('first_accepted', maps[update['evaluation_index']]), ('terminal', final)):
        name = row['evaluation_file']
        if Path(name).name != name:
            raise ValueError('NONLOCAL_SAVED_STATE')
        path = run/name
        if hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']:
            raise ValueError('SAVED_STATE_CHANGED')
        result.append(dict(role=role, evaluation=row['evaluation'], path=str(path), sha256=row['sha256']))
    return result
