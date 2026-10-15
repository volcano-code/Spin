"""Produce an immutable, explicitly partial timing report from an actual log."""
import hashlib
import json
from pathlib import Path
import sys

from make_it_stand.verification.faster_timing import summarize_solver_timing


def main(source, destination):
    raw = (source/'events.jsonl').read_bytes()
    events = [json.loads(row) for row in raw.splitlines()]
    report = summarize_solver_timing(events)
    terminals = [row['event'] for row in events if row.get('kind') == 'NEWTON_EVENT'
                 and row.get('event', {}).get('type') == 'attempt_terminal']
    report.update(source=str(source), events_sha256=hashlib.sha256(raw).hexdigest(),
        numerical_convergence_recorded=bool(terminals and terminals[-1].get('reason') == 'converged'),
        original_run_status_unchanged=True, performance_gate_passed=False)
    with destination.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps({k: report[k] for k in ('solver_completed', 'solver_wall_ns',
        'measured_algorithm_ns', 'validation_ns', 'recording_ns', 'unclassified_ns',
        'paper_comparable_seconds')}, indent=2))


if __name__ == '__main__':
    main(Path(sys.argv[1]), Path(sys.argv[2]))
