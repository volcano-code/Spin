"""Registered projective variant: same-input fixed gates precede derivatives."""
import json
from pathlib import Path
import sys

from verify_formal_fixed_states import main as verify_fixed
from probe_formal_fixed_derivatives import main as probe_derivatives

HERE = Path(__file__).resolve().parent


def main(output):
    output.mkdir(parents=True, exist_ok=False)
    frozen = json.loads((HERE / 'same-input-job-03/output/frozen-inputs.json').read_text())
    frozen['settings']['paper_projective_hyperboloid'] = True
    (output / 'registered-inputs.json').write_text(json.dumps(frozen, allow_nan=False))
    verify_fixed(output / 'fixed', frozen)
    report = json.loads((output / 'fixed/report.json').read_text())
    failures = []
    for row in report['rows']:
        relation = row.get('exact_relation', {})
        if (row.get('failure') or not row.get('cgal_derivative_defined')
                or not row.get('exact_cgal_valid') or not row.get('exact_cherchi_valid')
                or not relation.get('equal') or row.get('moment_error', float('inf')) > 1e-10
                or row.get('jacobian_error') is None or row['jacobian_error'] > 1e-6):
            failures.append(row['state'])
    if failures:
        (output / 'report.json').write_text(json.dumps(dict(
            status='FIXED_STATE_GATE_FAILED', failures=failures, core_admission=False)))
        raise SystemExit(1)
    probe_derivatives(output / 'derivatives', frozen)
    (output / 'report.json').write_text(json.dumps(dict(
        status='PROJECTIVE_NUMERICAL_CHECKS_RETURNED',
        fixed_states_passed=True, material_branch_certification=False, core_admission=False)))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
