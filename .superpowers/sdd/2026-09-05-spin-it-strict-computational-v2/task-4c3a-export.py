"""Persist the prescribed engineering candidate and actual unsuccessful run."""
from pathlib import Path
import sys

from test_faster_material_bridge import problem, COEFFICIENTS
from make_it_stand.spin import SpinEvaluatorConfig, SpinSolverConfig
from make_it_stand.verification.faster_material import (
    evaluate_and_seal_faster_candidate, solve_and_seal_faster)
from make_it_stand.verification.spin_bundle import save_spin_artifact, load_spin_artifact
from make_it_stand.io import write_json_atomic, write_text_atomic

root = Path(sys.argv[1])
p = problem()
evaluator = SpinEvaluatorConfig(derivative_method='DISCRETE_ELLIPSOID', icosphere_subdivisions=1)
solver = SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1', maximum_iterations=1,
                         maximum_restarts=2, maximum_function_evaluations_per_restart=1)
candidate = evaluate_and_seal_faster_candidate(p, COEFFICIENTS,
    evaluator_settings=evaluator, solver_settings=solver)
save_spin_artifact(candidate.artifact, root/'material')
reloaded = load_spin_artifact(root/'material', expected_problem_digest=p.input_digest)
assert reloaded.scientific_digest == candidate.artifact.scientific_digest
failed = solve_and_seal_faster(p, initial_quadrics=(COEFFICIENTS, COEFFICIENTS),
    evaluator_settings=evaluator, solver_settings=solver)
assert failed.artifact is None and failed.status != 'CONVERGED'
write_json_atomic(root/'candidate-run.json', candidate.identity)
write_json_atomic(root/'unsuccessful-solve.json', failed.identity)
summary = {
    'path': str(root.resolve()), 'problem_digest': p.input_digest,
    'candidate_status': candidate.status, 'candidate_scientific_digest': candidate.scientific_digest,
    'artifact_scientific_digest': reloaded.scientific_digest,
    'artifact_evidence_digest': reloaded.evidence_digest,
    'mass': reloaded.facts['physical']['mass'],
    'fixed_normalized_moments': candidate.identity['fixed_moments_input_normalized'],
    'maximum_normalized_moment_error': max(reloaded.facts['moment_comparison']['absolute_errors']),
    'native_solve_status': failed.status, 'native_solve_scientific_digest': failed.scientific_digest,
    'native_solve_event_hash': failed.identity['native_run']['event_hash'],
    'native_solve_restart_evaluations': failed.identity['native_run']['solve_restart_function_evaluations'],
}
write_json_atomic(root/'summary.json', summary)
write_text_atomic(root/'README.md', '''# Engineering material candidate — NOT_OPTIMIZED_CANDIDATE

The prescribed cube/shell fixture passed independent mesh/material validation.
It is not an optimized solution: target feasibility is UNKNOWN, project status
is PARTIAL_BLOCKED, and physical validation is NOT_IN_SCOPE.

`material/` is the existing exclusive sealed-artifact bundle, reloaded with the
separately recorded expected problem digest in `summary.json`. Its manifest
maps the original, fixed, design, retained, carved and complete material meshes
to NPY arrays. `candidate-run.json` is the separate immutable bridge scientific
identity as exported to JSON, not a new authenticated persistence format.
`unsuccessful-solve.json` retains the actual native events/restarts and failed
status, with no final material artifact. No historical output was overwritten.
''')
print(summary)
