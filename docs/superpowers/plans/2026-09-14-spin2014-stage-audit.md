# Spin2014 staged residual audit

Approved scope: the user's instruction to execute staged residual audit, Laplacian/refinement semantics inspection, and a fixed G6 control. This is numerical diagnosis, not a new solver or paper completion.

## Global Constraints

- Work only in the existing spin-it-strict-computational-v2 worktree; preserve all prior dirty code/results. No commits, cleanup, merges, or edits to the protected main checkout.
- No objective, rounding, initial state, resolution or tolerance changes. Gate remains 1e-6. Reconstructed teapot stays RECONSTRUCTED; failures are not infeasibility certificates.
- Python diagnostic algebra does not certify geometry or paper conformance. Missing exported geometry must be NOT_AVAILABLE, never passed.
- Tests first. All diagnostic arithmetic must reject nonfinite/invalid data. Fully exhaust and verify archived event streams before returning a result.
- Runs sequential, one thread, 24 GiB address-space cap, 300 seconds per G6 solve; preserve failures. Do not edit effective sources during a run.

## Task 1: Staged trace diagnostic

Add experiments/teapot_spin2014/stage_audit.py and tests/python/test_stage_audit.py; no native solver edits. Extend/reuse existing resolution_rounding_audit only if necessary, preserving semantics. Consume actual source-trace-archive-v1 header/pass/iteration records, bounded by largest record/pass (do not materialize full history). Inspect actual schema first.

Report each pass's post-optimization pre-transition material moments, post-split/merge material moments, constraint vectors at its actual final phi and input-fixed normalization, split/merge moment delta, final rounding delta and existing binary-check comparison. Last iterate beta must be taken from the recorded accepted terminal state, not guessed from arbitrary iteration input; check dimensions and consistency with solver termination and trace metadata. Distinguish pre-transition KKT evidence from terminal-partition KKT applicability. Unknown/missing states must be explicit, never fabricated.

Report actual fractional cell count/volume distribution and level summaries. Geometry export stage is NOT_AVAILABLE for nonadmitted results. No success-granting flags. Bind input archive SHA256 and checker/source digests, validate footer and sequence and reject corruption/truncation/inconsistent terminal state. CLI writes an exclusive new output and refuses overwrites. Test hand-derived asymmetric stencils, nonzero phi, density/scale, split preservation, changed terminal partition, rounding error, bad dimensions/NaN, invalid archive completion. Derive expected values independently.

## Task 2: Laplacian/refinement audit and fixed control

Controller reads complete native adjacency and functional path, checks beta polarity, 1/2 energy factor, symmetric graph weights, derivatives, refinement rebuilding and physical scale. Run existing relevant native tests and document precisely what is and is not established; do not assert a bug from resolution sensitivity alone. Run fixed online_trace_g6.json again using unchanged native binary and current provenance, then stage-audit new and prior archives. If new checker sources make old provenance stale, label historical analysis explicitly and do not rewrite old manifests. Preserve fresh source and binary hashes.

## Task 3: Review and report

Review Task 1 independently (spec and quality), address findings via original implementer, run final focused tests. Publish concise report with numerical stage table and next evidence-based hypothesis. Do not claim binary failure fixed unless measured. Broad review is scoped to this increment, not a retrospective approval of all dirty repository work.
