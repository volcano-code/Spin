# Controller handoff verification checklist

Task5b starts at `080b1f1d1f275cc873b1ef9317af0e0153e85d50`.
This is a checklist, not a record of passing tests.

1. Read implementer report and all relevant RED/GREEN evidence, including actual fixed boundary-release failure before audited behavior. Confirm all promised adversarial cases exist and were exercised; do not infer this from the first two passing core tests.
2. Match final source, native executable, Python extension, and CMake cache hashes against reported observations. Confirm settings/profile defaults and feature flags, including actual 2014 ON / reachability OFF linking. Full-suite counts must refer to the final implementation, not an earlier intermediate build.
3. Inspect scoped commit range from the recorded base and package it for fresh read-only task review. The review must return both specification and quality verdicts. Resolve external-evidence limitations separately; do not treat the reviewer as having re-run tests.
4. After sole build ownership returns, run focused public and native tests against final code with explicit new output logs. Independently recompute KKT components from exported state, gradient, Jacobian, constraints and multipliers; do not reuse the native residual evaluator as an independent oracle.
5. Check old default SLQP and end-to-end golden hashes, Task5a default formulation identity, and audited two-run scientific determinism. Times are excluded from scientific identity. Check final accepted state remains audited when iteration budget ends immediately after a step.
6. Check scientific scope: restoration failure/timeout/singular Jacobian are not infeasibility certificates; first-order KKT is not a local-minimum or physical stability certificate. A failed rounded material cannot inherit continuous-solve success.
7. Recheck dirty-main tracked diff hash `93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`; preserve historical outputs. No tag, release, merge or push.
8. If source-bound Faster Section6.6 results become stale from CMake changes, retain them unchanged and, after final code freeze, use the unchanged runner/spec for two new output directories. Do not reimplement or silently relabel the old proofs. Old DSO-bound material candidate is already stale and remains historical.
9. Record actual verification and remaining gaps in a separate final observation document. Overall reproduction remains `PARTIAL_BLOCKED`; this numerical increment cannot close full SLQP, source-surface sparse octree, original models/CAD, or original experimental matrices.
