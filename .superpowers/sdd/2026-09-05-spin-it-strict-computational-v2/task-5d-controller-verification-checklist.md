# Controller Task5d acceptance checklist

Final: bounded Task5d accepted after fixround1/5 and fresh independent review. The entries below preserve the planned checks and intermediate observations; final results and remaining limits are in task-5d-controller-verification.md and task-5d-independent-review.md. No full-paper completion grant.

## Intermediate observed evidence (not final acceptance)

- Read `build/task-5d-red-availability-01.log`: missing source-native integration failed before production implementation. This is an API-availability RED, not a numerical-behavior RED claim.
- Read first linked-build failure `build/task-5d-build-focused-01.log`: bundled Eigen lacks range iteration on this vector type. Worker reports indexed iteration fix and preserves the green second build separately; controller has not independently rebuilt this intermediate state.
- Read `build/task-5d-behavior-01.log`: native fixed fixture made an accepted step from (.5,.5) to (.25,.25), then four audited solves overall. The second solve reduced stationarity from about1.52662e-5 to9.29344e-10. Subsequent changed-partition solves reset history. Final fractions stay near.25 and status is BINARY_VALIDATION_FAILED, consistent with the independent derivation. Full partition levels and final residual still require final report/test verification.
- These intermediate runs do not replace final code-bound verification, a fresh task review, or original-paper experiments.
- Read `build/task-5d-green-geometry9-01.log`: actual source geometry depth9 has2,097,993 nodes; coarse optimization partition8leaves,48face keys,RSS351580KiB. This observes sparse graph storage for a coarse partition of a deep geometry tree, NOT complete level9 optimization. `/data` still1.4TiBfree; main tracked-diff SHA remains93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e.
- Actual freshness call on both `build/task-5c-refresh-faster-section-6-6-run-{1,2}/evidence.json` now returns STALE_RESULT: cmake_configuration content binding changed. Preserve them unchanged; regenerate two new Task5d runs only after code/review freeze. This does not imply the unchanged Section6.6 mathematics became false, or that new Faster functionality was implemented.

- Read full report and actual preserved RED/GREEN/final-suite logs. Do not equate compile-only API absence with numerical behavior RED. Fixed quarter-fill example must have real accepted native steps and failed binary outcome, not only mocked control flow.
- Verify new functional uses original outer moments and actual accepted design contributions; partial-cell binary output must list true accepted cuboids, not whole optimization bounding boxes. Independent tetrahedron checker must not call source moment columns.
- Compare fixed-case final normalized residual with exact1/131072 if all rounded fractions are0. Compare relaxed solver against its own actual KKT gate, not controller SciPy trajectory. Cases with insufficient convergence retain failure; no parameter/input tuning to force pass.
- Inspect actual partition-level progression, remap/merge counts, and history restart when partition changes. Empty/zero carve domain, resource budgets, nonfinite/invalid states cannot return ordinary valid empty results.
- Check post-round material moment disagreement/constraint failure means explicit BINARY_VALIDATION_FAILED, not infeasible or sealed. Successful cases limited to numerical/cuboid material contract.
- Check fixed-step gradient/constraint finite differences, face graph exhaustive small oracle, shifted-root query coverage and recursive mass preservation.
- Independently match source/binary/cache SHA and main tracked-diff identity; no historical experiment changes. Final native focused rerun after implementation releases build ownership. Full suite not duplicated without code changes.
- Fresh task reviewer over BASE777c8f7..HEAD package with brief/report. Preserve findings and resolve all Important/Critical through original implementer + scoped re-review. Carry remaining minors to next task/final whole-branch review.
- Refresh unchanged FasterSection6.6 proofs only after review/code settles, in new unique Task5d directories; leave prior scientific outputs intact. No newly implemented Faster capability inferred.
- Final bounded handoff: clear original geometry -> continuous solve -> binary moment/constraint distinction; no original cases, manifold seal, physics, CAD or full reproduction claim.
