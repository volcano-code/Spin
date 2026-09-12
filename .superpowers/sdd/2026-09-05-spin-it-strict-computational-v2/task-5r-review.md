# Task5r independent task review

Reviewer: `review_spin2014_source_cage_paper`, fresh read-only task-scoped
specification and quality seat. BASE `86234ff`, reviewed HEAD `a603876`.
Review package `review-86234ff..a603876.diff`; same `task-5r-brief.md` and
implementer report, including the numerical verification addendum.

## Initial verdict: NEEDS_FIXES

No Critical issue. One Important resource-contract violation:

- `SourceCagePaperFunctional.cpp:158–159` transforms dynamic Jacobians through
  Eigen GEMM, whose packed buffers are omitted from the declared derivative cap.
  Bundled Eigen threshold8 selects the GEMM path; GeneralMatrixMatrix.h:159–163
  allocates `kc*mc` and `kc*nc` buffers. In the actual static fixture B8,H26,N112,
  cap1464 is exceeded by raw-output/fill1128 plus packed-RHS1040 alone, even
  before packed LHS. Not-yet-allocated gradient/constraint outputs do not cover
  this peak. Fix by constant-size per-column transformation or an explicit
  correct scratch bound; strengthen the capacity regression.

Controller checked the production expressions, actual contract and reviewer's
bundled-Eigen dispatch evidence before assigning repair to the original worker.
No formula or tolerance change is needed for this resource defect.

Minor: inherited libigl CMake deprecation warning in `task-5r-red.log:2`.
Keep the original diagnostic; it does not invalidate the expected missing-header
RED and does not authorize dependency edits. Deferred as inherited noise.

## Checks and strengths

- Actual source snapshot, fixed frame, full handles, optional nonstatic angle,
  existing Givens objectives, sparse fill and scientific limits match the brief.
- All30 long-double five-point comparisons passed unchanged thresholds;
  ordinary double failures remained explicit. Controller independent checks
  reconstruct geometry and inertia without production moment/Givens formulas.
- All10 fixture/mode combinations, two identical exports and20 corruption
  rejections corroborated. Reviewer performed read-only SHA checks.
- Focused upstream checks resolved finite Givens outputs, immutable source
  snapshot and sparse const-reference Laplacian risks.
- Inspected build/suite evidence confirms36/36 scoped tests; reviewer did not
  rerun the suite. No full-project or full-paper completion verdict is granted.

## Fix round1

Original worker implemented `c79b9ed`, report `77e4f38`; original reviewer
performed scoped re-review using `review-a603876..77e4f38.diff`.

Final verdict: **All findings addressed, no new Critical/Important breakage**.

- Uncounted packing ADDRESSED: both production Jacobian paths call the fixed
  ten-vector helper (cpp164–165; header75–81), writing preallocated columns.
- Allocation regression ADDRESSED: actual former GEMM is rejected and the
  production-used helper passes scalar expectations under the guard at both
  sizes (test17–46). Retained RED/GREEN logs demonstrate the repair.
- Inherited deprecation warning documented and retained (report385–388);
  dependency changes were not required.
- Probe macro settings are target-local; registration retains2014/exact gate.
  Symbols inspected show no upstream evaluator implementation in the probe.
- Fresh37/37 suite, allLD5/complex-step assertions, retained double diagnostics,
  both controller10-case checks and20 negative checks corroborated.
- Reviewer read-only checks verified identical exports, all19 repair-manifest
  entries, and absence from the actual exact-OFF build graph. No tests rerun.
- No new breakage or out-of-scope observation. This closes Task5r only, not
  the full joint solver, geometry admission or original-paper reproduction.
