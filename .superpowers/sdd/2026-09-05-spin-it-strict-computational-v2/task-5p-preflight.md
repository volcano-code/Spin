# Task5p preflight — shared source material surfaces

Approved architecture is full strict M0–M8 plan; this is its next already-approved material geometry increment. Implementation brief task-5p-brief.md; code not yet written at preflight. Existing linked worktree verified branchspin-it-strict-computational-v2, gitdir differscommon and no superproject; startingHEAD4c6462d. NoAGENTS found in scopedsearch. ProtectedmaintrackedSHA93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e unchanged. SixFasterpaths preserved. Baseline actual sourcebinding1/1PASS0.97s. /data1.4TiBavailable.

| Shared interface | Disposition |
|---|---|
| Binding Full primitives -> boundary geometry | Keep all primitives and beta owner order; cannot use bounding adaptive leaves. |
| Coarse/fine coincident faces -> triangle consumers | Sparse dyadic planar overlay stores one atomicpatch and one trianglefan, referenced with opposite orientation. |
| Perpendicular planes -> hanging edge vertices | Global integer axis-line lookup splits all incident edges; face-only matching insufficient. |
| Canonical keys -> patch centers | HalfstepD+1 sufficientprecisionguard, exactFMA; reject unsupported frames rather than perturbing points. |
| Original outer -> new material vertices | Preserve outer rows/triangles, reindex only internal materialpoints explicitly. |
| Conforming primitives -> future moments/SLQP | Closed/shared indexing necessary but not selfintersection/wallthickness/manifold/optimization certificate. |
| Resource plan -> topology growth | Sparse planequadtrees and linepoint maps bounded before insertions; explicit caps, no finestrootarray. |
| Task tests -> scientific claims | Actual source/fullBBW synthetic tests, independently enumerated smallunitfaces and tetrahedral moments; no originalcase substitution. |

Ruling: use sparse dyadicface overlay and centerfan as IMPLEMENTATION_DEFAULT — deterministic shared topology without uniform root expansion — costifwrong: downstream triangulation/derivative adaptation and reruns, not author-exact claim. No new physicsororiginalassetassumption. Ruling: this task closes topology prerequisite, not jointsolver; force actualdeformedpoint tests so it is more than metadata, retain all finalgeometrychecks separately. PreviousTask5o exactcounterexample demonstratesneed. No deletion/merging/release authority inferred.
