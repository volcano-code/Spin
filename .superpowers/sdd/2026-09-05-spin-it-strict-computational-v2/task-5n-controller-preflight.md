# Task5n controller preflight and next source-binding boundary

Base99cfbf07b58967a8243108fee2f7590b2a021bd8. Existing linked worktree and branch verified, no submodule; main tracked diff hash93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e unchanged. Six unrelated Faster paths now exist, including new src/spin/QuadricMaterialExportGuard.h; all preserved. Approved prior Task5m is complete, not redispatched.

| Interface / requirement | Inspection and ruling |
|---|---|
| CMake registration -> dependency script | Pass invoking binary root explicitly; no implicit SOURCE_DIR/build cleanup. Missing destination fails before writes. |
| Per-run allocation -> repeated/concurrent runs | Allocate retained unique directories, never reuse or remove; allocation lock must release before expensive compile. |
| Test behavior -> project dependency assertion | Harness fixtures test routing/safety, but actual2014ON/FasterOFF linked build must still execute. No fake fixture results count as actual contract success. |
| New output paths -> evidence | Preserve logs and announce path before child work, so failures and interruption remain discoverable. Old interrupted log is immutable historical evidence. |
| CMake shared edits -> next native binding | One implementer/build owner at a time; source-binding implementation only after task5n review. Controller reads interfaces while task5n runs. |

Ruling: fix only the2014 destructive contract now, not unrelated Faster/reachability scripts. Cost if wrong: a future broad suite remains unsafe until those scripts are separately audited; no whole-suite safety claim. Scoped2014 runs remain the authorized verification target.

## Source-binding audit (read-only during Task5n)

SourceAdaptiveDomain leaves can aggregate multiple included Full source subcells; sourceBinaryCuboids(domain, ones, cap) traverses that actual conservative hierarchy and yields canonical boxes, not bounding-box substitutions. An immutable future binding can assign each emitted cell to exactly one adaptive ancestor leaf, preserve that ownership when weights are transferred, and verify primitive moment sums per leaf.

SourceCageBbw is an immutable actual-run object; use its real grid/field and check success. Its source digest is outer-as-design, whereas SourceSurfaceOctree may have a different design domain. Comparing those digest strings directly is incorrect: bind the exact original outer arrays and root frame, then retain both distinct upstream identities. A failed BBW or another model's field must not be accepted.

SourceCageGrid.embed provides checked canonical trilinear embeddings, with no snapping or fallback; cache sparse embeddings rather than constructing a dense global3N-by4H matrix. All handles remain present. Resource failures must be explicit and preserve resolution. This prerequisite alone does not establish deformed cell-face conformity, a final material partition, source-adaptive joint SLQP or full paper regularizers. Those gates remain required before successful joint optimization claims.

## Actual source binding preflight

Controller compiled a separate read-only probe against the unchanged FasterOFF archive after the worker confirmed no archive writes. The initial diagnostic compile used the wrong accessor name sourceHash; corrected to sourceSha256. This was a probe compile error, not a production regression or TDD RED. Final executable ran twice with identical output. Outer cube(.3,.7), design cube(.34,.66), sourceDepth5, initialOpt2/max5 produced 8 leaves and 496 Full primitives. Every primitive had exactly one adaptive ancestor owner. Independently accumulated separable box moments differed from each leaf's ten stored moments by at most 1.0842021724855044e-18. Total volume equaled (10/32)^3 =0.030517578125; full leaf bounding cubes would incorrectly yield0.125.

Synthetic regression only; no geometry optimization, full2014 convergence, or paper-case claim. Probe source is task-5n-controller-source-binding-probe.cpp. Executable and two retained logs are in build/task-5n-controller-zTtQ8D. Executable SHA2567a6c7cb564cd39bdbae66eb277b923f7baec0704566aacc58a12676a19d17b4b; each log SHA256e6c4d82552ff860baf46b9fc523980087db8cc4a557b6e7d4d16d032ef6971ea. Linked archive remained85a80297bbb7d82dd176df016f821757f11d3e6b1ff4012a304de9d35f9fdc24. Task5o brief records the scoped actual-material/full-BBW binding design; implementation waits for Task5n review.
