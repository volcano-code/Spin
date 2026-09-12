# Task5q independent verification — 2026-09-12

BASE `6f61946cbb2b7d3329321f6741c4667ad616042f`.
Actual worktree: `/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`.
Existing `/data/xvxinyao/make-it-balance-paper-sim` alias still resolves here.
This task did not perform the external move. See `task-5q-relocation-build.md`
for interrupted regeneration, dependency re-creation and containment measures.

## Independent computation

Controller standalone exporter invokes actual source-octree, adaptive domain,
full cage/grid/FEM/BBW, binding, conforming boundary and the new moment API.
Two fixtures: baseline cube and nonsymmetric box transformed by scale1.7 and
translation(-.2,.15,.6), densities0.8/2.3. Both retain all26handles and8beta
variables. SourceD5, actual496/147primitive geometry, output2965/912points.
Every handle uses a distinct deterministic scale/translation pattern; beta is
nonuniform. Nothing is initialized from a precomputed optimizer solution.

The checker imports NumPy/stdlib only, not project/native code. It reconstructs
each rest point's embedding by coordinate cell lookup and manual trilinear
coefficients, then applies the exported actual nodal BBW field and all handles.
It integrates signed tetrahedra, not the solver's supplement surface-normal
polynomials. Complex-step on this independent integral checks all112columns
(8beta and104handle parameters) including mixed second moments. Geometry and
moment values from value-only/evaluate paths must also agree exactly.

| Fixture | Max normalized moment error | Max normalized complex-step Jacobian error | Point motion discrepancy |
|---|---:|---:|---:|
| Baseline | 3.953e-16 | 3.717e-16 | 0 |
| Transformed nonsymmetric | 5.583e-16 | 2.941e-16 | 0 |

Normalization uses density and the fixed rest input outer diagonal to powers
3/4/5 by moment degree. Derivative gate combines absolute1e-10 of those fixed
scales with relative1e-7 of expected entries. Complex-step epsilon1e-30.
These engineering checks supplement, not replace, the native three-step FD
entrywise1e-6/near-zero absolute1e-11 checks. No paper-authored tolerances claimed.

Both pipelines ran twice. Full exports are byte-identical, not merely close.
Twelve deliberately corrupted exports (betaJac, handleJac, movedpoints, moments,
value-only moments and mode flag in each fixture) all fail the independent
checker. This is synthetic regression, not original-paper experiment admission.

## Commands / artifacts

Helpers are `task-5q-controller-{export.cpp,check.py,negative.py}` here.
Retained data: `build/task-5q-controller-DhEfuX/{export,run1.txt,run2.txt}`.
Actual rechecks of preserved exports write `check1.log`, `check2.log` and
`negative.log`; they do not regenerate or rewrite historical solver outputs.

Exporter compiled with `/usr/bin/c++ -std=c++17 -O2 -g -DNDEBUG -frounding-math`,
2014/exact definitions ON, include paths include/tests/libigl Eigen, assigned
core archive, fixed available Boost/HiGHS/GMP/MPFR/zlib. Address-space limit
25165824KiB; thread environment OMP/OPENBLAS/MKL=1. Python executable:
`/home/lin/anaconda3/bin/python3`.

Run `export > run1.txt`, repeat to run2.txt, then run check.py on each file,
`cmp` the files, and negative.py on run1.txt. Each observed command exited0.

## SHA-256 at independent test

| Artifact | SHA-256 |
|---|---|
| SourceCageMomentFunctional.h | 0a4fc174adca7d9189c84558f42fda6c35d51920d466fe497b9e84c16e1d0148 |
| SourceCageMomentFunctional.cpp | 753ebe6061b5a7ff1de3cd738d346d7d6d7e049eff46991ae3a024473c3dc6b2 |
| SurfaceTriangleMoments.h | 7a54977ae05c541b5e20d524de3c3a1c64856aa6db091166f64c3f49b5a51f84 |
| SurfaceRawMomentDifferential.cpp | fb6c2d0835f1b4630507836a9d83ad909221012415db7f7e4bf7f7cd4a4bd220 |
| Core archive | 52754814994bfa24dfd84d335e9cfad1c64b9e330cb94510d94254a36a98eb94 |
| Exporter source | 422d3b28555bf709cd0232e290a145c1f08de48b2a14f901025c37421e4d06c0 |
| Checker | 6f486609def5eff96ad2a63114120ea416958c89f39b636c949ca9120d1dfc35 |
| Negative checker | e1b482967c959e0a2dfcde93d9435e141397d7e3d8a78935e1a8efca55c22524 |
| Linked exporter | 97652217dfbdf5cb8ab44a80b8fcf96add0d3629a76f595ef64cc57810b8517c |
| run1.txt and run2.txt | 3c4610f96af31682ea7073d1d269f13292b8a7ce335ee51d191c9cc6933e4617 |

Final scoped suite and fresh review are recorded below. This layer
does not grant a joint SLQP, paper shape regularizer, geometric legality,
sealed deformed material, original cases or full reproduction. PARTIAL_BLOCKED.

## Post-relocation rebind

Canonical configuration/full build completed. Controller recompiled a separate
`export-canonical` against the final core and rebuilt HiGHS, from the canonical
worktree. Its binary hash is identical to the table above. Two new actual
`canonical-run{1,2}.txt` outputs compare byte-for-byte to one another **and** the
pre-canonical pair, with the same recorded SHA. Independent checking of the
new export reproduces all reported moment/Jacobian/motion errors. HiGHS archive
SHA: `323bb3bdd8d2dee87feb814abd07f99b4570c5b4d330f4bbb041dd662f736117`.
This directly checks final linked outputs rather than assuming path relocation
is scientifically irrelevant.

The previously failing source Python evidence test passed1/1 in23.85s after
canonical generated import paths; controller read its complete retained log.
No evidence validator or numerical acceptance threshold was weakened.

## Initial pre-review scoped suite

Implementation commit `e6fc2d07549a9b96ca699a845670b14381b0dc2a`;
implementer report commit `e2499187b5d7470fbe4eb454c8d12d26a4b92419`.
Controller inspected the retained canonical suite: **35/35 PASS**, 237.47s,
with no Not Run or skipped entries. Log `task-5q-canonical-safe-suite.log`
SHA-256 `279d8f67221c9bcc9d6c5534cdda5c94dc3c8f85c2d88dbc1424929237b04f13`.
The final linked core archive still matches the table above.

Main project's tracked-diff SHA remains
`93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`.
This is a scoped tracked-change check, not proof that all external/untracked
data is unchanged. Concurrent Faster/Plane/Polyhedral/FullRank/Truncation work
is outside this review and remains unstaged by this task.

## Compatibility fix: independent rebind

Review found one exception-category compatibility defect; see `task-5q-review.md`.
Repair implementation `47f969538d4e7bd985bf17b3d5b3808df1d6dbba` and report
`82d7ad1cf7cdb38609f1ddf9e6566c8080e71d3c`.
After the bounded production fix, controller compiled `export-fix` against the
new core, without reusing an old linked executable:

- Core SHA `1bdd1c8eb43e474359d6cbf6732f0c5faf18768fb7d4939ccaf4eff3306e86f2`.
- Exporter SHA `83b5b87a7d8b1c2666d140e440fc8f8d26123e228eea0cc6c3c5cf798b6336ad`.
- Actual `fix-run1.txt` and `fix-run2.txt` both SHA
  `3c4610f96af31682ea7073d1d269f13292b8a7ce335ee51d191c9cc6933e4617`.

The two new exports compare exactly with each other and the earlier outputs.
Both actual independent checks reproduce the reported normalized moment and
112-column derivative errors; all12 corrupted exports again fail. Actual logs
in `build/task-5q-controller-DhEfuX/` are `fix-check1.log`, `fix-check2.log` and
`fix-negative.log`. Commands exited0; exclusive output creation preserved old
evidence. Numeric output identity is observed, not assumed from a narrow diff.

The repair suite's source Python test correctly detected a concurrent external
source change. Controller parsed the actual retained
`build/task-5f-runner-test-gtpjptst/success/evidence.json`: begin source SHA
`95bdc72c6326bdfb551cb49ee2c8e4947b983e02e8c67b6e893a150c067e4d26`, end SHA
`33a5f7897ea8f33a41dc1f76dbd00ae94310ceb4e648a74b64eaa6d98a1ee60a`.
The only differing source entries were external
`experiments/m1_general_differential/CMakeLists.txt` and
`tests/polyhedral_moment_differential_regression.cpp`. The original
`STALE_RESULT` failure remains evidence; do not label that suite all-green or
weaken the validator. A separately identified retry is required for recovery.

The repair suite completed **34/35 PASS**, 228.47s, with that sole failure.
One separately recorded retry then passed **1/1**, 24.06s, log
`task-5q-fix-python-recovery.log` SHA
`dba0de5ccaa3327b8d8a292781bc9161f4c2ba9c6ec7e5b3389ef496dae57453`.
Controller read the full log and the retained freshness comparison. New actual
`build/task-5f-runner-test-9nio4som/success/evidence.json` has begin=end source SHA
`a29cf9639d9993ea279c34fe2b0c699efe002067251d3f1dc1a49a760d9aa5b2`, no differing
source files. All35 scoped tests therefore have passing executions against
the repaired owned implementation, but no single all-green repair-suite claim
is made, and no full-project/foreign-source build completeness is claimed.

Final changed kernel header SHA
`ab069c6de25befd41da88190a468528844c07383bafb5001dd562db9b92ffea9`;
whole-surface implementation SHA
`6b65b0dcde17904039588a1f8ffd0d223d6dd88ad62c2307875f4965b20e2989`.
The source functional/header and controller helper hashes are unchanged.

## Acceptance

Fresh task review required one compatibility repair. The original reviewer
approved the scoped repair through `82d7ad1`: all findings addressed, no new
Critical/Important breakage. Scope: actual-source/full-BBW analytic moment
layer only. Deferred canceled-motion hash-test polish and inherited build
deprecation warnings are not hidden completion requirements.

Scoped source/helper/document whitespace checks passed. Raw retained build/RED
logs contain tool-produced trailing spaces; their bytes and hashes are kept
unchanged, so an unfiltered log-inclusive whitespace check is not claimed clean.
Only explicit task records, helpers, logs and the public gap report are staged.
The cumulative ledger, builds and all foreign edits remain outside this commit.
Global status remains PARTIAL_BLOCKED; next algorithmic work is the real-source
joint objective/regularizer/SLQP and geometric precommit chain, followed by
remaining Faster and original-paper gates in the strict plan.
