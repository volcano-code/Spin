# Task5m controller verification

Base: f78eac7aed4a8d47af99f3412092dfda2e4f5e42, existing isolated branch spin-it-strict-computational-v2. The user approved continuing joint deformation. Skills used: executing-plans, subagent-driven-development, TDD and verification-before-completion. Fresh implementer /root/spin2014_cage_explicit_functional owns the scoped Faster-OFF build; controller only writes these independent check artifacts during implementation.

## Preflight and scope ruling

Ruling: first close the explicit Givens derivative gap in the existing mesh functional with actual full-cage BBW integration coverage, without prematurely identifying this adapter as the complete source-adaptive joint solver. Cost if this ordering is wrong: some mesh-adapter plumbing may need replacement by the later sparse source bridge. No original-case or geometry gate is lowered.

The existing SourceAdaptiveDomain leaf is NOT necessarily its entire bounding cuboid. A boundary leaf aggregates only included Full source subcells. A future adapter must carry those actual cells/material boundaries and cannot deform a bounding cube and retain the old leaf moment claim. Shared faces across different cell levels also need a conforming deformation representation; an arbitrary per-box triangle export is not an established final material partition.

| Interface | Current fact and next obligation |
|---|---|
| Givens -> cage moments | New separate evaluation must chain beta, explicit phi, all existing free handle columns; old eliminated evaluation remains compatible. |
| Source cage -> actual BBW -> mesh weights | Full source cage plus actual finer-grid solve and trilinear transfer can be exercised now; caller-provided meshes are not immutable source partition binding. |
| Regularizers | Local authenticated Spin2014 text section6.1 specifies absolute deformed-position bi-Laplacian plus displacement. Existing uniform graph is a stated discrete implementation, not full recovered Sorkine rotation-aware operator. |
| General audited SLQP -> joint functional | EqualityBoxProblem supports bounded beta prefix and unbounded trailing variables, but no production joint bridge is added in Task5m. |
| SLQP -> geometry rejection | Current audited line searches catch callback exceptions and reduce step, but generic INVALID_EVALUATION_OR_CALLBACK_EXCEPTION does not preserve geometric audit facts and catches resource exceptions too. Future joint bridge needs specific failures, resource termination, trial hashes and accepted geometry evidence. |
| Old cage solver -> geometry | Old checks each cell vs outer but does not establish the complete deformed cell partition; do not transfer its geometryValid status to the future strict path. |

## Independent checker design

`task-5m-controller-export.cpp` runs actual original-mesh cage/FEM/BBW, retains all 26 handles, passes actual trilinear fields to the new functional, and exports all scientific inputs plus five joint evaluations. This is SYNTHETIC_REGRESSION, unit density, not an original paper case or completed optimization.

`task-5m-controller-check.py` imports no project module. It independently rebuilds LBS from exported pivots/rotations/weights/state, integrates oriented tetrahedra, constructs the inertia matrix and Givens rotation, and derives every state-column derivative through complex-step arithmetic. It also checks central differences at three steps. Exported native geometry is compared with the independently deformed geometry; native moments and all objective/constraint derivatives are not assumed correct.

Checker self-test executed before reading native output: literal ten moments of tetrahedron (0,0,0),(2,0,0),(0,3,0),(0,0,4); winding reversal; ten complex-step scale derivatives at1.7. All passed. Exact command/result is present in controller tool history; no fabricated log file.

## Protected work

Main checkout tracked-diff SHA at preflight remains 93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e. This check is not an exhaustive untracked inventory. Existing unrelated Faster edits remain untouched. Build/strict-v2-spin is not rebuilt; its old whole-suite pass is historical, not current verification.

## Actual final results

Independent verification completed against the stable implementation archive while the worker ran scoped regression tests. Full scoped test outcome and independent task review are recorded below when available; overall PARTIAL_BLOCKED.

Output directory: `build/task-5m-controller-pBCwYb/`.

- Standalone C++ build exit0, empty compiler log. Linked only Faster-OFF `libmake_it_stand_core.a` and the pinned existing HiGHS library, with exact-geometry/2014 macros and `-frounding-math`. Virtual-address-space cap24GiB; standalone compilation was the fourth job alongside worker-j3 test builds, after its archive was declared stable.
- Two actual export invocations completed with valid JSON and empty stderr. `cmp export-1.json export-2.json` exit0: all exported scientific inputs, actual BBW hash and five evaluations are byte-identical. No timing or timestamps in this scientific export.
- Actual original nonsymmetric box cage has26 handles; two internal cuboids have separate beta values; several handles carry fixed nonidentity rotations and every handle has nonidentity scale/translation. State count107 for rotating modes,106 for StaticLowering.
- Independent checker exit0 for all five objectives. Max absolute discrepancies: outer geometry3.3306690738754696e-16; raw moments5.551115123125783e-17; raw-moment Jacobians2.2551405187698492e-17; joint objective/constraint derivatives4.787836793695988e-16; objective/constraint values7.216449660063518e-16.
- Independently recomputed central differences at1e-3,1e-4,1e-5 show the expected decreasing error before floating-point noise. For rotating modes the largest error decreases about1.36e-6 ->1.36e-8 ->1.5e-10; StaticLowering6.21e-9 ->6.27e-11 ->1.31e-11. This is engineering derivative verification, not Newton convergence evidence.
- `bindings.sha256` binds the header, implementation, actual linked archive, exporter executable, both helper sources, both actual exports and independent check log. All9 entries were recomputed and matched after the run. This scoped binding is not a claim that the entire M2 transitive runtime provenance gate is complete.

Invocation details (from the repository root):

```bash
ulimit -v 25165824
/usr/bin/c++ -std=c++17 -O2 -g -DNDEBUG -frounding-math -DMAKE_IT_STAND_HAS_SPIN_IT_2014=1 -DMAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS=1 -Iinclude -Itests -isystem libigl/external/eigen .superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5m-controller-export.cpp -o build/task-5m-controller-pBCwYb/export build/task-5c-no-faster-20260909/libmake_it_stand_core.a /home/lin/anaconda3/lib/libboost_system.a /usr/lib/x86_64-linux-gnu/libgmp.so /usr/lib/x86_64-linux-gnu/libmpfr.so /home/lin/anaconda3/lib/libboost_thread.a build/task-5c-no-faster-20260909/lib/libhighs.a /usr/lib/x86_64-linux-gnu/libz.so
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 build/task-5m-controller-pBCwYb/export
# Run twice into separate export-1.json/export-2.json files, then cmp them.
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 60 /home/lin/anaconda3/bin/python3 .superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5m-controller-check.py build/task-5m-controller-pBCwYb/export-1.json
sha256sum -c build/task-5m-controller-pBCwYb/bindings.sha256
```

### Test-oracle correction during implementation

The first post-implementation focused run passed actual source-BBW integration but failed the core legacy comparison. Diagnosis identified YoYo paper15.6941 vs legacy8.42843 at off-centroid state. Controller independently checked SpinMomentProblem.cpp centralInertia and the authenticated paper coordinate/constraint formulas: old central inertia and explicit raw-moment extension are equal only on required centroid equalities. Ruling: preserve the original off-manifold independent-value/derivative fixtures, and add a separate centered posed-geometry fixture solely for diagonal-phi legacy value equivalence. Full ambient gradients need not equal even on the manifold; do not require that false identity. No tolerance relaxation or production algorithm change follows from this oracle correction. Cost if wrong: legacy comparison coverage would need correction; direct independent mathematical tests remain unaffected.

### Scoped test harness side effect and containment

The31-test scoped CTest invocation passed the first30 tests, then entered existing `tests/cmake/assert_spin_it_2014_dependency.cmake`. That script unconditionally removes and recreates `build/spin-2014-disabled-contract` and `build/spin-2014-enabled-contract` under SOURCE_DIR, violating this task's intended single-build-directory restriction indirectly. The worker noticed during its nested build and terminated the validated CTest process group3680723 with TERM; result143, no children remained. Controller read the entire script and verified both paths/source settings point within this linked worktree. The removed data were generated dependency-build caches, not main checkout or experiment assets; rebuilding regenerates their products, but old cache contents were not backed up.

Ruling: do not rerun the destructive script. Permit only completion of the already configured enabled-contract build and its actual audited executable as a separately logged dependency check, using the same24GiB/-j3 limits and no further deletion. Preserve the interrupted run and count it honestly:30 actual passing CTest cases plus1 interrupted case, not31/31. Cost if wrong: the separate dependency result would need a fully isolated rerun; this does not grant full-suite or release success. Mandatory followup: replace this preexisting fixed-path destructive test workspace with isolated per-run retained scratch before another broad suite/release. User was explicitly informed of the removed/recreated build caches and interruption.

The permitted existing-build completion and actual audited executable subsequently exited0. Output includes `[PASS] audited actual KKT/restoration/analytic/Givens tests`, with Givens KKT6.3137399822221774e-11; this is an existing audited solver fixture, not the new joint optimization. No exact-OFF joint-functional regression execution is inferred from it. Raw logs are `build/task-5m-logs-9xUnJr/dependency-existing-{build,test}.log`.

Controller parsed the scoped raw log and verified30 unique passing test IDs; test42 only started, with no PASS or100%success statement. All16 file/hash entries from the worker report were recomputed and matched, as were all9 independent check bindings. Implementation commit ce733ac915f322d939e887f04dea40fc5e4a1f58 uses authorized volcano identity. Only the original five unrelated Quadric paths remain dirty; protected main tracked-diff SHA still matches preflight. Fresh reviewer /root/review_spin2014_cage_explicit_functional is reviewing the complete f78eac7..ce733ac scoped diff.

Final review returned SpecPASS/QualityApproved with no correction required. Exact returned scope and controller disposition are saved in task-5m-review.md. No source changes followed the reviewed commit. Actual output digest for both exports:53a4c12dd246e5b4c634ff6538eb851d06491eeb655b7ce685587029d27aa993. Independent check log:23523bee62976eacdb2eec8780616c7dad3061bb6525a71a8147831edca39da5. Export executable:644869042a769ab013d3db0243058942d1625f168877d2138189407f0239afeb. Full bounded task accepted; all original full-plan gates remain outstanding as documented. No merge, release, tag, main-tree integration, physical validation or original-case claim.
