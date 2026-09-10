# Task5i independent controller evidence

BASEa22de6b69b0dd0d550ad57157fbd9bac834cac60; worker owns production/builds. Controller baseline source-octree/cage-grid2/2 passed2.59s. Existing linked worktree verified, no AGENTS found in repository or ancestors, branch correct. Main tracked diffSHA93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e; this protects tracked-diff identity only, not an exhaustive untracked inventory. No main writes.

## Independent coordinate check

`task-5i-independent-planes.py` constructs exact Fraction(origin)+k*Fraction(side)/2**level and rounds once using Python binary64 conversion; it does not call native octree/grid or libm FMA. It independently reproduces the positive binade's7collapsed intervals while first interval and all coarse intervals are positive. Negative/aligned examples have0collapsed intervals; rejection there is explicitly conservative admission policy.

Two actual runs exit0 and compare byte-identically in `build/task-5i-controller-W6WCqr/planes-run{1,2}.json`. Sweep1400roots over signs, binades and scales;1020strict-bound admitted,380rejected; no admitted collapsed intervals. This finite sweep supplements, not replaces, the following sufficient-bound argument, and does not audit source-material geometry.

For a monotone root interval, maximum endpoint-magnitude ULP bounds the width of a round-to-nearest binary64 rounding cell throughout the interval, including a power-of-two transition. Two exact planes spaced strictly more than this width cannot round to the same value. The old equality case can put two planes at opposite tie boundaries of an even significand's rounding cell. Checking only the first finest interval misses interior ties. Conservative strict admission retains sparse processing without a root-wide enumeration and never changes input geometry.

Native actual-input RED/GREEN and final corrected bindings must be appended after worker execution and independently inspected; this mathematical probe alone does not grant the fix's acceptance.

```text
23e30e01d6855398a57ef82bc0354a258bdfd628d28c40742a395600b567b25f  task-5i-independent-planes.py
19549dce2510e3004ea65013283066e5d2a9efafb740e5c0d4d8a800f0cf65b0  build/task-5i-controller-W6WCqr/planes-run1.json
19549dce2510e3004ea65013283066e5d2a9efafb740e5c0d4d8a800f0cf65b0  build/task-5i-controller-W6WCqr/planes-run2.json
25b76f47ca180c1070c966ec269b57c823adc20f317ced34a43e0da3e2a8195d  libigl/include/igl/massmatrix.cpp
```

Separate read-only BBW audit: existing tetra massmatrix usesV/4 barycentric lumping, and current cage triangle faces plus fine-grid trilinear transfer require separate nodal and off-node boundary checks. Documented in docs/spin2014_bbw_contract_audit.md; no solver behavior or dependency changed. The interpolation counterexample is derived mathematically, not reported as an actual BBW solve.

## Actual native behavioral RED/GREEN inspected

Controller read raw `build/task-5i-evidence-20260910-0FF6Vb/red-ctest-01.log`: compiled test ran and failed1/1, positive root had10distinct planes/7collapsed intervals and reached late collapse (`early_numerical=0 late_collapse=1 completed=0`); two aligned cases completed. Thus the demonstrated upstream defect is inadequate early admission, not a claim that this tested collapsed material was successfully exported. The constructor's existing late cuboid check did catch this particular material.

Controller read `green-ctest-01.log`:1/1passed1.52s; all three cases nowearly_numerical=1 with no late collapse or completed unsupported tree. Positive actual collapse and conservative aligned rejection remain distinct. Ordinary source/config/partition printed identically before/after:

```text
d05f34b01ef35d581ec903f9e803deb12402faddf93a0b6f209a464d0bf4297a
fc5b062c289ae11c288f542f9bf193a0a43f656f81d42f9bc8ca1143d46b43ef
3c4531f7bab657f722244a74b82dad6bc08c63a23ad6747e0d0e54df8538565f
```

These are intermediate focused results; final whole-suite and review are pending and must not be inferred from them.

## Actual downstream grid compatibility

After main covering tests, while worker owned the separate FasterOFF build/test, controller compiled the unchanged diagnostic `task-5h-controller-probe.cpp` against the amended stable main archive (single compiler process, no CMake rebuild; coordinated with worker). Both new actual exports passed unchanged `task-5h-check-probe.py` independent geometry/field checks.1728cells/2197nodes/6084edges/218handles; maximum field error1.3877787807814457e-17. Both outputs byte-identical to each other AND the previous Task5h ordinary output. No old file overwritten. Compiler, both probe runs, checker and cmp exit0.

```text
c48ea709bc7c1aac919ec59d063f8d980acfcd25e9baa6509cd92d313c1c4580  build/strict-v2-spin/libmake_it_stand_core.a
54ec05f4851a404ceeb566400347b8a49374992b3284436acc6d7bf5277e4649  build/task-5i-controller-W6WCqr/grid-probe
b4193f3e510cc6f142765ec90e5e9c3f164d86e97e03cc8e9c65add3e7b1dd6e  build/task-5i-controller-W6WCqr/grid-run1.json (also run2)
```

Logs grid-compile.log/grid-check.log are in the same unique controller directory. Final source/archive binding recheck remains required after all-target verification.

## Final source and all-target verification

Final main all-target build exit0; full CTest201/201passed474.14s. Controller parsed201actual passing records, no duplicate/missing IDs, NotRun/Skipped/failure markers; raw capture and final summary agree. Main relevant integration8/8passed133.11s; FasterOFF8/8passed120.16s. FasterOFF cache explicitly2014ON/exactON/PythonON/FasterOFF, disconnected fetch settingsON. Final archive rechecked identical to the independent actual probe binding above. Main tracked diffSHA unchanged. This full result is fresh Task5i evidence, not Task5h's prior201pass.

```text
57618e3c3dc3a092dfd85c0a91d2c606419bb7827881644b0014415085e6a5cf  include/spin2014/SourceSurfaceOctree.h
96b842b65a3697e9c9455c385e57b403e2f503f6591b51272f20ff7f0df63b14  src/spin2014/SourceSurfaceOctree.cpp
562e584f850e57354d824d54a1aed891cd8f3fd572c3c6f6eb938c166e28226a  tests/spin_it_2014_source_surface_octree_regression.cpp
01ac9ff30de0812fc263dd420e906e4d5b148607b4eb55886c20bbc815415259  build/task-5i-evidence-20260910-0FF6Vb/full-build-01.log
b97b2c375046d4e028427c8195183d182f76243f84599c827b2b8b9311828b79  build/task-5i-evidence-20260910-0FF6Vb/full-ctest-01.log
```

Independent task review and actual fresh material replays follow separately; passing regression counts do not grant full BBW/original-case/physical verification.

## Fresh actual source-to-sealed-material replays

First controller launch omitted the CTest runtime's libstdc++ preload: actual Conda NumPy import selected an old runtime lackingGLIBCXX_3.4.32, so no optimization occurred. The failed attempt is preserved at `build/task-5i-controller-W6WCqr/material-run1/failure.json`; stderr appeared in the tool output, initial tee stdout log is empty. It is not counted as a successful scientific run. Controller read actual CTest ENVIRONMENT, then used the same resolved system libstdc++ path; no dependency install, source/fixture alteration or tolerance change.

Successful attempts are NEW directories `build/task-5i-controller-W6WCqr/material-run2` and `material-run3`. Each runner and subsequent independent checker exited0. Environment: `/home/lin/anaconda3/bin/python3`, PYTHONPATH=project/python:build/strict-v2-spin/python, LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6, OMP/OPENBLAS/MKL threads1; ulimit-v25165824; sequential timeout3600 runs. Each invocation uses unchanged `experiments/spin_it_reproduction/strict/run_source_spin2014_material.py`, exact fixed fixture/settings, main build directory, declared local libigl, main `_deps/highs-src`, installed pybind11/include. Declared roots are observed provenance, not complete compiler-discovered dependency closure or loaded-memory attestation.

Actual native status BINARY_MOMENT_VALIDATED; sealed=true, freshnessCURRENT. The independent `task-5i-check-material-replays.py` reloads/reintegrates both bundles and checks literal world moments, full inertia, outer bytes, actual1cavity,72cuboids,39ordered records and exact wall clearance25/64. Inputs/scientific/native-run bytes match pairwise. Source was frozen after final tests but still uncommitted while these runs were captured; committing the same effective source is contextual metadata only in the existing freshness checker, not a rewrite of begin/end observations.

```text
artifact scientific digest: 778f8dfb1dd0dc601152e2ffb78d5756efe8c42b1e61c1b3c9c8c03af9f0a97c
solver/material scientific digest: 780431ce4e1970da972488483a575b2b566a61fa072fdc6d09b65b5399c6e6a7
inputs.json: c415ea8f44faa23ee2c2e7e77ca3ed3ffbf8c8a77c56662893948097b1ba3978
scientific.json: 80cc6eb362b1f46cbcfc7e3094f69860cdb5f79dbf4b8dfdc05d6fa3e3d10227
native-run.json: fa943e91287b47254d38a8a2f0f0fb4aedfe990f681ab4d21cc30c22cdb50dd9
```

These match the prior ordinary synthetic fixture scientific outputs; provenance is fresh and separate. OverallPARTIAL_BLOCKED; solver_target_feasibilityUNKNOWN; geometry seal only, no original-case/BBW/CAD/physical-spin authority. No new OpenGL rendering or native viewer reload performed this turn; current mesh bundles are available, not a claimed desktop acceptance.

Postcommit check against `17f113ac84fcb5fa058f101b186f2accd3664de5` independently reloaded both existing bundles again (no optimization). Exit0, bothCURRENT; literal moments/inertia/cavity and byte comparisons still pass. Captured in `material-postcommit-check.log`. Source/header/test/core hashes match final-test bindings. This explicitly verifies that the contextual commit change did not stale these unchanged effective inputs.

```text
c3096f9dd1659d5dd7142c6ffe427b19ceb395abb8e81774be1d440ad69f95b6  build/strict-v2-spin/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so
5305322a8f58308c0c1c0782800a92b2c9cc0100b0649aefa56657b4e0a3f73d  build/task-5i-controller-W6WCqr/material-run1/failure.json
425d0984e333f259f8ff463f9eb04170490d400d70ff2715675ada8b9df84e5b  build/task-5i-controller-W6WCqr/material-run2/summary.json (also run3)
045829011a17cbc8d8718e771fccea44d293397fa34dc771c7d5d53214bd25f0  build/task-5i-controller-W6WCqr/material-check.log
```

Worker report also discloses an intermediate command-name error: nonexistent trailing target `make_it_stand_native` failed after five real targets built; actual `_make_it_stand` build corrected it. Controller inspected the preserved error in main-covering-build-01.log. This was not behavioral RED or a final build failure. Existing Boost global-placeholder deprecation note and prior shifted-root partition-coverage minor remain disclosed/deferred. No dependency refresh or hidden failed command.

## Task acceptance

Fresh independent reviewer `/root/review_spin2014_octree_precision` reports Spec compliant / qualityApproved, no Critical/Important/Minor findings. Controller resolved external CannotVerify items through the checks above. Implementation17f113a is accepted as a bounded numerical-admission correction; controller documentation/checkers are committed separately. No production changes followed the final201/201 run. Full reproduction remainsPARTIAL_BLOCKED and all unfulfilled algorithm/asset/experiment requirements remain present.
