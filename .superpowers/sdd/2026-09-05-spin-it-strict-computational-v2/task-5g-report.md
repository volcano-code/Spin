# Task5g implementation report

Status: DONE. Base `13afcf0f0036c38548c156eb7d0e7284ff670897`, branch `spin-it-strict-computational-v2`.

## Implemented scope and interface

Added native immutable `spin2014::SourceSurfaceCage` under the conjunction of 2014 and exact geometry, independently of Faster. Constructor takes source V/F and an explicit config with root origin/side, cage level (0..20), node, selected-cell, boundary-vertex/handle and triangle budgets. It copies V/F and config, and publishes const references to source, selected integer cell keys, indexed boundary V/F, boundary integer vertex keys, identity handle-to-boundary-vertex mapping, exact boundary audit, hashes and stored source-node count. `shared_ptr<const Data>` owns immutable results. Legacy `AutomaticCageBinding`, all original source-tree algorithms/hashes and existing fixtures are unchanged.

This is complete cage **geometry**, not BBW, a finer interpolation domain, joint optimization, physics, feasibility, an original-paper experiment or a sealed artifact. Analytic moment tests integrate the cage envelope solely as external validation; they never substitute for original material mass.

## Algorithm and scientific conventions

1. Validate finite representable resource budgets; construct the unchanged exact `SourceSurfaceOctree` with the original signed mesh as BOTH outer and design, depth equal to cage level, and explicit source-node guard.
2. Scan the stored sparse hierarchy. Boundary terminals contribute one cell. Full terminals contribute `2^(3*(cageLevel-level))` descendants. Count against remaining selected-cell budget in 64-bit arithmetic before any Full expansion or selected-key reserve. Depth<=20 bounds the shift by 60. Empty/internal nodes are never expanded for extraction.
3. Generate the accepted Full descendants and finest Boundary keys, sort lexicographically (x,y,z), and remove shared INTERNAL faces by integer-neighbor lookup. No root lattice scan or dense lattice allocation.
4. Visit all remaining quads in sorted-cell, -x/+x/-y/+y/-z/+z face order, with x-bit-first local corners and fixed corner0--corner2 diagonals. Count/check triangle and distinct vertex budgets before allocations. Shared integer corners have one indexed vertex. Canonical represented coordinates come exclusively from `SourceSurfaceOctree::cellBounds`, including its FMA dyadic planes.
5. Emit outward material-union triangles. Every retained corner, including intermediate coplanar corners, is a handle. All topological boundary components remain, including inward-oriented bounded cavities and disconnected islands. No hull, smoothing, filling, repair, culling, duplicate-vertex workaround or sampling occurs.
6. Run `auditExactTriangleBoundary(..., AnyNonzero)` independently on actual indexed V/F. AnyNonzero admits negative cavity shells; direction relative to material is fixed by exposed-cell face winding. Closure, edge/vertex manifoldness, degeneracy, consistent component orientation and exact self-intersections must all pass. Invalid products throw `SourceSurfaceCageGeometryError` carrying the immutable actual audit; no usable cage escapes.

Cell contact is explicitly `IMPLEMENTATION_DEFAULT`: CLOSED cells that touch original triangles count, including face/edge/corner zero-volume contact. This conservatively envelopes source material. The paper does not prescribe these exact tie rules. The convention and retention of all boundary components are in public documentation and config identity.

Resource limits use `length_error` with `RESOURCE_LIMIT`, separate from `INVALID_GEOMETRY` and never infeasibility. The exact auditor's signed-int face capacity (INT_MAX/2) is included in pre-extraction representability validation. No resource failure downgrades depth or emits fewer handles. Canonical SHA-256 encoding uses fixed eight-byte big-endian integers, length-prefixed strings, row/field order and binary64 bit patterns, never native padding or native endian. Source identity reuses the outer-as-design source-tree hash. Configuration identity includes source-oracle config identity, explicit root/depth, every budget and contact convention. Output identity binds source/config, selected keys, vertex keys, actual V/F and actual handle mapping.

## Actual TDD evidence and all failed attempts

These are verbatim relevant excerpts of actual tool command outputs, recorded here from tool sessions; separate raw RED stdout files were not generated. They are not reconstructed test runs.

All build/test invocations below ran in `/data/xvxinyao/make-it-balance-paper-sim`, with `ulimit -v 25165824` (24 GiB virtual-memory cap), `timeout 3600`, maximum three compile jobs and sequential heavy build/test commands.

### Initial build-system attempt (not a behavioral RED)

After adding the tests and CMake registration, before configuring the new target:

```text
timeout 3600 cmake --build build/strict-v2-spin --target spin-it-2014-source-surface-cage-regression -j3
gmake: *** No rule to make target 'spin-it-2014-source-surface-cage-regression'.  Stop.
exit 2
```

Configured the existing build with `cmake -S . -B build/strict-v2-spin`. Configure completed; the pinned libigl emits its preexisting CMake minimum-version deprecation warning. No dependency refresh/pin changes were made.

### Availability RED, before adding production API

```text
timeout 3600 cmake --build build/strict-v2-spin --target spin-it-2014-source-surface-cage-regression -j3
.../tests/spin_it_2014_source_surface_cage_regression.cpp:1:10: fatal error: spin2014/SourceSurfaceCage.h: No such file or directory
exit 2
```

This demonstrates absent API only. Then added a minimal compilable API scaffold that copied the source/config but returned empty extraction outputs, before implementing extraction.

### Behavioral RED, before extraction implementation

```text
timeout 3600 cmake --build build/strict-v2-spin --target spin-it-2014-source-surface-cage-regression -j3 && build/strict-v2-spin/spin-it-2014-source-surface-cage-regression
[100%] Built target spin-it-2014-source-surface-cage-regression
complete selected cell count
exit 1
```

The actual public API returned zero selected cells instead of the independent 216-cell expectation. This is the behavioral RED beyond compilation.

### First extraction run: independent test expectation correction

```text
timeout 3600 cmake --build build/strict-v2-spin --target spin-it-2014-source-surface-cage-regression -j3 && timeout 3600 build/strict-v2-spin/spin-it-2014-source-surface-cage-regression
[100%] Built target spin-it-2014-source-surface-cage-regression
every boundary corner is one handle, no sampling cap
level9 GEOMETRY: nodes=73 selected=1 handles=8 triangles=12
exit 1
```

The initial tunnel fixture expected 448 vertices/896 triangles, incorrectly retaining the two hole end patches. Controller independently flagged the error before GREEN. I verified the arithmetic: from an 8-cell-side cube, remove two 2x2 end-face patches (16 triangles) and add four 2x8 tunnel walls (128 triangles): `768-16+128=880`. Vertices: `386-2+8*7=440`. Only test expectations were corrected, with derivation comments. Production was not changed to fit this result.

### Initial geometry GREEN

The same build/run command after this test-only correction returned:

```text
[100%] Built target spin-it-2014-source-surface-cage-regression
level9 GEOMETRY: nodes=73 selected=1 handles=8 triangles=12
source surface cage regression passed; complete geometry only, no BBW/physics authority
exit 0
```

### Additional capacity regression RED and GREEN

Self-review found the downstream exact auditor requires at most INT_MAX/2 triangles. Added a test with an unrepresentable INT_MAX triangle budget before implementing that guard. The same focused build/run command returned:

```text
[100%] Built target spin-it-2014-source-surface-cage-regression
exhausted budget must produce RESOURCE_LIMIT, no partial cage
level9 GEOMETRY: nodes=73 selected=1 handles=8 triangles=12
exit 1
```

Added the representability guard, then ran:

```text
timeout 3600 cmake --build build/strict-v2-spin --target spin-it-2014-source-surface-cage-regression -j3 && timeout 3600 ctest --test-dir build/strict-v2-spin -R '^spin-it-2014-source-surface-cage$' --output-on-failure
1/1 Test #10: spin-it-2014-source-surface-cage ... Passed 2.15 sec
100% tests passed, 0 tests failed out of 1
Total Test time (real) = 2.16 sec
exit 0
```

Further self-review strengthened caller-config ownership coverage so the config actually passed to the constructor is mutated, rather than an unrelated local config; final suite rebuild includes this test adjustment.

## Coverage and independent expectations

- Source [.25,.75]^3, level3: indices1..6 on all axes, 216 cells, 218 distinct boundary handles, 432 triangles, extrema .125/.875; no duplicate faces/corners; every edge has two incident triangles; all ten signed-tetra moments agree with an independent cuboid formula (volume27/64).
- Level2:64 cells/98 handles/192 triangles and unit-root moments. Level0:1/8/12. Repeat inputs yield byte-identical V/F, keys/maps and hashes; changing source/root/depth/each budget changes appropriate identities. Copied source V/F and constructor config resist caller mutation.
- Genuine level9 tiny source strictly inside cell(100,100,100): exactly73 stored oracle nodes, 1 selected cell/8 handles/12 triangles under budgets73/1/8/12. Level20 sparse supported-depth fixture also passes. These are geometry extractions only.
- Bounded cavity:504 cells/412 handles/816 triangles, two boundary components with one negative sign, volume63/64; all ten moments independently match outer cuboid minus cavity.
- Through-hole:480 cells/440 handles/880 triangles, one genus1 boundary, volume15/16; moments match unit cube minus the open tunnel cuboid. Disconnected islands:2/16/24, two components and volume1/256.
- Valid disjoint original boxes whose selected cells touch only at an edge or vertex fail explicitly with immutable invalid exact boundary diagnostics.
- Aligned face/edge/corner contact box [.125,.25]^3 selects27 cells and 56 handles/108 triangles. Non-dyadic affine root(-.3,-.2,-.1), side1.7 retains216/218/432 and exact canonical FMA corner equality.
- Invalid/empty/open/nonfinite/outside-root source, invalid indices/root/depth; exhausted/zero node/cell/vertex/triangle budgets; unrepresentable node/cell/vertex/triangle budgets; Full-subtree selected-cell rejection all covered. No partial usable result, fallback or infeasible status is accepted.

## Final build, full suite and feature-off evidence

All-target main build and full suite completed successfully. Full command:

```text
ulimit -v 25165824; timeout 3600 cmake --build build/strict-v2-spin -j3 && timeout 3600 ctest --test-dir build/strict-v2-spin --output-on-failure -j1
```

CTest normally writes raw per-test stdout to `build/strict-v2-spin/Testing/Temporary/LastTest.log` (during execution `LastTest.log.tmp`). During this full run an existing nested contract test replaced that pathname and unlinked the parent CTest temporary log. The controller detected the collision and preserved the parent process's actual open regular-file descriptor in a new evidence path; the overwritten pathname is **not** evidence of the full run. This suite was run once after all production/test changes; focused iterations above did not repeatedly rerun the full suite.

Actual parent command output:

```text
[100%] Built target partition-manifest-regression
10/200 Test #10: spin-it-2014-source-surface-cage ... Passed 2.14 sec
200/200 Test #200: carving-no-candidates ... Passed 0.11 sec
100% tests passed, 0 tests failed out of 200
Total Test time (real) = 492.41 sec
exit 0
```

Immediately copied final main cost data to new `build/task-5g-evidence-vL6KVn/main-full-CTestCostData.txt`. Existing full-suite tests, including all source/cage and both genuine level9 source regressions, pass.

Actual parent-log capture: `build/task-5g-evidence-vL6KVn/controller-parent-ctest-full.log`, 313892 bytes, SHA-256 `23fcc77d6eb664c3c8281c013ee92f95fccb5baf61e25905dfe5ee97fbd5840c`; capture metadata in companion `controller-parent-ctest-full.json`. I independently checked its digest, 200 exact `Test Passed.` records, the new cage output and its final 200/200 test/end-testing marker. This is actual parent log content retained after exit, not a recreated summary. Controller's capture helper is not staged.

Controller traced the existing overwrite to `tests/cmake/assert_ctest_commands.cmake:8`, which queries `ctest --test-dir` on the current binary directory with `--show-only=json-v1`. This preexisting registration-query side effect is outside Task5g and was not changed. Future full runs should use explicit tee capture in unique log paths.

Faster-off configure/build/test completed successfully, using the existing feature-off build with its existing local pinned dependencies, network/dependency updates disabled:

```bash
set -o pipefail
ulimit -v 25165824
{ timeout 3600 cmake -S . -B build/task-5c-no-faster-20260909 -DFETCHCONTENT_FULLY_DISCONNECTED=ON -DFETCHCONTENT_UPDATES_DISCONNECTED=ON && timeout 3600 cmake --build build/task-5c-no-faster-20260909 --target spin-it-2014-source-surface-cage-regression -j3 && timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 -R '^spin-it-2014-source-surface-cage$' --output-on-failure; } 2>&1 | tee build/task-5g-evidence-vL6KVn/faster-off-configure-build-test.log
```

Actual result:

```text
[1/4] Building CXX object CMakeFiles/make_it_stand_core.dir/src/spin2014/SourceSurfaceCage.cpp.o
[4/4] Linking CXX executable spin-it-2014-source-surface-cage-regression
1/1 Test #10: spin-it-2014-source-surface-cage ... Passed 2.00 sec
100% tests passed, 0 tests failed out of 1
Total Test time (real) = 2.03 sec
exit 0
```

Confirmed cache flags in `build/task-5g-evidence-vL6KVn/faster-off-feature-flags.log`: 2014=ON, exact=ON, Faster=OFF, viewer=OFF, fully-disconnected=ON, updates-disconnected=ON. Copied `faster-off-LastTest.log` and `faster-off-CTestCostData.txt` immediately after the run into the same new evidence directory. No hidden Faster dependency. The focused test has no warnings/errors; feature-off configure emits existing third-party CMake deprecation warnings, and compilation succeeds without new warnings.

Captured log SHA-256:

```text
5478aaf9c530b7ecefa910c1e68898adfd836b81175d40baa6873f8056c45d2e  build/task-5g-evidence-vL6KVn/faster-off-configure-build-test.log
793361b5308b7a56282dc0c28d2e664e3603fdcd05aa2aa4f69d0823f04a415a  build/task-5g-evidence-vL6KVn/faster-off-LastTest.log
```

## Files and self-review

Scoped files: `include/spin2014/SourceSurfaceCage.h`, `src/spin2014/SourceSurfaceCage.cpp`, `tests/spin_it_2014_source_surface_cage_regression.cpp`, `CMakeLists.txt`, and this report. No legacy source changes, no main-checkout writes, no GUI, dependency updates, original-data changes or history deletion. Controller acceptance/preflight and BBW gap audit are not included in my commit.

Final production/test file SHA-256 values (computed with `sha256sum` after all changes):

```text
2ede99a842720d92daf10d0ce6ace6e1d9f688adf9cd900be695c14398e43a0a  include/spin2014/SourceSurfaceCage.h
413d092ce8067782f98bced3f9abb0dfa366323b36b4fd5362b8521951bf9084  src/spin2014/SourceSurfaceCage.cpp
bae5da24f1d8249422e6ec1a9a15fd6c099da80deef9dae8b2a2ea62eadf9970  tests/spin_it_2014_source_surface_cage_regression.cpp
b428d4f79b98a4c0f0705d910bce05fdb8f9106e6423a33897f918c464db5b0a  CMakeLists.txt
```

Read-only `git -C /data/xvxinyao/make-it-balance diff --binary | sha256sum` remains `93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`, exactly the preflight main-checkout tracked-diff hash. `git diff --check` is clean. The raw no-index diff inspection of the new source returned git's normal exit1 for a new file, not a test/build failure.

Reviewed the actual diff, exposed-face winding for all six axes, Full-descendant count shifts (maximum60 bits), remaining-budget subtraction before additions, canonical neighboring vertex identity, auditor capacity and cavity sign handling. Corrected the test-only tunnel count and caller-config ownership test; added the auditor capacity guard through RED/GREEN. No helpers/reviewers/subagents dispatched. Remaining limitation is the deliberate geometry-only scope; future finer-grid BBW/joint integration must not reuse this envelope as original mass or imply original experiment evidence.
