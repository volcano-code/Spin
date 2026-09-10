# Task5h controller verification

Bounded scope: sparse finer cage grid and trilinear field transfer, not BBW, joint deformation or original-case reproduction. BASE `59ed96c12e8f24582f14459056c05f74d9c15492`. This record does not grant publication authority or claim full build-dependency attestation.

## Preflight / independent actual output

Existing linked worktree `/data/xvxinyao/make-it-balance-paper-sim`, strict branch, no new worktree. Baseline source-cage/source-octree tests2/2 passed3.35s. Original complete-cage increment was not reimplemented.

After implementer confirmed the all-target build complete and main archive/source stable, controller compiled only a separate diagnostic executable, without CMake/project rebuild or mutation of main build files. This read-only linkage/check ran alongside the full suite with explicit coordination. Native probe calls actual `SourceCageGrid`, constructs fields from its actual nodes and exports all cell/edge/handle connectivity plus six query embeddings and values twice. Python checker constructs independent Cartesian expectations and brute-force containing-cell searches, not the production sparse lookup.

Actual outcome (exit0):

| Quantity | Checked outcome |
|---|---|
| Fine cells / nodes / edges / handles | 1728 / 2197 / 6084 / 218 |
| Complete connectivity and all handles | Match independent grid expectations |
| Closed lexicographic boundary ties | Match brute-force containing-cell oracle |
| nextafter outside envelope | Rejected, not snapped |
| Maximum multiaffine-field absolute error | 1.3877787807814457e-17 |
| Repeated actual export | Byte-identical |

Fields `[1,x,y,z,xy,xz,yz,xyz,2]` are diagnostic values, not optimized weight solutions. Preserving the constant2 additionally checks there is no field renormalization. This establishes interpolation precision, not BBW linear precision or KKT.

Retained helpers alongside this document: `task-5h-controller-probe.cpp`, `task-5h-check-probe.py`, and hand-derived `task-5h-controller-acceptance.md`. Initial checker invocation without file arguments correctly exited1 (`two actual exports required`); it is not counted as a successful native test. Actual successful two-file invocation and output are in `build/task-5h-controller-8bSTgg/check.log`.

Controller commands:

```sh
# All commands in worktree; ulimit -v25165824, compiler timeout120, each probe timeout60.
/usr/bin/c++ -std=c++17 -fPIC -O2 -g -DNDEBUG -frounding-math \
  -DMAKE_IT_STAND_HAS_SPIN_IT_2014=1 -DMAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS=1 \
  -Iinclude -isystem libigl/external/eigen \
  .superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5h-controller-probe.cpp \
  -o build/task-5h-controller-8bSTgg/probe build/strict-v2-spin/libmake_it_stand_core.a \
  /home/lin/anaconda3/lib/libboost_system.a /usr/lib/x86_64-linux-gnu/libgmp.so \
  /usr/lib/x86_64-linux-gnu/libmpfr.so /home/lin/anaconda3/lib/libboost_thread.a \
  build/strict-v2-spin/lib/libhighs.a /usr/lib/x86_64-linux-gnu/libz.so
build/task-5h-controller-8bSTgg/probe > build/task-5h-controller-8bSTgg/run1.json
build/task-5h-controller-8bSTgg/probe > build/task-5h-controller-8bSTgg/run2.json
/usr/bin/python3 .superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5h-check-probe.py build/task-5h-controller-8bSTgg/run1.json build/task-5h-controller-8bSTgg/run2.json
```

Actual shell used `set -euo pipefail`, compile/check stdout tee, and `noclobber` for export paths. Virtual-address-space bound is not an RSS/cgroup claim. Compile log is empty with exit0.

## Observed initial source/binary binding (pre-fix)

```text
3ad0e0989f486a8a6408bdee12b4948f8997a7dff86fe623cd47b7858584a950  include/spin2014/SourceCageGrid.h
0084902b1d78c24f7a3ae70a16f66123defc0213a4247363934cb07e71d5a88b  src/spin2014/SourceCageGrid.cpp
0a2dd5899d1a7f8cced537d638e5a5fedac702996da3a836208d653193d49d99  tests/spin_it_2014_source_cage_grid_regression.cpp
7858ba45ad02501072ab84caf8fcb03a1db09b7f992452f14ee3d4d442d8d1ef  CMakeLists.txt
1ffd25645562946405cf2db29c54ba6d93acb8aef91903596e351bcf0195d056  build/task-5h-controller-8bSTgg/probe
419273419afb3c33209545aeb1b39f1b19fccd5c4b14608095aa8c725a6ac514  build/strict-v2-spin/libmake_it_stand_core.a
7b105c58615c14a47776e447657d77b7f568bc522f43338b5e8f667dd95b4111  build/strict-v2-spin/spin-it-2014-source-cage-grid-regression
b4193f3e510cc6f142765ec90e5e9c3f164d86e97e03cc8e9c65add3e7b1dd6e  build/task-5h-controller-8bSTgg/run1.json (and run2)
```

## Initial integration verification (pre-fix)

Implementation commit `4c6f89b8ef5a72a47e4ca6ac872ecf22dd122990` contains only production/header/test/CMake/worker-report files. Controller reread the final report and independently inspected unique original integration logs and hashes. Main source/header/test/CMake/core-archive digests above were rechecked unchanged after both build configurations, so the earlier actual probe remains bound to the tested implementation.

- Main all-target build exit0; full registered CTest **201/201 passed**, 474.63s. Actual outer stdout contains201 passing test records and final100% summary, no NotRun/skip/failure. New source-cage-grid test0.94s.
- FasterOFF configure/build/focused CTest all exit0; **1/1 passed**, nine named groups,0.95s. Controller inspected exactON/2014ON/FasterOFF and both FETCHCONTENT disconnectedON cache flags.
- Full suite ran once after final production changes. No rebuild or repeated full suite for later controller documentation.
- Two genuine behavioral RED/GREEN cycles have original raw tee logs in `build/task-5h-evidence-20260910/`: empty scaffold (eight failures) and temporary-edge-set capacity guard (one new failure). Unavailable target before CMake regeneration is separately recorded as availability failure, not behavioral RED. Controller read the actual eight-failure RED capture.
- Unlike Task5g, outer tee capture avoided reliance on the default CTest log pathname; no descriptor rescue or reconstructed log needed.

```text
43b0f8e2b82366407da1f4fbd2b4380a12cebce1ad7bb942f0009e8b0495153d  build/task-5h-evidence-20260910/full-build-01.log
ce506d59c10047a772dd7a95f48a7547745e57cdd680171a3ee9e7165193a623  build/task-5h-evidence-20260910/full-ctest-01.log
99a685426208814fd2efcaad50b6e0fb23a0778d4d8d448a56d64029ad357a9f  build/task-5h-evidence-20260910/fasteroff-configure-01.log
b6c6db9dc0d5a412babece7a3d48efb32beee076437f650f875d98c512ce2388  build/task-5h-evidence-20260910/fasteroff-build-01.log
7d0a761e3399ac4db707296dbeb785d0431cf220ee4f0fddb2778179d60507b4  build/task-5h-evidence-20260910/fasteroff-ctest-01.log
7d1f2ca9da83aa7b1d17ce8da9eebd0b1fa4c56c06c7e16819901da91b3d6017  build/task-5h-controller-8bSTgg/check.log
```

Existing pinned dependency CMake deprecation warnings and a `cp -n` informational deprecation warning during log preservation are disclosed, not hidden. New compile/test output has no new compiler warnings or unexpected skips. No dependencies updated or downloaded for the build.

## Preservation / scope

Main tracked diffSHA remains `93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`; no main writes, old-output deletion, merge/push/tag/release or GUI changes. This hash only covers tracked differences; it is not falsely presented as a full untracked-tree audit. Source/build ownership has been released. `git diff --check` passes. Independent task review has its separate record; this verification document does not replace it.

The uniform finer cage-volume domain and strict represented-boundary tie rule are explicit implementation conventions, not recovered author settings. BBW energy assembly, cage-face boundary constraints, applicable solver/KKT protocol, joint beta/phi/handle optimization, deformed-material sealing and missing original experiments remain required. Overall status remains `PARTIAL_BLOCKED`; this task is an implemented/verified computational prerequisite, not a new original-paper case result.

## Fix1 verification and replacement bindings

Initial review found a real floating-coordinate counterexample despite ordinary tests passing. Corrective commit `da1bd34741c776275f420d274d4479f1291de00c` now requires fineStep strictly greater than the largest root-endpoint ULP. Actual upstream-backed RED proves7 coincident fine intervals for x-origin nextafter(2^48,0),side1,coarse3/fine4. Two additional tests document conservative rejection of valid aligned equality cases; these are not misreported as two further collapse defects. This guard can reject otherwise valid numerically marginal roots, explicitly documented; it does not alter input positions or claim physical infeasibility.

Fix covering tests: main4/4 (23.47s), FasterOFF4/4 (23.61s), including all12gridgroups, sourcecage and sourceoctree/level9. Original raw RED/GREEN logs are in`build/task-5h-fix1-evidence-20260910/`, with hashes in the committed fix appendix. Original reviewer scoped re-review: sole Important ADDRESSED, no new breakage, SpecPASS/qualityApproved. See task-5h-review.md.

Controller recompiled the exact same separate C++ probe against the amended stable archive, with output directory`build/task-5h-controller-fix1-GSRMx9/`, and reran both actual exports and independent Python checks. Exit0 throughout; output bytes match each other AND the ordinary pre-fix output, maximum field error still1.3877787807814457e-17. No prior outputs were overwritten. This establishes unchanged ordinary-domain behavior, not validity of old implementation provenance.

```text
7acc63d8862e611d208d77a9175f8cb34f3dbc60bd520b32f2d8dd40a947964f  include/spin2014/SourceCageGrid.h
b9e21d6468a2bd3ea9bc2303637f829e29c10fcaef4eff8a00d3ed0a71c6f9d9  src/spin2014/SourceCageGrid.cpp
0ca7e7b92141d906d2efe15ca8bb14e89533bebca6aa6f392f046032960a77e0  tests/spin_it_2014_source_cage_grid_regression.cpp
b85b3d04a7af9b86a811e6d02b9dd0c06fc0a7bbdd27520b88d468c1080ca6dd  build/strict-v2-spin/libmake_it_stand_core.a
6a5375fb5c5e0f5f2e55d6ed97046332336df52876aa1f17cdc2dfdcf0159bbd  build/strict-v2-spin/spin-it-2014-source-cage-grid-regression
80293a970393ba8acbb082510024ff98c8e553be4e34d9adb5d050339f7601cd  build/task-5h-controller-fix1-GSRMx9/probe
b4193f3e510cc6f142765ec90e5e9c3f164d86e97e03cc8e9c65add3e7b1dd6e  build/task-5h-controller-fix1-GSRMx9/run1.json (and run2)
```

Final current-source full-suite execution is recorded after actual completion below; initial201/201 remains pre-fix, not upgraded by a documentation change.

## Final post-fix all-target and whole-suite verification

Controller-owned run against `da1bd34741c776275f420d274d4479f1291de00c`: all-target build (`-j3`) and full CTest (`-j1`) both completed, combined shell exit0. **201/201 passed, 482.35 seconds**. The unique outer capture contains201 actual passing test records, no duplicate/missing test IDs, no NotRun/Skipped/CTest failure markers, and the final100% summary. This is the fresh post-fix run, not the earlier474.63-second run. All twelve finer-grid groups are also covered by the separately recorded main/FasterOFF4/4 runs.

Actual commands used `set -euo pipefail`, `ulimit -v 25165824`, individual3600-second timeouts and unique tee logs under `build/task-5h-final-AkTS5P/`. The limit is virtual address space, not a measured RSS/cgroup bound. Original logs remain intact. Header/source/test/CMake/core-archive/grid-executable hashes were rechecked after this final run and match the replacement bindings above. Controller evidence/documentation edits do not change those computational inputs.

```text
135f1cb7ecfd26313c934a6d9dd528aca53c90cf8bfb552f80848c2f20839a5d  build/task-5h-final-AkTS5P/build.log
40ed63dd6418ad3467eae19e0f39e1e47da9e5ef78267ea11b0fd87449f11c98  build/task-5h-final-AkTS5P/ctest.log
0fd922731366cb0314844a22ec4b2dadb1bdae5c52dcf2bb4907b1e15da231e7  build/task-5h-controller-fix1-GSRMx9/check.log
```

Main tracked-diff preservation SHA rechecked unchanged (`93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`). Task5h's bounded grid/interpolation implementation is verified; the strict reproduction remains `PARTIAL_BLOCKED`. No full-original-case, BBW optimality, joint-deformation, or new physical/visualization result is inferred from this test count.
