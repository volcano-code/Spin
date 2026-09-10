# Task5g controller verification — cage geometry only

Date: 2026-09-10. Implementation: `b61993772bd9429b08425353af61a969aba6ce71`, base `13afcf0f0036c38548c156eb7d0e7284ff670897`. No production changes occurred after the reported full test run. This document is evidence, not publication authority.

## Actual verification

The implementer completed the all-target build and registered suite: **200/200 passed, 492.41 s, exit 0**. Controller independently inspected the preserved actual parent CTest raw log, its new-cage output, 200 exact `Test Passed.` records and final end-testing marker. Faster-off configure/build/test also passed **1/1**, test 2.00 s, total 2.03 s. Cache flags confirm 2014/exact ON, Faster/viewer OFF and disconnected dependencies. Feature-off configuration has preexisting third-party CMake deprecation warnings; these were not hidden and do not establish pristine dependency output.

Independent controller C++ probe calls the actual immutable public API, exports its full indexed boundary twice, and does not call the native moment integrator. The Python checker imports the separate signed-tetrahedron checker, validates actual output against literal cuboid expectations and compares the two actual output byte streams:

| Measured quantity | Actual |
|---|---|
| Selected cells / handles / faces | 216 / 218 / 432 |
| Envelope volume | 0.421875 = 27/64 |
| Maximum absolute moment/inertia error | 1.113567084231137e-19 |
| Repeated actual JSON | Byte-identical |
| Probe compilation / checker | Exit 0 / exit 0 |

These are cage-envelope moments, **not original material mass**. Inputs are analytic regression fixtures, not original paper cases. Sparse level9/20 native tests concern geometry extraction, not level9/20 optimization.

Controller compile used `/usr/bin/c++`, C++17, O2/g/NDEBUG/fPIC/frounding-math, 2014/exact guards, `include`, project `libigl/external/eigen`, the actual main core archive and its existing Boost/GMP/MPFR/HiGHS/z dependencies. Exact probe sources are retained alongside this document. Calls were bounded by `ulimit -v 25165824` and 120/60-second timeouts (virtual-address-space limit, not an RSS/cgroup guarantee).

```sh
build/task-5g-evidence-vL6KVn/controller-cage-probe > build/task-5g-evidence-vL6KVn/probe-run1.json
build/task-5g-evidence-vL6KVn/controller-cage-probe > build/task-5g-evidence-vL6KVn/probe-run2.json
PYTHONPATH=python /usr/bin/python3 .superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5g-check-probe.py build/task-5g-evidence-vL6KVn/probe-run1.json build/task-5g-evidence-vL6KVn/probe-run2.json
```

Outputs were created with noclobber; full checker output and the empty successful compile log are retained under the same evidence directory. No optimization was rerun by this diagnostic.

## TDD and log limitations

Worker report preserves availability RED, behavioral scaffold RED, final focused GREEN and the later INT_MAX auditor-budget RED/GREEN. Early RED stdout survives in tool history/report excerpts, not standalone original raw log files. One tunnel test literal was corrected from 448/896 to independently derived 440/880; no production logic, source fixture or tolerance was altered to manufacture success.

The existing registration-query test (`tests/cmake/assert_ctest_commands.cmake:8`) invokes nested CTest `--show-only=json-v1` in the current binary directory. During the full run it overwrote the normal log pathname and unlinked the parent's temporary log. Controller detected this while the parent was alive, opened the parent's **regular file** descriptor read-only, held it until parent exit and copied those actual bytes exclusively into `controller-parent-ctest-full.log`. Companion JSON records the actual parent command, regular-file target, size and digest. No pipe was consumed and no log was recreated from a summary. The capture helper is retained. This preexisting test-infrastructure side effect is not fixed by Task5g; future runs must additionally tee to unique log paths.

## SHA-256 binding

All paths below relative to worktree; evidence directory is `build/task-5g-evidence-vL6KVn`.

```text
23fcc77d6eb664c3c8281c013ee92f95fccb5baf61e25905dfe5ee97fbd5840c  evidence/controller-parent-ctest-full.log
5478aaf9c530b7ecefa910c1e68898adfd836b81175d40baa6873f8056c45d2e  evidence/faster-off-configure-build-test.log
793361b5308b7a56282dc0c28d2e664e3603fdcd05aa2aa4f69d0823f04a415a  evidence/faster-off-LastTest.log
6936a5172b555c1e143c63bb967b27e0f172effea6001ac6b81a05ebc7d01954  evidence/controller-cage-probe
43558c1bc866ac47b3ed9c9904972df29ab04c6466148b8a4168c42792e3bf06  build/strict-v2-spin/libmake_it_stand_core.a
401eafb1be343e7fdd7a828a046751cdfdfc13cacf25c5a1e2f23bc7b88eb5e2  build/strict-v2-spin/spin-it-2014-source-surface-cage-regression
1c1b66da380d89afd1089cafcc22d07524bbbd44b279cbebf39df3870c8fd627  evidence/probe-run1.json (and run2)
a83fc163abc424cdf0ff13b79f4e8b06572340d4fe2c5bd445a379cdeed0e314  evidence/controller-independent-check.log
7b7aca64a5c11c0d4a330e8fcc19a5601310f971e154f9beb8e638e485f8fab7  actual cage scientific output identity
```

Production/test source digests are preserved in the committed worker report. This additional report changes no scientific source or fixture, and its own text is not included in the cage scientific identity.

## Preservation and unresolved scope

Dirty main checkout tracked diff SHA remains `93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`; no integration, push, tag, old-output deletion or GUI mutation. `git diff --check` is clean. Original main untracked work is not overwritten. Work remains on the existing isolated branch.

The closed-cell contact convention is explicitly `IMPLEMENTATION_DEFAULT`; it is not recovered author tie behavior. All boundary components and all their corners are retained; singular unions fail closed. Complete finer-grid BBW, interpolation, joint beta/phi/handle derivatives and optimization, final deformed-material sealing, original Dancing Couple inputs and experiments remain open. This task grants neither BBW convergence nor physical/feasibility/original-reproduction success. Strict overall status remains `PARTIAL_BLOCKED`.
