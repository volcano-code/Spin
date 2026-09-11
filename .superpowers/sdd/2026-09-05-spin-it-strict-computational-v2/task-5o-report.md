# Task5o implementation report: immutable actual source material/full BBW binding

## Scope and status

Task5o implements rest geometry/interpolation only. Global strict Spin-It 2014
status remains **PARTIAL_BLOCKED**. No deformed material boundary, hanging-face
conformity, joint SLQP, paper regularizer, material seal or original-case
acceptance is claimed. All numerical fixtures are SYNTHETIC_REGRESSION.

Worktree: `/data/xvxinyao/make-it-balance-paper-sim`.
Verified initial HEAD: `bf925db611195942a06d119862d61a36d383d49e`.
Assigned build: `build/task-5c-no-faster-20260909`, Release, 2014 ON, exact
geometry ON, Faster OFF. Builds used max `-j3`; CTest heavy work ran sequentially
with `ulimit -v 25165824`. No protected-main build or mutation, tag, push, merge,
deletion or subagent dispatch occurred.

## Implementation

- Added `SourceCageBinding` as a shared-private-const-Data value object. The
  constructor accepts only actual immutable source/BBW inputs and explicit
  resource settings. It requires BBW success and available passing final primal
  validity, retaining both protocols and their original raw/final assessments.
- Compared original outer vertex/triangle dimensions and scalar bytes, and
  root origin/side bytes including signed zero. Distinct design domains with an
  exactly matching outer are accepted; unequal upstream source digests are
  expected and both are preserved.
- Copied the adaptive partition and enumerated actual included Full cuboids
  with `sourceBinaryCuboids(domain, ones, maximumPrimitives)`. An ancestor hash
  map assigns exactly one original-order beta owner to every primitive in
  expected O(depth) lookups. Every adaptive leaf owns material.
- Retained original outer rows, followed by sorted unique material-corner
  lattice keys. Canonical source endpoints, exact corner-key deduplication and
  x-bit-first connectivity are retained without snapping or outer reindexing.
- Cached actual eight-node grid embeddings. Bounded `pointWeights(begin,count)`
  returns every handle from the actual final BBW field, without dense global
  field/LBS allocation, clamp, normalization or subsampling.
- Guarded malformed/overflow budgets, incremental primitive/point counts and
  output chunks. Resource failures use RESOURCE_LIMIT length_error; invalid
  inputs use invalid_argument; chunk-range failures use out_of_range. Allocation
  failure is translated and no partial binding escapes.
- Added versioned binary SHA identities covering upstream scientific inputs,
  partition order, budgets, actual final weights, original geometry/root bytes,
  derived point/primitive arrays and cached embeddings. Timings are excluded.
- Registered one native regression and gated all new compilation under both
  exact geometry and 2014. Added the rest-binding contract and remaining gates
  in `docs/spin2014_source_cage_binding.md`.

## TDD evidence and focused verification

The regression file was written before the new header/source. It specifies the
actual source/FEM/BBW fixture, owner moments, canonical geometry, bounded full
weights, immutability, protocol, identity and failure behavior.

RED command, before implementation:

```sh
ulimit -v 25165824
/usr/bin/c++ -std=c++17 -DMAKE_IT_STAND_HAS_SPIN_IT_2014=1 \
  -DMAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS=1 -Iinclude \
  -I/usr/include/eigen3 -fsyntax-only \
  tests/spin_it_2014_source_cage_binding_regression.cpp
```

Exit 1, expected missing-interface failure:

```text
tests/spin_it_2014_source_cage_binding_regression.cpp:1:10: fatal error:
spin2014/SourceCageBinding.h: No such file or directory
compilation terminated.
```

GREEN commands:

```sh
ulimit -v 25165824
cmake --build build/task-5c-no-faster-20260909 \
  --target spin-it-2014-source-cage-binding-regression -j3
ctest --test-dir build/task-5c-no-faster-20260909 \
  -R '^spin-it-2014-source-cage-binding$' --output-on-failure -j1
```

Initial build exit 0 and focused CTest 1/1 passed (0.53 s). Existing CMake/libigl
minimum-version deprecation warnings were visible during regeneration; they
were not suppressed. No compilation warnings or test failure occurred.
The final focused test expanded byte assertions (face order, one-ULP vertex
change, root signed zero, outer signed-zero mismatch and preserved matching
signed-zero arrays); test-only rebuild exit 0, focused 1/1 passed (0.92 s).
The implementation header/source/archive did not change after initial GREEN.

Main fixture: original outer `[.3,.7]^3`, design `[.34,.66]^3`, root `[0,1]^3`,
source maximum depth 5, initial optimization 2 / max 5, coarse cage 2, fine grid
3, actual FEM and coupled 26-handle BBW. Result: 8 adaptive leaves, 496 Full
primitives, 997 points, volume `0.030517578125`. The material matches literal
`[11/32,21/32]^3`, including ten independently integrated nonzero moments.
Nonuniform owner beta matches independent per-cuboid integration. The chunk
budget is 182 values (7 full rows); the 25,922-value global field is rejected.
Tests also cover split/merge snapshot isolation, changed design, actual failed
and invalid BBW, independent-normalized primal acceptance without coupled KKT,
all three budget identities, invalid/overflow ranges and transformed geometry.

## Safe scoped suite and exact-OFF check

The assigned safe regex was invoked exactly once after final focused GREEN:

```sh
ulimit -v 25165824
ctest --test-dir build/task-5c-no-faster-20260909 \
  -R '^(spin-it-2014-|spin-2014-|source-spin-material-python$)' \
  --output-on-failure -j1 \
  --output-log build/task-5c-no-faster-20260909/task-5o-safe-scoped-suite.log
```

Final result: **33/33 passed**, exit 0, total 223.40 seconds. The actual
dependency contract passed in 89.68 seconds and its isolation checks in
7.32 seconds. Suite log SHA-256:
`4baef868d43fe8ea47b15b3d0d38447818f3f02be647df76547f286cd0cc944c`.
The new registration increases this safe selection from 32 to 33 tests,
including the three previously unlabeled tests. No unrelated destructive
dependency-contract tests or broad main suite were run.

The retained actual dependency run is
`build/task-5c-no-faster-20260909/spin-2014-dependency-contract-runs/run-2edd7750b6e83d13`.
Its enabled CMakeCache has 2014 ON and exact geometry OFF; read-only inspection
of its generated `build.ninja` found no `SourceCageBinding` entry. This verifies
actual exact-OFF graph exclusion without creating an additional build.
The retained run's `contract.status` is `PASSED`; enabled graph SHA-256 is
`cec2bb7e4ce179404d0badcdb87cc1a4e9246c0ece7d4979387db6a4e33b49b4`.

## Scientific input and binary hashes

SHA-256 values at final implementation/test source:

| Item | SHA-256 |
| --- | --- |
| Initial assigned archive (controller-provided) | `85a80297bbb7d82dd176df016f821757f11d3e6b1ff4012a304de9d35f9fdc24` |
| New header | `993887d711560770368331d314156ceaf57fd6162d2161a865a22de143aad991` |
| New source | `69f09a52a0536bccd7903be8a592aad921a5588c01e7edbcf82d32bc6f5793fd` |
| Final regression source | `d92cf6ec45a16aeabbb76aa36c3faf34f81a65d3c631d1070b96d6d8f62293ae` |
| Documentation | `1a918967ed136cbf29cdcf38a72bcbaadbf4b75d4186ee7fab0923049e22c186` |
| Assigned `libmake_it_stand_core.a` | `4c7b75d1620190594c06e9bee10bb0abf37cf8e6f6edd91865b7f84142404d43` |
| Focused regression executable | `128306b8bf6b182f51533750acad3411cdbcc4d66022b0edc38ffbf5f181b560` |
| Main original source content identity | `7c4c4be1b8bd5c9817074dda710fc5b2cee258df2e499e2371c52417fa476cb9` |
| Main source partition identity | `6beb11aa4064bc447809986e92296a3c68c8476193a0f2d77f78c02a59b24fa6` |
| Main actual coupled BBW output identity | `c4fccbc0be3bca40485bfb169479b7ec4d35b49050e5e4ec3ef2b72be40da2b2` |
| Main 182-value-budget binding configuration | `9943e0f3bde071c82d49b6bd09179d0b272b22ec1140cc4c7cc7b7d5c2ccf835` |
| Main 182-value-budget binding output | `ac1698162cd63efb42eb63dad1e4bb94d8fec05e3f415d1e8a651e7f3dbfb469` |

## Independent controller evidence (attributed)

The controller independently compiled its exporter against the same header,
source and archive hashes above. Its reported probe covers both protocols and
two geometries, each run twice. Exports were byte-identical with SHA-256
`a26451a8e3f32644e86bb114922837338fa4edc72bf8e4b451ca5837a9b9b00a`.
Fixed fixture counts: 8 / 496 / 997 / 26; transformed nonsymmetric fixture:
8 / 147 / 324 / 26. Full-unit coverage, beta ownership, canonical corners,
closed-cell embedding and all 26 handle weights passed independently. Reported
maximum per-leaf moment error was `5.204e-18`, total ten-moment error `6.939e-18`,
weight difference zero. Six corrupted record types were rejected per case
(24 negative checks). Evidence is under `build/task-5o-controller-jOBBE6`.
These are controller results, separate from the native tests above.

## Self-review, ownership and commits

Self-review checked every brief requirement against the code and tests:
canonical shift arithmetic is bounded by source maximum depth 20; range checks
use subtraction before multiplication; actual immutable final fields preserve
handle order; hash iteration follows deterministic arrays rather than hash-map
order; corner insertion is incremental; no temporary dense global handle field
is formed. Count budgets bound containers, Eigen indices and relevant byte
allocations; they are not a combined resident-memory budget for retained
upstream objects. No implementation correctness concern remains identified.
`git diff --check` passed.

Owned files: `CMakeLists.txt`, new binding header/source, one native regression,
one documentation file, and this report. The following six preexisting Faster
paths remain unmodified by this task and excluded from commits:
`include/spin/QuadricMomentEvaluator.h`,
`src/spin/QuadricHollowingSolver.cpp`,
`src/spin/QuadricMomentEvaluator.cpp`,
`tests/quadric_moment_evaluator_regression.cpp`,
`src/spin/QuadricMaterialExportGuard.h`,
`tests/quadric_certified_empty_cases.h`.

Implementation commit: `dd94566929fc3beb8e5a843a8fa773f43f186ffe`,
`feat(spin2014): bind immutable source material to full BBW`. Author/committer: volcano
`<xuxinyao307@gmail.com>`. The report follows as a separate evidence commit so
that it can record the exact implementation commit without self-referential
hashes. Build ownership is returned to the controller after the passing suite
and scoped commits; no build remains active at handoff. Controller verification
files and its future ledger are excluded from these commits.
