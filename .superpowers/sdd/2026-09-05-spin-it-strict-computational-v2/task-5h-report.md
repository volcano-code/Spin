# Task5h implementation report

Base: `59ed96c12e8f24582f14459056c05f74d9c15492`, branch `spin-it-strict-computational-v2`.
Implementation identity: `volcano <xuxinyao307@gmail.com>`.
Initial implementation status at `4c6f89b` (pre-review/pre-fix): **DONE**. The original verification evidence below describes that commit, not the later review correction. The reviewer subsequently found an equal-ULP coordinate-collapse defect; see the appended Fix round 1 report for the corrected behavior and its covering evidence.

## Scope and API

Added immutable native `spin2014::SourceCageGrid` in `include/spin2014/SourceCageGrid.h` and `src/spin2014/SourceCageGrid.cpp`. Constructor takes an already validated `SourceSurfaceCage` and a copied `SourceCageGridConfig`. Immutable `shared_ptr<const Data>` owns the results and an immutable cage value; cage internals remain copy-owned immutable upstream data, so reassigning/destroying the caller's cage cannot affect the grid.

Configuration has `fineLevel` (default 4), `maximumCells` (1,000,000), `maximumNodes` (2,000,000), `maximumEdges` (6,000,000), and `maximumInterpolationValues` (10,000,000 scalar output values). The public API provides:

- `cage()`, `config()` by const reference.
- `cells()`, `nodeKeys()` as lexicographic integer `array<uint32_t,3>` vectors.
- `nodePositions()` as a finite Nnodes x 3 matrix with canonical represented coordinates.
- `cellNodes()` as `vector<array<int,8>>`, x-bit-first local corner order.
- `edges()` as sorted unique `array<int,2>` pairs of increasing node indices.
- `handleNodeIndices()[i]` mapping `cage().handleVertexIndices()[i]` to its exact fine node.
- `sourceSha256()`, `cageConfigSha256()`, `cageOutputSha256()`, `configSha256()`, `outputSha256()`.
- `embed(Vector3d)` returning `Embedding{cellIndex, cellKey, nodeIndices, coefficients}`.
- `interpolate(MatrixXd nodalField, MatrixXd queries)` returning Nqueries x K, from Nnodes x K and Nqueries x 3 inputs.

The interpolator supports arbitrary finite fields with K > 0. Empty queries explicitly return 0 x K, but still validate field finiteness. It does not clamp or normalize fields. It does not certify weight bounds, partition of unity, KKT, BBW convergence, or BBW linear precision. Tests of polynomial reconstruction and a positive partition fixture are **SYNTHETIC_REGRESSION** of generic transfer, not solved BBW.

No edits were made to `SourceSurfaceCage`, `SourceSurfaceOctree`, their identities, legacy dense deformation, graph-Laplacian BBW, or existing solver APIs. No solver, optimizer, dependency, Python/mesh bridge, visualization, geometry seal, feasibility, physics, or paper-case success is claimed. Controller documentation changes were left outside this commit.

## Domain, indexing, and represented-coordinate conventions

The explicit **IMPLEMENTATION_DEFAULT** is to uniformly subdivide every selected coarse cage cell to `fineLevel`, requiring `cageLevel < fineLevel <= 20` and retaining the same root. This preserves the union of the validated coarse cage cells and its holes/components. It is a deformation/weight computational domain, not original or final material, and not recovered author discretization metadata. No independent fine-level source classification removes coarse handles.

1. Validate budgets and finer coordinate representability. As in the upstream oracle, compare the finer step against the largest ULP at the root endpoints; numerically collapsed finer planes reject.
2. Count the exact refinement product against `maximumCells` by division before expansion/reserve. Per-coarse count is `uint64_t(1) << (3 * levelDifference)`, with shift at most 60. Enumerate only selected coarse descendants and lexicographically sort them.
3. Insert each cell's integer corners into a budget-checked ordered map. Assign node indices in lexicographic order. Position coordinate is exactly `fma(ldexp(double(integerCoordinate), -fineLevel), rootSide, rootOrigin[axis])`. Coarse integer keys scaled by the power-of-two refinement factor produce bit-identical handle positions.
4. Emit eight-node connectivity, then only the twelve axis edges of each actual cell into a budget-checked ordered set. Copy unique sorted edges to the public vector. Keep all coarse handles exactly once through direct integer-key mapping.
5. For embedding, binary-search canonical represented planes independently on each axis to find the first plane >= the coordinate. A non-plane point has one possible incident cell per axis; an exact shared plane has at most two. Enumerate at most eight candidates in lexicographic `(x,y,z)` order and binary-search the sparse active cell vector. Select the first active candidate. This is not a linear-root-index tie rule and does not scan every fine cell per query.
6. Outside-root, outside-envelope, and missing-cell/hole points throw `out_of_range`. Exact boundary points are accepted. There is no geometric tolerance, nearest-node fallback, or clamping of an outside point.
7. Fractions use the actual canonical represented lower/upper planes, with long-double subtraction/division and double coefficients. Tensor-product coefficients are finite and nonnegative, sum to one to roundoff, and give exact nodal identity. Shared-face/edge restrictions agree continuously. Newly split adaptive corners use the same embedding/transfer calls without rebuilding a weight field.

The construction uses O(Ncells + Nnodes + Nedges) storage, including temporary node/edge containers. Lookup is O(fineLevel + log Ncells), up to a fixed eight cell searches. Transfer uses the output matrix and a single fixed-size embedding; no Nnodes x Nnodes matrix or stored Nqueries x 8 intermediate is formed.

## Resource and identity contracts

All four configured limits must be positive and fit their container, integer, Eigen-index, and byte capacities. Nodes fit signed int connectivity, Eigen rows x 3, and position bytes. Cells fit key and connectivity vectors. Edges fit both the final vector and temporary ordered set. Output limits fit Eigen scalar counts and double bytes. New distinct nodes/edges are checked before insertion; actual final matrix/vector allocations follow those checks. Output rows are checked against `maximumInterpolationValues / K` before any output allocation, avoiding product overflow. Configured resource failures are `length_error` containing `RESOURCE_LIMIT`, never infeasibility, coarsening, culling, or partial results. Invalid dimensions/levels/nonfinite inputs throw `invalid_argument`. Nonfinite interpolation arithmetic results throw `overflow_error`; accumulation uses long double without changing user fields.

SHA-256 serialization follows the upstream fixed eight-byte big-endian integer, length-prefixed string, and IEEE binary64-bit conventions. It never hashes padding or native endian memory. Source and cage hash getters expose the exact upstream values. Grid config identity binds the versioned selected-coarse-union/closed-represented/lexicographic convention, cage configuration hash, level, and all four budgets. Grid output identity binds source, cage config/output, grid config, all cell/node keys, actual FMA positions in row/axis order, complete connectivity, sorted edges, and handle mapping. Repetition yields identical arrays and hashes; source changes affect output even if the selected cell domain is unchanged, and every budget affects config/output identity.

## Behavioral coverage

Nine named test groups exercise the real API and upstream exact source fixtures:

1. Unit-root original cube [.25,.75]^3: validated level3 cage 216 cells/218 handles; fine level4 exactly 1728 cells/2197 nodes/6084 edges, literal fine indices 2..13, envelope [.125,.875]^3. All connectivity and edge integer adjacencies and exact coarse-handle positions are checked.
2. Independent analytic fields `[1,x,y,z,xy,xz,yz,xyz]` at original surface vertices, all fine nodes, outer corners, coarse/fine shared points, non-grid (.3,.45,.6), and newly split (.28125,.34375,.40625), tolerance 1e-12. Positive arbitrary partition fixture preserves sums; constant 2 stays 2; empty query preserves K.
3. Literal interior coefficients [.064,.256,.016,.064,.096,.384,.024,.096], literal first incident face/edge/corner keys, exact split-point eight 1/8 coefficients, finite/nonnegative/sum checks, embedding index/connectivity agreement. Exact boundaries accept; nextafter outside rejects and nextafter inside accepts on all six envelope sides. Shared-plane transfer is continuous.
4. Analytic bounded cavity (4032 fine cells), through-hole (3840), and disconnected islands (16) retain their missing cells/components. Hole-boundary ties skip absent earlier candidates. Every edge belongs to actual cells; no orphan nodes.
5. Sparse coarse level9 tiny box in cell(100,100,100), fine10 under exact small limits yields 8 cells/27 nodes/54 edges/8 handles. Huge requested refinement rejects under the tiny cell limit. A genuine sparse coarse19/fine20 fixture likewise yields 8/27/54 with no root lattice expansion or optimization claim.
6. Root origin (-.3,-.2,-.1), side1.7: exact canonical FMA node positions/handle matches, analytic affine and multilinear transfer, represented shared-corner tie, nextafter outside/inside on all axes.
7. Invalid/repeated/>20 levels, zero/unrepresentable budgets, separate exact cell/node/edge exhaustion and interpolation-output exhaustion, row/column mismatch, zero field columns, NaN/infinite queries and fields (including empty queries) reject.
8. Temporary edge-set capacity budget exceeding its `max_size` rejects before construction; this caught and fixed a real self-review gap beyond the initial vector-capacity guard.
9. Caller cage/config reassignment leaves the grid unchanged; repeat arrays, positions, mappings and hashes agree; changed source, upstream cage budget, each grid budget and fine level bind the appropriate identities. Const-reference ownership API assertions are included.

## Actual RED/GREEN commands and evidence

All commands ran from `/data/xvxinyao/make-it-balance-paper-sim`. Every build/test/configure invocation used `set -o pipefail`, `ulimit -v 25165824`, `timeout 3600`, and original unique `tee` stdout/stderr paths. Builds were sequential and used at most `-j3`. No download or dependency update was requested.

Persistent evidence directory: `build/task-5h-evidence-20260910/`. Initial `/tmp/task-5h-*-20260910-*.log` files were copied unchanged here using `cp -n`; originals were retained. The copy tool printed an informational deprecation warning for `-n`, with no failed copy. Each named log below is an actual original stdout capture, not a reconstructed transcript. Final log SHA-256 values are listed below.

Focused build command B:

```bash
set -o pipefail
ulimit -v 25165824
timeout 3600 cmake --build build/strict-v2-spin --target spin-it-2014-source-cage-grid-regression -j3 2>&1 | tee /tmp/LOG
```

Focused test command T:

```bash
set -o pipefail
ulimit -v 25165824
timeout 3600 ctest --test-dir build/strict-v2-spin -R '^spin-it-2014-source-cage-grid$' -V 2>&1 | tee /tmp/LOG
```

Run sequence with exact LOG names:

| Command | LOG | Exit/result |
|---|---|---|
| B before regeneration | task-5h-red-build-20260910-01.log | 2: no rule to make target; not behavioral RED |
| `timeout 3600 cmake -S . -B build/strict-v2-spin` with the same shell/tee wrapper | task-5h-red-configure-20260910-01.log | 0; existing libigl CMake deprecation warning |
| B, compiled empty scaffold | task-5h-red-build-20260910-02.log | 0 |
| T, empty scaffold | task-5h-red-ctest-20260910-01.log | 8; all eight original groups fail |
| B, first implementation | task-5h-green-build-20260910-01.log | 0 |
| T, first implementation | task-5h-green-ctest-20260910-01.log | 0; 8/8 groups, 1/1 CTest pass |
| B, additional capacity test | task-5h-capacity-red-build-20260910-01.log | 0 |
| T, missing temporary-set capacity guard | task-5h-capacity-red-ctest-20260910-01.log | 8; only new capacity group fails |
| B, capacity guard added | task-5h-capacity-green-build-20260910-01.log | 0 |
| T, final production behavior | task-5h-capacity-green-ctest-20260910-01.log | 0; 9/9 groups, 1/1 CTest pass, 0.95 seconds |

Original behavioral RED includes:

```text
FAIL refine_selected_volume_without_losing_full_cells_or_handles: 12^3 cells, 13^3 nodes, 3*12*13^2 edges
FAIL multilinear_transfer_original_and_new_adaptive_vertices_without_resolve: query by field shape preserved
FAIL closed_represented_ties_and_nextafter_do_not_snap_outside: face chooses lexicographic first active incident cell
source cage grid SYNTHETIC_REGRESSION: 8 failed groups; no BBW or physics authority
```

The later capacity RED actually printed:

```text
FAIL temporary_container_capacity_budget_rejects_before_construction: resource failure must be explicit RESOURCE_LIMIT
source cage grid SYNTHETIC_REGRESSION: 1 failed groups; no BBW or physics authority
```

Final focused GREEN printed all nine named PASS lines and:

```text
SYNTHETIC_REGRESSION sparse fineLevel20: 8 cells / 27 nodes / 54 edges; no optimization claim
source cage grid SYNTHETIC_REGRESSION: 0 failed groups; no BBW or physics authority
100% tests passed, 0 tests failed out of 1
```

## Initial integration checks (pre-fix commit 4c6f89b)

After the final production/test/CMake changes, ran the main all-target build and full CTest **once each**, sequentially, with the same `set -o pipefail`, `ulimit -v 25165824`, `timeout 3600` shell guards:

```bash
timeout 3600 cmake --build build/strict-v2-spin -j3 2>&1 | tee build/task-5h-evidence-20260910/full-build-01.log
timeout 3600 ctest --test-dir build/strict-v2-spin --output-on-failure 2>&1 | tee build/task-5h-evidence-20260910/full-ctest-01.log
```

Build exit 0, all targets built. Full CTest exit 0, **201/201 tests passed**, no failures, total **474.63 seconds**. The new test passed in 0.94 seconds. The captured outer log independently contains 201 `Test ... Passed` records and the final 100% pass summary. This unique `tee` file preserves the run despite the existing nested registration test overwriting default `LastTest.log`; no rescue or reconstructed log was needed. Existing CMake minimum-version deprecation diagnostics were emitted by pinned build dependencies; no new compiler warnings occurred.

Then ran the existing FasterOFF build configuration sequentially (same guards):

```bash
timeout 3600 cmake -S . -B build/task-5c-no-faster-20260909 2>&1 | tee build/task-5h-evidence-20260910/fasteroff-configure-01.log
timeout 3600 cmake --build build/task-5c-no-faster-20260909 --target spin-it-2014-source-cage-grid-regression -j3 2>&1 | tee build/task-5h-evidence-20260910/fasteroff-build-01.log
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 -R '^spin-it-2014-source-cage-grid$' -V 2>&1 | tee build/task-5h-evidence-20260910/fasteroff-ctest-01.log
```

All three exits 0. Actual cache flags were checked: `MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS=ON`, `MAKE_IT_STAND_BUILD_SPIN_IT_2014=ON`, `MAKE_IT_STAND_BUILD_SPIN_IT_FASTER=OFF`. Focused CTest **1/1 passed**, all nine groups passed, total **0.95 seconds**. This proves the new component builds/runs without Faster. Configure emitted existing CMake deprecation warnings; compile/test output had no warnings/failures.

Final `git diff --check` passed. Final main tracked diff SHA-256 remained exactly `93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`. Only the scoped files below are committed; unrelated controller documentation changes and all prior outputs remain intact.

| Integration log (under build/task-5h-evidence-20260910) | SHA-256 |
|---|---|
| full-build-01.log | 43b0f8e2b82366407da1f4fbd2b4380a12cebce1ad7bb942f0009e8b0495153d |
| full-ctest-01.log | ce506d59c10047a772dd7a95f48a7547745e57cdd680171a3ee9e7165193a623 |
| fasteroff-configure-01.log | 99a685426208814fd2efcaad50b6e0fb23a0778d4d8d448a56d64029ad357a9f |
| fasteroff-build-01.log | b6c6db9dc0d5a412babece7a3d48efb32beee076437f650f875d98c512ce2388 |
| fasteroff-ctest-01.log | 7d0a761e3399ac4db707296dbeb785d0431cf220ee4f0fddb2778179d60507b4 |

### Initial independent controller output check (pre-fix)

After the final main all-target build, the controller compiled a separate diagnostic executable with `/usr/bin/c++` against the stable core archive, into `build/task-5h-controller-8bSTgg/`, without touching project sources or the CMake build. I read the actual resulting `check.log`: two independent checks passed, with literal actual counts `[1728,2197,6084,218]`, complete actual connectivity/handle-map validation, and maximum six-query multilinear field error `1.3877787807814457e-17`. `run1.json` and `run2.json` are byte-identical. This remains SYNTHETIC_GRID_INTERPOLATION_NOT_BBW evidence.

| Controller evidence | SHA-256 |
|---|---|
| compile.log (empty, no diagnostics) | e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855 |
| run1.json and run2.json | b4193f3e510cc6f142765ec90e5e9c3f164d86e97e03cc8e9c65add3e7b1dd6e |
| check.log | 7d1f2ca9da83aa7b1d17ce8da9eebd0b1fa4c56c06c7e16819901da91b3d6017 |

### Focused original-log hashes

All paths below are relative to `build/task-5h-evidence-20260910/`; retained `/tmp/` originals have the same filenames/content.

| Log | SHA-256 |
|---|---|
| task-5h-red-build-20260910-01.log | b70196df98bbabae1cc66c26251fbaf86544940818916ed7005b55eae088c028 |
| task-5h-red-configure-20260910-01.log | 5750ab71c1cdccc8d399616d7bc5ff7eee03dca0a3dcdda8d01e346b0cefac8f |
| task-5h-red-build-20260910-02.log | df7e5fbe5e4f1ce48532a71e6ead1264385be309c13c2bad7028651644ef56a0 |
| task-5h-red-ctest-20260910-01.log | d1bae45864f615e11c36612fe540c416c45f1be017242f655c1dc263885c81b7 |
| task-5h-green-build-20260910-01.log | df7e5fbe5e4f1ce48532a71e6ead1264385be309c13c2bad7028651644ef56a0 |
| task-5h-green-ctest-20260910-01.log | 8c78715b2dfcb292a5d373c9fb1fb4376ad90f501171319bae113cfa605a26f4 |
| task-5h-capacity-red-build-20260910-01.log | 0ea92c9c1318cbe2a0534622ab018a30c2e775d93e97729d1fe7e7d2673883bf |
| task-5h-capacity-red-ctest-20260910-01.log | 771b45c4482ef03bbb40f8e5b8226ee7512f3a65594cd1fe50aa282dbfc48610 |
| task-5h-capacity-green-build-20260910-01.log | ddea5ded42588fc1bb9e187012f16c2e0e4de8c245d2932458a71c6ef902fbb4 |
| task-5h-capacity-green-ctest-20260910-01.log | 7894f72cc3fcaf5d06db27dc7b465d8cbab69ecf4f17b7512bf730f28055ac70 |

### Initial implementation artifact hashes (pre-fix)

| Artifact | SHA-256 |
|---|---|
| include/spin2014/SourceCageGrid.h | 3ad0e0989f486a8a6408bdee12b4948f8997a7dff86fe623cd47b7858584a950 |
| src/spin2014/SourceCageGrid.cpp | 0084902b1d78c24f7a3ae70a16f66123defc0213a4247363934cb07e71d5a88b |
| tests/spin_it_2014_source_cage_grid_regression.cpp | 0a2dd5899d1a7f8cced537d638e5a5fedac702996da3a836208d653193d49d99 |
| build/strict-v2-spin/libmake_it_stand_core.a | 419273419afb3c33209545aeb1b39f1b19fccd5c4b14608095aa8c725a6ac514 |
| build/strict-v2-spin/spin-it-2014-source-cage-grid-regression | 7b105c58615c14a47776e447657d77b7f568bc522f43338b5e8f667dd95b4111 |
| build/task-5c-no-faster-20260909/spin-it-2014-source-cage-grid-regression | e5edfb2ffaba84db41c952e659464b49cc75fac9d619728c43f8d96f84d6ec06 |

## Initial self-review and remaining scope (pre-fix)

Read the new header/source/tests and CMake diff. Confirmed selected-volume refinement, strict plane ordering, bounded candidate enumeration, lexicographic rather than linear-root tie priority, integer handle mapping, actual represented interval fractions, per-allocation capacity checks, and complete canonical hashing. Self-review found the temporary edge-set capacity gap; it received its own actual RED/GREEN cycle. No known correctness concern remains in the requested increment.

The new API/header/source/test are gated by 2014 AND exact geometry, outside Faster. The legacy deformation API and graph BBW are unchanged. Test code and API documentation explicitly avoid granting geometry seal, physical feasibility, solved-weight/KKT properties, optimization, or paper reproduction success. BBW discretization/solver and subsequent deformation integration remain separate work.

No helper/reviewer subagents were used. Implementation, configure and project build/test ownership remained with this implementer; the controller separately derived analytic acceptance and prepared an external diagnostic probe without modifying the production build. Main checkout was never written; tracked main diff SHA-256 matched `93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e` before work and during verification. Final preservation check is recorded with the integration results.

Scoped commit files: this report, `CMakeLists.txt`, `include/spin2014/SourceCageGrid.h`, `src/spin2014/SourceCageGrid.cpp`, and `tests/spin_it_2014_source_cage_grid_regression.cpp`. Controller notes/brief/progress and documentation audit changes are excluded; raw evidence remains under the isolated build evidence directory.

## Fix round 1 — collapse-proof finer-coordinate guard

Fix base: `4c6f89b8ef5a72a47e4ca6ac872ecf22dd122990`. The review's sole Important finding was verified before changing production. No upstream source-oracle code, legacy API, CMake registration, or controller documentation is changed in this correction.

### Actual reproduction and cause

Added three real source/cage/grid tests before implementing the fix. Each constructs an actual exact-audited `SourceSurfaceCage` with 216 selected cells, outside the catch that expects grid rejection. Thus an upstream fixture failure cannot count as grid RED. Each independently checks that all nine coarse x-planes are strictly increasing. Unit side, coarse level3/fine level4:

- Positive binade crossing: x-origin `nextafter(2^48,0)`. The actual upstream cage succeeds; seven adjacent canonical FMA fine-plane pairs coincide. The current grid constructor nevertheless returns successfully, causing the behavioral test to fail.
- Negative binade crossing: x-origin `-2^48-.5`. Actual upstream cage succeeds; zero fine intervals coincide. The equal-ULP grid nevertheless must reject under the new explicitly conservative support boundary.
- Positive aligned root: x-origin `2^48`. Actual upstream cage succeeds; zero fine intervals coincide. This explicitly tests/document the otherwise valid aligned equal-ULP inputs deliberately rejected by the guard.

The old test `step < spacing` allowed equality. Under nearest-even rounding, the exact preimage of a floating value can include both midpoint endpoints. Two exact FMA plane values separated by exactly one ULP can therefore round to the same value, notably when the origin has a half-ULP phase inherited from the smaller binade. This breaks distinct-node positions, positive cell widths and per-node interpolation identity.

The minimal production correction requires `step > maximumEndpointULP` (rejects `step <= spacing`). For either sign, throughout the root interval, the width of any rounding preimage is at most the largest endpoint ULP; at a binade boundary the width is no larger than that bound. Exact neighboring dyadic planes separated by strictly more than that width cannot round together. The power-of-two-scaled root step is exact unless subnormal rounding occurs; even then a rounded step strictly above the representable ULP bound implies the exact step is above it, by monotonic rounding. This gives the sign-general strict-increase condition required by represented-plane binary search. The public header explicitly documents conservative rejection of equality, even for otherwise valid aligned roots. Rejection remains `invalid_argument` with a numerical-unsupported diagnostic; no recentering, coarsening, relocation, clamping or changed FMA coordinates.

All accepted ordinary inputs retain the same domain, arrays, mapping, interpolation and identity serialization. No version-string change is made: this repairs the constructor's enforcement of its existing distinct represented-plane contract. Numerically marginal equal-ULP inputs now reject instead of publishing a grid. The upstream oracle's related guard is deliberately not edited; it remains separate parked work.

### Actual fix RED/GREEN evidence

Every command ran in `/data/xvxinyao/make-it-balance-paper-sim` with `set -o pipefail`, `ulimit -v 25165824`, `timeout 3600`, at most three compile jobs, sequential project builds/tests, and a unique original `tee` log in `build/task-5h-fix1-evidence-20260910/`. No full 201-test suite was rerun during the fix iteration: its earlier 201/201 result is pre-fix evidence only, explicitly labeled above.

RED commands, after adding the tests and before changing the production guard:

```bash
timeout 3600 cmake --build build/strict-v2-spin --target spin-it-2014-source-cage-grid-regression -j3 2>&1 | tee build/task-5h-fix1-evidence-20260910/red-build-01.log
timeout 3600 ctest --test-dir build/strict-v2-spin -R '^spin-it-2014-source-cage-grid$' -V 2>&1 | tee build/task-5h-fix1-evidence-20260910/red-ctest-01.log
```

Build exit0. CTest exit8 with three failed groups (all nine existing groups passed). Actual captured output:

```text
binade fixture 0: upstream valid 216 cells; independently collapsed fine intervals=7
FAIL positive_binade_crossing_must_not_create_zero_width_fine_cells: equal-ULP finer root must reject explicitly before publishing grid
binade fixture 1: upstream valid 216 cells; independently collapsed fine intervals=0
FAIL negative_binade_crossing_equal_ULP_requires_explicit_rejection: equal-ULP finer root must reject explicitly before publishing grid
binade fixture 2: upstream valid 216 cells; independently collapsed fine intervals=0
FAIL aligned_equal_ULP_root_documents_conservative_rejection: equal-ULP finer root must reject explicitly before publishing grid
source cage grid SYNTHETIC_REGRESSION: 3 failed groups; no BBW or physics authority
```

After the strict guard/header documentation change, covering GREEN commands:

```bash
timeout 3600 cmake --build build/strict-v2-spin --target spin-it-2014-source-cage-grid-regression spin-it-2014-source-surface-cage-regression spin-it-2014-source-surface-octree-regression -j3 2>&1 | tee build/task-5h-fix1-evidence-20260910/main-green-build-01.log
timeout 3600 ctest --test-dir build/strict-v2-spin -R '^spin-it-2014-source-(cage-grid|surface-cage|surface-octree(-level9)?)$' -V 2>&1 | tee build/task-5h-fix1-evidence-20260910/main-green-ctest-01.log
```

Both exit0: **4/4 CTests passed**, including source cage, source octree, and source octree level9; total23.47 seconds. The grid regression passes all **12 groups** including the three new binade cases, in1.20 seconds. The old upstream regression prints `faster=1 exact=1`; its level9 result remains 2,097,993 stored nodes and the same partition hash `652de7d6b99c64445640902e891a98f85f76c4a8b5a93b6ef34471954a5bc137`.

FasterOFF covering commands, after the main covering run finished:

```bash
timeout 3600 cmake --build build/task-5c-no-faster-20260909 --target spin-it-2014-source-cage-grid-regression spin-it-2014-source-surface-cage-regression spin-it-2014-source-surface-octree-regression -j3 2>&1 | tee build/task-5h-fix1-evidence-20260910/fasteroff-green-build-01.log
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 -R '^spin-it-2014-source-(cage-grid|surface-cage|surface-octree(-level9)?)$' -V 2>&1 | tee build/task-5h-fix1-evidence-20260910/fasteroff-green-ctest-01.log
```

Both exit0: **4/4 CTests passed**, grid **12/12 groups** passed in1.20 seconds, total23.61 seconds. Upstream octree output explicitly reports `faster=0 exact=1`; the same level9 partition hash is preserved. No configure step or dependency update was needed. Build and test outputs contain no new diagnostics/warnings.

### Fix evidence hashes and final review

Original, retained log paths are under `build/task-5h-fix1-evidence-20260910/`:

| Log | SHA-256 |
|---|---|
| red-build-01.log | 0ea92c9c1318cbe2a0534622ab018a30c2e775d93e97729d1fe7e7d2673883bf |
| red-ctest-01.log | b3c410a64a25877badc5c348d7dff2181aba5efce9896c4c5d413cc4100f08e1 |
| main-green-build-01.log | 48f3043a9fad67f5c5a55824916f68901c43e354bcb4fe12dbf222c595bb4a14 |
| main-green-ctest-01.log | 840ecefa00c042940ab70582a108f2840071ee6d30e1a3f9b94ea53bf2bf00cb |
| fasteroff-green-build-01.log | f88c350b0a6484d214576c134f6e44f2783d8c6208d628e8a11f725823432baa |
| fasteroff-green-ctest-01.log | cfb50af28280454a745231e94e94b4b4b9312bb1215d42b728f006cf0ab415f0 |

Corrected artifact hashes (supersede the pre-fix artifact hashes above):

| Artifact | SHA-256 |
|---|---|
| include/spin2014/SourceCageGrid.h | 7acc63d8862e611d208d77a9175f8cb34f3dbc60bd520b32f2d8dd40a947964f |
| src/spin2014/SourceCageGrid.cpp | b9e21d6468a2bd3ea9bc2303637f829e29c10fcaef4eff8a00d3ed0a71c6f9d9 |
| tests/spin_it_2014_source_cage_grid_regression.cpp | 0ca7e7b92141d906d2efe15ca8bb14e89533bebca6aa6f392f046032960a77e0 |
| build/strict-v2-spin/libmake_it_stand_core.a | b85b3d04a7af9b86a811e6d02b9dd0c06fc0a7bbdd27520b88d468c1080ca6dd |

Fix status: **DONE**. Self-review confirms the single stricter guard plus its diagnostic/docs and three regressions address the reported collapse cause. The conservative equal-ULP rejection is deliberate and visible in the API; no accepted grid is silently modified. No solver/scientific authority changes. The upstream numerical-guard issue remains outside this bounded fix. Fresh full-suite verification/re-review and regeneration of the controller probe are left to the controller's next explicitly scheduled stage; no pre-fix artifact is presented as post-fix evidence.

`git diff --check` passed; main tracked diff SHA-256 is still `93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`. Only the grid header/source, regression test, and this appended report are committed in fix1; controller documentation and all previous outputs/logs remain untouched. Source/build/test ownership is released with the fix handoff.
