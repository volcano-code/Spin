# Task5j source cage FEM and complete nodal constraints

Date: 2026-09-10. Base: `011f1677766369bed65f67d21455771f8f9919f2`.
Worktree: `/data/xvxinyao/make-it-balance-paper-sim`, existing
`spin-it-strict-computational-v2` branch. Implementation commit:
`22cd391e7836d5e81a133a294dd051d6443e7bbc` (four explicitly scoped files).
Final validation: all-target build passed; full suite 202/202 passed in 480.33s;
FasterOFF focused suite 3/3 passed in 5.12s. No merge, push, tag, release, GUI, helper,
legacy solver mutation, source reclassification or whole-plan redispatch.

## Interface and implementation

New `include/spin2014/SourceCageFem.h` and `src/spin2014/SourceCageFem.cpp`
implement one immutable native module, compiled only with 2014+exact and
independently of Faster. `grid()` retains the copy-owned immutable upstream
object; config and all output getters return const references.

- `tetrahedra()` and `tetrahedronVolumes()`: six conforming Freudenthal tets per
  actual fine cell, x-bit-first local diagonal0--7, explicitly positive fixed
  orientation. The split is an `IMPLEMENTATION_DEFAULT`, not recovered author
  meshing metadata. Every element uses actual represented positions.
- `stiffness()`: row-major sparse linear-tetrahedron FEM K; `lumpedMass()`:
  positive diagonal computational M. Deterministic long-double element
  arithmetic/accumulation checks finite values and binary64 overflow/underflow.
- Controller's pre-implementation amended mass contract is followed:
  `SourceCageFemMassRule::RESTRICTED_GRID_VORONOI` is the default, each actual
  represented brick contributing V/8 per corner; `BARYCENTRIC_TET` is the
  explicit Vtet/4 alternative. Restricted Voronoi clips nearest-node dual
  volumes to the actual selected-cell union, including missing cells. These
  masses are not physical density, material moments or certificates.
- `boundaryConstraints()` returns sorted `BoundaryRow{nodeIndex, entries}`;
  entries are sorted `BoundaryEntry{handleIndex,weight}` with at most three
  nonzeros. Coarse triangle integer face boxes enumerate actual fine boundary
  nodes; long-double barycentrics use represented physical coordinates. Exact
  zeros are omitted, no negative weight is accepted, and there is no clamp,
  normalization, nearest fallback or tolerance-based geometric snapping.
- Duplicate rows keep the first deterministic row after verifying consistent
  sparse coefficients; partition and local position reconstruction are also
  checked. The declared validation tolerance is128*binary64 epsilon. All
  exposed fine-cell faces receive an independent completeness pass. No
  boundary-node-by-all-handles or node-by-all-triangles dense structure exists.
- `applyQ(Nnodes x K)` applies Q=K^T M^-1 K without constructing Q. Symmetric
  bit-identical K entries permit row access for both K and K^T. `energy()` is
  0.5*sum_columns(x^TQx), evaluated by the nonnegative squared-gradient sum;
  `applyQ()` is its gradient. Both accept finite fields with K>0, preserve
  arbitrary field values, and validate shapes and arithmetic.
- Limits independently cover tets, all16 potential assembly coefficients per
  tet, unique sparse entries, boundary nonzeros, and two Nnodes*K evaluation
  arrays. Assembly streams into a bounded unique-entry map rather than
  allocating16*Ntet triplets. Integer/container/byte validation precedes
  allocation; resource violations throw `length_error` containing
  `RESOURCE_LIMIT`, invalid configuration/fields throw `invalid_argument`, and
  unsupported arithmetic throws `overflow_error`. Storage follows actual
  sparse domain size, including level20 with only8cells/27nodes.
- `sourceSha256`, `gridConfigSha256`, `gridOutputSha256` preserve upstream
  identities. Config hashing binds grid config, mass/tet/energy conventions,
  numerical tolerance and every budget. Output hashing additionally binds
  upstream output/source and every tet, volume, M value, sparse K coefficient
  and boundary row. Evaluated fields change results, not immutable identity.

The complete constraints are NODAL conditions. Existing trilinear transfer is
unchanged and can differ between nodes from coarse triangle functions. The
native analytic example measures0.125 transfer against true triangular
weight0 on the diagonal. Neither a BBW solve nor convergence/KKT, feasibility,
geometry seal, physical simulation or original-case reproduction is claimed.

## Test-first and focused evidence

All logs are retained under the fresh unique directory
`build/task-5j-evidence-20260910-YGUc1N/`. Commands below ran from the worktree
with `set -o pipefail`, `ulimit -v 25165824`, `/usr/bin/time -p`, and
`timeout 3600`, with combined stdout/stderr piped to the named log. Builds used
`-j3`; heavy commands were sequential.

1. Test/CMake registration was written before the production header/source.
   The first attempt against the stale existing Makefile had no new target
   (`red-missing-api-build.log`, exit2,0.00s); this is administrative evidence,
   not the required API RED.
2. `cmake -S . -B build/strict-v2-spin`, `red-configure.log`, exit0,0.97s.
3. `cmake --build build/strict-v2-spin --target spin-it-2014-source-cage-fem-regression -j3`,
   `red-missing-header-build.log`, exit2,0.23s. Actual error:
   `fatal error: spin2014/SourceCageFem.h: No such file or directory`.
   This is the genuine missing-API RED, before implementation existed.
4. Both mass-rule tests and other behaviors were added before production
   implementation. Initial implementation compiled successfully
   (`green-configure-01.log`,0.98s; `green-build-01.log`,26.15s).
5. `ctest --test-dir build/strict-v2-spin -R '^spin-it-2014-source-cage-fem$' -V`,
   `focused-01.log`:1/1 passed,9groups,1.80s wall.
6. Self-review strengthened represented-coordinate coverage with a large
   offset/non-dyadic root producing materially unequal actual fine intervals,
   and changed the test reconstruction check to translation-stable local
   differences. No production fix was needed. `shifted-build-02.log`:5.04s;
   `shifted-focused-02.log`:1/1 passed,10groups,1.82s wall.
7. An upstream-valid subnormal-volume grid demonstrates explicit rejection
   of unrepresentable FEM mass contributions under both rules.
   `underflow-build-03.log`:5.06s; `underflow-focused-03.log`:1/1 passed,
   11groups,1.92s wall. These later tests are additional coverage, not a claim
   of an additional pre-implementation behavioral RED.

The eleven native groups substantively cover the1728cell/2197node fixture,
10368tets,432coarse triangles/218handles,866boundary/1331interior nodes;
positive orientations/volumes, conforming face multiplicity and exact union
volume to roundoff; independent literal brick K and both M formulas; natural
boundary affine energy and interior patch consistency; dense-small Q and
energy, three-step finite differences and multiple columns; translation/scale
covariance and actual unequal represented intervals; midpoint edge/diagonal
and three-corner face rows, all handles one-hot and triangle support; both mass
rules on cavity/tunnel/disconnected domains and component nullspaces;
level20 sparse operation, individual caps/overflow/shape/nonfinite failures;
deterministic hashes and exact equality of boundary/K under the two mass rules.

The independent brick formula starts from the reference right tetrahedron
gradients(-1,-1,-1),(1,0,0),(0,1,0),(0,0,1), volume1/6, mass1/24, and K=Gram/6;
its mapped path-simplex sum yields a literal cube edge stencil independently
of production connectivity/assembly. The eight barycentric cube corner masses
are1/4 at0and7 and1/12 elsewhere. Voronoi dual volumes are separately checked
with Cartesian neighbor-midplanes, including the unequal-interval fixture.

## FasterOFF and final validation

Existing dedicated build profile `build/task-5c-no-faster-20260909` was retained:
Release,2014+exact ON, Faster OFF, Python ON, reachability OFF, viewer OFF.

```sh
cmake -S . -B build/task-5c-no-faster-20260909
cmake --build build/task-5c-no-faster-20260909 --target spin-it-2014-source-cage-fem-regression spin-it-2014-source-cage-grid-regression spin-it-2014-source-surface-cage-regression -j3
ctest --test-dir build/task-5c-no-faster-20260909 -R '^spin-it-2014-source-(cage-fem|cage-grid|surface-cage)$' -V
cmake --build build/strict-v2-spin -j3
ctest --test-dir build/strict-v2-spin --output-on-failure
```

- `faster-off-configure.log`:exit0,0.77s.
- `faster-off-build.log`:exit0,3.12s.
- `faster-off-ctest.log`:3/3 actual tests passed, no skips/NotRun,5.12s wall.
  FEM11groups passed in1.93s, grid1.17s, surface cage2.02s.
- `final-all-target-build.log`:exit0,27.74s, including the Python module and
  every previously configured target. Main profile RelWithDebInfo,
  2014+exact+Faster+reachability+Python ON, viewer OFF.
- `final-full-ctest.log`:exit 0, all 202/202 actual tests passed, no skips or
  NotRun, 480.33s wall. Previous 201 tests plus the new
  `spin-it-2014-source-cage-fem` test, labeled `synthetic-regression`. The final
  main FEM test passed in 1.95s; grid 1.19s; surface cage 2.01s.

Existing CMake/libigl minimum-version deprecation warnings remain in configure
logs. Boost1.82/CGAL are unchanged. The new source and focused/all-target build
logs emitted no compiler warnings; no dependency update was attempted.

## Independent controller cross-check and source bindings

Controller independently derived an exact Fraction oracle and compiled its
standalone public-API probe against the stable main archive, without modifying
the shared build. Its small27node fixture checks all729 Kentries, both M,
Qx/energy,48tets and all26boundary rows; both rules passed, and repeated JSON
exports were byte-identical. These are controller-owned evidence under
`build/task-5j-controller-kb6lIC/`, not substitute native expected values.
For source cube `[.3,.4]^3`, cage level 2, fine level 3, its actual API hashes
are:

| Mass rule | Configuration SHA-256 | Output SHA-256 |
| --- | --- | --- |
| RESTRICTED_GRID_VORONOI | `c1bfbf2d3b462a9226712b0ab9f553c1df4be2f79493bdca407b067acec04e1f` | `9e81be494dac2c163608bf9fb054abd8e43e46c86794bda7896bc47cad98866c` |
| BARYCENTRIC_TET | `1d66a79a9adad9aa0d43ae5d7065ef0f0f80d8e0c3a347608a07709fa2c68a72` | `d75a3a0a094e6205a4458586be57699af03a64bb7e2fbb28a9fb0a639a9ae6a0` |

Final source SHA-256 bindings:

```text
2206eb4eda4aab8bf348a5202e4c606467a318d6f9d6b3dc8a1249bb3e55f531  src/spin2014/SourceCageFem.cpp
029ecfe50bdf2b8cfaac24ad57c8d52fa37bf83d6b00078c3e4f8d6bf1e7f09f  include/spin2014/SourceCageFem.h
668ac3630b620ca06d7c1c3e1887a76b098c52899662e4af5f6f422ed8862442  tests/spin_it_2014_source_cage_fem_regression.cpp
464af02e50757e72b128c3afc06f4a3988a2f378c48abf3d29d70315790547a7  CMakeLists.txt
```

Final binary SHA-256 bindings (main archive unchanged since controller probe):

```text
dff5252cf596ce549f58653a54cb759e7b7dec836ce542d814529caa5a228e4c  build/strict-v2-spin/libmake_it_stand_core.a
f7480bfd1c05aa7724b5b4403092b9fb10a6532ec6ecb5bcf505f81f1d54f977  build/strict-v2-spin/spin-it-2014-source-cage-fem-regression
c8e9eab651d80faa4dd1b8165771030509ddc3f348f8e08c0a24185e63151883  build/task-5c-no-faster-20260909/spin-it-2014-source-cage-fem-regression
```

Final log SHA-256 bindings:

```text
10fb8927d58abc670cb87d24df90bf651a42127d6b77f4146337bf94f429365f  red-missing-header-build.log
91969eee8d2b645bcaaafa9dc82305acfbcd23c016a52a20cc023e16dc492964  faster-off-ctest.log
5f86b08941931db1657fa8a0814c3804b22fdd06a5ce502895eb22619ca6e66d  final-all-target-build.log
f0d53b5bbdc90fe5754921cc53417a3fbbb4665d195e61949e3d80e8d7cb4db5  final-full-ctest.log
```

All four logs are in `build/task-5j-evidence-20260910-YGUc1N/`.

## Self-review and limitations

Reviewed the amended brief against implementation and substantive assertions:
all6tets/cell and positive orientation, both explicit mass rules, symmetric
sparse stiffness, full nodal boundaries, resource checks before allocation,
immutable deterministic identities, sparse high-level domain, no dense Q or
nodal-by-handle matrix, no source/grid/legacy semantic changes. Small dense
matrices occur only in the independent27node test reference. The enum rejects
unknown values; shared-boundary duplicates validate rather than overwrite.
Mass-rule differences preserve K and boundary coefficients exactly.

The actual interval test could have been too weak at ordinary offsets; the
added2^45-scale shifted1.7root fixture directly distinguishes actual unequal
intervals from nominal spacing. Native underflow coverage constructs the
valid upstream grid outside the expected FEM rejection catch, so an upstream
failure cannot falsely satisfy the test. No behavioral failure was concealed.

No subagents/reviewers were spawned. TDD guided the genuine missing-header RED;
verification-before-completion requires the final suite evidence before the
scoped implementation commit. The controller retains its docs and independent
checker ownership. No cumulative progress file or historical evidence was
edited. No author-exact meshing/clipping, successful BBW solve, material or
original-case authority follows from this prerequisite.
