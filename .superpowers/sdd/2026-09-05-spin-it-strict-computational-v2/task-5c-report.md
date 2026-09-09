# Task 5c implementer report — original-surface boundary subtrees

Base: `c07663bba4bb8b82b3f8cde9e31e4eac302dafe0`.
Status: bounded Task 5c foundation implemented and verified; controller's fresh
independent review remains pending. This does not complete Task 5.
Worktree: `/data/xvxinyao/make-it-balance-paper-sim`, existing isolated branch
`spin-it-strict-computational-v2`. Sole production/configure/build writer; no
subagents or reviewer were dispatched. Controller owns independent notes and
subsequent review. No main-checkout writes, merge, push, tag, or historical
experiment edits. Overall strict reproduction remains `PARTIAL_BLOCKED`.

## Implementation and scope

Added public `integrateSurfaceRawMoments` next to the existing differential API.
Both call a private shared Algorithm 1 polynomial accumulator; a null Jacobian
destination bypasses allocation and every Algorithm 2 derivative operation. The
differentiator preserves its validation, signed-component, density, ordering and
Jacobian semantics. Surface moments refer to the unapproximated original
polyhedron, evaluated with binary64 polynomial arithmetic, not rational moments.

Added `SourceSurfaceOctree.h/.cpp`, a copied-source immutable PIMPL value object
containing root/depth/node configuration, copied original outer/design V/F,
surface moments at density one, stored hierarchical cells, approximate carveable
moments, unresolved volume, counts by node kind and level, and three SHA-256s.
No occupancy input or reconstruction, lattice adapter, dense export, Python API,
solver/cage changes, repair, or optimizer was added. No private helper file was
needed: pinned CGAL templates compile within the focused implementation file.

The existing exact material-region constructor performs the exact boundary audit
and signed nesting validation; an exact relation establishes design subset outer.
Nef is used only for legality, never moments. Original signed disconnected
components, reversed cavities and positive islands within cavities are supported.
Invalid or incomplete geometry fails construction.

Cell classification uses original design triangles. Inclusive original-coordinate
bounding boxes and inherited intersecting triangle lists form a sparse broad
phase. The actual pinned CGAL API is `do_intersect(Triangle_3<Gmpq>, Bbox_3)`;
this revision has no triangle/Iso_cuboid overload. Its implementation in
`Intersections_3/include/CGAL/internal/Intersections_3/Bbox_3_Triangle_3_do_intersect.h`
and the corresponding plane helper was read: box endpoints become exact kernel
coordinates before plane and SAT axis operations; strict separating comparisons
preserve closed touching. Thus enclosed triangles are not lost. For boundary-free
cells only, AABB-backed exact `Side_of_triangle_mesh` parity establishes material
membership, following prior signed nesting validation.

Intersecting cells recurse through all eight deterministic children. At finest
depth, closed boundary-touch cells contribute zero carveable moments and their
whole volume as uncertainty (`IMPLEMENTATION_DEFAULT`). Parent contributions sum
actual children. Full/empty terminal descendants have analytic nonallocating
queries. This is conservative geometry discretization; volume uncertainty is not
a ten-moment relative-error certificate, cut-cell integration, recovered author
sampling rule, continuous infeasibility proof, or completed adaptive optimizer.

Finite roots, containment, depths 0..20, dyadic collapse and volume/moment overflow
or underflow are checked. Node growth checks the configured budget before root or
eight-child allocation, caps reserve growth, and throws `RESOURCE_LIMIT` when
exceeded; no valid partial tree is returned. The immutable configuration defaults
to 5,000,000 nodes. Copied source identities use fixed-width big-endian field
encoding and binary64 bits; partition identity binds actual tree scientific fields.
Per controller clarification, the **requested maximumNodes guard is included in
configuration identity** because it controls construction success/failure.
Measured runtime/RSS and timestamps are excluded. Public algorithm, units,
provenance and limitations are in `docs/spin_it_2014_source_surface_octree.md`.

## Tests and TDD evidence

All log paths below are relative to this report's directory. The preexisting
controller baseline is preserved at `build/task-5c-baseline.log` (3/3 existing
adaptive-domain, audited solver and surface-differential tests). No full baseline
rerun was substituted for the requested final full suite.

- `task-5c-value-red-01.log`: newly written analytic/value/allocation test failed
  compiling because the public integrator was absent. This is API absence RED,
  not a runtime numerical failure.
- `task-5c-value-green-build-01.log` and `task-5c-value-green-test-01.log`:
  value, allocation and existing differential tests passed 3/3. Differential
  finite-difference maximum error remained `4.45924e-12`.
- `task-5c-octree-red-01.log`: requested hierarchy/geometry/query tests failed
  compiling because `SourceSurfaceOctree.h` did not exist.
- `task-5c-octree-green-build-01.log`: retained first implementation compile
  failure. The pinned Eigen lacks vector range iteration and pinned CGAL lacks
  the attempted triangle/Iso_cuboid overload. Root-cause inspection led to indexed
  moment iteration and the pinned exact triangle/Bbox overload; no dependency
  change or substitute predicate.
- `task-5c-octree-green-build-02.log` and `task-5c-octree-green-test-01.log`:
  initial complete geometry regression passed (1.27 seconds, 24,568 KiB RSS).
- `task-5c-guard-red-build-01.log`, `task-5c-guard-red-test-01.log`: added the
  clarified budget identity regression, which failed with
  `configured node guard must bind scientific identity`. Added the input guard to
  config hashing; its GREEN build is `task-5c-guard-green-build-01.log`.
- `task-5c-focused-green-01.log`: geometry including depth 9 and value/differential
  passed; allocation harness failed (4/5 overall) under the newly inherited
  24 GiB hard address-space limit. The child tried to set its hard limit to
  infinity. Fixed the test harness to preserve the inherited hard limit using
  `getrlimit`; the same 32 MiB additional allowance and million-vertex input remain.
  The failed attempt is retained. No numerical tolerances or geometry were retuned.
- `task-5c-focused-build-02.log`, `task-5c-focused-green-02.log`: 4/4 focused
  value/allocation/differential/geometry tests passed after the harness correction,
  including added depth-20 and non-default root coverage. Depth-9 was already
  green and was exercised again in the independent feature build and full suite.

Independent moment expectations use closed-form cuboid integrals and simplex
barycentric expectations `E(lambda_i)=1/4`,
`E(lambda_i lambda_j)=(1+delta_ij)/20`. Fixtures include a translated asymmetric
tetrahedron with nonzero cross moments, translated anisotropic box, reversed inner
cavity, density/translation covariance, and malformed/nonfinite/overflow inputs.
Differential agreement is a compatibility check only.

The allocation regression forks real subprocesses after allocating one million
input vertices. A measured current process address-space size plus 32 MiB is
imposed with `RLIMIT_AS`. Value-only succeeds; the dense 10x3N differential's
240 MB Jacobian throws `bad_alloc`. This checks a real allocation/resource
consequence rather than a claimed allocation counter. The intentionally
unreferenced extra vertices are valid for the signed-chain integrator contract;
the geometry oracle separately rejects unreferenced mesh vertices.

Geometry tests cover all ten analytic dyadic moments at levels 5 and 9, original
volume 1/8, conservative enclosure, root moments, every stored parent/child and
query contribution, lineage/order, virtual terminal descendants, copied source
V/F and configuration mutation, repeated hashes, changed sources/guards/depth,
reversed nested cavity `[7/16,9/16]^3`, nested material islands, disconnected small
islands, an asymmetric non-axis-aligned tetrahedron, and tiny enclosed cavity/island
triangles that defeat sample-only classification. Bad open/flipped/misnested/
self-intersecting/unreferenced/outside/NaN geometry, design in an outer cavity,
invalid depths/coordinates/roots, numerical collapse and resource budget failures
are rejected. Depth zero and a sparse actual depth-20 case are exercised.

## Commands and resource evidence

Configure and every native build/test ran in the specified worktree. Shells for
the geometry build/test batches used `ulimit -v 25165824` (24 GiB per-process
address-space cap); tests used `timeout 3600`. Compilation used at most `-j3`,
and test/experiment batches were sequential. `/usr/bin/time -v` records actual
elapsed time and maximum RSS; no aggregate cgroup-memory claim is made. Available
space was 1.4 TiB, exceeding the 100 GiB start/20 GiB stop limits. No files were
deleted to make space.

The commands below describe the preserved exact build/test invocations; each
was redirected to its unique log above or below. Initial RED commands deliberately
returned nonzero, while the later corrected checks must pass.

```sh
cmake -S . -B build/strict-v2-spin
cmake --build build/strict-v2-spin --target spin-it-2014-surface-value-regression -j3
cmake --build build/strict-v2-spin --target spin-it-2014-surface-value-regression spin-it-2014-surface-differential-regression -j3
ctest --test-dir build/strict-v2-spin -R '^spin-it-2014-surface-(value.*|differential)$' -V
cmake --build build/strict-v2-spin --target spin-it-2014-source-surface-octree-regression -j3
timeout 3600 /usr/bin/time -v build/strict-v2-spin/spin-it-2014-source-surface-octree-regression
timeout 3600 /usr/bin/time -v ctest --test-dir build/strict-v2-spin -R '^spin-it-2014-(source-surface-octree.*|surface-value.*|surface-differential)$' -V
cmake --build build/strict-v2-spin --target spin-it-2014-source-surface-octree-regression spin-it-2014-surface-value-regression -j3
timeout 3600 /usr/bin/time -v ctest --test-dir build/strict-v2-spin -R '^spin-it-2014-(source-surface-octree$|surface-value.*|surface-differential)$' -V
cmake -S . -B build/task-5c-no-faster-20260909 -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DMAKE_IT_STAND_BUILD_SPIN_IT_2014=ON -DMAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS=ON -DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=OFF -DMAKE_IT_STAND_BUILD_REACHABILITY=OFF -DMAKE_IT_STAND_BUILD_PYTHON=OFF -DMAKE_IT_STAND_BUILD_VIEWER=OFF -DFETCHCONTENT_SOURCE_DIR_HIGHS=/data/xvxinyao/make-it-balance-paper-sim/build/strict-v2-spin/_deps/highs-src
/usr/bin/time -v cmake --build build/task-5c-no-faster-20260909 --target spin-it-2014-source-surface-octree-regression spin-it-2014-surface-value-regression spin-it-2014-surface-differential-regression -j3
timeout 3600 /usr/bin/time -v ctest --test-dir build/task-5c-no-faster-20260909 -j1 -R '^spin-it-2014-(source-surface-octree.*|surface-value.*|surface-differential)$' -V
/usr/bin/time -v cmake --build build/strict-v2-spin -j3
timeout 3600 /usr/bin/time -v ctest --test-dir build/strict-v2-spin -j1 --output-on-failure
/usr/bin/time -v cmake --build build/spin-2014-enabled-contract --target spin-it-2014-surface-value-regression spin-it-2014-surface-differential-regression -j3
/usr/bin/time -v ctest --test-dir build/spin-2014-enabled-contract -j1 -R '^spin-it-2014-surface-(value.*|differential)$' -V
sha256sum -c .superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5c-primary-isolation-sha256.log
git diff --check
git diff c07663b -- experiments
```

Feature isolation uses the existing pinned CGAL and existing HiGHS source checkout,
with no dependency download/upgrade. Its actual linked runtime prints
`faster=0 exact=1`, and passes all five selected native CTests, including depth 9.
The isolation build took 92.69 seconds and maximum compiler RSS 1,265,336 KiB
(`task-5c-isolation-build-01.log`). The first isolation test batch passed 5/5 in
21.17 seconds (`task-5c-isolation-test-01.log`). A final test-only rebuild and run
retained in `task-5c-isolation-final-build.log` and
`task-5c-isolation-final-test.log` includes precise moment/level logging and the
explicit saved-query-after-buffer-mutation assertion.

The real depth-9 tree contains 2,097,993 nodes: 262,249 internal, 786,448 unresolved
boundary, 506,360 full and 542,936 empty. Actual maximum visited level is 9.
Each public cell is 144 bytes, totaling 302,110,992 bytes of stored scientific
cells (vector capacity and other allocations are additional). Initial primary
geometry construction took 20.0084 seconds, peak RSS 349,060 KiB. Final isolation
batch peak RSS is 349,724 KiB and elapsed time 21.30 seconds. The independent
analytic closed-boundary moments pass normalized error <=1e-10. The accepted
volume is `2048383/16777216`; unresolved volume is `49153/8388608`.
Primary and Faster-OFF native builds both produce partition SHA-256
`652de7d6b99c64445640902e891a98f85f76c4a8b5a93b6ef34471954a5bc137`.

The final main all-target build passed in 26.53 seconds, maximum RSS 1,756,776 KiB
(`task-5c-final-all-build.log`). The final main full registered suite passed
**194/194, zero failed/skipped/Not Run**, in 340.56 seconds, maximum RSS
4,253,124 KiB (`task-5c-final-full-ctest.log`). It was run once after the focused
iteration. No source/configuration/test changes followed this full run.
Its actual exact-OFF `spin-2014-dependency-contract` native fixture passed in
94.81 seconds. In addition, the **changed surface value and differential paths**
were linked and run in that real 2014 ON / exact OFF / Faster OFF configuration:
3/3 native tests passed, including the memory-limited value allocation regression.
The focused exact-OFF build took 1.43 seconds (321,716 KiB peak compiler RSS) and
tests took 0.01 seconds (18,492 KiB RSS), recorded in `task-5c-exact-off-build.log`
and `task-5c-exact-off-test.log`. Actual cache flags, link command, dynamic
dependencies and artifact hashes are in `task-5c-exact-off-identity.log`.
The exact-OFF value executable SHA-256 is
`3937b70f2ef32d72d4dfd04318b27f7b03775835671af138b2e790cb60fd9d83`, differential
executable `36b8379cbdfc191f6258a9891e1a379fe40b6310f93f6f281345067b8530dca6`,
core archive `99b19de105654c62d6373891f23622788918158ac89cd8f216f64a32f1b2818d`,
cache `72fcd9470ddc0bed149523aeeaf10ce8389d446a6cf4fe838c28433f0891a61d`.

`task-5c-primary-isolation-sha256.log` binds every scoped source/header/test/doc,
the primary and isolation executables, core archives, primary native Python
extension, caches, and linked GMP/MPFR libraries. Actual link commands, dynamic
dependencies and feature flags are in `task-5c-linked-feature-identity.log`.

| Final artifact | SHA-256 |
|---|---|
| Primary source octree executable | `ddd076c9c97b256744631a14ed7a67bd93390e1fc5e26bdaf7b4e9cd3a813b1b` |
| Primary value executable | `fc15ff6cefbba4eeada92bd55e52adf5dcb8d4eb39f9bfb58642a8f0cf7f80f4` |
| Primary differential executable | `151d1bc7af69b740bd24635264c12f280be0eacae8577583001ad3b431fb6e47` |
| Primary core archive | `98c2e229331219341efd6c5ff0ca0fe86ebf8933ec661426b32b42cab318bc5c` |
| Primary Python native extension | `f0c8e1f02f41a75de28809832b3a6ea5773d6595615482b738806c03910c1d02` |
| Primary CMake cache | `2c3f5a72c82d0ba5d8b5d0cd8869d13d89ea8a21e0daa405e88242e1778e517c` |
| Faster-OFF source octree executable | `ef21bb37bed16ce5ededaa08a4ff59ebf3560ff8989b16a32f6bf7b21e96d93a` |
| Faster-OFF value executable | `3937b70f2ef32d72d4dfd04318b27f7b03775835671af138b2e790cb60fd9d83` |
| Faster-OFF core archive | `70fe78d883571c09144f4035d59ffc91b9b5b1b31d83bd22c46adbfef6df835c` |
| Faster-OFF CMake cache | `e57a7b85de312d5a27cb6979c88b67c1b2806d43dcb8f884a4f711c2650b01e2` |

Primary flags remain Faster ON, 2014 ON, exact ON, Python ON, reachability ON,
paper Python tests ON, viewer OFF and PyBullet OFF, with compiler observation OFF.
Isolation uses Faster OFF, 2014 ON, exact ON, Python/reachability/viewer OFF.
Pinned dependency warnings (libigl's CMake deprecation and Boost bind-placeholder
pragma notes) are retained in raw configure/compile logs. No dependency was changed
to remove warnings. Final incremental main all-target output has no warning/error
diagnostic. Whitespace and historical tracked experiment diff checks are clean.
After the full suite and exact-OFF check, all 22 recorded source/primary/isolation
and linked-library identities passed `sha256sum -c`
(`task-5c-final-sha256-check.log`). Historical tracked experiment diff remains
empty (`task-5c-historical-experiments-check.log`).

## Scoped files and implementer self-review

Production/configuration: `CMakeLists.txt`,
`include/spin2014/SurfaceRawMomentDifferential.h`,
`src/spin2014/SurfaceRawMomentDifferential.cpp`,
new `include/spin2014/SourceSurfaceOctree.h` and
`src/spin2014/SourceSurfaceOctree.cpp`.
Tests: new `tests/source_surface_test_utils.h`,
`tests/spin_it_2014_surface_value_regression.cpp`,
`tests/spin_it_2014_source_surface_octree_regression.cpp`.
Documentation: new `docs/spin_it_2014_source_surface_octree.md`, this report, and
the explicitly listed Task 5c raw evidence. Controller's brief/math notes are not
staged by the implementer.

Self-review read the complete production/header/test/configuration diff and new
files, checked signed material validation before parity, closed-box exactness,
candidate inheritance for enclosed triangles, numerical/root failure behavior,
preallocation budget arithmetic, vector reallocation-safe lineage, actual child
sums and virtual descendants, copied-source/config const ownership, canonical
hash field coverage, independent ten-moment expectations, and feature guards.
The unchanged differential finite-difference regression protects shared-kernel
compatibility. Realistic wrong-boundary, sample-only, cavity-as-full, altered
source/config, missing-budget and missing-Jacobian-bypass changes are caught by
observable tests. Test-only process-limit correction is explained above; no
geometric tolerances or inputs were weakened. No independent reviewer was
self-dispatched; controller review follows this scoped handoff.

Remaining limits: explicit conservative sampling policy, binary64 moment arithmetic
and conservative coordinate-spacing rejection; node budget does not bound the
separate exact input-mesh audit complexity; no original-paper experiment or
complete Task 5, solver, cage, feasibility or release claim follows.
