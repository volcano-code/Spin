# Task5i implementation report

Status: **DONE**. Base `a22de6b69b0dd0d550ad57157fbd9bac834cac60`, branch
`spin-it-strict-computational-v2`. This is a bounded upstream source-octree
numerical-admission correction. It is not a solver change, an infeasibility
certificate, a BBW result, or an original-paper reproduction case.

## Root cause and actual native RED

`SourceSurfaceOctree.cpp::validateConfig` computed the finest dyadic step and
the largest binary64 spacing at the represented root endpoints, but rejected
only `step < spacing`. At equality, the validator then checked only finest cell
`(0,0,0)`. That first interval can be positive even though later FMA planes
coincide: a rounded value's exact preimage can include both midpoint endpoints
under ties-to-even, and crossing a binade changes the surrounding spacings.
Therefore one local cuboid check does not establish strict ordering of every
implicit dyadic plane.

Before any production edit, the regression was extended with separate closed,
strictly nested outer and design boxes. For each root, the identical geometry
was first constructed successfully at depth 3, where the step is strictly wider
than the maximum endpoint ULP. Positive outer/design masses and actual depth-3
traversal establish that these are valid native source inputs rather than empty,
open, invalid, or duplicated-formula fixtures.

At depth 4, the test independently evaluates all 16 consecutive x intervals
with literal FMA fractions and records plane/order counts:

| x origin | represented planes | ordered intervals | collapsed intervals | pre-fix result |
|---|---:|---:|---:|---|
| `nextafter(2^48,0)` | 10/17 | 9/16 | 7 | late `overflow_error` containing `collapsed coordinates` |
| `-2^48-0.5` | 17/17 | 16/16 | 0 | constructor completed |
| `2^48` | 17/17 | 16/16 | 0 | constructor completed |

The positive case's first interval is explicitly asserted positive. Its seven
collapses are the pairs `(2,3),(4,5),...,(14,15)`. The late exception arises
during actual tree traversal after `validateConfig`, containment, exact material
region audits, subset relation, and surface moment integration; it is not an
early malformed-geometry rejection. The negative and aligned cases document
the conservative equal-ULP boundary: they are valid alignments, not additional
collapse defects.

Actual RED command, from repository root:

```bash
set -o pipefail
ulimit -v 25165824
timeout 3600 cmake --build build/strict-v2-spin --target spin-it-2014-source-surface-octree-regression -j3 2>&1 | tee build/task-5i-evidence-20260910-0FF6Vb/red-build-01.log
timeout 3600 ctest --test-dir build/strict-v2-spin -R '^spin-it-2014-source-surface-octree$' -V 2>&1 | tee build/task-5i-evidence-20260910-0FF6Vb/red-ctest-01.log
```

Build exit 0. CTest exit 8, 0/1 passed. The output contains the literal three
rows above and fails only at `equal-ULP source roots reject early as unsupported
numerical input`. It therefore distinguishes the actual late collapse from the
two accepted conservative-boundary inputs.

## Minimal correction

The admission now requires a finite, positive finest step **strictly greater**
than the maximum root-endpoint ULP on every axis (`step <= spacing` rejects).
This is the same sufficient fail-closed bound used by `SourceCageGrid`. The
diagnostic is `source octree finest coordinate spacing is numerically
unsupported` and remains `std::invalid_argument`; neither the code nor the test
uses an infeasibility label.

The public header documents that equal-ULP aligned inputs can be conservatively
rejected and that the implementation does not recenter, edit geometry, or
coarsen the request. There is no root-lattice scan, tolerance change, geometry
repair, output identity-version change, or solver behavior change.

The regression additionally pins the pre-fix ordinary default identities,
derived from the actual RED run before production changed:

```text
source    d05f34b01ef35d581ec903f9e803deb12402faddf93a0b6f209a464d0bf4297a
config    fc5b062c289ae11c288f542f9bf193a0a43f656f81d42f9bc8ca1143d46b43ef
partition 3c4531f7bab657f722244a74b82dad6bc08c63a23ad6747e0d0e54df8538565f
```

## GREEN and integration evidence

Focused GREEN repeated the RED commands after the single correction. Build
exit 0; CTest exit 0, 1/1 passed in 1.52 seconds. All three cases report
`early_numerical=1 late_collapse=0 completed=0`, while the ordinary identities
remain exactly the literals above.

Relevant main-configuration coverage used:

```bash
timeout 3600 cmake --build build/strict-v2-spin --target \
  spin-it-2014-source-cage-grid-regression \
  spin-it-2014-source-surface-cage-regression \
  spin-it-2014-source-material-regression \
  spin-it-2014-source-adaptive-regression \
  spin-it-2014-source-surface-octree-regression -j3
timeout 3600 cmake --build build/strict-v2-spin --target _make_it_stand -j3
timeout 3600 ctest --test-dir build/strict-v2-spin \
  -R '^(spin-it-2014-source-(cage-grid|surface-cage|material|adaptive(-geometry9)?|surface-octree(-level9)?)|source-spin-material-python)$' -V
```

The C++ targets and Python extension built. The covering CTest exited 0 with
8/8 tests passed in 133.11 seconds: grid, surface cage, Python material bridge,
native material, adaptive, adaptive geometry9, source octree, and source octree
level9. The depth-9 tree remains 2,097,993 nodes with partition hash
`652de7d6b99c64445640902e891a98f85f76c4a8b5a93b6ef34471954a5bc137`.

One preliminary combined covering-build command appended a nonexistent target
name `make_it_stand_native` after the five real C++ targets. Those five targets
had already built successfully, then Make exited 2 for the unknown final target.
This was a command-name error, not a behavioral RED. The correct `_make_it_stand`
target was immediately built separately with exit 0, and both the covering and
full tests below use that extension. The failed command is retained, not hidden,
in `main-covering-build-01.log`.

The existing FasterOFF directory was discovered as
`build/task-5c-no-faster-20260909`. Its actual cache records exact geometry ON,
Python ON, Spin-It 2014 ON, Faster OFF, viewer OFF, and both FetchContent
disconnected flags ON. Without configure, download, or dependency update, the
same five C++ targets plus `_make_it_stand` built with exit 0 and the same eight
CTest entries passed 8/8 in 120.16 seconds. Output explicitly records `faster=0`
for source/adaptive tests and preserves both ordinary and depth-9 hashes.

All-target final verification, after production/test changes were frozen:

```bash
set -o pipefail
ulimit -v 25165824
timeout 3600 cmake --build build/strict-v2-spin -j3 2>&1 | tee build/task-5i-evidence-20260910-0FF6Vb/full-build-01.log
timeout 3600 ctest --test-dir build/strict-v2-spin --output-on-failure -j1 2>&1 | tee build/task-5i-evidence-20260910-0FF6Vb/full-ctest-01.log
```

Both commands exited 0. The all-target build completed. Full CTest passed
**201/201**, zero failures, in **474.14 seconds**. The outer tee contains exactly
201 `Test ... Passed` records and the final 100% summary, independently of the
nested registration test's writes to the default CTest temporary log. Commands
were sequential; builds used at most `-j3`; the full test was serial. Every
build/test command used the 24-GiB virtual-memory cap and 3600-second timeout.
No downloads, dependency updates, file deletions, or main-checkout writes were
performed. Recompiling the touched CGAL translation unit emitted only the
existing Boost global-placeholder pragma note.

### Evidence hashes

All owned logs are retained in the fresh unique directory
`build/task-5i-evidence-20260910-0FF6Vb/`:

| Log | Exit/result | SHA-256 |
|---|---|---|
| `red-build-01.log` | 0 | `5c433bc8bd7a5da2be7481dd339f0aff0e1d7923e45ab8e71a5e54dcf3c180a8` |
| `red-ctest-01.log` | 8; behavioral RED, 0/1 | `0d40d486efad22494c50dd73e4e99fcdfbb3ace75d7bfe8b10c63716bba9e773` |
| `green-build-01.log` | 0 | `cc69a371ed812f15923848822a0d5e14b2a6191bb1da626a1b06a53437f7889e` |
| `green-ctest-01.log` | 0; 1/1 | `4399dce22430d1671f775122f3fb4d3653514f615a7eca252cf84a72bfb12370` |
| `main-covering-build-01.log` | 2; five real targets built, bad trailing target | `552eace13b3d1a5954834f8801b786712f1c7ff122384de975ba96f30ea8f4e2` |
| `main-covering-python-build-01.log` | 0; correction | `6c964b75ef64517322539320708bd09d18eab06fe5ed48a941b2facadfa89ce4` |
| `main-covering-ctest-01.log` | 0; 8/8 | `ff8d1280e9432e2d56f36ffa913baf38b18e558fef02a66e53f5b56230b02cc8` |
| `fasteroff-cache-flags.log` | 0 | `b4e93b823c4b9f341cb120a40e4d0ed5d2578908bc30a503e5cb9949350a0dfe` |
| `fasteroff-covering-build-01.log` | 0 | `03e96fdb6b9a0531dd58b7206e046db435dc7b5230adac317a764e725fdb743d` |
| `fasteroff-covering-ctest-01.log` | 0; 8/8 | `a6e30e4bf11ae39f7fe61eb15691b8bd0491e1e95a355b2bfbe7fe0c9b8325cf` |
| `full-build-01.log` | 0; all targets | `01ac9ff30de0812fc263dd420e906e4d5b148607b4eb55886c20bbc815415259` |
| `full-ctest-01.log` | 0; 201/201 | `b97b2c375046d4e028427c8195183d182f76243f84599c827b2b8b9311828b79` |

As supplementary controller-owned read-only verification, the existing
standalone Task5h grid probe was compiled against the frozen final main archive
without a project rebuild. Two actual exports are byte-identical to each other
and to Task5h ordinary output: 1,728 cells, 2,197 nodes, 6,084 edges, 218
handles, maximum field error `1.3877787807814457e-17`, JSON SHA-256
`b4193f3e510cc6f142765ec90e5e9c3f164d86e97e03cc8e9c65add3e7b1dd6e`.
The checked log is `build/task-5i-controller-W6WCqr/grid-check.log`, SHA-256
`c5543dc9282b329b32a6640b6b32403a89bb186381a3d40a48db317bf69d8242`.
This is `SYNTHETIC_GRID_INTERPOLATION_NOT_BBW`, not paper or BBW authority.

## Files, hashes, and self-review

Scoped implementation files:

| Artifact | SHA-256 |
|---|---|
| `include/spin2014/SourceSurfaceOctree.h` | `57618e3c3dc3a092dfd85c0a91d2c606419bb7827881644b0014415085e6a5cf` |
| `src/spin2014/SourceSurfaceOctree.cpp` | `96b842b65a3697e9c9455c385e57b403e2f503f6591b51272f20ff7f0df63b14` |
| `tests/spin_it_2014_source_surface_octree_regression.cpp` | `562e584f850e57354d824d54a1aed891cd8f3fd572c3c6f6eb938c166e28226a` |
| main `libmake_it_stand_core.a` | `c48ea709bc7c1aac919ec59d063f8d980acfcd25e9baa6509cd92d313c1c4580` |
| FasterOFF `libmake_it_stand_core.a` | `05b3b65e05d1c6a732bcd45ef8f1d9b042986df43a4417234405702eca7a70dd` |

Self-review re-read the complete scoped diff and checked the mutation that
restores `<`: the native test would reproduce either late collapse or completed
construction and fail. The expectations use literal plane fractions/counts and
pre-fix output hashes rather than production helpers. Exact native construction,
not a mock, establishes closed geometry and consumer-visible exception type.
The strict comparison is axis-local, finite/positive guarded, and makes no
accepted-input arithmetic or identity change. `git diff --check` passed.

The dirty main checkout's tracked-diff SHA-256 remains exactly
`93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`,
matching preflight. Concurrent controller edits to
`docs/spin2014_bbw_contract_audit.md` and
`docs/spin_it_2014_source_surface_octree.md` are explicitly outside this commit,
as are controller plane scripts/notes, probes, logs, and material replays.

## Scope limits and residual risks

- Equality is deliberately fail-closed. Negative and aligned equal-ULP roots
  with all planes distinct are rejected; this is documented numerical scope,
  not a claim that they collapse and not physical infeasibility.
- The largest-endpoint-ULP test is a sufficient bound, not an exhaustive scan
  of all implicit planes. It is intentionally conservative and resource-bounded.
- Existing accepted configurations keep their geometry, moments, tree traversal,
  hashes, tolerances, and solver behavior. Only the unsupported admission domain
  narrows.
- No failed construction or solve is treated as an infeasibility certificate.
  No original asset/paper case was reconstructed or run.
- No BBW implementation, FEM discretization decision, deformation integration,
  optimizer change, sealing claim, or release claim is included.
- The previously parked Task5c minor observation remains deferred: the shifted
  non-dyadic-root test comment mentions query partition arithmetic while that
  section directly asserts enclosure/config ownership. Expanding that unrelated
  coverage is outside this numerical-bound correction.

Only the source-octree public header, implementation, native regression, and
this report are intended for the Task5i commit.
