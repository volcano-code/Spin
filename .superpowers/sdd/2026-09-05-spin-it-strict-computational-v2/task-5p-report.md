# Task5p implementer report — shared source material boundary

Approved architectural slice of Task5/M3; base
`4c6462d1ec3d3510f2963beff64a5b8f496e0d58`. Global strict reproduction remains
**PARTIAL_BLOCKED**. This implementation closes a necessary shared topology
condition; it does not deliver a joint solver or legal deformed artifact.

## Ownership and implementation

Owned files: `include/spin2014/SourceConformingBoundary.h`,
`src/spin2014/SourceConformingBoundary.cpp`,
`tests/spin_it_2014_source_conforming_boundary_regression.cpp`,
`docs/spin2014_source_conforming_boundary.md`, this report, and the two scoped
CMake hunks adding the source/header and native test. No helper source needed.
Public API was proposed to and accepted by the controller before implementation.
Patch plane is explicitly source-D, while material keys are D+1.

The immutable private-shared-const-Data object retains the actual binding and
all 496 baseline Full primitives/eight adaptive beta owners. Sparse per-plane
quadtrees inherit side ownership and emit atomic square patches. Shared global
axis-line point dictionaries split edges at perpendicular-plane hanging points.
Each patch has one shared center fan; both owners reference the identical
triangle indices in opposite orientations. Primitive boundaries are verified
closed and consistently oriented before object acceptance.

Original outer rows/topology stay intact; internal binding corners are
explicitly reindexed into sorted D+1 material keys. Every original material
corner must survive. Canonical FMA coordinates, strict endpoint-ULP spacing,
finite positive triangle area, checked sparse growth and derived resource
counts, and versioned identities are implemented. Resource failures use
`RESOURCE_LIMIT` length errors. Empty carveable bindings emit no material.

Only assigned Release build `build/task-5c-no-faster-20260909` was used:
2014 ON, exact geometry ON, Faster OFF. Build commands used `-j3` at most and
`ulimit -v 25165824`; suite ran `-j1`. No subagents were spawned. The six
protected Faster files and concurrently appearing PlaneMoment files/directory
were not edited. No main-project build, cache deletion, merge, push, or tag.

## TDD evidence

The regression and test registration were written before the public header or
implementation. This exact command produced missing-interface RED:

```sh
ulimit -v 25165824
cmake --build build/task-5c-no-faster-20260909 \
  --target spin-it-2014-source-conforming-boundary-regression -j3 \
  > build/task-5c-no-faster-20260909/task-5p-red.log 2>&1
```

Exit 1, expected diagnostic:
`fatal error: spin2014/SourceConformingBoundary.h: No such file or directory`.
This is the exact planned interface, not a misspelled include. No production
SourceConformingBoundary file existed when RED was observed.

After implementation the same build command (log `task-5p-build.log`) exited
0. Focused GREEN command:

```sh
ulimit -v 25165824
ctest --test-dir build/task-5c-no-faster-20260909 \
  -R '^spin-it-2014-source-conforming-boundary$' --output-on-failure \
  --output-log build/task-5c-no-faster-20260909/task-5p-focused.log
```

Initial result 1/1 PASS, 0.44 s. Additional edge tests required no production
change; their target build exited 0 (`task-5p-edge-build.log`) and final focused
run passed 1/1 in 0.61 s (`task-5p-edge-focused.log`). Existing CMake/libigl
configuration notices are retained in the initial build log; compiler and
focused test output contains no new warning or failure.

## Actual behavioral checks

- Baseline outer cube (.3,.7), design cube (.34,.66), source D5, adaptive
  initial 2/max 5, actual 26-handle BBW coarse 2/fine 3: eight leaves and
  496 Full primitives retained. Output: 1,968 patches, 2,965 points including
  original outer rows, and 7,944 triangles.
- Included unit cells independently enumerate primitive face coverage and side
  orientation. Every primitive edge occurs twice in opposite directions.
  Every rest triangle points outward. Signed tetrahedral integration matches
  all ten analytic box moments at normalized tolerance 1e-10. Total mass is
  literal `(10/32)^3`, and nonuniform adaptive beta totals match binding moments.
- Actual coarse/fine patches and inserted hanging edge points are required.
  Shared IDs/opposite references cancel under the injective map
  `(x,y,z) -> (x,y,z+xy/16)`. A separate analytic naive coarse/fine center
  counterexample differs by exactly 1/64; it is not a source model result.
- Actual final BBW weights are interpolated in chunks of 11 rows using all
  26 handles, with nonuniform scales and translations. Fixed topology remains
  closed and cancels shared triangles. Reintegrated primitive moments agree
  with independently summed external/shared surfaces after motion.
- Repeated constructions have identical coordinate bytes, keys, patch fields,
  triangles, references and hashes. Each of four budget changes alters both
  identities. Caller split/merge followed by new bindings preserves the union
  and correctly updates owners while the original snapshot stays immutable.
- Translated/scaled nonsymmetric design box and a small true source-D9 fixture
  pass the same independent coverage/moment checks. The latter has lattice
  level 10 and fewer than 10,000 points; no global 512^3 array or level-9
  optimization is used.
- All budgets reject zero, size_t maximum and cap 1 with RESOURCE_LIMIT.
  Exact point/patch/triangle caps succeed, and one-under caps fail. An empty
  actual source binding emits no material and preserves the outer rows.
  An actual upstream-accepted binding at origin `(2^40,0,0)`, side `2^-6`,
  source D5 has D-step 2 ULP but D+1-step 1 ULP; boundary construction rejects
  it with the explicit D+1 spacing diagnostic. No tolerance was relaxed.

Baseline configuration identity:
`f9c191745d8158702fa8bb12113a4b36bd5ba477e1bb27fcdccd461ccdae2fab`.
Baseline output identity:
`00e3bcdbfca9b992f963a221fb3cee7e8ead67fe560cab0112cfce241e21871e`.

## Safe suite and dependency graph

The safe suite is invoked exactly once after final focused GREEN:

```sh
ulimit -v 25165824
ctest --test-dir build/task-5c-no-faster-20260909 \
  -R '^(spin-it-2014-|spin-2014-|source-spin-material-python$)' \
  --output-on-failure -j1 \
  --output-log build/task-5c-no-faster-20260909/task-5p-safe-scoped-suite.log
```

Result: **34/34 PASS**, exit 0, total 226.63 s. Dependency contract passed in
89.96 s and isolation in 7.30 s. The retained actual dependency build is
`build/task-5c-no-faster-20260909/spin-2014-dependency-contract-runs/run-d38760ca539db268`.
Its enabled CMake cache has 2014 ON, exact geometry OFF, Faster OFF, and its
generated `enabled/build.ninja` has no SourceConformingBoundary source or target.
This checks actual exact-OFF graph exclusion without another build.
Retained enabled `build.ninja` SHA-256:
`3169055488c94400f30661100789cf66d926ffee9c3452c1db14dd57a4b77637`.

## Independent controller evidence

Controller independently compiled its exporter against the unchanged archive
below. It reports three actual pipelines, each repeated twice, passing
unit-face ownership/closure, analytic and signed-tetrahedral moments,
nonuniform beta weighting, and manually reconstructed full BBW motion.
Maximum normalized moment error 4.163e-12; motion discrepancy 4.441e-16.
Exports were byte-identical with SHA-256
`0e91494a9dc711d77bf2b160072181772a43b688ca611a84dab71b2e99214c44`.
Controller-owned evidence: `build/task-5p-controller-V23Tpb`.
Controller also reports 18 deliberate corruption checks passed (six for each
fixture): missing references, reversed triangles, wrong patch owners, moved
moments, moved point, and duplicated key were detected.
This paragraph attributes controller evidence; it is not implementer rerun
evidence and does not replace the scheduled fresh review.

## Artifact and log SHA-256

Paths below are relative to this worktree. Log and binary paths use prefix
`build/task-5c-no-faster-20260909/`.

| Artifact | SHA-256 |
|---|---|
| `CMakeLists.txt` (only six owned added lines) | `a0015aa9ceb76e476cadfe0e99729287778a9e74c648169d14044b14fce5cae4` |
| `include/spin2014/SourceConformingBoundary.h` | `c4d99b94cd12614f829793e648679a0f6312110815d1b5c7e75a974080349655` |
| `src/spin2014/SourceConformingBoundary.cpp` | `808efe1e5b9f6eaa0b724256cacdfc30c47a67c57a6743fbc9fd6a27b543fc02` |
| `tests/spin_it_2014_source_conforming_boundary_regression.cpp` | `38507f08e44d816d74545cb5491e5e3c76023c3cda4c9499bb05140420bb657a` |
| `docs/spin2014_source_conforming_boundary.md` | `b2ebdf7eb5c23a30b771108d9167c67850048d2cfc866fe18485c5c6979eece7` |
| `libmake_it_stand_core.a` | `ba312a72f57108e4fe9ae9f39bc9627ddcff2fab1e8368c053d8ecbd378d4361` |
| `spin-it-2014-source-conforming-boundary-regression` | `65ad18825aa63ea2a0bd5e7330675e1ce696bd4f374bedbced397de09fbd1f9c` |
| `task-5p-red.log` | `849d35e5049af4022af752a2f2e78aa97ea66475e0a95ea2831525d2a792d0d9` |
| `task-5p-build.log` | `2e8b391461e33f7b583228d9714bbc742bc971be056b987404466fca8546ff0c` |
| `task-5p-focused.log` | `a4536f879563e24f4ff780891068399e712069986f4f162559a4e32a7b63b653` |
| `task-5p-edge-build.log` | `fbfb4ccf5cf6d70d58581a78502f57c4632c10c02e94bca8ab277bdb7c383bfb` |
| `task-5p-edge-focused.log` | `55e3862b3875151c091723a0401a0c773aa043ef70b63e9a5058f49900e3a2e9` |
| `task-5p-safe-scoped-suite.log` | `1a1a9b561753ab44301a0766923c683f5cb47c839965e8ee4a6bcc6b90a80c30` |

## Self-review and limitations

Reviewed ownership signs, cyclic axes, all-plane line queries, missing-quadrant
inheritance, canonical ordering, outer preservation, primitive ID retention,
checked insertion/count arithmetic, immutable identity fields, exception paths,
and the complete scoped diff. `git diff --check` passed. Added explicit
outward, repeated-reference, exact-cap, empty-source and D+1 frame tests during
self-review; all passed without production fixes. No unresolved implementation
concern identified. The geometry unit fits one 238-line implementation file.

Shared indexed closure is not global material manifoldness at edge-only
contacts, self-intersection prevention, shell-thickness enforcement, a material
seal, or a joint solver. Full source beta/phi/all-handle SLQP, paper
regularizers, analytic handle derivatives, geometric precommit rejection,
final exact material partition/seal and original-paper assets remain required.
All motion checks here are SYNTHETIC_REGRESSION. Resource caps bound counts,
not exact combined RSS. No scientific completion flags were changed.
