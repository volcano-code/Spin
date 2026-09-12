# Task5p independent verification — 2026-09-12

Scope: synthetic source-material shared boundary, not original-paper experiments,
joint SLQP, deformed geometry safety, physical validation or complete reproduction.
BASE: `4c6462d1ec3d3510f2963beff64a5b8f496e0d58`.

## Executed independent checks

The standalone exporter uses the actual source surface octree, adaptive domain,
complete cage, FEM and coupled BBW solve, immutable binding and conforming boundary.
No reconstructed project Python classes or replacement surface inputs are used.
The Python checker imports only NumPy/stdlib; its integer unit-face coverage,
directed-edge closure, signed-tetrahedron moments and manual eight-node BBW
interpolation do not call project checker/integration functions.

| Fixture | Source depth | Primitives | Points | Patches | Triangles | Coarse/fine patches | Patches with hanging-edge subdivisions |
|---|---:|---:|---:|---:|---:|---:|---:|
| Baseline | 5 | 496 | 2965 | 1968 | 7944 | 384 | 12 |
| Nonsymmetric, scale 1.7 and translation | 5 | 147 | 912 | 588 | 2412 | 120 | 48 |
| Small genuine depth-9 source | 9 | 57 | 360 | 228 | 912 | 24 | 0 |

The last fixture tests depth-9 geometry, **not** a level-9 optimization claim.
All 26 actual handles are used for the first two fixtures; the small third cage
has all its 8 handles, not a reduced subset. Motion applies distinct scales and
translations to every handle. Independent interpolation checks every exported
point. Rest and moved primitive moments, exterior-union moments and nonuniform
beta-weighted moments are all checked. Rest moments also use analytic boxes.

Maximum observed relative ten-moment error: `4.163e-12` (threshold `1e-10`).
Maximum observed full-BBW motion discrepancy: `4.441e-16` (threshold `1e-12`).
Every fixture was run twice in the same process configuration; complete text
exports compare byte-for-byte. Eighteen deliberate record corruptions are
rejected (six categories per fixture): missing face reference, reversed triangle,
invalid ownership, incorrect moments, incorrect moved coordinates, duplicated key.

## Reproducible commands and retained files

Working directory: `/data/xvxinyao/make-it-balance-paper-sim`.
Evidence directory: `build/task-5p-controller-V23Tpb`.
Helper prefix: `.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5p-controller-`.

Compile used `/usr/bin/c++ -std=c++17 -O2 -g -DNDEBUG -frounding-math`,
2014/exact feature definitions ON, includes `include`, `tests`,
`libigl/external/eigen`; linked assigned Faster-OFF core archive, static Boost
system/thread, static HiGHS and system GMP/MPFR/zlib. Address-space limit was
`ulimit -v 25165824` (not a measured RSS assertion).

For each run:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 build/task-5p-controller-V23Tpb/export > build/task-5p-controller-V23Tpb/run1.txt
/home/lin/anaconda3/bin/python3 .superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5p-controller-check.py build/task-5p-controller-V23Tpb/run1.txt
```

Repeat with `run2.txt`, then `cmp` the two files. Run the negative checker with
the first file as its sole argument. All commands exited zero, including both
independent checks, comparison and the final negative-checker revision.

## SHA-256 bindings

| Input / output | SHA-256 |
|---|---|
| SourceConformingBoundary.h | c4d99b94cd12614f829793e648679a0f6312110815d1b5c7e75a974080349655 |
| SourceConformingBoundary.cpp | 808efe1e5b9f6eaa0b724256cacdfc30c47a67c57a6743fbc9fd6a27b543fc02 |
| Core archive | ba312a72f57108e4fe9ae9f39bc9627ddcff2fab1e8368c053d8ecbd378d4361 |
| Exporter source | 13331ee979a3464f6b71101f12875e38746e6a12ebbf3582add2835827ecfce3 |
| Independent checker | f55fa0fc9046ba604554df8b6a9f79e50cf94316a99f75973d446f58f1ebfb9b |
| Negative checker | fae760bcb0f4cbbb6642931695f8f86ba62eda1e63e3bb572dc4383a98cf869f |
| Linked exporter | 0bc879bda99551cc4ed11a61a9012997fbc3ac7842d29d8c149e0bc0bc742285 |
| run1.txt and run2.txt | 0e91494a9dc711d77bf2b160072181772a43b688ca611a84dab71b2e99214c44 |

Main worktree tracked-diff SHA remained
`93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`.
This is not an exhaustive untracked preservation assertion. Existing Faster and
new external plane-family paths remain excluded from this work and its claims.

Implementation committed as `ca9f976ba8a27349f8a24eecdf846d67ba382d87` using the
authorized Git identity. Controller read the complete actual scoped suite log:
34/34 tests passed in 226.63 seconds, including the new regression and actual
dependency compilation/isolation checks. This is not the whole-project suite.
All 13 hashes in the implementer report were independently rechecked and match.
Actual retained exact-OFF build graph contains no new boundary source/target;
its SHA-256 is `3169055488c94400f30661100789cf66d926ffee9c3452c1db14dd57a4b77637`.
No rebuild or deletion was needed to verify the existing run.

Fresh read-only review approved spec compliance and task quality, without
Critical/Important findings; `task-5p-review.md` records the two deferred minor
maintenance items and scope. Nothing grants final release: `PARTIAL_BLOCKED`.
