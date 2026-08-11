# R0 baseline and acceptance record

Audit date: 2026-08-11 (Asia/Shanghai)

## Frozen starting point

- Branch: `reproduction-v1`
- Upstream baseline: `71e08361ca255a1cf6b89de6e4c92ef33024f432`
- Upstream remote at audit: `https://github.com/helloxinw/make-it-balance.git`
- Pinned libigl gitlink: `f06fb13014aadff7eb1d7423ac018620441a3657`
- All Phase 1.1–1.6 changes were initially uncommitted.

## Pre-R0 diagnostic observation

Before the reproducible target split, the legacy libigl headers required a
temporary compiler-wide `-include cstdint` workaround to build the regression
target. With that workaround, the one-step and full-run smoke paths contained
no NaN/Inf and produced CSV output. The full current-heuristic run nevertheless
showed the following non-conformances:

| Metric | Initial | Final/observed |
|---|---:|---:|
| Legacy projected error | 0.032699 | 0.252500 |
| Geometric projected error | 0.028996 | 0.811397 |
| Incremental CoM vs voxel-proxy CoM discrepancy | — | 0.800743 |
| Carved/carveable ratio | 0 | 80.293% |
| Carved/entire-grid ratio | 0 | 14.766% |
| Candidate exhaustion | — | step 17 |
| Actual stop | — | step 21, maximum iterations |

These values are a known-failure baseline. They are not target numbers and are
not evidence of successful balancing. They justify the R2 mass/CoM work and R4
paper-faithful best-prefix carving.

## R0 acceptance gate

R0 is accepted only when all of the following are recorded after the new files
are committed:

- [x] libigl submodule is exactly at the pinned commit and clean;
- [x] headless preset configures from an empty build directory;
- [x] `make_it_stand_core` and `carving-regression` build;
- [x] headless CTest passes without global forced includes;
- [x] viewer preset configures from a separate empty build directory;
- [x] `make_it_stand_core` and `make-it-stand` viewer build;
- [x] README, manifest, provenance and conformance status agree;
- [x] `git diff --check` is clean;
- [x] committed checkout leaves `git status --porcelain` clean after builds.

GUI interaction is not automated in R0. A successful viewer compile is therefore
recorded separately from an eventual interactive behavior regression.

## R0 acceptance evidence

The following commands were executed with the locked system toolchain:

```text
/usr/bin/cmake --preset headless                         configured
/usr/bin/cmake --build --preset headless --parallel     13/13 built
/usr/bin/ctest --preset headless                         1/1 passed, 0.13 s
/usr/bin/cmake --preset viewer                           configured
/usr/bin/cmake --build --preset viewer --parallel       41/41 built
```

The headless cache records all OpenGL/GLFW/ImGui options as `OFF`; the viewer
cache records them as `ON`. Both compile-command databases use `-std=c++17`.
The libigl working tree is clean at the pinned gitlink. The full 21-step
diagnostic was rerun without a forced-include workaround and reproduced the
known-failure values above, including candidate exhaustion at step 17 and the
current maximum-iteration stop at step 21.
