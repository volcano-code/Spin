# R1 diagnostic-foundation acceptance

Acceptance date: 2026-08-12 (Asia/Shanghai)

## Scope

R1 closes the diagnostic/control-flow foundation inherited from Phase 1.6. It
adds explicit coverage for zero steps, one ordinary step, a complete Bunny run
and an initially empty candidate set. It also makes candidate exhaustion an
actual production solver stop reason.

R1 deliberately does not repair:

- the repeated-square legacy projected error;
- inconsistent mass/first-moment/CoM updates;
- fixed batches of at most 1000 removed voxels;
- the whole-grid empty-ratio denominator;
- paper best-prefix selection or temporary-energy-increase handling;
- support polygons, safe regions, multiple targets or deformation.

Those unresolved numerical and algorithmic issues remain assigned to R2–R4.

## Test matrix

| CTest | Expected behavior | Result |
|---|---|---|
| `carving-initial-state` | Executes zero steps; iteration, carved count and history remain zero; CSV contains only its header. | Passed |
| `carving-one-step` | Executes one ordinary positive-candidate step and verifies all state/count/finite invariants. | Passed |
| `carving-full-known-baseline` | Runs production `carveAll()` and stops at the first empty candidate search, iteration 17. | Passed |
| `carving-no-candidates` | Removes all test-fixture carveability before construction; one search removes nothing and stops with `no_candidates`. | Passed |

Every populated CSV row is checked for exactly 26 columns. All 25 serialized
numeric fields are parsed and mapped back to their corresponding in-memory
diagnostic values, while the stop token must match exactly. CSV row counts must
equal history size, with a single additional header row.
The two terminal scenarios also call `carveAll()` a second time and require
iteration, history, carved count and CoM to remain unchanged.

## Acceptance evidence

```text
/usr/bin/cmake --preset headless                         configured
/usr/bin/cmake --build --preset headless --parallel     built
/usr/bin/ctest --preset headless --output-on-failure    4/4 passed
/usr/bin/ctest --preset headless --repeat until-fail:10 4/4 passed ten times
Debug build with assertions enabled                            4/4 passed
Fresh viewer build                                            41/41 built
```

Two independently generated full-run CSV files had the same SHA-256:

```text
3cc4d09ff0dc310d636751a71351418989be6406e64b5fbf7cb0e3a53f810f44
```

The full live R1 trace has 17 diagnostic rows plus one header. The zero-step,
one-step and no-candidate files contain 1, 2 and 2 total lines respectively.

## Frozen known-failure result

The R1 full run still reaches the same incorrect numerical state as R0:

| Characterization metric | R1 value |
|---|---:|
| Initial geometric projected error | 0.0289961646 |
| Final geometric projected error | 0.8113967646 |
| Final incremental/voxel-proxy CoM discrepancy | 0.8007426348 |
| Total carved voxels | 14028 |
| First empty candidate set | iteration 17 |
| R1 stop reason | `no_candidates` |

Thus R1 fixes empty-loop termination but does not improve balance. The exact
pre-R1 21-step trace remains in
`tests/baselines/bunny_r0_known_failure.csv` with SHA-256:

```text
3f9108649fc329a136a5371597335a99275021b6ddd9a3984198ce2b8e8c1856
```

R2 is expected to invalidate the frozen floating-point characterization when
it replaces the mass and center-of-mass implementation. At that point the
characterization assertion must be deliberately replaced by analytic and
incremental/full-recompute correctness tests, not relaxed to retain bad output.
