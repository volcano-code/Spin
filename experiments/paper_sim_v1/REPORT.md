# Make-It-Stand paper-simulation implementation status

Last audited: 2026-09-01 (Asia/Shanghai)

## Current claim

The closeout infrastructure is implemented and all registered pre-reachability
regression tests passed.  Five current-source full-resolution `RELEASED_BETA_PROFILE`
diagnostic runs now exist, but the branch is **not** yet
`VERIFIED_PAPER_SIM`: two cases miss the 5-degree safe region, the paper-profile
lambda/Table 1 matrix and Figure 14 are incomplete, and the current Spheres
collision proxies disagree with the static oracle.  Spin-It Faster therefore
is not evidence for a `VERIFIED_PAPER_SIM` release.

The later user-approved Spin-It track now has a separate carving-reachability
gate.  It does not upgrade the paper-simulation claim or waive any blocker in
this report, but it may unblock Spin-It Faster independently once its own build,
test, provenance, and determinism checks pass.

## Alternative-pose reachability result

The official-cache Spheres 45-degree pose, official-cache Dinosaur head-contact
pose, and reconstructed-asset Toy head-down pose were run twice against the
fixed-region continuous carving relaxation.  All three are
`PROVEN_UNREACHABLE`, with byte-identical repeat artifacts and valid Farkas
certificates.  Their stationarity infinity norms are respectively
`7.782e-10`, `6.106e-10`, and `0`, within the `1e-9` validation threshold.

This resolves the prior plane-carver question without claiming general
physical impossibility.  In particular, Dinosaur's earlier successful run is
consistent with this result because that workflow allowed deformation and a
changed contact geometry.  The auditable report is in
`experiments/spin_it_reproduction/reachability/results/REPORT.md`.

## Evidence reset

Every pre-existing experiment manifest in this tree was audited and marked
`STALE_RESULT` because it had no execution fingerprint.  In particular, the
old Armadillo diagnostic was produced before the released-beta one-carving-pass
semantic correction and is not evidence about the current optimizer.

Historical observations such as beta spheres/dinosaur convergence and the
5/18 spheres drop result remain useful debugging records only.  They cannot be
used by the release gate or reported as current reproduction results.

## Implemented gates

- Exact source/executable/configuration/runtime fingerprints and automatic
  `STALE_RESULT` classification.
- Separate `PAPER_PROFILE` and `RELEASED_BETA_PROFILE` with a schema-v2 trace
  that does not conflate paper iterations, beta calls, carving, deformation,
  or committed events.
- Strict pre-commit deformation checks for self-intersection, flipped or
  degenerate faces, with line-search reduction or rollback.
- Post-flatten `SealedBalanceArtifact` with source/final outer and cavity
  topology audits, deterministic structured-cavity refill, collapsed support
  patch remeshing, exact shell clearance/containment checks, and recomputed
  mass, center of mass, and center-frame inertia.
- Fixed-support quasi-static tilt oracle and separately named reseeded-pose
  PyBullet protocol with contact points, normals, distances, and forces.
- Full-resolution lambda search, five-case Table 1 orchestration, official
  N=80 Figure 14 orchestration, and a fail-closed final release gate.
- Explicit horse and rocking-toy reconstruction products.  Rocking toy uses a
  deterministic farthest-surface handle policy and can never be upgraded to
  `EXACT_CACHE` without the missing author files.

## Fresh test evidence

After reconfiguration and rebuild on 2026-09-01:

| Preset | Result |
|---|---:|
| `headless` | 116/116 passed |
| `paper-python` | 117/117 passed |
| `python` | 122/122 passed |
| `physics-python` | 117/117 passed |

No test was `Not Run` or skipped.  The Python extension suite includes the
official read-only archive gate and took 40.83 seconds in this run.  Presets
were run sequentially because their regression executables intentionally use
some fixed diagnostic paths under the system temporary directory; launching
multiple presets concurrently can create a cross-build-tree file race.

## Current full-resolution beta diagnostics

| Case | Evidence | Result | Events / first safe | Sealed state |
|---|---|---|---:|---|
| Spheres N=80 | `EXACT_CACHE` | stable | 80 / 79 | `SAFE`, residual 0 |
| Dinosaur N=100 | `EXACT_CACHE` | stable | 100 / 99 | `SAFE`, residual `6.43e-7` |
| Armadillo1 N=80 | `EXACT_CACHE` | `mu_exhausted_unstable` | 145 / none | `OUTSIDE`, residual `1.051e-3` |
| Horse N=100 | `RECONSTRUCTED_ASSET` | `mu_exhausted_unstable` | 97 / none | `OUTSIDE`, residual `2.158e-2` |
| Toyrocking N=50 | `RECONSTRUCTED_ASSET` | stable | 2 / 1 | `SAFE`, residual 0 |

All five final artifacts pass source/final outer, cavity, wall-clearance and
containment audits.  Spheres refilled 2 cavity cells, Dinosaur 8, Armadillo1
6, Horse 12 and Toy 6.  Dinosaur's archived export projection collapsed 2,370
support triangles; deterministic planar remeshing replaced them by 474
triangles and reduced 6,564 self-intersection pairs to zero.  The current
fixed-support oracle accepts Spheres, Dinosaur and Toy and rejects Armadillo1
and Horse.

The current Spheres PyBullet discriminator contains 18 unique scenes and 36
simulations.  Default V-HACD fails the 45-degree 2-degree reseeded drop;
robust V-HACD fails nominal and six perturbation directions.  Both repeats are
identical.  This is an active collision-model blocker, so the full 90-scene / 
180-simulation matrix was not run.

## Remaining release blockers

1. Run current-source paper-profile lambda reconstruction and two-repeat Table
   1 matrices; the present runs are beta diagnostics, not Table 1 evidence.
2. Resolve Armadillo1 and Horse so all cases have safe-region residual at most
   their normalized bounding-box diagonal times `1e-6`.
3. Recover the required Figure 14 qualitative ordering with two deterministic
   runs per ablation.
4. Pass fixed-support 0°, eight-direction 2°, and 5° tilt checks for all cases.
5. Resolve the demonstrated Spheres collision-proxy disagreement, then pass
   default and robust V-HACD for five cases, nominal plus eight 2°
   perturbations, twice: 90 unique case-scenes and 180 total simulations.
6. Do not relax angle thresholds or friction to hide collision-model errors.
7. Configure the user's Git name and email, commit the final evidence, obtain a
   clean tree, rerun the release gate, and only then create
   `make-it-stand-paper-sim-v1`.

Unavailable Toy caches, some paper parameters, and ancillary experiment assets
mean that the final permitted statement remains `VERIFIED_PAPER_SIM`, not a
perfect or bitwise reproduction of the authors' original experiments.
