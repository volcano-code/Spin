# Task5w implementation report

Status: bounded native transaction implemented and focused regression GREEN.
One requested negative remains explicitly unattained: actual passing input
geometry followed by rejected output geometry. No acceptance was mocked and no
clearance threshold was tuned to roundoff. Controller owns independent review,
graph evidence, broad verification and progress updates.

Base: `60b31e0f3d8568b045a0428b94c54e973c3228a2`.
Repository: `/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`.

## Implementation and exact API

```cpp
SourceCageAdaptiveRebindAudit rebindSourceCageAdaptive(
    const SourceCagePaperFunctional&,
    const Eigen::VectorXd&,
    const SourceCageAdaptiveRebindSettings& = {});

// SourceCageAdaptiveRebindSettings:
// fractionalEpsilon=.05 (IMPLEMENTATION_DEFAULT)
// geometry: SourceCageGeometrySettings
// budget.maximumOwnedBytes=256*1024*1024

const SourceCageAdaptiveRebindFacts& SourceCageAdaptiveRebindAudit::facts() const;
std::shared_ptr<const SourceCageAdaptiveCheckpoint>
    SourceCageAdaptiveRebindAudit::checkpoint() const;
const SourceCagePaperFunctional& SourceCageAdaptiveCheckpoint::functional() const;
const Eigen::VectorXd& SourceCageAdaptiveCheckpoint::state() const;
```

The immutable audit has shared const backing storage. A successful privately
constructed checkpoint owns a const functional and remapped state; every failed
transaction has a null checkpoint. No caller-provided replacement partition,
BBW, audit or passing flag exists. Copied facts cannot authorize a checkpoint.
Moving an audit moves ownership, not the referenced data; moved-from accessors
throw and previously obtained references require another live owner.

Eight stage slots distinguish preflight, input geometry, input full material
evaluation, local split, local merge, immutable rebuild, output geometry and
output full material evaluation. The original `SourceAdaptiveDomain` is copied,
so its original hidden leaf budget is retained. Split then merge are the existing
native methods; the legacy fill-only solver is untouched. The entire immutable
binding/boundary/moment/paper chain is reconstructed from retained source/full
BBW/configs/density/frame. Fresh full evaluations establish material hashes and
ten-moments and validate derivative caps before publication. Remapping copies
beta from the local transition and phi/each original handle parameter using old
and new functional offsets, including static no-angle layouts.

Facts retain ordered old/new leaves and beta, split/merge maps and REST geometric
moments, actual geometry audits, material world/spin moments and full-evaluation
hashes, before/after difference, and deterministic state/functional/partition/
binding/config/scientific digests. Geometry audit hashes bind independent
value-only evaluations. Timing and numerical history are absent. Nothing claims
convergence, KKT, feasibility, stable spin, a seal or paper acceptance.

Owned evidence/state/evaluation estimates use checked addition/multiplication,
up to eight new leaves/cuboids per old leaf/cuboid, source-depth-bounded merge
chains, capacity-growth allowances, retained checker records/strings and
evaluation dimensions. Upstream construction/checker caps remain unchanged.
The additional cap is conservative, separate from domain/rebuild/checker scratch,
and not a native RSS guarantee. Tiny rejecting caps still permit a small fixed
failure diagnostic. Diagnostic allocation failures may propagate. External
24 GiB / 60 minute guards remain required; no budgets were silently raised.

## TDD and verification evidence

Used `superpowers:test-driven-development` (including writing-good-tests),
`superpowers:systematic-debugging` for a test-assumption failure, and
`superpowers:verification-before-completion`. No subagents or reviewer dispatch.

All build/run commands used this shell prefix:

```sh
ulimit -v 25165824
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
```

Only `build/task-5c-no-faster-20260909` was used; its existing Release cache has
2014/exact geometry ON, Faster OFF and `LIBIGL_SKIP_DOWNLOAD=ON`.

RED, before the production header or implementation existed:

```sh
timeout 3600 cmake --build build/task-5c-no-faster-20260909 \
  --target spin-it-2014-source-cage-adaptive-rebind-regression -j3
timeout 3600 build/task-5c-no-faster-20260909/spin-it-2014-source-cage-adaptive-rebind-regression
```

The test compiled and ran, exiting **1**, with:

```text
FAIL: actual source adaptive rebind cannot publish a freshly checked immutable checkpoint
```

Build output is in `task-5w-red-build.log`; the RED runtime output was captured
in the tool transcript and is quoted above. This was an executable assertion
failure for the missing feature, not a compiler failure.

Initial production build succeeded. Two early integration runs failed at the
same config-retention test because the test assumed upstream config hashes were
only config-field hashes. Inspection confirmed the moment/paper hashes include
dimensions and boundary config hashes include binding output. The test was
corrected to compare every retained config field; production needed no change.
These attempts remain in `task-5w-green-build.log`, `task-5w-green-run.log`,
`task-5w-green-build-2.log`, and `task-5w-green-run-2.log`.
`task-5w-green-build-3.log` / `task-5w-green-run-3.log` show the first complete
GREEN run (exit 0).

Final expanded regression:

```sh
timeout 3600 cmake --build build/task-5c-no-faster-20260909 \
  --target spin-it-2014-source-cage-adaptive-rebind-regression -j3
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 \
  -R '^spin-it-2014-source-cage-adaptive-rebind$' --output-on-failure -V
```

Both commands exit **0**. `task-5w-final-build.log` has no compiler errors or
warnings. `task-5w-final-ctest.log`: **1/1 passed, 0 failed, 17.93 seconds**.
Scoped `git diff --check` also exited 0. Controller will perform the affected
broad suite after reviewing the final commit; worker did not duplicate it.

Coverage includes:

- Actual source outer cube `.3125..6875`, design `.36..515`, source depth 4,
  optimization 3→4, actual coupled 26-handle full BBW, one Full cuboid→eight.
  Hand-checked leaf lattice coordinates, primitive ownership and inherited beta;
  byte-identical full handle state and untouched originals; new points and binding.
- Equal `.01` siblings recursively merge to level 0 with four transition records;
  unequal `.01/.02` siblings never merge. No-change rebuild/check path, twelve
  independent cube adjacency edges and analytic fill penalty `4.5e-6`.
- Static and non-static offsets; all binding/boundary/moment/paper config fields,
  density/frame and shared exact old source/BBW retained.
- Analytic rest and globally affine all-ten-moment checks; independent signed
  tetra reintegration of non-affine before and after. REST transition moments
  checked independently without density. Non-affine difference is recorded,
  not constrained to zero.
- Actual input insufficient-clearance rejection; byte budget exhaustion;
  inherited leaf cap failure at Split; binding cap failure at Rebuild; checker
  exhaustion at OutputGeometry; derivative cap failures at InputMaterial and
  OutputMaterial; actual finite-input arithmetic Unknown at InputGeometry.
  Every failure has no checkpoint and leaves the original intact.
- Dimension, NaN beta/phi/handle fields, and epsilon rejection. Deterministic
  repeated hashes and changed state/settings binding. Const/constructor API
  restrictions and safe shared-evidence move behavior.
- Checkpoint consumption by actual guarded SLQP, new state dimension, initial
  `historyBefore==0`, independent initial and final geometry checks and final
  material hash. Observed numerical status is **RESTORATION_STALLED**, while
  final represented geometry checks pass and feasibility remains **Unknown**.

Non-affine world moment delta (after minus before):

```text
-7.4993997123407041e-11 -3.1112459830673345e-11 -3.4897203182726599e-11
-2.9348870556056283e-11 -1.2718031489056258e-11 -1.5724765239921723e-11
-1.1315934300704100e-11 -9.3754854574301660e-11  8.6623555550779940e-12
 3.5822202842927808e-11
```

Guarded consumer scientific hash:
`5b19b0302ab39c0b47ce9f8770f3efa5d5ce32d73e7d0746df7eebe1522cd084`.

## Bounded output-rejection search and limitations

Optional `--probe-clearance` runs actual transactions with translations
`amplitude*cos((h+1)*(axis+1)*mode)`, all scales 1, original source/BBW unchanged.
Modes 1–12 use amplitudes `.02,.04,.06` in groups of four; modes 13–44 use `.1`.
The two sequential observed batches are retained in
`task-5w-clearance-probe.log` and `task-5w-clearance-probe-2.log`, with corresponding
`task-5w-negative-build.log` and `task-5w-negative-build-2.log`. Both probes exited
0. All 12 mild modes passed both checks. Of 32 stronger modes, 21 passed both
checks and 11 genuinely rejected input geometry. No output geometry rejection
was found; tiny clearance differences were not used to tune a rejection.
The final optional probe combines these same two bounded batches. These are
synthetic adversarial fixtures, not paper cases or optimizer validation.

The unattained output-only geometry-rejection negative stays open. Output-stage
resource/material-cap failures demonstrate nonpublication, not geometric
rejection. Earlier exact-contact and positive primitive-overlap source negatives
remain open and unchanged. This increment adds no adaptive multi-pass loop,
binary validation, Python binding, new artifact seal or paper profile. Native
RSS/CGAL allocator behavior remains guarded externally.

## Self-review and task-owned file list

Reviewed the single checkpoint-constructor path and all failures: no state or
functional escapes before both actual gates and full material evaluations.
Reviewed static offsets, local domain/beta ownership, unchanged config/cap
forwarding, full-BBW retention, separate REST/deformed moments, deterministic
hash inputs, moved-reference lifetime, and absence of numerical history.
Checked explicit scopes against foreign dirty Faster files and controller
artifacts; no changes to them, `progress.md`, dependencies, or legacy solver.

Files for the implementation commit:

- `include/spin2014/SourceCageAdaptiveRebind.h`
- `src/spin2014/SourceCageAdaptiveRebind.cpp`
- `CMakeLists.txt` (source/header/test registration only)
- `tests/spin_it_2014_source_cage_adaptive_rebind_regression.cpp`
- `docs/spin2014_source_cage_adaptive_rebind.md`
- `.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5w-report.md`

Runtime/build logs are retained locally, not staged by the worker. Controller
owns independent evidence artifacts and any later report supplements.
