# Task5w controller verification — bounded adaptive rebind

Implementation: `434be7d0a9a323048e058f18bb913cf9ee1a7bb5`.
Base: `60b31e0f3d8568b045a0428b94c54e973c3228a2`.
Working directory: `/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`.
Scope: **SYNTHETIC_REGRESSION**, overall **PARTIAL_BLOCKED**.

## Outcome

The bounded immutable split/merge + full-BBW/material rebind transaction is
implemented, reviewed and verified on actual small source geometries. It
reconstructs the immutable dependency chain, remaps all state coordinates,
rechecks geometry/material and publishes only a complete checkpoint. It is not
the full adaptive optimization loop, binary validation, sealing or either
paper's original experiment matrix.

## Independent final verification

The controller used a separately fixed synthetic source: outer cube
`.3125..6875`, design `.355..52`, source depth4, optimization initial3 or4,
maximum4, actual coupled BBW with all26 handles, density **1.17**. These values
differ from the worker fixture and must not be mixed into one protocol.

Six cases: rest/affine/non-affine/static splits, equal-sibling merge, unequal
near-bound no-merge. Four splits actually change one Full primitive into eight
and add sampled geometry points, not only ownership labels.

- Each case ran twice against the frozen core; both complete JSONL files are
  byte-identical, SHA-256:
  `f11a9ed09e4f471809e002da319a470d1bb28dc1e7de0f6171612cdf4d8400d2`.
- Independent Python closed-cell lookup and eight-node interpolation reconstruct
  all deformed vertices from exported actual grid/BBW/handle data. No native
  material/checker integration function is called by the checker.
- Separate exact rational signed-origin tetrahedral integration verifies ten
  material moments, all beta derivative columns, primitive ownership and rest
  leaf contributions. Maximum normalized material error: **1.8000784988681468e-15**.
- Raw state comparison verifies unchanged phi/all104 handle parameters across
  old/new beta layouts; rest/full affine moment conservation is checked.
- Non-affine subdivision changes the represented moments: maximum absolute
  component change **7.003670560168374e-11**. This is recorded, not hidden by
  forcing deformed conservation from rest-domain facts.
- The static checkpoint is fed to the actual guarded solver. Initial history
  is empty, remapped state is used exactly, and final geometry is checked anew.
  Status: **MAXIMUM_ITERATIONS**, not convergence. Worker's distinct consumer
  reports RESTORATION_STALLED; neither is promoted to feasible/stable spin.
- Sixty actual-export corruptions are rejected by the independent checker:
  stale/nonfinite moments, replacement BBW identity, invalid handle layout/count,
  wrong density/point/owner/derivative/rest moment. These are checker-negative
  tests, not native geometry rejection experiments.

Final native runs: **13.22s / 13.23s**, maximum RSS **28,484 / 27,812 KiB**, exit0.
Complete resource logs and raw JSONL are retained as task-5w-controller-run1/2
and task-5w-controller-resource1/2. Earlier preview logs remain clearly preview;
they are not substituted for these final runs.

## Build and tests

All-target build succeeded in the existing Release, 2014/exact ON, Faster OFF,
skip-download configuration, with three compile/link jobs. Worker focused test
passed1/1 after executable missing-feature RED; its earlier failed test
assumptions and all negative-search attempts remain preserved.

Controller freshly ran:

```sh
ulimit -v 25165824
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout 3600 cmake --build build/task-5c-no-faster-20260909 -j3
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 \
  --output-on-failure -j1 \
  -R '^(spin-it-2014-|spin-2014-|source-spin-material-python$)'
```

**44/44 passed, 0 failed, 369.08s**. No skipped, disabled or Not Run result is
counted as a pass. This is the related2014/dependency/native/Python group, not
a claim that every registered project test or any Faster experiment ran.
The initial controller baseline was3/3. Heavy native runs were sequential;
compiler address-space guard24GiB and per-command timeout3600s were inherited.

## Binding, staleness and replay

Actual core SHA:
`108cd55607df2a9d022413314e5814de01be85889481494ab4adcfd5829cb357`.
Actual final controller executable SHA:
`0ce558d79bdb4e55a063ed1ce0c84efb38eab2dde9b37abe1b013d98ed7ba140`.
Actual native regression executable SHA:
`5c5ff233e50afcc78a24ceae2463641e5055f799ea7574e36501446d5d1fb3e4`.

task-5w-controller-inputs.sha256 contains **3,436** entries: actual Ninja
compiler-discovered archive source/header inputs, CMake/cache, compiler,
Python, core/dependency libraries, controller/checker/helper source and runner.
The replay verifies this manifest before compiling/running. It does not include
report output, avoiding circular source/report hashes. Runtime ELF libraries
identified with ldd are additionally recorded in task-5w-controller-runtime.sha256
and independently checked; that supplemental runtime inventory is not silently
described as part of the replay's primary input manifest.

The old Task5v replay now correctly stops before execution with **STALE_RESULT**
because CMake and core hashes changed. Old logs are unchanged; historical results
remain evidence for their old source, not newly blessed current results.

Replay on this local environment:

```sh
sha256sum --check .superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5w-controller-runtime.sha256
bash .superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5w-controller-replay.sh
```

Full replay stdout is task-5w-controller-final-replay.log. It verifies all primary
inputs, rebuilds an independent exporter, runs twice, compares outputs, checks
independent geometry/BBW/moments and executes all60 corruptions. It does not
grant a paper or release status.

## Build Graph and independent review

Build Graph skill was explicitly requested and used. Unavailable MCP tools were
replaced by the already installed matching CLI, not a new installation. The
graph plus direct source inspection identified the old fill-only adaptive
endpoint and immutable actual-material chain as the integration boundary.

Initial index410files/5569nodes/59303edges. Its diff-based updater excludes
untracked files; after the normal source commit, it parsed3changed files.
Final index **413files / 5624nodes / 60014edges**, branch/commit matched434be7d.
The new exact-qualified function query returns54extracted callees. C++ member
resolution/import coverage remains incomplete: graph absence is not proof of
no dependency. Raw commands/statistics and interpretation are preserved in
task-5w-build-graph.md and task-5w-graph-*.

Independent task reviewer: spec compliant, quality **Approved**, no Critical
or Important finding. Minor: compressed transaction-stage lines at
SourceCageAdaptiveRebind.cpp:168; retained with prior readability findings.
Review is scoped, not the whole-branch release audit.

## Preservation and remaining work

Protected main tracked-diff SHA was checked before and after:
`93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`.
Only the isolated branch and explicit task-owned files were committed. Foreign
Faster tracked/untracked work remains excluded from changes/staging and no
Faster validation is claimed. No dependencies installed, main integration,
merge, push, tag, original logs deletion or asset substitution occurred.

Open: reliable passing-input/rejected-output geometry negative. All44 bounded
probe modes are retained (33both-pass, 11input rejections); resource/material-cap
failures and roundoff-level clearance differences do not close it. Old actual
contact/positive-overlap negatives remain open. Next substantive implementation:
full adaptive guarded multi-pass orchestration with partition/history records,
then binary material recomputation, geometry validation and sealing. Continuous
BBW-map legality, original author assets/parameters, full case matrices, Faster
gaps and physical-scope exclusions remain unchanged.

## Evidence packaging

Raw configure/build logs preserve original Boost trailing spaces and dependency
warnings. The verbatim parent-task extract preserves its trailing blank line.
The unrestricted staged whitespace check reports only those three archived
files (task-5w-green-build.log, task-5w-red-build.log, task-5w-parent-task.md).
All other staged files receive the normal strict whitespace check; the complete
stage is additionally checked with only those whitespace categories disabled.
Scientific outputs are not reformatted. Artifact copy emitted the system's
`cp -n` portability warning; copied source/destination JSONL bytes were compared.
The packaging warnings are not omitted test failures or scientific acceptance.
