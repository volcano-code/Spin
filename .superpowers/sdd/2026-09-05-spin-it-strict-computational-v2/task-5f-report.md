# Task5f implementation report

Base: `d664ab6812b8bbb476b7fb77272671df62464262`, branch
`spin-it-strict-computational-v2`, repository `/data/xvxinyao/make-it-balance-paper-sim`.
Scope is the bounded source Spin2014 actual solve to independent material seal
and synthetic runner. Project status remains `PARTIAL_BLOCKED`.

## Implemented

- Added focused `bindings/SourceSpinMaterialBindings.{h,cpp}`, conditionally
  compiled/registered for 2014 plus exact geometry, independently of Faster.
  Registration reuses existing lossless mesh readers and scalar ingestion.
  The entrypoint parses complete settings, constructs a copied source octree,
  and calls the existing private-terminal material wrapper exactly once.
  Neither native nor Python run results have a public constructor accepting
  status, moments, hashes or an arbitrary successful result dictionary.
- Added frozen source, adaptive/audited SLQP and boundary-cap Python settings.
  Every native control is bound and validated; defaults remain
  `IMPLEMENTATION_DEFAULT`, formulation is `PAPER_GIVENS`, and the only supported
  algorithm is `AUDITED_KKT_RESTORATION_V1`. Actual parsed native settings are
  returned as owned diagnostic data, including the fixed audited HiGHS controls.
- Exposed source/config/partition identities and source statistics, every pass's
  original/final partition and state, split/merge transitions, complete ordered
  SLQP audit records through `writeSlqpAuditJson`, beta/phi, binary cuboids,
  numerical checks, exact boundary/partition/clearance facts and scoped resource
  facts. Nonfinite unavailable values are null with enclosing state/status;
  rejected initial evaluations do not publish default zero objectives/moments
  as measurements. Arrays and nested diagnostics are copy-on-read.
- Added `verification/source_spin_material.py`. The facade reauthenticates the
  locked problem, requires one design region/uniform density/native +z and
  origin zero, checks objective/mode/support compatibility, and retains original
  model-coordinate surface bytes. Boundary wall/density come from the problem.
  It independently checks original design/cavity boundary noncontact, constructs
  retained design from original design plus reversed signed cavity, and invokes
  the existing complete exact partition and independent Python tetra seal.
  The claim comes from the actual native source-subtraction moments, transformed
  by full affine first/second moment formulas and R L^3/4/5 scaling. No geometry
  normalization roundtrip or fake full attestation was introduced.
- Added a reproducible fixed-fixture synthetic runner, with exclusive output
  creation, actual pre-solve NPY/config inputs, existing begin/finish/freshness
  observations, complete scientific/native trace outputs, independently saved
  and reloaded geometry bundle, and retained failures. Its path-independent
  scientific identity is distinct from actual path/time-dependent evidence.
  A generic trust-the-run-JSON loader was intentionally not added.
- Added user-facing documentation of entrypoints, frames, controls, supported
  modes, status separation, provenance and geometry-only bundle persistence.

The native geometry budget is explicitly scoped to the source-material chain;
it does not count the subsequent independent Python seal's native audits, total
process memory, CGAL allocator internals or triangle-distance storage.
Build/source attestation stays `NOT_VERIFIED`, physical validation
`NOT_IN_SCOPE`, original admission `BLOCKED_UNKNOWN`, and target feasibility
`UNKNOWN`.

## TDD and iterative evidence

All paths below are relative to repository `build/`. Commands used the selected
`/home/lin/anaconda3/bin/python3`, `PYTHONPATH=python:build/strict-v2-spin/python:tests/python`
and `LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6`. Native work was bounded by
`ulimit -v 25165824` and `timeout 3600`, compile parallelism at most 3.

1. Before production changes, `python3 -m unittest test_source_spin_material_bridge -v`
   failed its actual native availability assertion: `missing actual native
   source material solve entrypoint`. Retained `task-5f-red-missing-01.log`.
   Expanded contract scaffolding then retained the same assertion plus missing
   Python-module setup error in `task-5f-red-contract-01.log`. This establishes
   availability RED, not separate numerical-behavior RED for every later assertion.
2. First focused build succeeded (`task-5f-build-01.log`). The first test run
   exposed a test harness mistake: class attribute `run` shadowed unittest's
   `run` method. The log is retained (`task-5f-focused-01.log`); the fixture
   attribute was renamed `solved`. The actual 8-test bridge suite then passed
   in 16.357 s (`task-5f-focused-02.log`).
3. Runner missing-file assertion failed before its implementation in
   `task-5f-red-runner-01.log`. The first implemented runner attempt failed the
   existing strict evidence schema because a settings tuple had not been
   canonicalized to JSON arrays (`task-5f-focused-03.log`). The attempt's inputs
   and `failure.json` remain under `task-5f-runner-test-5_eu2g7h/success/`.
   Corrected canonical input data and used the existing seed/protocol schema.
4. Actual missing bridge observation property was captured in
   `task-5f-red-observations-01.log`. Added before/after observations separate
   from scientific identity. The resulting full 10-test suite passed in
   33.197 s (`task-5f-focused-04.log`). Expanded unsupported/ownership boundary
   checks passed in `task-5f-focused-05.log` (8 tests, 23.024 s).
5. A direct native run lacked actual parsed settings diagnostics; the missing
   field was observed in `task-5f-red-native-settings-01.log`. Added immutable
   native settings serialization and scoped-budget labels. Build succeeded in
   `task-5f-build-02.log`; all 10 focused tests passed in 35.345 s in
   `task-5f-focused-06.log`.

All failed logs and runner directories were retained. Later code changes make
older source observations stale as expected; they are historical attempts, not
the final demonstration evidence. Successful runner tests deliberately tamper
their own persisted inputs at the end to verify stale detection; those test
directories should not be presented as current demonstration bundles.

## Covered behavior

Real native regressions cover repeat actual successful solves, source fixture
byte ordering, unchanged golden input/trace/output hashes, all 39 trace records
(0..37 then 0), 72 finest cuboids, original output bytes, all ten hand-derived
normalized moments, complete central inertia, 503/64 mass, 5/8 clearance and one
new cavity. Density two doubles model moments but leaves normalized moments
unchanged. The independent translated analytic check exercises all affine cross
terms without presenting a new solved target.

The exact previous quarter-fill fixture stays `BINARY_VALIDATION_FAILED` with
no artifact. All five empty-source objectives stay `RELAXATION_UNKNOWN`, with
null invalid KKT measures and no artifact. Cap 4 and nextafter(5/8) wall requests
retain `BINARY_MOMENT_VALIDATED` numerically but return no geometry seal.
Tests reject unsupported multiple design regions, unequal densities, arbitrary
axis/origin/support, zero axis, buoyancy, deformation scope, invalid enums and
settings, bool scalars/counts, lossy integer coordinates/settings, float faces,
extended precision, unknown settings and caller construction of result objects.
Native and facade array/nested-copy ownership are tested against mutation.

Actual saved source artifacts are reloaded through the existing mesh parser,
exact checker and independent integration; wrong expected problem digest and
tampered mesh bytes reject. Existing generic `test_spin_bundle` coverage is
reused for tampered claims and stale native/checker observations, since the
source result uses that same unchanged geometry seal persistence contract.
Actual runner tests retain success and cap failure, confirm original event
outcomes, reject overwrite, and re-read changed inputs as `STALE_RESULT`.

## Independent Faster-OFF verification

The separate existing `build/task-5c-no-faster-20260909` cache was configured with
Python ON and Faster OFF, preserving its 2014 ON/exact ON settings. Configuration
and actual extension build passed (`task-5f-no-faster-configure-01.log`,
`task-5f-no-faster-build-01.log`). The build took 6.50 s, maximum reported process
RSS 824728 KiB (not aggregate concurrent RSS).

An actual import asserted the loaded module from that cache exposes the source
entrypoint and does not expose `solve_spin_it_faster`; the printed result was
`source True faster False` (`task-5f-no-faster-loaded-01.log`).
`ctest --test-dir build/task-5c-no-faster-20260909 -R '^source-spin-material-python$' -j1 -V`
passed 1/1 registered test with all 10 native-backed Python tests in 24.109 s,
CTest 24.20 s, no skipped/absent-extension success (`task-5f-no-faster-python-01.log`).

## Final main verification

Main all-target build passed, exit 0, with `cmake --build build/strict-v2-spin -j3`
(`task-5f-all-target-build-01.log`). The final all-target invocation was incremental;
the changed extension TUs had already been compiled by the two focused builds.
The final registered main `ctest --test-dir build/strict-v2-spin -j1 -V` passed
**198/198 in 445.61 s**, exit **0**, no failed/NotRun/Skipped tests. Output was
captured directly to `task-5f-full-ctest-verbose-01.log`, not copied from the shared
LastTest file: **2,418 lines / 267,121 bytes**, SHA-256
`0dead466ab01f868a1377450b23366353e7cae77ccee87a0e54e616fc3977409`.
The measured maximum reported process RSS was **4,253,316 KiB**, not aggregate
concurrent RSS. This includes the new source bridge (10 Python tests, 35.596 s),
unchanged source native regressions and depth-9 cases, generic exact/seal/bundle
checks, Faster bridge, old lattice/deformation APIs and feature dependency contracts.
No source/test edits, builds or reconfiguration followed this final full suite.

Final source/binding/checker/test/cache/actual module file hashes are recorded in
`task-5f-source-build-checker-hashes-01.log`; all 15 entries matched after the
suite (`task-5f-source-build-checker-recheck-01.log`). The original solver/material/source files and existing native regression
are unchanged from base; no native algorithms, thresholds, geometry math or
golden values were retuned.

## Changed files and self-review

Scoped changes: `CMakeLists.txt`, narrow registration in `bindings/pybind_module.cpp`,
the new binding TU/header, `python/make_it_stand/verification/source_spin_material.py`,
`experiments/spin_it_reproduction/strict/run_source_spin2014_material.py`,
`tests/python/test_source_spin_material_bridge.py`,
`docs/spin_it_2014_source_material_bridge.md`, and this report.

Self-review inspected strict input ingestion, conditional linkage, owned result
lifetimes/GIL release, all mapped config fields, null/availability diagnostics,
complete native traces, unchanged byte geometry, normalization cross terms,
independent partition/integration, failure propagation, exclusive persistence,
and separation of scientific/provenance identities. No unaddressed functional
issue was found within this bounded increment. Native diagnostic serialization
is intentionally kept in one focused TU; legacy solver/core APIs are untouched.
Historical deferred CMake/Boost/readability/Task5b diagnostics remain out of scope.

Controller brief/math/review/ledger/global progress and original outputs are not
staged. No helper/reviewer agents, GUI edits, paper/asset searches, main writes,
pushes, tags, releases or merges were performed. Two final frozen-source
demonstration runs and independent review belong to the controller after handoff.

The runner tests used the actual installed pybind11 include tree as their declared
dependency root. Controller demonstrations can declare the actual `libigl/` and
main `_deps/highs-src/` trees; a read-only symlink scan found none under either
tree, and the HiGHS source path resolves to itself. Dependency coverage remains
`DECLARED_ROOTS_ONLY`. No extra demonstration was run during the final suite.

Production/configure/build ownership is released with the scoped commit at
handoff; no background build or test process remains.
