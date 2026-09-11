# Task5k controller verification — bounded increment accepted

Worktree: `/data/xvxinyao/make-it-balance-paper-sim`.
Base: `29faa6423d96bccaf7b80a16a502c4686715884a`.
Scope: actual source-cage BBW numerical increment only; `PARTIAL_BLOCKED`
globally. No original-paper case, material-seal, GUI or physics grant.

## Independent reference and first actual native check

Controller-owned scripts in this directory:

- `task-5k-reference-bbw.py`: literal rational brick/triangle construction,
  independent NumPy/SciPy QP and original-KKT active-set polish.
- `task-5k-controller-probe.cpp`: actual native C++ API export.
- `task-5k-check-native.py`: reconstructs geometry/operator independently and
  checks native weights; does not call project numerical code.

Output directory: `build/task-5k-controller-Nmogqv/`.
Reference inputs are fixed analytical regression cubes, not original assets.
NumPy 2.1.3, SciPy 1.15.3, `/home/lin/anaconda3/bin/python3`.
Threads OMP/OpenBLAS/MKL fixed to 1; shell virtual-address limit 24 GiB;
each native probe timeout 300 seconds. This is not an exact RSS cap.

The first probe was linked only after the implementer released the stable
archive; source/header/archive/probe hashes are in
`probe-build-inputs-01.sha256`. The implementer resumed self-review afterward,
so these are **pre-final** results until a final source binding/replay occurs.

Eight configurations (small/medium × two mass rules × two protocols), each
executed twice, independently passed. Every node and handle identity matched;
all actual weights were checked, not only the result's success flag. JSON
scientific outputs including actual backend counts were identical per repeat.

- Maximum raw-weight difference from independent polished/literal reference:
  `1.4224732503009818e-15`.
- Maximum independently reconstructed raw stationarity across the eight:
  `2.979387569990166e-16`.
- Medium normalized-independent final coupled stationarity:
  `0.014648954198848494`; correctly **not** granted coupled optimality.
- Small/medium coupled actual backend calls: 1 each.
- Small/medium independent actual backend calls: 8/26 each.

Evidence: `native-v1-*.json`, `check-v1-*.json`, and per-run stderr files.
Seven checker mutations all rejected: wrong weights despite success=true,
substituted node, duplicate handle, false normalized coupled-optimality claim,
dropped backend call, NaN, and inconsistent energy. See
`checker-negative-v1.log`. These are checker regression tests, not a release
attestation or proof that an arbitrary submitted JSON is an authentic run.

## Adapter decision and upstream limitation

The first actual native runs exposed duplicate sparse `start_` initialization;
the implementer diagnosed and corrected it, preserving failed logs.

The pinned HiGHS bounded-only start also returned `Optimal` with original-FEM
KKT about `0.005138` on an independent column. Controller inspected
`highs/qpsolver/a_quass.cpp:154` and `feasibility_bounded.hpp`: its zero-row
branch forms a dense Cholesky buffer using the sparse supplied pattern,
truncates entries below `1e-4`, and can return interior status immediately.

Approved adapter representation: add one explicit redundant `0=0` row to
**every** independent QP before solving, selecting the general start. It
changes neither objective nor feasible set; no retry/fallback or dependency
edit is permitted. Its identity, resource footprint and rationale must be
reported. Original KKT/finite gates remain mandatory.

A non-dyadic covariance fixture (scale 1.7, root origin [-0.3,-0.2,-0.1])
failed the upstream FEM boundary assembly before BBW. The failure must remain
reproducible and documented as an open upstream limitation, not invalid
geometry or physical infeasibility. Exact dyadic covariance is a separate
positive numerical fixture, not a silent replacement paper experiment.

## Final-source independent native verification

After trace/resource self-review, the implementer froze the header and
source. The controller rebuilt `bbw-probe-final` against that archive and
executed the same eight configurations twice again. All eight independent
checks passed, including `passModelStatus == 0` and actual `runInvoked` for
every expected backend call. The numerical maxima above were unchanged.

Evidence: `native-final-*.json`, `check-final-*.json`, per-run stderr files,
and `checker-negative-final.log` (7/7 mutations rejected). These supersede
the preliminary `v1` outputs as current task verification; old files remain.
`task-5k-negative-checks.py` makes the checker mutations rerunnable.

`probe-final-bindings.sha256` binds the final header, production source,
native test source, actual core archive, probe executable/source, checker and
reference. All eight hashes were rechecked successfully after execution.
Registered-test discovery was separately saved as `registered-tests.json`.
The protected main tracked diff still hashed to
`93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`.
This preservation check covers tracked changes, not an exhaustive historical
inventory of all untracked main-directory files.

An additional actual partial-run probe (`task-5k-partial-failure.cpp`) used the
unchanged medium input and `iterationLimit=30` for independent columns. Under
both mass rules, the first column validated in 26 actual iterations, and the
second ended at the backend-reported 31 with `ITERATION_LIMIT`. The one-count
backend overshoot is preserved, not rewritten to the configured 30. Exactly
two calls were retained; full raw/final assessments were unavailable, final
weight storage empty, and feasibility UNKNOWN. Both checks passed; evidence
is `partial-failure-final.log` and executable `partial-failure-probe`.

## Task acceptance

Implementation `6ba2aa1`, report `1016bbb`. Fresh independent reviewer returned
Spec PASS / quality Approved, no Critical or Important issues. Review and its
limited authority are retained in `task-5k-review.md`; only existing dependency
CMake warnings were deferred. All 11 implementer source/binary/dependency
bindings and every raw-log hash were independently rechecked after review
dispatch. No production change or duplicate full-suite rerun was needed.

Only controller-owned evidence/docs will be committed; externally created
Faster edits must remain untouched. This accepts Task5k, not the full plan.

## Full-suite verification and concurrent external changes

Observed `faster-off-ctest.log`: 4/4 passed, 6.74 seconds. Main all-target build
completed, and `full-ctest.log` recorded 203/203 passed, 499.12 seconds. The
controller parsed all 203 actual unique passing names and compared them with
the complete registered-test JSON: exact set equality, no Not Run or skips.

During that run four unrelated Faster paths changed outside this task:
`include/spin/QuadricMomentEvaluator.h`, `src/spin/QuadricMomentEvaluator.cpp`,
`tests/quadric_moment_evaluator_regression.cpp`, and new
`tests/quadric_certified_empty_cases.h`. The controller and BBW implementer did
not edit them. Observed modification times were 18:12–18:14 +0800; their first
observed hashes are preserved in `external-faster-observation.sha256`. The
user was asked non-blockingly whether another task owns this work.

These changes are preserved, excluded from scoped BBW commits, and not rebuilt
or tested by this controller. Therefore **203/203 is evidence for the actual
tested binary checkpoint, not a full-current-working-tree validation claim**.
Frozen BBW source/header/archive bindings still match the actual independent
probe. A clean-worktree/full-project release gate cannot be granted here.
