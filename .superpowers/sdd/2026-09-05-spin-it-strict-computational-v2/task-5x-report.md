# Task5x implementation report

Scope: actual-source guarded adaptive multi-pass orchestration on BASE
`d81bc808806fbff24fd40fad4c6d5e16e5dc3dc0`, in the canonical isolated checkout
`/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`.

Status: DONE. Final focused native verification passed, 1/1 CTest, 35.49 seconds.
Scoped commit: `97f7e12c4a8152779f5ad0691325e261e63b769e`,
`feat(spin2014): orchestrate guarded actual-source adaptive passes`, authored and
committed as `volcano <xuxinyao307@gmail.com>`. Owned paths are clean after commit.

## Implementation

Added `solveSourceCageAdaptive` using only the existing actual immutable source,
fresh full-state guarded SLQP, and private successful adaptive rebind checkpoint.
The public audit retains each returned child audit, source and starting state,
dimensions, observed fresh history, raw numerical/checker evidence and scientific
hash links. Last successfully solved state and changed pending checkpoint are
separate immutable owned objects. Copying facts cannot manufacture either one.
Changed partitions are detected by actual partition hash, not beta count.

Only numerical Converged plus the actual final geometry and full material gate
allows rebind. Initial rejection is classified from the first actual checker
attempt because upstream finalStatus is intentionally Unavailable when no final
check ran. Unchanged completion additionally requires complete identity and no
refinable inclusive-epsilon leaf. Maximum-depth fractional leaves can remain;
there is no binary, feasible-witness, physics, geometry-seal or paper-case grant.

Before work, checked cumulative nonrefunded reservations cover fixed pass
metadata, state/wrapper payload, the complete existing guarded known-buffer
envelope, and each rebind's full configured owned cap. A separate source-pair
construction slot cap bounds retained immutable rebuilt graphs by both
maximumPasses+1 and maximumRetainedSourcePairs. These graphs are not falsely
treated as fully shared or included in the known-buffer byte count. Their
upstream per-object caps and mandatory external process bound remain in force.
This exact bounded-resource design was raised to and approved by the controller.

Nullable pre-admission pass fields are documented: resource failure before
snapshot admission leaves null source/solve and empty initialState, with proposed
input hashes/dimensions only. No attempted numerical work is asserted there.

## TDD and preserved logs

All logs below are raw preserved bytes in this report's directory. The initial
logs were created in the canonical repository root and moved with `mv -n` to
this SDD directory; no earlier evidence was overwritten or reconstructed.

Initial RED, before any production header/source existed:

```sh
g++ -std=c++17 -Iinclude \
  tests/spin_it_2014_source_cage_adaptive_solver_regression.cpp \
  -o build/task-5c-no-faster-20260909/task-5x-red
build/task-5c-no-faster-20260909/task-5x-red > task-5x-red.log 2>&1
```

The compiled native test exited 1 with:
`FAIL: actual guarded adaptive solver cannot retain solved and pending endpoints across fresh passes`.
The pre-header branch follows the existing native regression pattern and signals
the missing feature without a compiler error. `task-5x-red.log` SHA256:
`4f4d3d5b8012cde32cd4acbf9a2aa01cfd0bee71b15e5e6b8389653d51f2162c`.

Focused build command throughout (at most 3 jobs):

```sh
ulimit -v 25165824
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cmake --build build/task-5c-no-faster-20260909 \
  --target spin-it-2014-source-cage-adaptive-solver-regression -j3
```

Build flags remained exact geometry ON, 2014 ON, Faster OFF, and
LIBIGL_SKIP_DOWNLOAD ON. No dependencies were added/downloaded. Controller and
worker explicitly transferred the heavy native slot; worker native runs were
sequential under the stated ulimit, thread limits and `timeout 3600`.

Early real integration runs and their failures are retained:

- `task-5x-native-first.log`: the actual changed/unchanged two-pass path passed;
  aggregate exhaustion assertion failed because the test compared maxPasses=1
  metadata reservations with maxPasses=3. The test now uses the same-settings
  first-pass reservation. This was a test setup correction, not a relaxed cap.
- `task-5x-native-second.log` and `task-5x-native-third.log`: initial actual
  geometry rejection was incorrectly summarized as unknown. Third run records
  finalStatus=Unavailable and first actual attempt=RejectedGeometry. Production
  mapping was corrected using that actual attempt; upstream evidence unchanged.
- `task-5x-native-fourth.log`: all then-current default regression assertions
  passed, including nontrivial nonconvergence and geometry/resource rebind
  failures. Positive cases remained stationary known-solution validation.

The controller's integration check identified an epsilon predicate mismatch.
Added real exact-lower/upper-epsilon regression before changing production:

```sh
timeout 3600 build/task-5c-no-faster-20260909/\
spin-it-2014-source-cage-adaptive-solver-regression --epsilon \
  > task-5x-epsilon-red.log 2>&1
```

Exit 1: `inclusive exact-epsilon leaves refine and remain counted at maximum depth`.
Then changed the descriptive counting predicate from strict comparisons to the
domain's actual inclusive comparisons. `task-5x-epsilon-red.log` SHA256:
`0ebf1358edbba01626b30c8c53998e222ba0b1133da279e0e4f7c88bcb102454`.

The final fixture also uses phi=.23 and nonidentity parameters on all 26 handles;
its two equal transverse principal values make this a stationary known solution
without tolerance relaxation. These tests validate orchestration, not nontrivial
optimizer convergence. Nonzero objective weights remain separately tested as a
bounded failure case.

## Owned changed files

1. `include/spin2014/SourceCageAdaptiveSolver.h`
2. `src/spin2014/SourceCageAdaptiveSolver.cpp`
3. `tests/spin_it_2014_source_cage_adaptive_solver_regression.cpp`
4. `docs/spin2014_source_cage_adaptive_solver.md`
5. `CMakeLists.txt` — only source/header registration and one native CTest target.

No legacy fill solver, Faster source, foreign tracked/untracked work, controller
artifacts, public remaining-gap document, ledger, paper assets, seals or Python
facades were modified. No subagents or reviewers were spawned. Controller owns
independent review, full all-target build and related regression suite.

## Self-review

Checked actual partition identity chaining, preservation of full offsets,
immutable ownership after move/lifetime, separation of last solved and pending
state, all failure paths retaining prior endpoints and returned audits, and
checked arithmetic before reservation/allocation. Reviewed scoped CMake diff and
ran `git diff --check`. All finite-input geometry and numerical failures remain
honest dispositions; no solver infeasibility is inferred.

Material limitations: conservative known-buffer reservations are not complete
object-graph bytes or RSS; distinct immutable graphs have separate count/upstream
caps and the external 24 GiB guard. Allocation failures while allocating or
finalizing diagnostic evidence may propagate, as in upstream interfaces. Positive
fixtures are zero-objective stationary validations, not paper convergence. No
broader pipeline completion is asserted.

## Final GREEN evidence

Final scoped build exited 0 with clean compilation/linking output in
`task-5x-build-final-verification.log`. Then ran:

```sh
ulimit -v 25165824
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 \
  -R '^spin-it-2014-source-cage-adaptive-solver$' --output-on-failure \
  > .superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5x-green.log 2>&1
```

Exit 0. `100% tests passed, 0 tests failed out of 1`; target time 35.49 seconds,
total 35.50 seconds. Preserved actual test stdout in
`task-5x-green-native-output.log` (copied without overwrite from the fresh CTest
LastTest.log). Output includes `KNOWN_SOLUTION_VALIDATION inclusive epsilon
endpoints PASS`, `KNOWN_SOLUTION_VALIDATION multiple passes PASS`, the expected
actual initial rejection facts, and `adaptive solver regression PASS`.

Final `git diff --check` exited 0. The staging area was checked empty before
staging only the five owned files. Native slot was released to the controller
after this pass. The controller independently reported a three-pass fixture with
two partition-hash changes at unchanged beta dimension; its native evidence and
arithmetic checker are controller-owned, not substituted for this worker run.
