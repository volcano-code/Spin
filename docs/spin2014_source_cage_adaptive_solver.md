# Actual-source guarded adaptive solver

`solveSourceCageAdaptive` is bounded multi-pass orchestration of the existing
actual-source `SourceCagePaperFunctional`, `solveSourceCageGuardedSlqp`, and
`rebindSourceCageAdaptive`. It does not replace the legacy fill-only
`SourceAdaptiveSpinSolver`. Every pass invokes a fresh guarded SLQP solve over
all beta values, optional phi, and every original cage handle. A changed
partition is recognized by its actual partition hash, including changes that
leave the number of beta values unchanged.

Only raw numerical `Converged`, together with the guarded audit's freshly
completed final represented-geometry check and full material evaluation, can
enter rebind. Limits, stalls, initial rejection and unknown outcomes stop the
loop without declaring infeasibility. The first actual guard attempt identifies
an initial geometry rejection even though no final geometry status exists.
Every returned guarded and rebind audit remains owned in its pass, including
callbacks, KKT, trial records, geometry records, material hashes and raw numerical
termination. No user callback, copied scalar facts or caller-provided boolean
can authorize progress.

The successful private rebind checkpoint supplies the entire next state and
immutable functional. Rebind preserves phi and every handle using the actual
offsets; the loop transfers no KKT or geometry grant. `freshHistory` is only a
description derived from the first actual numerical trace row's history counts.
The underlying trace remains available for inspection.

`lastSolved()` owns the most recent source/state whose own solve passed these
gates. `pendingCheckpoint()` separately owns a changed, successfully rebound
checkpoint which has not yet passed its own fresh solve. A maximum-pass exit or
failure on the next pass preserves both. Once that pending state successfully
solves, it becomes the solved endpoint; a subsequent failed rebind cannot
publish its candidate. An audit and its endpoints expose only const access.
Their shared ownership supports use after the calling source or audit dies;
moved-from audit accessors throw. Previously obtained references remain valid
while their shared storage has an owner.

An unchanged rebind completes only if the old/output partition, full state and
functional identities match and no refinable fractional leaf remains. The
fractional predicate is the domain's inclusive `[epsilon, 1-epsilon]` predicate.
Leaves at maximum optimization depth can still be fractional. Accordingly,
`RelaxedAdaptiveComplete` is **not binary validation, a feasible witness, spin
stability, a geometry seal, or original-paper case completion**. Fractional
counts describe the last successfully solved endpoint when available, not an
unsolved pending checkpoint.

## Audit and resource accounting

The facts retain settings, numbered passes, old/solved/new state hashes,
old/new partition and functional hashes, both child scientific digests, starting
dimensions, observed fresh history, cumulative reservations and the precise
terminal reason. A streaming scientific digest binds this data and endpoint
presence. It includes configured numerical time limits as settings, but no
observed elapsed time. Configured limits and implementation defaults are not
claimed to be recovered paper parameters.
If admission fails before a pass snapshot is allocated, that pass can have a
null `source`, empty `initialState` and null `solve`; its dimensions and hashes
describe proposed input only. An empty solve pointer never denotes attempted
numerical work or a numerical termination.

`maximumKnownBufferBytes` defaults to 8 GiB. Checked addition and multiplication
admit these conservative reservations before the corresponding work:

- All pass metadata storage: twice `maximumPasses * (sizeof(pass) + 4096)`, plus
  the fixed facts object and 8192 bytes for control/diagnostic/hash overhead.
- Per pass, a full starting-state snapshot and 4096 bytes for immutable wrapper,
  shared-control and endpoint overhead. Solved state storage aliases the guarded
  result's final vector and is not copied into a second endpoint vector.
- Each guarded solve's entire `planSourceCageGuardedSlqpResources().requiredBytes`
  envelope, covering known source arrays, callback and numerical trace tapes,
  KKT/trials, known solver workspaces, retained guard records, bounded strings
  and its conservatively planned checker scratch.
- Before each rebind, its full configured `maximumOwnedBytes` plus 4096 bytes.
  The rebind's own checked plan independently bounds its leaf/state snapshots,
  transition evidence, retained geometry/material evidence and checkpoint state
  payload. Reserving the full cap avoids relying on an after-the-fact estimate.

Reservations are cumulative and never refunded, including scratch which has
already been released and failed transaction reservations. This deliberately
conservative known-buffer bound can reject work that would fit in memory. A
small fixed failure-diagnostic object can exist even below the requested cap.
Allocation failure while constructing/finalizing diagnostics may propagate.

The immutable rebuilt domain, binding, boundary, moment and paper object graphs
are **not included as complete graph bytes** in this known-buffer envelope.
Graphs for different partitions are not assumed to be shared. A separate
`maximumRetainedSourcePairs` cap (default 16), also bounded by checked
`maximumPasses + 1`, admits the initial source pair and one possible retained
graph per actual rebind before construction. Unchanged successful rebuilds also
consume slots. `admittedSourcePairs` counts admitted construction slots, including
failed rebinds; `retainedSourcePairs` counts the initial graph and successful
checkpoint graphs actually retained. Existing upstream per-object caps remain
unchanged. Caller-created copies of facts/endpoints are outside this run's
accounting. This is not an RSS guarantee or a bound on CGAL/HiGHS opaque native
allocators. Mandatory external process guards remain:

```sh
ulimit -v 25165824
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 \
  -R '^spin-it-2014-source-cage-adaptive-solver$' --output-on-failure
```

## Focused regression scope

The native regression uses actual source cubes, a coupled full 26-handle BBW
solve, and real guarded geometry/material checks. Positive multi-pass and exact
epsilon fixtures are explicitly `KNOWN_SOLUTION_VALIDATION`: zero objective
weights and a stationary principal frame isolate orchestration behavior. They
do not demonstrate nontrivial optimization or original-paper convergence.

The regression covers changed then unchanged passes, full dimension/state/hash
links, fresh histories, phi/all-handle preservation, depth-limited fractional
completion, both exact epsilon boundaries, maximum passes with distinct pending
state, deterministic hashes, aggregate reservation exhaustion, checked arithmetic
overflow, graph-slot admission, later adapter dimension failure, actual
nontrivial nonconvergence, initial geometry rejection, failed geometry/resource
rebinds, and immutable move/lifetime behavior. This increment supplies no binary
rounding/validation, sealing, physics validation, Python facade or paper assets.
