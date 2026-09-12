# Actual-source joint SLQP adapter

`solveSourceCageSlqp` integrates the immutable `SourceCagePaperFunctional` with
the existing `EqualityBoxProblem` / `solveSlqpLbfgs` outer algorithm. It requires
`AuditedKktRestorationV1`. The state is `[beta, optional phi, all original handles
(scale, tx, ty, tz)]`. Only beta has box bounds. Static mode has no angle: its
first scale remains a scale, and all static phi diagnostics are exactly zero.
Rotational phi is the raw coordinate at B; it is never wrapped. The source
functional's strictly positive scale domain is checked through numerical trial
evaluation, without clamping scale or omitting a handle.

Every actual callback evaluates the source paper functional and moves its joint
gradient, constraints, and Jacobian into the generic evaluation. Its objective
is paper objective plus fill, shape, and displacement penalties. Generic mass
moments are in the functional's fixed spin frame. The returned final source
evaluation separately retains world and spin moments, all original/material
deformed points, and all derivatives. It is evaluated explicitly at the final
accepted raw state, never selected from a rejected trial. Frame origin does not
follow the moving center of mass.

## Numerical status and downstream admission

`outcome` says whether the adapter completed, encountered `ResourceLimit`, or
had an ordinary `InitialEvaluationFailure`. `numericalResultAvailable` is false
if the first actual source callback fails before the generic engine returns.
That result preserves its attempted-state record and source failure digest,
without manufacturing a generic state, trace, KKT witness or termination.
Whenever available, `numerical.status` preserves the underlying termination.
A resource-limited result is not a successful solve even if that underlying
status is `Converged`. `finalEvaluationAvailable` is false if the resource latch
or exhausted source-call budget prevents the final accepted-state evaluation.

Every returned result has geometry `NotEvaluated` and feasibility `Unknown`.
Convergence means only the generic `FIRST_ORDER_KKT_TOLERANCES_MET` test. It is
not a minimum, second-order condition, continuous feasibility certificate,
original-case result, recovered author trajectory, or original-paper full-SLQP
grant. This type has no route to a physical/sealed artifact release gate.
Inversion, intersection, overlap, wall-thickness, final partition and binary
material checks remain downstream work. Adaptive partition/BBW transitions
require a new immutable functional and fresh solve/history; this adapter does
not orchestrate those transitions.

## Deterministic provenance

The call tape preallocates full state columns and fixed-size records. Only its
first `recordCount` entries are events; unused storage is zero/`Unused`.
`callbackCount` counts generic callback invocations, including invalid or
resource-latched attempts. `sourceEvaluationCount` counts actual source calls,
including the explicitly marked `FinalAcceptedState` evaluation. Each event has
its sequence, purpose, actual-call flag, full-state digest, output digest or
failure category and reason digest. Displayed reasons may be truncated to 255
characters; SHA-256 covers the entire exception message without copying it to
an unbounded string. Invalid numerical trials can be backtracked and do not
certify infeasibility.

The configuration identity binds source/paper/shape weights, density, upstream
source/BBW/frame identities, full initial state, all solver/LP controls and all
adapter budgets. Source identities transitively bind upstream budgets. The
result binds that identity, call tape, corrected generic trace hash, raw active
membership/counts, final state, outcome, and final geometry/moments output digest
when available. SHA streams fixed-endian integers and IEEE double bits. Generic
trace SHA is recomputed after phi correction. No clock values enter scientific
identity. These digests are internal scientific provenance, not executable or
source-build attestation and not a manual approval.

## Explicit implementation defaults and simultaneous memory

All controls in the public settings structs and the accounting version are
`IMPLEMENTATION_DEFAULT`, not recovered author values. Solver defaults are 120
iterations, 24 search steps, memory 8, equality tolerance 1e-8, first-order
tolerance 1e-7, active tolerance 1e-10, merit penalty 100, reduction .5 and Armijo
1e-4. Audit defaults are trust bound 1, LP iterations 10000, LP time 10 seconds,
and feasibility tolerance 1e-9. HiGHS uses one thread, parallel off, seed zero,
simplex, and presolve on. Adapter caps default to 100000 state values, 10000
actual source calls, 1000 trace rows, 100000 call records and 512 MiB known-buffer
budget. Caps reject zero and unrepresentable values; they never shrink the
problem or increase silently. The public state cap cannot exceed `INT_MAX/6`,
conservatively covering every objective because the unchanged generic trace
serializer visits K×N matrix coefficients using an `int` counter and K can be
six. Preflight additionally checks the actual K×N coefficient count with checked
multiplication, as well as generic/native state loops and HiGHS N+2K indexing.

`planSourceCageSlqpResources` checks additions, products, bytes and matrix index
capacity before size-dependent adapter allocation. Let N=B+a+4H, K=2/5/6,
I=iterations, L=search steps, M=history capacity, P=all boundary points,
O=original outer points, T=I+1, A=2+3IL and D=(K+5)N+2K. The deliberately
conservative envelope sums these simultaneous categories:

| Category | Accounted storage |
| --- | --- |
| Call tape | A full N-vectors and A fixed records; final-purpose slot included |
| Source outputs | Two complete 3P geometry, 10(B+4H) world derivatives, (11+K)N joint/spin derivatives and K constraints, plus object/string overhead |
| Source workspace | Split state, beta fill vector, two H weight rows, 6O shape arrays and 4H shape handle gradient, plus fixed scratch |
| Generic trace | (3T+2) complete rows: capacity growth up to 2T, finish-return full-result copy T, local rows; D doubles, N membership bytes, bounded strings, and twice capacity for up to 3L trial records per row |
| History | 2MN vectors, M reciprocal curvatures and deque/container overhead |
| Solver workspace | 32N+12KN+16K²+8M doubles plus fixed overhead; current/candidate/result evaluations, directions/KKT, free/inverse constraint arrays, Schur/COD and known Eigen GEMM packing |
| Explicit LP arrays | Four times 16 max(2N(K+1), K(N+2)) + 24 max(2N,K) + 24 max(K+1,N+2K) bytes, plus fixed scratch; value/index growth, row scratch and known model copy |
| Generic hash workspace | Eight times the full hexfloat serialization bound, including stringstream capacity and str() copies; 32 bytes per number, bounded scalar/status/trial text; fixed additional SHA workspace |
| Result | Initial state, result object, fixed overhead; generic final evaluation/state duplication is included in solver workspace |

The exact checked formula and `sizeof` object sizes are exposed through the
resource breakdown. Adapter result return uses moves/return-value elision, so
it does not copy the entire tape. No N×N/full Hessian, P×H weight matrix, or
duplicate mapped world/Jacobian is created. The existing generic hash function
does materialize a full serialization, which is explicitly budgeted; adapter
hashes stream. Existing upstream immutable geometry/BBW/operator storage and
construction caps remain separate. Source evaluation scratch is additionally
counted here, without weakening upstream derivative/shape caps.

This is a conservative bound on known native and adapter buffers, not exact
process memory. Allocator bookkeeping and opaque HiGHS factorization/internal
allocations are not falsely claimed exact. Run under the process ceiling of
24 GiB, LP iteration/time controls, and one OMP/OpenBLAS/MKL thread. Numeric
termination can depend on a real LP time limit although elapsed times are not
hashed.

Resource exhaustion latches before subsequent source work. The preallocated
tape still records later generic callbacks, whose exceptions may be swallowed
by line search, and the final-purpose attempt. Preflight resource failures before
a source call throw a clear `RESOURCE_LIMIT` length error; invalid initial
state/settings throw an invalid-argument exception. A failing first source call
returns typed evidence with `numericalResultAvailable=false`. Resource exceptions
inside recoverable trials preserve the engine's returned state and trace.
Callback bookkeeping allocation failures mark a fixed resource reason/status
without allocating. Finalization recovers missing state/reason digests and also
streams each retained raw state into the tape identity. A returned resource
result cannot have `Unused` or incomplete attempted records: if recovery itself
cannot finish, the adapter throws a RESOURCE_LIMIT exception instead.
An actual allocator failure outside callback recovery can only throw a
`RESOURCE_LIMIT` length error: the unchanged generic engine does not expose an
in-progress trace during its own allocation failure. No success is reported.

## Native regression

`spin-it-2014-source-cage-slqp` uses actual source octree, full cage, grid, FEM and
BBW on independently specified small synthetic boxes. They are validation
inputs, never original paper cases. All five objectives use nonstationary full
112/113-coordinate states with nonzero fill/shape/displacement weights and fixed
six-iteration settings declared before observations. Tests inspect accepted
steps, beta/handle participation, complete output/state/trace, scalar independent
KKT algebra, invalid/corrupted traces, repeated identities, resource boundaries
and callback exhaustion. A test-only global allocator sweeps real one-shot
65-byte allocation failures, including callback state hashing, and checks that
returned resource evidence is complete. No production test hook is present.
Neither the tests nor this integration assert that all
modes should converge within that small budget.
