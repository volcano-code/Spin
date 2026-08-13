# R7 alternating balance-optimization acceptance record

Specification audit date: 2026-08-12 (Asia/Shanghai)

Acceptance status: **accepted at the documented R7 status ceiling**

The sixteen focused R7 tests pass on the current working tree. This accepts
the deformed-cell mass sidecar, fixed-support target scheduler, transactional
deformed carving, all-filled occupancy reset, and production block-coordinate
controller at `VERIFIED_UNIT`, with the local overhanging-prism controller
fixture at `VERIFIED_REGRESSION`. It is not a claim of full-paper
reproduction, agreement with an official asset, fabrication readiness, or
physical standing success. No R7 result may be labelled `VERIFIED_PAPER`.

## Authority and phase boundary

The normative source is Prévost et al., *Make It Stand: Balancing Shapes for
3D Fabrication*, especially Sections 3.2--5.3 and Figure 3. R7 composes the
already accepted R2--R6 components into the paper's outer alternation. The
official beta archive remains a read-only behavioral oracle and is not copied
or linked into this repository. Where the printed paper and inspected beta
differ, the paper governs this clean-room implementation.

The optimizer state is represented explicitly as

```text
(alpha, q, M_O(q), M_I(alpha,q), m, Q, c, target context, eta, mu),
```

where `alpha` is binary voxel occupancy, `q` is the unlocked-handle parameter
vector, `m` and `Q` are signed mass and first moment, `c=Q/m`, `eta` is the
line-search step, and `mu` is the CoM/shape tradeoff. The selected support,
suspension attachments, rest lattice, BBW field, handle binding, fixed handle
rotations, density, and mesh-specific shape scale remain fixed throughout one
run. A standing target may be recomputed from the current CoM and the same
fixed support; recomputing a target does not rebuild that support.

Within this boundary R7 implements:

- all-filled occupancy reset at the start of every carving block;
- exact mass and centroid integration for synchronously deformed hexahedral
  cells;
- target refresh after the full-solid reset and after every committed carving
  or accepted deformation update;
- sequential carving steps whose next step consumes the newly refreshed
  target set;
- a fixed-`alpha` deformation block using the R6 analytic objective and paper
  line-search state machine;
- transactional state snapshots, rollback, context-generation records, and
  explicit terminal/guard statuses; and
- a returned final snapshot containing occupancy, deformed geometry, mass,
  first moment, CoM, parameters, target state, and both block traces.

R8 remains responsible for deformed minimum-wall verification, final mesh
construction and repair, watertightness/connectivity validation, slicing, and
manufacturing checks. R9 adds the Python/pybind11 experiment interface; R10
remains responsible for official case assets, paper-scale experiments,
parameter sweeps, and paper-table or paper-figure comparison.

## All-filled reset and per-step target schedule

The first target snapshot is constructed from the CoM of the all-filled outer
solid at the supplied initial `q`. Every later carving block begins from a
value copy of the immutable rest/topology grid and calls `refillAll()` before
planning any removal. Reset changes only fill flags: rest corners, voxel
indices, signed-distance metadata, carve eligibility, and shell metadata are
not rewritten. The trace records refill count, pre/post hashes, Hamming
distance, target generation, and geometry epoch so that a new block cannot
silently inherit the preceding carved pattern.

At the reset geometry, `BalanceTargetSchedule` refreshes the current target
snapshot. For a standing objective it projects the current CoM into the safe
region derived from the fixed support and requested toppling angle. For a
suspension objective it retains the fixed attachment and gravity while
reevaluating downstream direction and angular tolerance. Specifications retain
their original order, and termination uses all-of rather than any-of
semantics. `OutsideSupport`, `FeasibleButUnsafe`, an unattainable requested
margin, an upstream suspension state, and a zero attachment offset are all
non-terminal.

One carving iteration follows this sequence:

```text
freeze current target context
plan and transactionally commit one best carving prefix
record same-context energy before and after the committed prefix
refresh targets from the committed c
record the new generation and refreshed-context energy
stop if all constraints are reached; otherwise update the next-step targets
```

`CarvingSolver::updateBalanceObjectives()` changes only the target collection.
It does not change `alpha`, `m`, `Q`, `c`, counters, historical diagnostics, or
terminal state; empty input and reentrant in-step updates are rejected with
rollback. Thus the next carving step uses the fresh target without reopening a
terminal solver or retroactively reweighting old records.

The two energy values inside one frozen carving context are a valid descent
pair. A target refresh may change the objective itself. The trace therefore
marks that generation boundary as `objectiveContextChanged` and records the
refreshed energy separately; it does not claim descent between pre-refresh and
post-refresh objective values.

## Exact deformed hexahedral mass field

The rest-grid `Voxel.center` and axis-aligned `voxelSize` are invalid mass
oracles after LBS deformation. `DeformedVoxelMassField` instead binds every
active cell to its eight synchronously deformed global corners. It integrates
the cell as an outward-oriented closed boundary of twelve triangles. For each
oriented triangle and a local reference point it accumulates signed
tetrahedral mass and first moment, then derives the exact piecewise-linear
volume centroid as `Q_cell/m_cell`. This centroid is generally not the average
of the eight corners under non-affine distortion.

The field stores one positive `(m_cell,Q_cell,c_cell)` tuple per active cell.
It rejects non-finite coordinates, inverted orientation, collapse, mass below
the scale-aware numerical floor, a user-specified positive mass floor, and
out-of-range access. During deformed carving, candidate plane scores use
`c_cell`, and prefix simulation/subtraction uses the same `m_cell` and
`Q_cell`; the solver never falls back to rest centers or rest box volumes.

Identity deformation agrees with the legacy axis-aligned result. Translation,
rotation, and uniform scaling obey the expected mass/moment covariance. For a
distorted three-cell cavity, the sum of cell contributions agrees with an
independent inward cavity-surface integral. After a production carving block,
the committed outer-minus-cell state is also compared with the synchronized
outer-plus-cavity surface oracle before the block snapshot is returned.

## Closed oriented triangle-chain contract

Voxel cavities can contain closed components that touch along an edge. Such a
surface may be a valid integration chain even though more than two triangles
meet at that indexed edge and the connectivity is not a strict triangle
2-manifold. Requiring exactly two incident faces would incorrectly reject this
legitimate algebraic boundary.

The mass-differential path therefore uses a closed oriented triangle-chain
predicate: every directed use of an edge must be paired by an opposite use,
so the chain has zero boundary. Open or directionally unbalanced chains are
rejected. The supporting `surface-mass-triangle-chain-analytic` regression in
the complete headless suite checks mass, first moment, and CoM gradients on an
edge-touching closed chain and separately confirms that it is not a strict
manifold.

This is the mathematical closure condition needed by the divergence-theorem
integral and analytic differential. It is not a claim that the geometry is
self-intersection-free, vertex-manifold, watertight for a slicer, or directly
printable; those geometric and manufacturing predicates remain R8 work.

## Paper-governed block-coordinate controller

`AlternatingBalanceOptimizer` alternates two distinct transaction domains.
The carving block holds `q` and the deformed lattice fixed, starts from the
all-filled occupancy, and commits complete-prefix removals one step at a time.
The fixed-`alpha` deformation block then holds occupancy and cavity
connectivity fixed while varying only `q`. No deformation trial calls
`refillAll()` or invokes `CarvingSolver`.

For each proposed deformation state the controller constructs handles from the
immutable binding, evaluates synchronized outer/lattice/cavity geometry, and
constructs a `DeformedVoxelMassField`. A non-positive scale, non-positive
material mass, inverted cell, or collapsed cell is translated into an
`InfeasibleDeformationState` and treated as an expected trial rejection.
Unrelated internal errors propagate rather than being recast as ordinary line
search failure.

An accepted trial atomically advances `q`, geometry, `m`, `Q`, `c`, and the
geometry epoch, then refreshes the target schedule. A rejected trial leaves all
of those values and `alpha` unchanged; only the step/`mu` state may advance
according to the R6 state machine. Every deformation record stores parameter
vectors before and after, occupancy hashes and Hamming distance, mass/CoM after
the decision, target generations, and explicit zero refill/carving-call
counts. Accepted target refreshes are context boundaries and are not reported
as same-objective descent.

If all target predicates are satisfied after a committed carving step or
accepted deformation, the controller returns the corresponding stable status.
Otherwise, deformation gradient convergence transfers control to the next
all-filled carving block. The controller has distinct outcomes for no
deformation degrees of freedom, unstable `mu` exhaustion, no-progress stall,
per-block carving/deformation guards, and the outer block guard. These are
explicit operational statuses, not proofs of nonlinear convergence or global
optimality.

## Paper 1% fixed-alpha profile versus the beta

The paper-governed line-search constants remain:

```text
initial step             1.0
step reduction ratio     0.8
minimum step              0.4
required decrease          1%
mu decrement              0.05
```

For fixed `alpha`, current target context, and current `mu`, a proposed
`q_trial=q-eta*gradient` is accepted only when

```text
E(q_trial; alpha, targets, mu) <= 0.99 * E(q; alpha, targets, mu).
```

The exact 1% boundary is accepted. A sub-1% decrease, increase, or physically
infeasible state rolls back `q` and geometry. Failed steps follow the
deterministic sequence `1,.8,.64,.512,.4096,.4`; after the failed minimum-step
trial, `mu` decreases by `0.05` and `eta` resets to `1`. A `mu` transition
reweights the objective, so its before/after values are explicitly marked
`objectiveReweighted` and are not compared as a descent pair.

The inspected beta follows a different outer profile: it may refill and
recarve within a deformation trial and uses a different decrease threshold for
schedule changes. R7 follows the printed paper's block-coordinate profile.
The discriminator regression observes zero refill and zero carving calls in a
paper-style fixed-`alpha` trial, while the local beta-profile oracle records
one of each. This difference is intentional and documented rather than hidden
as implementation equivalence.

## Production-code integration fixture

The integration fixture is a closed overhanging prism at voxel resolution 18,
with a fixed bottom support and one unlocked vertically blended user handle.
It is deliberately asymmetric so that carving, deformation, target refresh,
and a second distinct carving block are all exercised. It is a local
production-code fixture, not an official paper model and not a manufactured
artifact.

The observed Carving--Deformation--Carving trace is:

```text
pattern A carved voxels              304
pattern B reset/refilled voxels      304
A-to-B occupancy Hamming distance     48
fixed-alpha deformation trials        27
```

Thus block B demonstrably starts all-filled rather than inheriting A, yet
recarves to a distinct occupancy after deformation. All 27 deformation trial
records report zero refill calls, zero carving calls, and unchanged `alpha`.
The companion transaction run observes at least one accepted and one rejected
trial, a committed target refresh, fixed support nodes, a non-terminal
`FeasibleButUnsafe` state, and all-of behavior for a mixed two-constraint
schedule. Repeating the end-to-end run from identical inputs produces a
bit-for-bit equal trace and final value snapshot under the tested comparison.

For context, the three-block end-to-end fixture reports 27 deformation trials,
final density-one mass `0.163856`, and final CoM approximately
`(0.218563, 0.426800, 0.500315)`. These values characterize this fixture only;
they are not calibrated physical measurements or paper-table results.

## Registered R7 test matrix

The values below are from the current RelWithDebInfo focused run on
2026-08-12. The R7 label contains exactly sixteen registered tests.

| CTest | Contract | Current result |
|---|---|---|
| `deformed-voxel-mass-identity-similarity` | Exact identity-cell mass/centroid and translation, rotation, and uniform-scale covariance of deformed cell mass properties. | Passed; mass `0.1472`, scaled mass `0.362167`. |
| `deformed-voxel-mass-cavity-sum` | Non-affinely distorted cell centroids differ from corner averages; three-cell mass, first moment, and CoM agree with an independent cavity-surface integral. | Passed; mass `0.134706`, centroid delta `8.18934e-4`. |
| `deformed-voxel-mass-invalid` | Reject non-finite, inverted, collapsed, and out-of-range deformed-cell inputs without partial state. | Passed. |
| `deformed-carving-identity` | The deformed-mass solver at identity reproduces the legacy axis-aligned committed prefix, fill state, plane, mass, and CoM. | Passed; 18 voxels carved. |
| `deformed-carving-mass-field` | Candidate ordering, prefix mass subtraction, and full recomputation use exact distorted-cell mass properties rather than rest boxes. | Passed; 18 voxels carved, remaining mass `4.34603`. |
| `deformed-carving-transaction` | A pre-commit exception rolls back staged fill flags, mass, first moment, CoM, counters, diagnostics, and output arguments. | Passed. |
| `deformed-carving-objective-update` | Between-step target replacement preserves committed material state and history, rejects empty/reentrant updates, preserves terminal state, and drives the next step with the new target. | Passed; refreshed next prefix length 18. |
| `alternating-target-refresh` | Standing target refresh uses fixed support, changes generation only when objective coordinates change, and distinguishes context changes from state-only changes. | Passed; two objective generations. |
| `alternating-stability-terminal` | Closed standing/suspension boundaries, unattainable margin, upstream and zero-offset suspension cases, and exact terminal predicates. | Passed. |
| `alternating-multi-target-terminal` | Ordered standing/suspension mixtures, summed objective energy, and all-of completion for partial and complete multi-target states. | Passed; partial energy `0.03125`, three ordered constraints exercised. |
| `alternating-refill-reset` | `refillAll()` is counted and idempotent, preserves immutable metadata, and permits pattern B to differ from pattern A. | Passed; analytic patterns of 3 and 2 empty cells. |
| `alternating-fixed-alpha-transaction` | Exact 1% acceptance, sub-1% rollback, and immutable occupancy during the paper line-search trial. | Passed. |
| `alternating-paper-beta-profile` | Paper deformation trials perform no refill/recarve; a separate beta-profile oracle discriminates the inspected behavior. | Passed; paper `0/0`, beta oracle `1/1` refill/carve calls. |
| `alternating-controller-refill-block` | Production controller executes Carving--Deformation--Carving with reset counts, hashes, epochs, fixed-alpha trace, and a distinct second occupancy. | Passed; A `304`, refill `304`, Hamming `48`, deformation trials `27`. |
| `alternating-controller-transaction` | Stable early exit, unsafe non-terminal state, accepted commit, rejected rollback, target refresh, fixed support, and multi-target all-of behavior. | Passed; accepted and rejected paths both observed. |
| `alternating-controller-end-to-end` | Multi-block synchronized final geometry/mass/CoM, surface oracle, support invariance, distinct patterns, and deterministic repeated runs. | Passed; three blocks, 27 deformation trials. |

The production fixture rows are `VERIFIED_REGRESSION` evidence for this local
geometry and configuration. They do not promote the other rows beyond their
unit/semantic contracts.

## Numerical and claim boundaries

The following claims are permitted after the R7 gates:

- the tested deformed-cell sidecar integrates positive mass and first moment
  from synchronized hexahedral corners and rejects inverted/collapsed states;
- deformed carving uses those cell properties and can refresh targets safely
  between committed steps;
- every tested carving block starts from all-filled occupancy, while every
  deformation trial retains exactly fixed `alpha`;
- same-context carving and accepted deformation records obey their tested
  non-increase/1% rules, with target and `mu` changes identified as objective
  context boundaries;
- the fixed support and ordered target specifications remain unchanged across
  the tested alternating run; and
- the production controller returns a synchronized, auditable state on the
  local overhanging-prism fixture.

These are `VERIFIED_UNIT` claims with one local integration fixture at
`VERIFIED_REGRESSION`. The following claims are not permitted:

- `VERIFIED_PAPER`, official-beta equivalence, or matching a paper asset,
  figure, table, trajectory, or final metric;
- convergence, global optimality, or successful stability for arbitrary
  geometry, initialization, support, or target combinations;
- preservation of model-space `t_min` after deformation;
- universal absence of floating components, self-intersection, inverted final
  faces, or topological defects;
- final watertightness, slicer compatibility, structural strength,
  printability, or fabricated standing success; or
- R10 experimental, robustness, resolution, timing, or official-asset claims.

## Final local verification gate

Observed on 2026-08-12:

```text
RelWithDebInfo headless complete CTest                 55/55 passed
R7 focused CTest                                      16/16 passed
R7 semantic/code audit                 0 blocker, 0 major findings
Independent Debug complete CTest                      55/55 passed
RelWithDebInfo viewer configure/build                 passed
Recursive clean-clone headless/viewer gate                    passed
```

The R7 tests were also repeated with CTest `until-fail:20`; all sixteen
completed twenty consecutive passes. The implementation is commit `eaf2c99`;
the conformance documents and README are commit `f638ccd`. A recursive,
no-hardlink clone of `f638ccd` independently configured and built both the
headless and viewer presets, passed all 55 headless tests, and ended with an
empty superproject status. Its libigl submodule was clean at
`f06fb13014aadff7eb1d7423ac018620441a3657`, and the historical R0 trace kept
SHA-256 `3f9108649fc329a136a5371597335a99275021b6ddd9a3984198ce2b8e8c1856`.
Even a completely passing repository gate cannot promote R7 above the
numerical and claim boundaries stated here.
