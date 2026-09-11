# Source-cage BBW: independent verification contract

Task5k adds an actual weight solve to the previously verified cage/grid/FEM
chain. The checks described here concern synthetic numerical fixtures, not
the original Spin-It models, joint deformation, physical material, or a full
paper-reproduction release. The global status remains `PARTIAL_BLOCKED`.

## Two optimization protocols, not one interchangeable result

With the prescribed boundary rows eliminated, the computational FEM objective
is `0.5 * sum_h w_h^T K^T M^-1 K w_h`, with `0 <= w <= 1`.

- `COUPLED_PARTITION` includes `sum_h w_ih = 1` at every free node in the solve.
  It must not normalize or clamp the result afterward.
- `INDEPENDENT_NORMALIZED` first solves the separate bounded problem for each
  handle, then explicitly normalizes each row. Raw optimality and final
  feasibility are different checks. The final normalized field must be
  separately assessed against the coupled problem; its optimality cannot be
  inherited from the raw solves.

Both alternatives are described by the [original BBW paper, Sections 3.1 and
3.4](https://igl.ethz.ch/projects/bbw/bounded-biharmonic-weights-siggraph-2011-compressed-jacobson-et-al.pdf).
Our existing pinned HiGHS backend is an independent implementation choice,
not the original MOSEK program or evidence of an identical optimizer path.

## Independent reference construction

The controller reference imports NumPy/SciPy, not project bindings or native
FEM routines. It constructs the stiffness matrix from a literal rational
brick stencil and the complete boundary rows from explicit piecewise-linear
triangle formulas. It integrates both restricted-grid Voronoi and
barycentric-tetrahedral lumping conventions separately.

The source cube `[0.3, 0.4]^3`, with cage level 2 and fine level 3, yields the
represented cage domain `[0.25, 0.5]^3`: 8 bricks, 27 nodes, 8 handles and one
free node. The exact free-node weight is `5/16` for each diagonal-endpoint
handle and `1/16` for each of the other six handles. The full objective is
`638/3` for Voronoi lumping and `659/3` for barycentric lumping.

The second fixed source cube `[0.3, 0.7]^3` yields the represented domain
`[0.25, 0.75]^3`: 64 bricks, 125 nodes, 26 handles and 27 free nodes. All 702
coupled variables are retained. These cubes are analytical regression inputs,
never stand-ins for original-paper animals or people.

The reference starts from uniform weights and uses SLSQP. On the medium
fixture, SLSQP alone left normalized stationarity around `2e-8`–`3e-8`.
That does **not** relax the native `1e-8` gate. A separately recorded polishing
step solves the original reduced KKT system on the active set found by SLSQP,
then checks the original bounds and multiplier signs. Pre-polish traces and
results are retained. The polished reference's stationarity is below `2e-15`.

The rectangular medium fixture exposes a useful distinction: the independent
raw weights have a row-sum residual about `0.09862`. Normalization changes a
weight by about `0.03238`, but the normalized field's coupled projected
stationarity remains about `0.01465`. Partition of unity alone therefore does
not certify coupled optimality.

Both mass conventions happen to produce the same weights on these regular
rectangular fixtures, with different boundary-constant energies. This is not
a general claim about other geometries or mass conventions.

## Native acceptance and remaining boundaries

The separate native probe exports actual weights, complete handle identities,
backend call counts, native assessments and scientific hashes. The controller
checker reconstructs `Q` independently, matches every represented node and
handle, and recomputes boundary, bound, partition, stationarity,
complementarity and projected-gradient residuals. Native success flags are
not the numerical oracle. Scientific repeats exclude wall-clock observations.

The fixture reference and checker are retained under
`.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5k-*`.
The final-source probe has executed all eight configurations twice. The
independent checker measured maximum raw stationarity `2.98e-16` and maximum
weight difference from the independent reference `1.42e-15`; each pair's
scientific outputs matched exactly. Seven malformed-result mutations were
rejected. Full paths, source/binary bindings and limitations are recorded in
`task-5k-controller-verification.md`. The frozen build passed 203/203 registered
tests and Faster-OFF 4/4; fresh independent code review approved this bounded
increment with no Critical/Important findings. Concurrent external Faster
source edits are excluded: these results do not certify the whole current
working tree, nor a clean release.

During implementation, independent checks caught the pinned HiGHS bounded-only
startup reporting `Optimal` despite a residual about `0.005138`. The adapter
now represents every independent problem with a redundant `0=0` row from the
outset, selecting its general startup without altering objective or feasible
set. This is explicit, hashed and not an automatic retry or dependency update.

One non-dyadic transform (scale 1.7 and origin [-0.3,-0.2,-0.1]) exposed an
upstream FEM incomplete-boundary error before BBW. It remains a documented
open limitation, not a proof of invalid geometry or physical infeasibility.
Positive dyadic covariance tests do not resolve that separate failure.

Complete nodal boundary constraints do not make off-node trilinear transfer
identical to the coarse cage's piecewise-triangular boundary function.
Computational FEM masses are not physical material masses. Actual BBW weights
alone do not complete joint handle/filling optimization, its shape derivatives
and safeguards, original-case evidence, sealed material or independent physics.

## Task5l follow-up (2026-09-11)

The specific non-dyadic FEM boundary error above is now fixed in `dd5428c`.
Direct opposite-edge areas avoid the `1-w1-w2` cancellation without widening
tolerances or altering points. The original1.7-scale input now passes all four
BBW combinations and deterministic repeats. Independent exact-rational boundary
checks and scoped review passed; see [boundary-fix evidence](spin2014_fem_boundary_fix.md).
This supersedes that particular open blocker, not the remaining joint-shape,
original-case, physics or whole-current-tree verification gaps.
