# R4 paper-governed fixed-geometry carving acceptance record

Specification audit date: 2026-08-12 (Asia/Shanghai)

Acceptance status: **accepted at the documented R4 status ceiling**

## Authority and scope

The normative source is Prévost et al., *Make It Stand: Balancing Shapes for
3D Fabrication*, Section 4 and Equation (8). The official beta archive remains
a read-only behavioral oracle. Where it differs from the paper, the paper
governs.

R4 implements the fixed-outer-surface inner-carving stage:

- one analytic cutting plane for every target/gravity pair;
- the paper multi-target score, summing independently positive plane values;
- deterministic descending score order with lexicographic grid-key ties;
- simulation of prefix zero and every valid ordered deletion prefix;
- selection and transactional commit of the lowest computed true Equation
  (5) CoM-energy prefix, including after temporary increases;
- repeated plane reconstruction from the new CoM;
- explicit no-candidate and too-few-removals termination;
- model-space `t_min` eligibility under a conservative voxel-AABB clearance
  bound for the fixed outer mesh; and
- six-neighbor component reporting for remaining material and carved void.

R4 does not implement BBW/LBS deformation, Laplacian preservation, gradients,
line search, target refresh, voxel refill between deformation alternations,
final-mesh manufacturing checks, or paper result reproduction.

## Cutting planes and candidate scores

For objective `k`, let

```text
P_k = I - ghat_k ghat_k^T
n_k = P_k (c_0 - c*_k)
d_ik = (r_i - c*_k) dot n_k
```

The local implementation follows the printed Equation (8): `n_k` is not
normalized. This distinction does not change a one-target order, but it can
change relative weighting for multiple targets. The inspected beta normalizes
each residual; this is recorded as a paper/beta discrepancy rather than copied
silently.

An already-satisfied target has `n_k=0` and contributes neither a plane nor a
score. A candidate score is

```text
s_i = sum_k max(0, d_ik).
```

A voxel enters the queue when `s_i>0`; unlike the beta two-target branch, not
every per-target distance must be positive. Candidates are ordered by
descending computed `s_i`, then ascending `(x,y,z)`.

## Complete true-energy prefix scan

Let the current mass properties be `(m_0,Q_0)` and the ordered removable boxes
be `v_1..v_N`. R4 evaluates prefix zero and every prefix whose remaining mass
is positive and finite:

```text
m_j = m_0 - sum_(l<=j) m_l
Q_j = Q_0 - sum_(l<=j) Q_l
c_j = Q_j / m_j
E_j = 0.5 * sum_k ||P_k(c_j-c*_k)||^2.
```

An energy lower by the documented scale-aware floating-point tolerance wins;
ties retain the shorter prefix.
Simulation never mutates the grid. Only after the complete scan succeeds is
the selected prefix committed. The inherited `MAX_CARVE=1000` batch cap and
finite-plane `signed_distance` path have been removed.

The `global-prefix-transient` fixture has prefix energies

```text
0.5 -> 12.625 -> 0.0
```

and therefore proves that an intermediate increase does not stop the scan.
`global-prefix-middle` proves that the optimum need not be the full queue;
`global-prefix-zero` proves that a non-improving queue performs no mutation;
and `global-prefix-over-1000` proves that 1001 deletions remain eligible.

## Stopping semantics

The paper says to stop when “too few voxels were removed” but gives no number;
the beta does not recover a missing threshold. R4 therefore exposes
`minimumRemovalsToContinue`, requires it to be at least one, and defaults to
one. A positive improving prefix smaller than a user-specified threshold is
committed and then terminates. A zero best prefix terminates without mutation.

The 100-iteration limit is an explicit local run guard, not a paper convergence
constant. Its terminal token remains distinct from the paper-style progress
stop. R4 deliberately does not use an inherited one-voxel CoM tolerance:
standing feasibility belongs to the support/safe-region objective, while a
fixed scalar tolerance is not specified by Section 4. The old public enum
tokens for target tolerance and maximum empty ratio remain source-compatible
but are no longer emitted by the R4 solver.

## Wall-thickness contract

For signed distance `s<0` at an axis-aligned voxel center and box half-diagonal
`h`, the signed-distance function's 1-Lipschitz property gives the conservative
clearance lower bound

```text
clearance_lower_bound = -s - h.
```

The complete voxel is carve-eligible only when this bound is at least the
finite nonnegative model-space length `t_min`. A physical-unit interpretation
requires the input mesh coordinates to be calibrated to that unit. This
protects the fixed outer surface, not the later deformed result. The analytic
unit-box test evaluates the same model-space `t_min` at longest-axis
resolutions 24 and 48 and checks every eligible voxel with an independent
box-to-surface clearance oracle.

## Topology boundary

The paper states that ordered plane carving avoids floating components, but it
does not define voxel adjacency or prove a general theorem for arbitrary
nonconvex inputs. R4 consequently reports six-neighbor material and carved-void
component counts after every committed step. The component oracle is tested on
an analytic connected box and carved interior voxels. It is diagnostic evidence,
not a universal connectivity guarantee, and the number of void components is
not required to be one.

## Verification matrix

The R4-focused CTest cases are:

```text
cutting-plane-analytic
multi-plane-score-analytic
carving-order-deterministic
global-prefix-transient
global-prefix-middle
global-prefix-zero
global-prefix-shortest-tie
global-prefix-scale-invariant
global-prefix-over-1000
minimum-wall-thickness-analytic
voxel-topology-report
```

Production integration is additionally covered by
`multi-target-solver-boundary`, `carving-one-step`, and
`carving-full-mass-consistency`. The latter independently recomputes mass,
first moment and CoM from voxel state after every step, checks accepted energy
is non-increasing, round-trips the 44-column CSV, and repeats the full run to
verify deterministic state and diagnostics.

The step CSV is an audit trace, not a self-contained experiment manifest: the
objective collection, grid resolution, `t_min`, and progress threshold must be
recorded by the calling experiment configuration. Reproducible paper-case
configuration and result artifacts remain R10 work.

## Status ceiling

Passing R4 can advance paper Section 4 cutting-plane, ordering, true-energy
prefix, repeated planning and arbitrary-`K` carving semantics to
`VERIFIED_UNIT`/`VERIFIED_REGRESSION`. The model-space `t_min` and topology rows
are verified only at their documented voxel-fixture scope. R4 cannot establish
`VERIFIED_PAPER`, successful fabricated balance, final printability, deformation
quality, or complete reproduction of *Make It Stand*.

## Final verification gate

The pinned R0 toolchain completed the following local gates on 2026-08-12:

```text
RelWithDebInfo headless CTest                         25/25 passed
R4-labeled tests, each repeated until-fail:10        14/14 passed
Debug headless CTest                                 25/25 passed
RelWithDebInfo viewer build                          passed
git diff --check                                     passed
Pinned libigl status                                 clean
Historical R0 trace SHA-256                          unchanged
```

The committed tree at `989f82fe3332d4e718bef099bbdb0bf5923d3cd5` was then
checked through a recursive clean clone. Its headless suite passed 25/25, its
viewer target built, the parent and pinned submodule remained clean, and the
historical trace hash remained unchanged.
