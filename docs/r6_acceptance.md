# R6 fixed-alpha deformation-optimization acceptance record

Specification audit date: 2026-08-12 (Asia/Shanghai)

Acceptance status: **accepted at the documented R6 status ceiling**

The five focused R6 tests pass on the current working tree. This accepts the
implicit-similarity shape energy, analytic differential chain, and fixed-alpha
optimization machinery at `VERIFIED_UNIT`, with one local end-to-end fixture at
`VERIFIED_REGRESSION`. It is not a claim of full-paper reproduction, successful
carving/deformation alternation, balance of a paper asset, or fabrication
readiness. In particular, no R6 result may be labelled `VERIFIED_PAPER`.

## Authority and phase boundary

The normative source is Prévost et al., *Make It Stand: Balancing Shapes for
3D Fabrication*, especially Section 5, Equations (11)--(14), and the Appendix
mass derivatives. The official beta archive remains a read-only behavioral
oracle and is not copied or linked into this repository. Where the paper and
the inspected beta differ, the paper governs this clean-room implementation.

R6 freezes the voxel occupancy `alpha`, cavity connectivity, support, target
set, BBW field, fixed handle rotations, rest geometry, and mesh-specific
`lambda`. It varies only the R5 free handle vector

```text
q = [t_x, t_y, t_z, s] for every unlocked user handle.
```

Within that frozen state R6 implements:

- the paper's implicit-similarity Laplacian editing energy on `M_O`;
- analytic signed surface mass, first-moment, and center-of-mass derivatives
  for the outer boundary and inward-oriented cavity boundary;
- the complete multi-target CoM and shape-energy gradient with respect to `q`;
- the fixed-alpha total objective from Equation (11);
- a transactional gradient trial with the paper's 1% decrease rule; and
- the paper step-size and `mu` schedule as an explicit deterministic state
  machine.

R6 deliberately does not implement the outer Section 3.2/Figure 3 alternation,
carving after deformation, voxel refill, target refresh, or the paper's
"stability achieved" stopping predicate. Those are R7 responsibilities. R8
remains responsible for deformed minimum-wall checks, final mesh construction,
watertightness/connectivity checks, slicing, and manufacturing validation.

## Implicit-similarity Laplacian energy

`LaplacianEditingEnergy` is not the ordinary graph-smoothing penalty
`||L*x||^2`. For every outer-surface vertex `i`, it forms the local set
`{i} union N_i` and the uniform umbrella differential operator `L_i`. The
rest positions in that neighborhood define a seven-parameter linear
similarity design `A_i`: one uniform scale, three infinitesimal rotations, and
three translations. Applying the same scale/rotation model to the rest
differential gives `K_i`. The local residual is

```text
U_i = K_i * pinv(A_i) - L_i.
```

Stacking the `U_i` blocks gives a sparse residual operator `M`. R6 stores the
unscaled energy and gradient as

```text
E_shape_unscaled(x) = 0.5 * ||M*x||^2,
grad_x E_shape_unscaled = M^T*M*x,
Q = M^T*M.
```

The paper's mesh-specific scale is applied outside this component:

```text
E_shape(x) = lambda * E_shape_unscaled(x).
```

This is the local implicit-rotation construction underlying the paper's
positive-semidefinite `M_Lap`, with `Q` corresponding to the quadratic matrix
in Equation (12). Translation and global uniform scaling lie in the tested
nullspace. As in the paper's linear variational model, this is only a local
linear treatment of rotation; R6 does not claim exact invariance under
arbitrary finite rotations.

Rest positions are centered before constructing each pseudoinverse as a
conditioning improvement. Translation is already represented in `A_i`, so
centering leaves the intended residual unchanged up to roundoff. A full-rank
SVD pseudoinverse is used instead of an unchecked normal-equation inverse;
closed-manifold, neighborhood-valence, rank, and finite-coefficient contracts
are enforced. These are numerical hardening choices, not a new energy.

## Analytic surface mass and CoM differential

`SurfaceMassDifferential` follows the same signed-boundary convention as R2:
the outer surface is outward and positive, while the carved cavity is inward
and negative. For a triangle with coordinates `a`, `b`, and `c` relative to a
reference point `r`, its signed tetrahedral contribution is

```text
m_f = (rho / 6) * a dot (b cross c),
u_f = r + (a + b + c) / 4,
Q_f = m_f * u_f.
```

The analytic vertex derivatives are

```text
d m_f / d a = (rho / 6) * (b cross c),
d m_f / d b = (rho / 6) * (c cross a),
d m_f / d c = (rho / 6) * (a cross b),

d Q_f / d v = u_f * (d m_f / d v)^T + (m_f / 4) * I.
```

After summing outer and inner contributions, `m` must remain finite and
strictly positive. With first moment `Q` and center of mass `c=Q/m`, the
quotient rule gives

```text
d c / d v = (d Q / d v - c * (d m / d v)^T) / m.
```

The implementation checks the expected outer/inner orientation and rejects a
non-positive material mass as an infeasible deformation state. It does not
silently continue with an absolute volume or a NaN center of mass.

## Free-parameter and total-energy gradient

R5 Equation (10) becomes an affine map from only the unlocked parameters:

```text
x_O(q) = b_O + J_O*q,
x_I(q) = b_I + J_I*q.
```

Locked handles and fixed rotations are absorbed into the offsets and
Jacobians. The CoM chain rule is therefore

```text
d c / d q = (d c / d x_O) J_O + (d c / d x_I) J_I.
```

For target `k`, let `P_k = I - ghat_k*ghat_k^T` be projection perpendicular
to its normalized gravity. R6 sums the independently projected target terms,
as required by the paper's multi-base extension:

```text
E_CoM(c) = sum_k 0.5 * ||P_k (c - c*_k)||^2,
grad_c E_CoM = sum_k P_k (c - c*_k).
```

For fixed occupancy and tradeoff `mu`, the evaluated objective is

```text
E(q;mu) = (1-mu) E_CoM(c(q)) + mu E_shape(x_O(q)),

grad_q E =
    (1-mu) (d c / d q)^T grad_c E_CoM
  + mu lambda J_O^T M^T M x_O.
```

This implements the analytic chain in Equations (11)--(14). The objective
returns the total energy, both component energies, signed material mass, CoM,
and analytic free-parameter gradient from one evaluation. Invalid parameter
dimensions and internal numerical failures propagate as errors. Only expected
finite-but-physically-infeasible trials, such as a non-positive scale or
non-positive material mass, are eligible for transactional rejection.

## Compact cavity memory design

The full voxel lattice may contain `(R+1)^3` corners, while the inner mass
integral depends only on corners referenced by the current cavity faces.
Building and retaining a free-parameter Jacobian for every global corner would
make the fixed-alpha objective unnecessarily large, especially at production
resolutions.

`DeformationObjective` therefore extracts the current global cavity faces,
sorts and deduplicates their referenced global corner indices, gathers only
those rest positions and BBW rows, and deterministically remaps the faces into
a compact local vertex array. `FreeParameterLbsModel` then stores only the
affine offset and free-parameter Jacobian for this compact boundary. Empty
cavities are represented by empty inner arrays and zero inner contribution.

This compaction changes storage, not geometry or integration semantics. The
objective-gradient fixture independently evaluates R5's complete synchronized
geometry and requires the compact objective's mass and CoM to match it within
`1e-11`. R6 does not claim that the complete R5 grid snapshot itself has
production-optimal memory use; it only avoids duplicating that full-grid
Jacobian inside the objective.

## Paper 1% transactional line search

The default state-machine constants follow Section 5.3:

```text
initial step             1.0
step reduction ratio    0.8
minimum step            0.4
required decrease       1%
mu decrement            0.05
```

For current parameters `q`, gradient `d`, step `eta`, and fixed current `mu`,
the trial is `q_trial=q-eta*d`. It is committed only when

```text
E(q_trial;mu) <= 0.99 * E(q;mu).
```

The exact 1% boundary is accepted. Every insufficient decrease, including a
smaller positive decrease, rolls back `q` before changing the state. A
physically infeasible trial also rolls back; unrelated internal evaluation
errors are not hidden as ordinary line-search failures.

The paper gives a ratio and a minimum but does not specify how floating-point
steps land on that minimum. R6 uses the deterministic sequence

```text
1.0, 0.8, 0.64, 0.512, 0.4096, 0.4.
```

After a failed `0.4` trial, it preserves `q`, reduces `mu` by `0.05`, and
resets the step to `1.0`. Floating-point endpoint residue is normalized so
that exhaustion terminates at exactly `mu=0`. The inspected beta uses a 3%
threshold; the paper's printed 1% rule is normative here, and the test
explicitly distinguishes the two by accepting a 2% decrease.

### `mu` reweight semantics

Changing `mu` changes the objective rather than the geometry. On a
`MuReduced` or `MuExhausted` transition, `energyBefore` is evaluated at the old
`mu`, while `energyAfter` is a fresh evaluation of the unchanged `q` at the
new `mu`. Such values are not a valid descent pair and must not be compared as
if one optimization step had reduced a single fixed objective.

The transition and every optimizer trace record therefore set
`objectiveReweighted=true`. Accepted trials and step reductions keep it false.
At fixed `mu`, accepted trials enforce the 1% decrease and rejected trials
preserve both parameters and energy. End-to-end evaluation also compares the
final state against the initial parameters under the same final `mu`.

## Registered R6 test matrix

The numerical values below are from the current RelWithDebInfo focused run on
2026-08-12. Gradient errors use centered finite differences with step `1e-6`
and a required relative error no greater than `1e-5`.

| CTest | Contract | Current result |
|---|---|---|
| `laplacian-editing-analytic` | Closed box; residual/operator dimensions; `Q=M^T M` symmetry, PSD quadratic identity, rest/translation/uniform-scale nullspace, positive shear energy, and analytic vertex gradient. | Passed; shear energy `0.162273`, gradient relative error `1.16168e-10`. |
| `surface-mass-gradient-analytic` | Outward outer box plus inward cavity; analytic signed mass, first moment, and outer/inner CoM Jacobians against independent surface integration and finite differences. | Passed; mass error `3.78574e-10`, moment error `5.39107e-10`, outer/inner CoM errors `6.37963e-10` / `5.20851e-09`. |
| `deformation-objective-gradient-analytic` | Two gravity targets, compact cavity, affine free-handle LBS, complete Equation (11)--(14) energy/gradient, and equality with the full synchronized-geometry mass/CoM. | Passed; total energy `0.113555`, gradient relative error `2.50748e-10`. |
| `paper-line-search-state-machine` | Exact 1% acceptance, 0.9% rollback, 2% paper-versus-beta discriminator, infeasible-trial rollback, internal-error propagation, exact step sequence, `mu` decrement/reset and exact exhaustion. | Passed; six failed step attempts before reweight, final `mu=0`. |
| `deformation-optimizer-end-to-end` | Fixed-alpha optimization trace, fixed support, unchanged voxel fill, positive scale/mass, fresh final reevaluation, per-record reweight semantics and final same-`mu` improvement. | Passed; 51 iterations, initial energy at `mu=0.20` `0.0361284`, initial-parameter energy at final `mu=0` `0.0451605`, final energy `4.84989e-05`, final `mu=0`. |

The end-to-end case is a small analytic box fixture. Its
`VERIFIED_REGRESSION` label means that the R6 components work together on this
local fixed-alpha case; it does not constitute a paper-asset or fabrication
regression.

## Numerical and claim boundaries

The following claims are permitted after the R6 gates:

- the local implicit-similarity Laplacian operator has the tested quadratic,
  nullspace, and gradient properties;
- signed surface mass, first moment, CoM, and complete free-parameter gradients
  agree with finite differences on the analytic fixtures;
- fixed-alpha total-energy trials are transactional and follow the paper's 1%
  schedule rather than the beta's different threshold;
- `mu` transitions are explicitly identified as objective reweighting; and
- the compact cavity evaluation agrees with the complete synchronized geometry
  on the tested fixture.

These are `VERIFIED_UNIT` claims, with only the small fixed-alpha integration
case at `VERIFIED_REGRESSION`. The following claims are not permitted:

- `VERIFIED_PAPER`, equality to the official beta, or matching a paper figure;
- convergence or global optimality of the nonlinear deformation problem;
- completion of the paper's alternating carving/deformation algorithm;
- successful target refresh, refill, or post-deformation grid update;
- achievement of the paper's stability predicate on a model;
- preservation of `t_min` under deformation; or
- final watertightness, connectivity, structural strength, printability, or
  fabricated standing success.

## Final verification gate

The focused working-tree run observed on 2026-08-12 is:

```text
R6 focused tests                                      5/5 passed
Maximum observed focused gradient relative error     5.20851e-09
Paper acceptance boundary                            exactly 1%
Deterministic failed-step sequence                    6/6 passed
Fixed-alpha end-to-end trace                          passed
```

Before the acceptance commit is treated as a release gate, the parent R6
workflow must additionally record a complete RelWithDebInfo suite, an
independent Debug suite, repeated focused tests, viewer build, `git diff
--check`, pinned-submodule cleanliness, unchanged historical R0 trace, and a
recursive clean-clone build/test. Passing those repository gates still cannot
promote R6 above its documented unit/local-regression ceiling.
