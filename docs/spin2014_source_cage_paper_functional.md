# Actual-source joint paper functional

`SourceCagePaperFunctional` owns an immutable `SourceCageMomentFunctional`
snapshot. It evaluates the actual source/full-BBW deformed material, transforms
its density-weighted raw moments into an explicit fixed spin frame, and calls
the existing `Spin2014GivensProblem`. It accepts no replacement mass, weights,
or Jacobian. This bounded integration is **SYNTHETIC_REGRESSION** evidence, not
an original-paper success or completion of the shape optimizer.

## Coordinates and frame

For B actual source partition leaves and H original cage handles, the state is
`[beta B, optional phi, handle 4H]`. Each handle keeps its upstream order
`[scale, tx, ty, tz]`, identity fixed rotation, and strictly positive scale.
Static lowering has N=B+4H and no angle coordinate; the other four objectives
have N=B+1+4H. `angleOffset()` returns `size_t::max()` for static, B otherwise.
`handleOffset()` is B or B+1. `initialState(beta, phi)` uses actual initial handle
parameters. Static split-state calls reject every nonzero phi, including when
constructing an initial state. All modes reject nonfinite phi, invalid beta,
nonfinite coordinates, invalid scales, and wrong dimensions without clamping.

`SourceCagePaperFrame` specifies world origin o and proper orthonormal rotation R:
`q = R (p-o)`. +z is the chosen axis, and o is the selected support/axis point.
Defaults o=0 and R=I are **IMPLEMENTATION_DEFAULT**, not recovered paper poses.
The frame never follows the current center of mass. R is checked for finite
entries, `R R^T = I` and determinant +1 with absolute tolerance 1e-10; scaling,
shear and reflection fail. Phi is the independent transverse Givens angle.

The exact constant 10x10 map uses native order
`[M,LX,LY,LZ,QXX,QYY,QZZ,QYZ,QXZ,QXY]`:

```
M' = M
L' = R (L-M o)
Q' = R (Q-o L^T-L o^T+M o o^T) R^T
```

Mixed tensor entries are stored once, not doubled polynomial coefficients.
There is no additional density factor. The same map multiplies both actual
upstream Jacobians. Its exposed accessor permits auditing against independently
transformed points and signed tetrahedral integration.

## Objectives, constraints and fill energy

The existing Givens evaluator preserves all five objective distinctions and its
inertia/center weight semantics. Raw equality ordering is unchanged:

| Objective mode | Equalities |
| --- | --- |
| Static | LX, LY |
| Top (all three variants) | LX, LY, QXZ, QYZ, transverse diagonalization |
| Yo-yo | LX, LY, LZ, QXZ, QYZ, transverse diagonalization |

Constraints have no result-dependent normalization. All moment and direct phi
derivatives are chained to every joint coordinate. The spin raw Jacobian has
an exactly zero phi column in nonstatic modes. No full joint Hessian is exposed:
source moment second derivatives are unavailable; later solver work uses L-BFGS.

The independent `fillWeight` gammaL contributes
`0.5 gammaL beta^T L beta`, with gradient `gammaL L beta`. L is the immutable
actual source partition's sparse unit face-adjacency Laplacian. Its equivalent
edge energy is `0.5 gammaL sum_edges (beta_i-beta_j)^2`. Handle and angle
derivatives are zero. This discretization is **IMPLEMENTATION_DEFAULT**, not
the full recovered paper Laplacian. No shape/displacement energy is implied.

Partition changes or BBW recomputation require construction of a new functional.
Adaptive split/rebind, joint SLQP, shape regularization, inversion/intersection,
overlap, thickness and geometry commit gates are outside this component. A finite
objective or positive mass does not establish feasibility or physical success.

## Results, identity and resources

Results retain the actual world evaluator output including deformed points,
spin moments, paper objective, fill penalty, total objective, equality vector,
joint objective gradient, joint constraint Jacobian and spin raw Jacobian.
`evaluateValues` returns the same geometry, moments, objective and constraints
with absent global derivative arrays, including absent upstream Jacobians.
The reused eleven-coordinate Givens implementation still uses a fixed-size
local derivative workspace, independent of B, H and N.

The default positive `maximumJointDerivativeValues` is 10,000,000 doubles.
The derivative path checks the combined `(10+K+1) N + B` capacity before
allocation: retained spin raw Jacobian, constraint Jacobian, objective gradient,
and B-value sparse multiply scratch. Index, byte, sum and product arithmetic
are checked. Upstream moment budgets remain independently enforced. No P×H,
P×4H, or dense B×B array is created. The value-only path does not require the
joint derivative budget and works with a cap too small for derivatives.

Invalid configuration/state throws `invalid_argument`; invalid capacities or
allocation failures throw `length_error` containing `RESOURCE_LIMIT`;
nonfinite frame-transform/fill/chain arithmetic throws `overflow_error`.
Existing upstream/Givens exception categories remain unchanged. Nonfinite
results are never returned.

Deterministic SHA256 configuration identity binds source configuration and
upstream identity, objective, all independent weights, frame, budget and layout.
Each output identity additionally binds the actual upstream evaluator output,
full beta/phi/handle state, value/derivative mode, every scalar, equality and
joint derivative result. Upstream output binds the actual moved geometry.
No timestamp or success boolean enters the digest.

## Verification scope

The native regression uses both actual Task5q fixtures: baseline and translated,
scaled nonsymmetric source/design, source depth 5, cage level 2, fine grid 3,
eight beta coordinates and all 26 handles. Nonuniform beta, distinct positive
scales/translations, nonzero frame origin, a nontrivial proper rotation and
nonzero nonstatic phi exercise all five objectives. Expectations use transformed
points, signed tetrahedral moments and explicitly rotated inertia tensors,
independently of production Givens and frame mathematics. Finite differences at
3e-4, 1e-4 and 3e-5 are recorded alongside an independent complex-step oracle.
The accuracy gate additionally requires five-point centered differences of an
independent long-double reference at all three original steps. The reference
uses at least 64 mantissa bits from actual BBW interpolation and handle motion
through point transforms, integration and objective evaluation. Every entry
must meet 1e-6 relative error when its analytic magnitude exceeds 1e-7, or the
predeclared 2e-9 absolute rule otherwise. No fixture or tolerance is changed.
The task report preserves the raw ordinary double finite-subtraction
discrepancies separately; these must not be called all-green double FD.
