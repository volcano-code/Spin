# Faster discrete Boolean-to-Newton integration — 2026-09-24

## Scope

This is a synthetic method-validation regression, not a Faster author case.
It verifies that `solveQuadricHollowing` can use the same discrete full-rank
Boolean-to-moment evaluator for its residual and Jacobian during a nontrivial
Equation-14 solve, then export a nonempty material and cavity boundary.

The fixture uses the existing full-rank hyperbolic yo-yo quadric and a frozen
off-target initial coefficient vector. Its design domain is `[-0.8, 0.8]^3`,
which avoids the separately preserved boundary-tangency negative control at
`z=1`. Fixed-shell moments are set before the solve from the fixed synthetic
target. The chosen geometry is method-validation input only; it is not an
author asset or a paper-case reconstruction.

## Result

The positive discrete-profile solve converged in five accepted trace
iterations with normalized equation residual `4.09926e-16`. It returned a
nonempty material boundary and cavity boundary, reported derivative profile
`discrete-fullrank-quadric-boolean-transport-v2`, and moved away from its
off-target initial quadric. A second independent invocation with the same
frozen request and start produced identical output, trace, event, and geometry
hashes; the regression asserts all four equalities.

The same regression preserves a negative control using the original `[-1,1]^3`
domain. There, a line-search trial reaches a near-contact at `z≈1`; that
trial's Boolean output contains two degenerate faces (minimum squared
double-area `2.97585e-23`) and is rejected. The frozen evaluation budget is
then exhausted, so the solver returns `RESOURCE_LIMIT_REACHED` without
accepting or exporting the invalid trial. The case is deliberately retained
as a fail-closed geometry blocker, not classified as mathematical
infeasibility.

## Verification and limits

- `spin-faster-end-to-end`: **1/1 passed** after rebuilding
  `spin-it-faster-end-to-end`.
- The complete registered CTest suite was not run in this batch.
- This closes one method-level integration gap only. It does not establish
  Newton convergence for arbitrary starts, all quadric types, topology
  transitions, author inputs, or Faster Table 1 / Figures 10–13.
- The regular default evaluator remains backward-compatible; this regression
  explicitly selects `DiscreteFullRankQuadric` and does not change defaults.
