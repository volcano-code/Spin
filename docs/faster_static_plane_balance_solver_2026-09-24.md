# Faster §4.2 static maximum-mass plane solver — 2026-09-24

## Scope

Added a dedicated numerical path for the paper's 3D static maximum-mass
problem. It uses the derived vertical interface
`beta(x) = 1 + lambda_1*x_1 + lambda_2*x_2 > 0` and solves the two total
horizontal first-moment equations with the true `10 x 4` plane-family moment
differential. It does not pad this derivative to ten dimensions or route it
through the general Equation-14 quadric Newton solver.

The optimum flag means that the paper's plane threshold/KKT conditions hold
within the declared numerical constraint tolerance. It is not an
interval-arithmetic proof and does not certify infeasibility on failure.

## Verification

The registered analytic regression uses a box design domain
`[-1,2] x [-1,1] x [-1,1]`, zero fixed-shell moments, and the prescribed target
at the origin. Starting from `lambda=(-2,0)`, the solver recovers
`lambda=(-1,0)`, which fills `[-1,1] x [-1,1] x [-1,1]`: analytic variable
mass is 8 and both horizontal first moments are zero. The test confirms that
the solve accepts nontrivial steps, retains rejected line-search evaluations,
and withholds the optimality flag when the iteration budget is zero.

The convergence residual is normalized by the fixed input scale
`materialDensity * (outer-target-frame bounding-box side)^4`. A 1/1000-scale
copy is included to ensure small raw moments cannot be mistaken for a solved
constraint system.

Verification command:

```text
cmake --build build/strict-algorithms-exact-current-20260922 --target plane-static-balance-regression plane-family-regression -j2
ctest --test-dir build/strict-algorithms-exact-current-20260922 -R '^(spin-faster-static-plane-balance|spin-faster-plane-family)$' --output-on-failure
```

Result: **2/2 passed**. The existing plane-family regression independently
checks analytic moments, derivatives, frame transforms, contact rejection,
and multiple finite-difference step sizes.

## Not established

- This does not exercise the original vase or any author CAD input.
- Lower-center-of-mass and Pareto static objectives remain outside this new
  solver path.
- A line-search, singular-Jacobian, or iteration-limit result is not an
  infeasibility certificate.
- Full rank-deficient Faster optimization, Cherchi/CGAL paper-backend
  equivalence, and the full original case matrix remain open.
