# Faster §4.3 buoyant plane solver — 2026-09-24

## Implemented path

Added `solveBuoyantStabilityByPlane`, a dedicated three-variable method for
the prescribed-waterline problem in §4.3. It parameterizes the interface as

`beta(x) = a + bx*x + by*y - z > 0`,

so the solid lies below a nonvertical plane. The unknowns are `(a,bx,by)`;
the fourth coefficient is fixed to `-1`. It calls the native plane Boolean
and moment evaluator, whose `10 x 4` coefficient derivative is reduced to the
three equations

* total mass minus `liquidDensity * displacedLiquidVolume`;
* horizontal first moment minus `centerOfBuoyancy.x * total mass`;
* horizontal first moment minus `centerOfBuoyancy.y * total mass`.

Each equation is normalized using fixed input density and design-domain
length scales. The offset parameter is internally scaled by the same fixed
length; the public interface and trace retain unnormalized coefficients.
Newton steps use a rank-checked reduced Jacobian and residual-decreasing
backtracking. The solver retains accepted and rejected evaluations and fails
closed on singularity, unavailable derivatives, invalid geometry, or budget
exhaustion. Nonconvergence is not classified as infeasibility.

On convergence, the result reports the independently computed center-of-mass
height, supplied-waterline metacenter, and their margin. This is only a
numerical prescribed-waterline check. It is not a dynamic fluid simulation,
not an interval proof, and not a proof that arbitrary non-planar feasible
solutions are absent. In particular, `equilibriumConstraintsVerified` means
that all three equations were met and the local reduced plane Jacobian had
full rank. The solver also converts the restricted plane coefficients to the
published equality multipliers and evaluates the finite-precision Eq. 7
conditions using the shared independent auditor: unit dual-monomial norm,
equality residuals, positive `mu`, and stationarity. The result keeps this in a
separate `eq7NumericalConditionsVerified` field; the audit's explicit proof
scope remains numerical Eq. 7/equality conditions only, not a proof of the
theorem or of global optimality. This does not reproduce an author example.

## Regression evidence

`spin-faster-buoyant-plane-balance` uses a closed, asymmetric rectangular
design domain with a known half-volume solution. Starting from a displaced
and tilted plane, it recovers the horizontal interface, checks the exact
mass and horizontal first moments, and tests both positive and negative
metacenter margins. A zero-iteration run must remain non-converged and must
not report stability. The test also checks the analytically recovered
multipliers (`lambda = [0, 0, 1/12]`, `mu = 1/6` for the unit-density
half-box), bounded Eq. 7 residuals, the restricted audit scope token, uniform
unit scaling, and that an iteration-limited result cannot carry an Eq. 7 pass.

Verified command:

```text
cmake --build build/algorithm-faster-current-20260924 \
  --target buoyant-plane-balance-regression plane-family-regression \
  plane-static-balance-regression hydrostatic-geometry-regression \
  spin-moment-problem-regression faster-theorem-conditions-regression -j4
ctest --test-dir build/algorithm-faster-current-20260924 \
  -R '^(spin-faster-plane-family|spin-faster-static-plane-balance|spin-faster-moment-problems|spin-faster-hydrostatic-geometry|spin-faster-buoyant-plane-balance|spin-faster-theorem-conditions)$' \
  --output-on-failure
PYTHONPATH=python:build/strict-current-tree-20260921/python python -m pytest -q \
  tests/python/test_spin_it_strict_audit.py::IndependentCheckerIntegrationTests::test_faster_buoyant_application_contract_repeats_verify_only_method_scope
```

Result: **6/6 focused CTest checks and 1/1 repeated evidence-audit integration
test passed**. This is a synthetic algorithm test. No Faster Table 1 asset,
author CAD model, or paper experiment is claimed as reproduced.
