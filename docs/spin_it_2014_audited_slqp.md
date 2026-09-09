# Opt-in first-order KKT and restoration profile

`SlqpAlgorithm::AuditedKktRestorationV1` is a bounded numerical reconstruction,
with **IMPLEMENTATION_DEFAULT** controls. Section 5.4 of the public 2014 paper
describes SLQP, an active set and L-BFGS, but does not specify these controls.
This profile does not establish original author trajectories, exact Byrd
subproblems, second-order optimality, a local minimum, or a global physical
certificate. A stationary maximum also has zero first-order KKT residual.

The default remains `LegacyV3`; existing default diagnostics and scientific
hashes retain their original serialization. At the spin and adaptive APIs,
the audited profile requires `PaperGivens`. Dedicated dual-material and cage
consumers remain unsupported. The general `EqualityBoxProblem` interface
accepts a differentiable evaluation callback, a variable count, and a bounded
prefix count. Both real `AdaptiveFillFunctional` adapters and analytic problems
consume that interface through the same outer SLQP engine. Variables in the
prefix lie in [0,1]; all remaining coordinates are unbounded. In particular,
phi is not wrapped or restricted to [0,1]. Callback hashes are trace identities,
not an attestation of arbitrary callback code or all possible inputs.

```cpp
spin2014::SlqpLbfgsSettings settings;
settings.algorithm = spin2014::SlqpAlgorithm::AuditedKktRestorationV1;
auto result = spin2014::solveSlqpLbfgsPaper(functional, beta, phi, settings);
```

The Python request appends `slqp_algorithm="AUDITED_KKT_RESTORATION_V1"` and
the controls `slqp_restoration_trust_bound`, `slqp_lp_iteration_limit`,
`slqp_lp_time_limit`, and `slqp_lp_feasibility_tolerance`. Use it with
`formulation="PAPER_GIVENS"`. Existing positional arguments keep their meanings.
The result's diagnostics export the profile and all selected audit controls.

## KKT witness and independent residuals

The convention is L = f + lambdaᵀc − alphaᵀx + etaᵀ(x−1).
`evaluateKkt` independently checks dimensions and finiteness, then reports
equality primal, bound primal, stationarity, dual sign, and complementarity
infinity residuals. `full` is their maximum. A report's `valid` means its
numerical residual can be evaluated; it does **not** mean its KKT tolerances
pass. Invalid dimensions, nonfinite values, derived overflow, and bound duals
on unbounded coordinates cannot produce a valid zero residual. In JSON,
nonfinite diagnostic quantities are `null`, with a false validity/status flag.

The deterministic witness LP minimizes t >= 0 over free equality multipliers.
Writing r = grad f + Jᵀlambda, interior/unbounded coordinates impose
−t <= r <= t; lower faces impose r >= −t; upper faces impose r <= t.
Lower/upper faces also impose x*r <= t and −(1−x)*r <= t respectively.
The latter rows prevent near-face multipliers from hiding complementarity.
Generated alpha=max(r,0) on lower faces and eta=max(−r,0) on upper faces;
all other bound multipliers are exactly zero. This supports empty equalities
and rank-deficient Jacobians without selecting one arbitrary least-squares
witness. The independent residual evaluator recomputes all components after
the LP. A failed witness LP returns numerical UNKNOWN and no convergence.

Convergence requires equality primal <= `equalityTolerance` and all remaining
components <= `projectedGradientTolerance`, without a relaxed hidden tolerance.
The legacy projected-gradient field remains a separate diagnostic.

## Restoration, history and evidence

When an infeasible state has an unusable/inaccurate equality direction, or its
merit line search fails, the restoration LP minimizes sum(u+v), subject to
J*d+c=u−v and u,v>=0. Beta step bounds intersect [−x,1−x] with [−Delta,Delta];
unbounded coordinates have only the finite trust **step** bound. The optimization
active mask does not freeze restoration variables. No target or binary value
is changed. With p = ||c||₁ − ||J*d+c||₁ > 0, restoration accepts only a finite
evaluation with ||c(x+a*d)||₁ < ||c(x)||₁ and <= ||c(x)||₁ − armijo*a*p.
It uses the existing backtrack reduction and maximum trial count. An accepted
restoration establishes reduction, then returns to objective SLQP. Failure is
RESTORATION_STALLED or numerical UNKNOWN; it is never a nonlinear/continuous
infeasibility certificate.

The trace records actual x, gradient, c, J, lambda, alpha, eta, all KKT
components, full active membership (0 free, 1 lower, 2 upper), membership hash,
history reset reason and before/after sizes, accepted/rejected curvature pairs,
objective/restoration trials and reasons, LP status, predicted residual,
restoration direction, accepted step, and a final-state marker. The actual
accepted state is audited even at a maximum-iteration exit. An invalid initial
callback evaluation is a rejected problem; invalid or dimension-changing trial
evaluations are retained as rejected trials. Equality dimensions stay fixed
throughout a solve. The callback must supply correct analytic derivatives.

History is cleared before applying a metric to changed membership, after
restoration, and on every new solve/adaptive partition. Each curvature pair
uses the same current witness lambda at both endpoint gradients. Zero,
negative and nonfinite curvature pairs are rejected. Adaptive transition
counts use membership identity, including changes with unchanged counts.
Binary validation retains its previous rule and may still fail after a
converged continuous relaxation.

## Exact defaults and remaining obligations

New controls: Delta=1, LP iteration limit=10000, LP time limit=10 seconds,
LP primal/dual feasibility tolerance=1e-9. Positive finite trust/time limits
and a positive integer iteration limit are required; LP tolerance must be
in [1e-10,1e-4]. The pinned HiGHS revision is unchanged:
`7df0786de3088c832297e5ed821db236d8fab281`. Both LPs use simplex, presolve on,
one thread, parallel off, seed 0, output off, checked option-return values,
and explicit time/iteration caps. Sparse assembly is O(n*m), without n×n
dense matrices. These fixed choices are included in audited trace identity;
elapsed execution time is not. A cap may make a solve return UNKNOWN.

Reused defaults: maximum iterations=120, maximum line-search trials=24,
L-BFGS memory=8, equality tolerance=1e-8, projected-gradient tolerance=1e-7,
active tolerance=1e-10, merit penalty=100, backtrack reduction=.5, Armijo=1e-4.
The audited profile rejects active tolerances >= .5. Existing numerical cutoffs
remain: unusable direction norm <=1e-15, step cutoff 1e-14, predicted merit
floor 1e-14, decomposition threshold 1e-12, curvature cutoff
1e-12*||s||*||y||, and inverse scale clamping to [1e-6,1e6].

The numerical implementation still needs broader globalization/subproblem
analysis and high-resolution/source-geometry validation. Synthetic tiny and
Givens integration successes do not complete Task5 or the whole 2014 method.
