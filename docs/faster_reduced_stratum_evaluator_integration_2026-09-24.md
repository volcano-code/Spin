# Faster reduced-stratum evaluator integration — 2026-09-24

## Change

The unified `evaluateQuadricMaterialMoments` API now accepts the explicit
`DiscreteReducedStratum` derivative profile. It routes the exact material
Boolean and its birth-face provenance through the existing fixed-stratum
derivative for planes, simple-spectrum elliptic cylinders, hyperbolic
cylinders, and simple-spectrum paraboloids. The result carries the true
`10 x d` moment Jacobian and `10 x d` normalized coefficient basis in a separate
`stratumDifferential` field, including source provenance and geometric
conditioning diagnostics.

The ordinary `momentJacobian` and `derivativeDefined` fields retain their
ambient/full-rank meaning. The reduced path does not populate or pad the
ambient matrix. `solveQuadricHollowing` explicitly rejects this evaluator
profile because the current Eq. 14 implementation requires the full
10-dimensional Jacobian; this prevents a reduced derivative from being
silently treated as a full-rank Newton input. Python exposes the same profile
and separately returns the reduced Jacobian, coefficient basis, and status.

## Verification

- The new evaluator regression was first observed failing to compile because
  the public profile and result field did not exist.
- A second regression failed at runtime when the reduced plane profile reached
  the Eq. 14 solver without an explicit guard; it confirmed the solver could
  accept the request and later proceed with an invalid zero ambient Jacobian.
- After implementation, the targeted CTest group passed **5/5**:
  `spin-faster-source-motion`, `spin-faster-discrete-moments`,
  `spin-faster-plane-family`, `spin-faster-quadric-moments`, and
  `spin-faster-end-to-end`.
- The Python material bridge passed **1/1**, including a public call returning
  the reduced plane, elliptic-cylinder, and paraboloid Jacobians and bases
  with the ambient derivative explicitly unavailable. A non-symmetric fixed
  design-domain offset keeps the synthetic Boolean fixtures away from exact
  input-vertex contacts; those contacts remain correctly rejected by the
  production differential.
- The standalone `quadric-moment-evaluator-regression` also passed, including
  its certified-empty and reduced-solver-rejection controls. Its unified plane
  evaluator output was checked against independent central differences at
  three steps (`1e-4`, `3e-5`, and `1e-5`), with the normalized error required
  to decrease and remain below `1e-6`. The same unified interface was then
  exercised on an actual clipped elliptic cylinder (`10 x 7`) and a paraboloid
  (`10 x 8`); both Jacobian and coefficient basis matched fresh calls to the
  independent reduced differential routine exactly.
- The same unified interface now covers a clipped hyperbolic cylinder (`10 x
  7`). Its multi-step Boolean finite-difference check uses a local retraction
  onto the rank-two cylinder stratum, preserving both the zero quadratic
  eigenvalue and the cylinder condition `b·n = 0`. A straight coefficient-space
  perturbation is not valid for this test: a tangent direction that rotates
  the null axis lifts the zero eigenvalue at second order and changes the
  surface type, so the cap certificate correctly rejects that off-stratum
  sample. The measured normalized errors are:

  | Step | Relative error |
  |---:|---:|
  | `1e-4` | `4.04395e-8` |
  | `3e-5` | `3.63816e-9` |
  | `1e-5` | `3.98140e-10` |

- The focused CTest group
  `spin-faster-(quadric-triangulator|source-motion|discrete-moments|stratum-tangent|quadric-moments)`
  was rebuilt and passed **5/5** after adding this regression.

## Limits

This connects existing fixed-stratum differentials to the material evaluator;
it does **not** implement a stratum-constrained optimizer or extend the
Equation-14 solver to rank-deficient quadrics. It also does not solve
repeated-eigenvalue, cone, or cross-stratum transition cases, nor establish
author-case reproduction. A fixed-stratum derivative
still fails closed when source-feature topology, chart, or Boolean provenance
is not locally stable.
