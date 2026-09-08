# Opt-in Spin-It 2014 Givens formulation

`SpinIt2014Request(formulation="PAPER_GIVENS", initial_phi=1.7, ...)`
selects a joint fill/angle solve. The native equivalent sets
`AdaptiveVoxelSpinRequest::formulation = Spin2014Formulation::PaperGivens`
and `initialPhi`. The default remains `LEGACY_ANGLE_ELIMINATED`.

The new evaluator implements the paper's raw moment equations, including the
explicit Givens constraint. Its 11 coordinates are
`(M,Lx,Ly,Lz,Qxx,Qyy,Qzz,Qyz,Qxz,Qxy,phi)`. For yo-yos, the inertia uses the
axis origin without a center-of-mass shift. Tops subtract only `Lz²/M` from
the two transverse diagonal components. Constraints use raw first and cross
moments. These functions differ from the existing full-COM, angle-eliminated
functions away from the equality manifold. Objective weights and the five
objective choices retain their existing meaning.

The angle is finite and unbounded in radians. It participates in the same
SLQP step and L-BFGS history as fill values; only fill values are bounded in
`[0,1]` or receive Laplacian regularization. No angle wrapping, diagonalization
postprocess, or angle reset occurs between adaptive passes. Static balance
has no angle unknown; static and legacy requests require `initial_phi=0`.
The standalone unified moment evaluator accepts finite static angles with
zero angle derivatives, while the solve interfaces reject irrelevant angles.

Explicit results expose `formulation`, `initial_phi`, and `final_phi` through
the Python facade. `diagnostics["adaptive_passes"]` includes each pass's
initial/final angle and `slqp_trace` entries containing the current angle,
objective, equality residual, and accepted step size. A row describes the
state before its recorded step. The final result can therefore differ from
the last row when the iteration limit stops immediately after an accepted
step. Scientific identities bind the selected formulation and actual angles,
including this last accepted state. Legacy hash preimages and JSON fields
remain unchanged.

Binary rounding still uses the existing `beta >= 0.5` rule. Validation uses
the selected paper equations at the final optimized angle; a residual above
the requested binary tolerance returns `BINARY_VALIDATION_FAILED`. Failed
optimization is not evidence of infeasibility.

The native/Python cage consumer and the dedicated dual-material product/export
wrapper reject this profile. Direct adaptive density weighting is supported.
This increment does not restore the full paper SLQP/restoration method,
surface-source octree, level-9 cases, cage optimization, dual digital cuts,
or original experiments. Overall reproduction remains `PARTIAL_BLOCKED`.

Evidence and implementation details are in
[`task-5a-report.md`](../.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5a-report.md).
