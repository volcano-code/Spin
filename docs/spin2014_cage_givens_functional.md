# Explicit Givens cage functional prerequisite

`CageSpinDeformationFunctional::evaluatePaper(beta, phi, handleParameters)`
evaluates the explicit paper-angle extension using `Spin2014GivensProblem`.
`evaluatePaperJoint(state)` accepts the assembled layout
`[beta, phi if applicable, existing free handle parameters]`. `betaCount()`,
`paperHasAngle()`, `paperHandleOffset()` and `jointParameterCount()` expose the
layout. StaticLowering has no angle column and requires `phi == 0`; rotating
modes accept any finite angle in radians.

The existing `evaluate(beta, handleParameters)` remains the eliminated-angle
API. Both paths use one geometry/moment chain and one regularizer evaluation.
The explicit path never evaluates the legacy objective as a prerequisite.
Fixed rotations remain in the existing LBS basis. Each unlocked user handle
contributes `[tx, ty, tz, uniformScale]` in handle order; locked handles
contribute geometry but no parameters. No handles are sampled or truncated.

For raw moments in order `[mass, x, y, z, xx, yy, zz, yz, xz, xy]`,

```
s(beta,h) = s_outer(h) - sum_k beta_k s_cell_k(h)
ds/dbeta_k = -s_cell_k(h)
ds/dh = ds_outer/dh - sum_k beta_k ds_cell_k/dh
ds/dphi = 0
```

The outer mesh and all void cell meshes are deformed before integrating these
moments. Objective and equality derivatives are chained through the full raw
moment matrices and the independent explicit angle column.

`CageSpinPaperEvaluation.deformation` contains the decomposed energies, total
fill/handle gradients, equality values and their fill/handle Jacobians, raw
fill/handle moment Jacobians, final material moments, and actual deformed
outer/cell vertices. The wrapper supplies `angleGradient`, `jointGradient`,
`jointConstraintJacobian`, and `rawMomentJacobian` with the assembled layout.
All computed output must be finite; malformed state, including beta outside
`[0,1]`, wrong sizes, nonfinite values and a nonzero static angle, throws.
This is a first-derivative interface, not a joint Hessian or a solver status.

All five objective selections and the independent inertia/center weights are
preserved. YoYo uses the two raw transverse ratios and centroid x/y/z,
xz/yz and Givens off-diagonal equalities. Top objectives use transverse
moments reduced by `s_z^2 / mass`, centroid x/y, xz/yz and off-diagonal
equalities. TopInertiaOnly adds no center energy; TopLoweringOnly adds the
weighted squared centroid height; TopMassAndLowering adds weighted `s_z^2`.
StaticLowering uses squared centroid height and centroid x/y equalities.

The legacy evaluator centers the full inertia tensor off the equality
manifold, whereas the paper form uses its centroid equalities. Therefore
their values agree at the angle diagonalizing transverse inertia **on the
required centroid manifold** (x/y zero for tops, x/y/z zero for YoYo).
Arbitrary off-manifold values and ambient gradients need not agree, even
when the Givens off-diagonal term is zero. Tests retain separate off-manifold
value/derivative checks and centered legacy-value comparisons.

## Regularization and scientific limits

The local authenticated paper text, `build/strict-paper-sources/spin2014.txt`
lines 360–385 (Section 6.1), gives
`gamma_D,L ||L_D v'||^2 + gamma_D,v ||v'-v||^2`. The implementation retains
these factors and evaluates `0.5 gamma_L beta^T L beta` for fill smoothing.
The supplied symmetric fill Laplacian is preserved. The outer operator is
the existing uniform undirected mesh graph Laplacian: an
`IMPLEMENTATION_DEFAULT`, not Sorkine's rotation-aware operator and not a
recovered author discretization. This increment does not upgrade the legacy
regularizers to exact ones.

This request is still a dense legacy mesh adapter at **unit material density**
for outer and void meshes. It does not introduce independently configurable
positive densities. The source integration regression obtains weights from
an actual `SourceSurfaceCage -> SourceCageGrid -> SourceCageFem -> SourceCageBbw`
coupled solve and trilinear transfer for the original mesh and two internal
mesh cells. Every cage handle is retained (26 handles, 104 free handle
parameters, 107 joint columns in the rotating integration fixture).

That regression is a mesh-based integration check. Native source-adaptive
partition ownership/binding, audited SLQP with geometry rejection, final
material integration, geometry sealing, physical validation and original
author cases remain mandatory subsequent work. There is no source-adaptive
solver, geometry certificate, pass-flag shortcut or reproduction claim here.
The overall reproduction status remains `PARTIAL_BLOCKED`.

## Verification scope

The core functional regression does not depend on exact geometry or Faster.
It independently reconstructs LBS and uniform tetrahedral moments, checks
all five objective values and every equality, and differences beta, phi and
every free handle column at `1e-3`, `1e-4`, `1e-5`. It includes nonidentity
fixed rotation, a locked handle, unequal beta, nonzero regularizers, malformed
state rejection, and the existing legacy regression.

The separate source-BBW joint regression is registered only with exact
geometry and Spin-It 2014 enabled. It checks actual solve evidence, all
handles, independent signed-tetra moment integration and Givens expressions,
decomposed regularizers and all joint derivatives with the same three steps.
Neither test imposes source-adaptive ownership or final geometry authority.
