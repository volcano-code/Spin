# Source-bound linearized-rotation shape regularization

This component implements the Section 6.1 shape and displacement terms of
Spin-It (2014), using the cited Sorkine et al. (2004), Section 3, Equations
2, 5, 6, and 8–12. It is a bounded mathematical/source-chain increment, not a
joint optimizer, geometry-legality certificate, or original-case reproduction.
The original source surface and all actual cage handles participate.

Source authority: original `SpinItPaper.pdf` SHA256
`3add70de03df9b0e173028d55c775e8520fa89a6b21d9b269f926ccbe8a76e4e`;
its existing text extract SHA256
`5cf5ab047189575df48b01796de73c72ceadc602a3380b74f33d2dda81459161`.
The [official Sorkine PDF](https://igl.ethz.ch/projects/Laplacian-mesh-processing/Laplacian-mesh-editing/laplacian-mesh-editing.pdf)
has SHA256 `246498341881005d40a27309c5aca0e3d6efd2de5d29a85f95bbcdcc978dfab4`;
the locally retained text extract has SHA256
`703ddc9364c235d6eb440920c04b04cd28f079588d258514a6bc460eb34b7095`.
Uniform one-ring weights are published in Sorkine. Applying that uniform
variant to Spin-It is a documented formula reconstruction, not recovered
author-2014 experiment metadata.

## Fixed mathematical operator

For each original rest vertex `i`, `Ni` is its unique, sorted triangular
one-ring and `Li v = vi - mean(vj, j in Ni)`. The local order is `i` followed
by `Ni`; all local vertices enter the equal-weight least-squares fit.
For `p=(x,y,z)`, the coefficient order is `(s,hx,hy,hz,tx,ty,tz)` and

```text
A(p) = [ x   0   z  -y   1 0 0
         y  -z   0   x   0 1 0
         z   y  -x   0   0 0 1 ].
```

The differential block `Di` uses the first four columns evaluated at
`delta_i = Li vRest` and zero translation columns. The three residual rows
are `Ri = Di Ai^+ - Li`, scattered into xyz-interleaved global coordinates.
`LD` stacks those rows. The sign is fixed, even though energy is sign-invariant.

Construction uses `qj=(vj-vi)/ri`, with `ri=max_j ||vj-vi||`, a rest-only
algebraic preconditioner. Both `Ai` and `Di` are built from `q` and its
differential. The centered, normalized linear-transform coefficients absorb
`ri`; translation absorbs the shift. Thus `Di Ai^+` is unchanged as an action
on original-unit deformed positions: `LD` is dimensionless and energy remains
in squared length units. No deformed state, beta, phi, density, or moving frame
enters this construction.

Each local `m=3*(degree+1)` by 7 system uses Eigen's column-pivoted QR
preconditioned Jacobi SVD with thin U (`m×7`) and V (`7×7`). All seven singular
directions are required: `sigma_min > relativeRankTolerance * sigma_max`.
The hashed, configurable `1e-12` tolerance is `IMPLEMENTATION_DEFAULT`.
Nonfinite tolerance or values outside `(0,1)` fail. There is no damping,
pinning, dropped ring, silent rank truncation, or normal-equation inversion.
Planar noncollinear patches are valid; numerically deficient patches fail.

For current original surface positions `vPrime`,

```text
Eshape    = gammaD,L * ||LD vPrime||²
Edisplace = gammaD,v * ||vPrime-vRest||²
gradient  = 2 gammaD,L LDᵀ(LD vPrime)
          + 2 gammaD,v (vPrime-vRest).
```

There is no half factor or density. Translation, uniform scale, and fields
`(sI + [h]cross)v + t` are null modes of shape alone. Displacement still
penalizes movement. A genuine finite rotation generally has nonzero shape
energy, as expected for the paper's small-rotation model. No ARAP/polar,
cotangent, frozen-differential, cage-surface, or material-cell substitute is used.

## Actual source chain and composition

`LinearizedRotationLaplacian` is immutable, mathematically evaluable, and exposes
copied rest geometry, CSR arrays, required workspace, and operator/config hashes.
Its direct vertex API is not a source-validity certificate.

`SourceCageShapeFunctional` accepts only the actual immutable `SourceCageBinding`;
there is no arbitrary-weight constructor. It uses cached original-rest grid
embeddings and every column of the binding's actual final BBW field. The row
is neither normalized nor clamped again. In original handle order,

```text
dvi / duh = wih * (viRest-centerh)
dvi / dth = wih * I.
```

The shape wrapper accepts a view of original deformed outer rows and returns
weighted energies, residual, vertex gradient, and all handle derivatives.
The supplied vertex view remains a mathematical state; the wrapper does not
certify arbitrary externally supplied positions. The paper functional obtains
that view directly from its actual moment evaluation and adds both energies
and all handle derivatives. It does not recompute or mutate geometry. No
material-corner rows enter the regularizer; beta and phi additions are exactly
zero, and moments/constraints plus their Jacobians remain unchanged.

Paper weights default to zero. When both are zero, no shape operator is built
or evaluated (`shapeFunctional()==nullptr`). Displacement-only has an explicit
enabled source-shape identity and the real operator. Value-only returns identical
energies and no derivative arrays, and works when only derivative capacity is
insufficient. Both weights and all shape settings are validated and hashed,
including when disabled. Source geometry, source/BBW identities, frame, operator,
state, weights, settings, mode, and raw outputs bind the resulting identities.
Immutable copies repeat exactly; changing source or actual BBW identity remains
distinguishable even at zero motion. No acceptance flags are introduced.

## Storage and fail-closed contracts

Let `V` be vertices, `F` faces, `d` maximum unique degree,
`S=sum_i(1+degree_i)`, and `m=3*(d+1)`. Block CSR stores `V+1` size_t offsets,
`S` size_t vertex indices, and all `9S` doubles (including exact zeros).
`coefficientCount()` reports blocks; `sparseValueCount()` reports scalar
contributions/stored values, not a thresholded numerical-nonzero count.
Within each block, coefficients are row-major; row and column components are
xyz. Only forward and transpose sparse residual actions are evaluated.
No dense `3V×3V`, Gram, global `V×H`/`V×4H`, or full joint Hessian is formed.

Positive defaults:

| Limit | Default |
| --- | ---: |
| Original vertices | 1,000,000 |
| Unique one-ring neighbors | 4,096 |
| Stored scalar contributions | 30,000,000 |
| Construction double-equivalent workspace | 10,000,000 |
| Evaluation workspace doubles | 10,000,000 |
| All-handle gradient doubles | 1,000,000 |
| Full weight-row doubles | 1,000,000 |

All integer/index/byte products are checked before the respective allocations.
Initial `6F` directed edge pairs are sorted/uniqued; their retained capacity is
still charged at `6F`. Before that allocation the geometry/edge byte bound is
checked; after the unique degrees are known, the full simultaneous construction
bound is checked before allocating copied geometry, CSR, or local SVD buffers:

```text
retainedBytes = 3V*sizeof(double) + 3F*sizeof(int)
              + (V+1+S)*sizeof(size_t) + 9S*sizeof(double)
edgeBytes     = 6F*sizeof(pair<size_t,size_t>)
constructionDoubles = ceil((retainedBytes+edgeBytes)/sizeof(double))
                    + 8*m*7 + 2048.
```

The eight local `m×7` allowances cover A, thin U, the scaled work matrix,
column-pivoted QR storage, and packing/scratch conservatively. The fixed 2048
doubles cover seven-dimensional matrices/vectors, object/index metadata and
bounded sort/local stack overhead. The installed Eigen implementation allocates
only the more-rows-than-columns preconditioner for these systems; its QR
Householder thin-U workspace has seven entries. It never requests full U.
The more permissive sparse cap never overrides the construction cap.

Operator evaluation needs exactly `3V` doubles for returned residual in
value-only mode or `6V` for residual plus returned vertex gradient. Views do
not copy the deformed matrix; contractions have only fixed small scratch. The
wrapper moves these arrays and additionally retains `4H` gradient doubles
and a simultaneous `H` weight row under its separate explicit caps. Value-only
needs neither. The paper uses moves as well, so no duplicate V-sized arrays
or hidden global GEMM packing enter this chain.

Caller inputs and upstream source/binding/BBW/conforming/moment state and
outputs are separately bounded by their existing owners and are excluded from
these additional component budgets. Constant control/hash/string overhead is
independent of V/H and included in the fixed construction allowance or outside
the variable evaluation-array budget. Resource failure, allocation failure,
zero budgets, and arithmetic size/cap overflow throw `length_error` with
`RESOURCE_LIMIT`. Invalid geometry/config/state throws `invalid_argument`;
nonfinite arithmetic throws `overflow_error`. They never become infeasibility
or silent zero energy. Repeated face indices, invalid indices, zero-area faces,
isolated vertices, malformed/nonfinite positions, and failed local rank are
rejected. Broader manifold/self-intersection/mesh-legality checks are excluded.

The regression uses independent raw-coordinate long-double QR fits on deformed
points, mathematical nullspace/covariance tests, three fixed centered FD steps,
actual source/BBW fixtures across five paper modes, all joint columns, cap
boundaries and source identity changes. These are `SYNTHETIC_REGRESSION`, not
original experiment or paper-completion evidence.
