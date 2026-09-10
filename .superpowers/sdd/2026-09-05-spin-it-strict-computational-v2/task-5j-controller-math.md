# Independent Task5j mathematical reference (not a solver result)

Current native implementation is in progress. This file establishes expected quantities independently before linking new code. Small analytic meshes are SYNTHETIC_REGRESSION only.

## Exact rational element assembly

`task-5j-check-fem.py` imports only Python standard library. It derives the eight cells in `[1/4,1/2]^3`, with grid spacing1/8, from literal integer Cartesian indices. Each monotone permutation path through a brick is a tetrahedron. Instead of native element gradient code, it solves four affine interpolation equations using Fraction Gauss-Jordan elimination; gradients are the last three coefficients of each exact affine shape function. It assembles exact rational stiffness and barycentric mass, then computes dense Q for this27-node independent fixture only. Production must remain sparse/matrix-free.

48tets, each volume1/3072, total1/64;26boundarynodes, onlynode13interior. Tetrahedral unordered-face multiplicities are exactly1or2, constant stiffness action exactly0. Field is fixed a priori by `x_i=((7*i)%11-5)/8`, unrelated to an optimized solution.

| Rule | Corner0 M | Center M | Exact E=0.5x^TQx |
|---|---|---|---|
| BARYCENTRIC_TET | 1/2048 | 1/512 | 70507/288 |
| RESTRICTED_GRID_VORONOI | 1/4096 | 1/512 | 37627/144 |

Voronoi reference is independently integrated via the clipped midpoint dual boxes: each coordinate width is1/16 at exterior endpoint,1/8 at middle. It does not call or copy production's per-cell accumulation. Both masses sum1/64 and share the same stiffness. Differences in energy are expected, not tolerance noise.

## Why per-brick V/8 is exact for this restricted grid Voronoi definition

For any point in a represented fine brick, every global grid-node coordinate lies outside or on that brick's interval along each axis. The nearest available node can be selected among the brick's eight corners, all of which exist by SourceCageGrid construction. Separability of squared Euclidean distance makes its nearest-corner region the product of three half-intervals. Each has volumeV/8. Sum these regions only over represented bricks to clip the global nearest-node cells to the selected domain. Holes and disconnected components do not add phantom volume. The proof uses actual represented interval lengths, so does not require exactly equal floating-point widths.

This definition is explicit and independently verifiable; it is not evidence of the author's unrecovered mesh and boundary-clipping settings. Freudenthal six-tet meshing remains IMPLEMENTATION_DEFAULT. Barycentric lumping remains a separate comparison, not called Voronoi.

## Boundary reference

For each cube face with normalized two-dimensional coordinates(u,v) and diagonal00–11, independent weights at corners(00,11,10,01) are `(1-max(u,v), min(u,v), max(u-v,0), max(v-u,0))`. Multiple incident faces agree on edge and vertex rows. The script checks all26boundarynodes and maps actual handle coordinates rather than assuming handle order. Interior is unconstrained. Off-grid trilinear transfer of these nodal weights need not equal this triangle function; the existing square-diagonal counterexample remains applicable.

## Primary source reread

[BBW2011](https://igl.ethz.ch/projects/bbw/bounded-biharmonic-weights-siggraph-2011-compressed-jacobson-et-al.pdf), Sec3.1 Eq3–6 and Sec3.4 Eq8, inspected2026-09-10: constrained energy includes inverse lumped mass; text specifically names Voronoi volumes, and both coupled and normalized-independent optimization variants exist. This subtask implements matrices/constraints, not either optimizer. The chosen grid construction is not the paper's recovered TetGen input. Future acceptance must retain both provenance and the applicable raw/final optimality checks.

## Pending actual checks

Once the implementer releases a stable main archive, controller will export its public API to JSON and compare native tetrahedra, all mass/stiffness entries, action, energy and complete boundary rows with this reference under both rules. Two fresh runs must match bytewise. Analytic reference execution alone proves nothing about native correctness.
