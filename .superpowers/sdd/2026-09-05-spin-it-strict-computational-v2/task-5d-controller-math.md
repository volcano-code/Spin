# Task5d independent mathematics and oracle

No native solver/integrator or production helper is called in the following checks. This is a fixed analytic regression, not an original-paper case or proof of general optimizer convergence. Inputs were fixed in the brief before implementation.

## Fixed continuous-success / binary-failure example

Outer box volume1; CoM=(1/8192,0,1/2). Root side2, geometry depth5 gives h=1/16. Conservative accepted design is [3/16,5/16] × [-1/16,1/16] × [11/16,13/16], volume1/512 and x-first-moment1/2048. Initial optimization level2 splits it only at y=0, into two leaves.

The y equality forces equal beta; x equality then forces `beta=(1/8192)/(1/2048)=1/4`. The prescribed initial beta1/2 is not feasible, so a genuine successful optimization must move.

At final optimization level5 there are8 accepted h-cuboids, each volume1/4096, centers x={7/32,9/32}, y={-1/32,1/32}, z={23/32,25/32}. Their graph has12 edges (cube graph). Let X,Y,Z be these center vectors, v=1/4096 and b the8 fractions. Independently,

`m=1-v*sum(b)`, `ell=(1/2-v*dot(Z,b))/m`, `f=ell^2+0.5*b^T*L*b`.

`grad=2*ell*v*(ell-Z)/m+L*b`.

Equalities reduce to `dot(X,b)=1/2`, `dot(Y,b)=0`. No mass/first-moment approximation is needed to derive these expressions.

Controller ran independent SciPy SLSQP on these literal8variables/12edges, starting b=.5, ftol1e-14/maxiter300. This is a cross-check, not the native engine's acceptance gate. In5iterations it reported f=.24987788473640274; fractions (x-fastest) approximately:

`.25000047910038525, .24999284403005860, .25000047910038640, .24999284403005778, .25000811035373340, .25000047528340480, .25000811035373266, .25000047528340513`.

Equality infinity residual after division by cell volume3.469446951953614e-18. Independent least-squares equality multiplier gives stationarity infinity residual9.811565264677274e-11. Every variable remains fractional and below.5, so binary rounding keeps all material. Raw binary first-x residual is then1/8192, and input length4 scaling gives1/131072≈7.62939453125e-6, above prescribed1e-6. Thus a binary failure is expected even with successful continuous relaxation. Native values need meet its own tolerances, not these exact SLSQP trajectory numbers.

## Sparse Laplacian contract

For undirected unique face pairs E with unit weights, `b^T*L*b=sum_(i,j in E)(b_i-b_j)^2`. Hence symmetry,zero row sums,PSD and constant nullspace are independently testable from pairs. Actual corners/edges do not contribute pairs. Graph is dyadic bounding-cell adjacency and need not prove physical material connectivity.

## Merge / acceptance

If disjoint child regions have exactly equal fraction b, their removed moments sum to b*parentMoments; retaining b under merge preserves all10moments. Snapping distinct values to a common near-bound value generally does not. Coverage must be proved from keys/source emptiness, not a small floating moment difference. Even mass-preserving merge can change the unit-weight graph, so old KKT cannot certify a new partition without another solve.

## Independent binary checker

Solver uses published Algorithm1 surface plus cuboid moment columns. The existing `computeSurfaceMassMoments` implementation instead accumulates signed tetrahedra referenced near each mesh, first moments with /4 and second moments `(sum sum^T + sum vertex vertex^T)/20`. Streaming one12triangle cuboid at a time provides an independent integration route without building a purported manifold union mesh. No seal/manifold/physical success follows from moment agreement alone.
