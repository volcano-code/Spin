# Spin-It Faster planar polygon/disk moment integration — 2026-09-21

## Scope and status

Implemented `integratePlanarDiskIntersection` as a density-free six-moment
integrator for the intersection of oriented polygon loops and a circular disk.
It accumulates clipped straight boundary segments and the circular arcs that
lie in the polygonal region, using Green's theorem in a disk-centered frame.
The moments are translated back to model coordinates after integration.

This is a component of the paper's 2D fidget-spinner path, not a complete 2D
Boolean representation. It does not certify that input loops are simple,
properly nested, or nonintersecting; callers must provide valid loops with
material on their left.

The companion `PlanarFidgetCircleProblem` parameterizes the removable disk by
its center `c` and radius `|c|`, so its boundary passes through the target
origin. It combines geometry-integrated removed moments with caller-supplied
fixed-material moments, then solves for zero remaining in-plane first moment.
The dogleg root solver uses a centered finite-difference Jacobian, a scaled
two-variable trust region, explicit iteration/function-evaluation limits and
non-success statuses for exhausted budgets or singular models. Before each
Jacobian evaluation it classifies the polygon--disk crossing signature. At an
edge/vertex contact, or when no reduced centered stencil remains in the same
signature, it returns `TOPOLOGY_TRANSITION` rather than reusing an invalid
derivative. A regular stencil is reduced geometrically when necessary; this
allows convergence towards a valid limiting root without evaluating a stencil
across the event. These numerical controls are implementation defaults, not
recovered MATLAB `fsolve` settings or a claim of identical optimizer behavior.
A converged root is local, not a global optimality certificate.

The solver regression uses a synthetic square design region and analytically
balanced fixed patch. Its known root is a unit-radius disk centered at `(1,0)`;
the dogleg reaches it from an ordinary initial point and from the degenerate
zero-radius initial point. It also checks an already-solved input and an
exhausted function-evaluation budget. This is `SYNTHETIC_REGRESSION`, not the
author Flower model or Figure/Table experiment.

## Regression evidence

The C++ regression uses independently derived analytic values for:

- full disks with nonzero center and non-unit radius;
- an upper semicircle clipped by a square;
- a polygon wholly inside a disk;
- a disk wholly inside a polygon hole;
- a disk intersecting both the outer loop and a hole;
- disjoint and tangent regions;
- topology classification at a tangent contact and away from a transition;
- invalid negative radius.

Build: `build/strict-current-tree-20260921`, target
`planar-disk-intersection-regression`. The four registered planar regression
groups (dogleg, disk intersection, planar moment problem, and planar raw
moments) passed in focused CTest runs. This does not imply the full CTest suite
passed.

The initial red build failed because the registered source and header did not
exist. After implementation, the focused tests passed. No scientific solver
parameters, geometry inputs, or historical run evidence were changed.

## Remaining Faster 2D work

1. The finite-difference solver now rejects clipping-topology events instead of
   differentiating through them. It deliberately does **not** claim an analytic
   Jacobian across such an event; a formal generalized-derivative treatment is
   still outside this component. The independent planar-to-3D extrusion moment
   and fidget-objective check is covered by `spin-faster-planar-extrusion`: it
   uses a translated ring and an inward-oriented 3D rectangular cavity, so all
   planar raw moments including `Qxy` are compared against an independently
   integrated closed 3D material.
2. Add source-bound fixed bearing-ring and upper/lower layer contributions for the
   actual fidget-spinner setup.
3. Obtain or certify the author's original input, pose, and design-domain
   parameters before counting the paper experiment as reproduced.
