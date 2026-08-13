# R3 support, target-construction, and multi-target acceptance record

Specification audit date: 2026-08-12 (Asia/Shanghai)

Acceptance status: **accepted at the documented R3 status ceiling**

The implementation, focused tests, R1/R2 non-regressions and local viewer
build described below completed on 2026-08-12. The clean recursive-clone gates
were executed after the implementation commit, as recorded below. This is R3
unit/integration evidence; it is not paper-level result reproduction.

## Authority and scope

The normative source is Prévost et al., *Make It Stand: Balancing Shapes for
3D Fabrication*, ACM TOG 32(4), 2013, especially Section 3, Figure 2,
Equation 5, and Section 5.3. As established in `provenance.md`, the authors'
beta archive is a read-only behavioral oracle; where it differs from the
paper, the paper governs this clean-room implementation.

R3 is responsible for:

- extracting a standing support polygon from mesh contacts relative to a
  gravity direction;
- constructing the toppling-angle safe region for a current center of mass;
- projecting the current CoM onto that region to obtain a standing target;
- constructing a suspension target from an attachment point;
- representing any finite, non-empty set of target/gravity pairs;
- evaluating the CoM part of paper Equation 5 for any number of targets; and
- exposing explicit classifications and failure states for degenerate or
  infeasible input.

R3 does not implement paper-faithful multi-target plane carving, deformation
handles, BBW/LBS deformation, the Laplacian energy, alternating optimization,
cord/mesh collision, mesh flattening, or physical fabrication validation.
Those phase boundaries are detailed below.

## Coordinate convention

Let `g` be a finite, nonzero gravity vector and

```text
ghat = g / ||g||
P_g  = I - ghat*ghat^T
```

where `ghat` points in the direction of gravity. Let `(u,v)` be an
orthonormal basis of the plane perpendicular to `ghat`. The support-height
coordinate is

```text
s(x) = -ghat dot x.
```

For the conventional `ghat=(0,-1,0)`, this is the world-space height `y`.
The lowest support height of mesh vertices `x_i` is

```text
h0 = min_i s(x_i).
```

All geometric predicates must be invariant to the magnitude of the supplied
gravity vector. Basis sign and ordering may differ, but changing the chosen
orthonormal basis must not change lifted 3D polygons, classifications,
targets, or energies.

## Standing support polygon

The paper defines the base of support as the convex hull of the mesh points
touching the ground, determined by the lowest mesh vertices relative to the
gravity direction. R3 therefore defines a contact set from vertices satisfying

```text
s(x_i) <= h0 + tau_contact,
```

using an inclusive, scale-aware comparison. With no requested flattening
band, `tau_contact` is only the numerical contact tolerance. With an explicit
nonnegative flattening threshold, the band additionally selects nearby
vertices that would be flattened onto the support plane. R3 computes the
support geometry but does not mutate or flatten the output mesh.

The selected points are projected into `(u,v)` coordinates and their unique
2D convex hull is the support polygon `S`. Duplicate points, collinear points
along an edge, vertex order, and triangle order must not change `S`.

R3 requires a nondegenerate support polygon with at least three non-collinear
projected contact positions and positive finite area. Empty, point, or segment
contacts are explicit construction failures rather than silently fabricated
bases. The support boundary counts as standing-feasible, in accordance with
the paper statement that a CoM projection at the boundary is sufficient for
static feasibility.

The selected support is immutable during a balancing optimization. R3 may
recompute a target from a new CoM height, but it must not silently recompute a
different contact hull from a deformed mesh. Turning contact vertices/voxels
into fixed deformation handles is deferred to R5--R7.

## Safe region and toppling angle

Let `q=P_g(c)` denote the current CoM's position in the support plane, expressed
in `(u,v)` coordinates, and let its perpendicular height above that plane be

```text
h = s(c) - h0.
```

Standing-target construction requires a finite, strictly positive `h`. For a
finite requested toppling angle

```text
0 <= gamma < pi/2,
```

the stability margin is

```text
r = h * tan(gamma).
```

The safe region is the inward Euclidean offset of the convex support polygon:

```text
S_safe = S eroded by a radius-r disk
       = {x in S : distance(x, boundary(S)) >= r}.
```

Equivalently, every supporting half-plane of `S` is shifted inward by `r`.
This form avoids a dependency on the beta archive's CGAL straight-skeleton
implementation and provides analytic polygon fixtures.

The current state is classified using inclusive, tolerance-aware predicates:

- `outside`: `q` is outside `S`;
- `feasible`: `q` is in `S` but outside `S_safe`;
- `safe`: `q` is in `S_safe`.

When `q` is in `S`, the corresponding toppling angle is

```text
theta = atan(distance(q, boundary(S)) / h).
```

The safe boundary is included: `theta=gamma` is safe. A support-boundary
projection is classified as safe at `gamma=0`; for every positive `gamma` it
is standing-feasible but not safe.

For `r=0`, `S_safe` must equal `S`; no artificial positive offset is allowed.
At the polygon inradius, the exact offset can collapse to a point or segment.
The result type must represent that lower-dimensional safe set, or report a
specific documented `degenerate_safe_region` state while preserving the
attainable target locus. If `r` exceeds the inradius, the requested safe region
is empty and construction must report `infeasible_safe_region`. It must not
silently claim that the requested toppling angle has been attained.

## Standing target construction

The desired tangential target is the Euclidean closest-point projection

```text
q_star = argmin_{x in S_safe} ||q-x||.
```

If `q` is already in the safe region, `q_star=q`; its projected CoM energy is
therefore zero. Outside a convex polygon, projection may lie on an edge or a
vertex and must be deterministic at equal-distance features.

Paper Figure 2 says the current CoM height is used to define the safe region.
Only the target component perpendicular to gravity affects Equations 2 and 5.
R3 should therefore store the 2D target with the support plane and expose a
canonical 3D target lifted onto that plane. Lifting it instead to the current
CoM height is energy-equivalent, but implementations and CSV/debug output must
use one documented canonical convention rather than mixing both.

Section 5.3 allows the target to be recomputed after each carving/deformation
iteration using the current rather than starting CoM. R3 target construction
is consequently a pure update from a frozen support polygon, current `c`, and
`gamma`; the scheduling of those updates belongs to the later alternating
optimizer.

## Suspension target

For suspension mode, the user supplies a finite attachment point `a`, desired
orientation, and nonzero gravity direction. The paper selects

```text
c_star = a
```

and equilibrium minimizes

```text
0.5 * ||P_g(c-a)||^2.
```

Thus all points on the gravity line through `a` have zero projected energy.
The paper assigns responsibility to the user for choosing a configuration in
which the cord does not intersect the mesh and the attained equilibrium is
stable. A diagnostic may check that the CoM is below the attachment point in
the direction of gravity, but that check must remain distinct from Equation 2
and must not alter its value. R3 does not claim cord clearance or suspended
stability from projected energy alone.

## Multiple targets and Equation 5

For a finite, non-empty collection of targets

```text
T = {(c_star_k, g_k)} for k=1..K,
```

R3 evaluates the CoM term

```text
E_multi(c,T)
    = sum_k 0.5 * ||P_{g_k}(c-c_star_k)||^2.
```

Each target has its own independently normalized gravity direction. The paper
uses an unweighted sum for any number of bases: it is not divided by `K`.
Duplicating a target therefore doubles that target's contribution. A standing
target and a suspension target share the same computational target/gravity
representation and may appear in the same collection.

Paper Equation 5 is the complete later optimization energy

```text
(1-mu) * E_multi(c,T) + mu * E_M(M_O).
```

R3 implements and verifies `E_multi`. Weighting by `mu` and evaluating the
shape-preserving deformation term `E_M` are R6 responsibilities. For multiple
standing bases, the paper notes that a solution requires a non-empty
intersection of the gravity-extruded support polygons. R3 may expose this as a
feasibility diagnostic; it is not evidence that the later optimizer finds a
solution.

## Expected implementation contract

Exact C++ names may differ, but the public behavior should separate four
concepts:

1. an immutable 2D convex `SupportPolygon` plus its 3D plane frame;
2. a `SafeRegion` result that can distinguish polygon, segment, point, empty,
   and invalid input where relevant;
3. a normalized `BalanceTarget` containing `c_star`, `g`, its construction
   mode, and optional support metadata; and
4. a `BalanceObjective` or equivalent non-empty vector of targets with
   per-target projected residuals and their summed energy.

Required operations are equivalent to:

```text
extractSupportPolygon(vertices, gravity, contact_threshold)
constructStandingTarget(support, current_com, gamma)
constructSuspensionTarget(attachment_point, gravity)
classifyStandingState(support, current_com, gamma)
projectedResidual(current_com, target)
multiTargetProjectedEnergy(current_com, targets)
```

All constructors must reject non-finite coordinates, non-finite thresholds or
angles, negative contact thresholds, invalid angle ranges, zero/non-finite
gravity, and empty target sets. Errors must be typed/status-bearing or throw a
documented standard exception; they must not propagate NaNs.

The R2 single-target energy helper remains a compatibility oracle. A one-item
R3 target set must agree with it value-for-value within the accepted R2
tolerances. The existing single-target carving constructor may remain as a
compatibility overload backed by a one-item target set. Full `K>1` carving is
not required until R4, but production code must not label a single-target
carve as a completed multi-target optimization.

## Paper/beta differences and ambiguities

The following beta behavior is informative but not normative:

1. `src/core/Support.cpp:67-105` forms a contact band and convex hull, but its
   strict threshold makes a zero threshold select no vertices. R3 uses an
   inclusive scale-aware contact predicate.
2. `src/core/Support.cpp:19-43` uses a strict bounded-side polygon test, which
   excludes the support boundary. The paper explicitly treats the boundary as
   statically feasible, so R3 includes it.
3. `src/core/Support.cpp:108-140` implements the `h*tan(gamma)` inset, but when
   the current projection is already safe it still projects to the nearest
   safe-region edge. R3 performs closest-point projection to the closed region
   and preserves a point already inside it.
4. The beta clamps the safe offset to at least `1e-4`. The paper specifies no
   such length-scale-dependent clamp; R3 uses `gamma=0 => S_safe=S`.
5. When the requested inset reaches/exceeds the inradius, the beta falls back
   to a straight-skeleton point or segment. R3 distinguishes an attainable
   degenerate locus from an empty, unattainable safe region and reports the
   latter explicitly.
6. Figure 2 discusses target height using the current CoM, while the beta
   stores the target on the support plane. These are identical under `P_g`.
   R3 adopts one canonical lifted representation and records that convention.
7. The paper does not prescribe the numerical contact tolerance, whether the
   flattening threshold comparison is strict, or the beta's use of
   `h0+threshold` as its display plane. R3 keeps the physical support plane at
   `h0`, makes the band inclusive, and records all user-provided tolerances.
8. `src/MIS/config.h:68-70` limits the beta to one or two objectives. The paper
   permits any number, so R3 uses a vector-sized abstraction.
9. `src/MIS/optimizer.cpp:270-289` correctly exposes the beta's sum of
   per-target energies. Its two-target carving path at
   `src/MIS/optimizer.cpp:70-134` requires both plane distances to be positive;
   the paper's later carving rule sums all positive contributions. That
   discrepancy is deferred to R4 rather than copied into R3 geometry.
10. The beta configuration path for initializing a suspended second objective
    calls its support-point setter without forwarding the objective index
    (`src/MIS/config.cpp:146-155`). R3 must not inherit this indexing behavior.

The beta defaults `contact threshold=0.02` and `gamma=5 degrees` in
`src/MIS_inc.h:26-27`. They are application defaults, not paper-mandated
constants or scale-independent acceptance criteria.

## Analytic and integration test matrix

| CTest | Minimum contract | Verification status |
|---|---|---|
| `support-polygon-analytic` | A rectangular contact base with duplicate/interior/collinear points; expected four-vertex convex hull and area; contact-band inclusion/exclusion; input-order, ordinary and large translation, rotation, gravity-scale, and plane-basis invariance; rejection of degenerate contacts and invalid values. | Passed |
| `safe-region-analytic` | For `S=[-2,2]^2`, `h=2`, and `gamma=atan(1/2)`, require `S_safe=[-1,1]^2`; outside, feasible, and safe states; inclusive boundaries; `gamma=0`; height validation and toppling-angle checks. | Passed |
| `standing-target-analytic` | Closest point/edge/vertex targets, identity for safe points, triangle/rectangle offsets, point/line maximum-inset loci, explicit over-inradius failure, and rigid-transform covariance. | Passed |
| `suspension-target-analytic` | A selected attachment and gravity; zero energy anywhere on its gravity line; known transverse `0.5*d^2`; gravity-scale invariance; invalid-input rejection; no cord-clearance claim. | Passed |
| `multi-target-energy-analytic` | `K=1` equality with R2; `K=2..4`; permutation invariance; duplicate-target contribution doubles rather than averages; independent gravity scaling; empty/invalid target rejection. | Passed |
| `multi-target-solver-boundary` | Production `CarvingSolver` consumes three objectives for aggregate queries, rejects an empty collection, and rejects incomplete `K>1` carving before any state mutation. | Passed |

Recommended numeric checks use scale-aware absolute-plus-relative tolerances,
with exact tolerance values frozen in the test source and reported here after
the implementation's conditioning is measured. Tests must compare geometry in
3D or through basis-invariant quantities rather than assuming one arbitrary
`(u,v)` orientation.

## Final integration and non-regression gate

The final gate used the pinned R0 toolchain and completed on 2026-08-12:

```text
RelWithDebInfo headless build and CTest                 14/14 passed
RelWithDebInfo full suite repeated until-fail:10        10/10 suites passed
Independent Debug build with assertions enabled        14/14 passed
Viewer configure/build                                 passed
Committed-tree recursive-clone headless build/CTest     14/14 passed
Committed-tree recursive-clone viewer build             passed
Parent repository and pinned libigl clean-state checks  passed
```

The six R3 tests passed without relaxing the R2 tolerances. The historical R0
defect trace remained byte-identical at SHA-256
`3f9108649fc329a136a5371597335a99275021b6ddd9a3984198ce2b8e8c1856`.
No official beta source or restricted beta asset was copied into this tree.

The current Chebyshev-inset construction enumerates active edge triples and
has worst-case `O(H^4)` work for a support hull with `H` edges. That is an
explicit R3 performance boundary: the verified fixtures and normalized sample
meshes have small hulls, while a future scalable implementation should use a
fixed-dimensional linear program. It does not change the accepted geometric
semantics or their unit evidence.

## R3/R4 and later-phase boundary

R3 establishes the geometry and objective semantics needed by later phases;
it does not establish a paper-faithful carving optimizer.

Specifically, R4 remains responsible for:

- defining one cutting plane per target;
- scoring a voxel by the sum of its positive per-target plane distances;
- sorting those scores deterministically;
- simulating every ordered deletion and selecting the globally best prefix;
- handling temporary energy increases without premature termination;
- applying the paper's too-few-removals rule; and
- validating cavity/material connectivity and wall-thickness behavior.

R5--R7 remain responsible for freezing support voxels as deformation handles,
BBW/LBS, Laplacian preservation, gradients, line search, and alternating
carving/deformation. R8 provides coherent output meshes and fabrication checks,
R9 provides the Python experiment interface, and R10 remains responsible for
paper experiment reproduction.

Consequently, successful R3 unit tests may advance support construction,
standing/suspension target construction, and Equation 5 CoM evaluation to at
most `VERIFIED_UNIT`. They cannot establish `VERIFIED_PAPER`, successful
multi-base carving, static stability of a manufactured object, suspension
safety, printability, or complete reproduction of *Make It Stand*.
