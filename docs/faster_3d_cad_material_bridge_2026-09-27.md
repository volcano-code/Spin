# Faster 3D quadric-to-CAD material bridge (2026-09-27)

## Delivered scope

Added `export_spin_it_faster_quadric_material_step`, which accepts only a
`FasterMaterialResult` bound to a revalidated `SealedSpinArtifact`. It recovers
the locked Faster coefficients, type, sign and frame, transforms the quadric
into the artifact's normalized model frame, and sends the original editable
design domain plus every nonempty fixed region to the pinned OCCT bridge. The
bridge constructs exact analytic ellipsoid, finite-carrier elliptic-cylinder,
and one-/two-sheet hyperboloid Booleans; fixed regions are imported as closed
faceted solids. It writes/reloads both B-rep and STEP and reports per-region
and aggregate ten-term mass moments.

The Python handoff checks result identity, artifact freshness, region labels,
normalized densities, exact source-region moment recomputation, CAD/STEP
round-trip identities and moments, and refuses to overwrite existing output.
The sidecar binds the Faster/artifact digests, CAD input digest, OCCT executable
SHA-256, and output-file SHA-256 values.

## Verification evidence

The integration fixture is explicitly `SYNTHETIC` / `METHOD_VALIDATION_ONLY`:
it evaluates and seals a small ellipsoid-carved design region with six fixed
shell regions, exports the exact analytical OCCT result, then checks the
reloaded entity moments against the independently integrated sealed material
at the existing `1e-6` normalized threshold. It is not an optimized solver
result, an author asset, or a paper-case reproduction.

Verification run on 2026-09-27: OCCT 7.9.3 bridge rebuilt successfully;
`test_occt_step_roundtrip` passed **20/20** (including the new sealed-material
handoff), `py_compile` passed for the edited Python modules, and `git diff
--check` reported no whitespace errors in the touched paths.
After the final CAD JSON evidence-key correction, the sealed 3D integration
test was rebuilt and rerun independently: **1/1 passed**.

Follow-up verification on 2026-09-27: sign-invariant cylinder triangulation
and analytic cylinder material export were added. The native triangulator
regression and Faster material bridge CTest passed **2/2**; the complete OCCT
round-trip file passed **21/21**, including both sealed ellipsoid and sealed
elliptic-cylinder paths. The cylinder fixture's CAD-to-sealed normalized
ten-moment error is checked against `1e-6`. Python compilation and
`git diff --check` also passed.

## Scope limits

- Current typed Boolean support includes sealed plane candidates and
  converged prescribed-waterline plane results, plus ellipsoids, elliptic
  cylinders, one-/two-sheet hyperboloids, and hyperbolic paraboloids. Elliptic
  paraboloids and other unsupported types fail explicitly; unbounded carriers
  are bounded beyond the design-domain projection.
- Each fixed input mesh must be a single closed shell that OCCT can import as a
  solid. Arbitrary nested/multi-shell fixed entities are not yet supported.
- The CAD B-rep is analytic for the supported quadric; the sealed solver
  output is a triangulated Boolean. If their independent moments differ beyond
  `1e-6`, export is rejected rather than relabeled as validated.
- This is an OCCT integration/method check. It does not replace the Faster
  paper's Boolean backend, prove the full set of quadric types, or satisfy
  author-CAD/table/figure reproduction requirements.

## Reduced-stratum fail-closed regression (2026-09-27)

The candidate-evaluation entry permits `DISCRETE_REDUCED_STRATUM` because a
geometry-only candidate does not claim an ambient optimizer Jacobian. It still
has to pass the same sealed-material partition audit before it can be exported.
An anisotropic elliptic-cylinder fixture (`q=x^2+2y^2-r^2`) initially failed:
`q` and `-q` selected the radial eigenvectors in different signed-eigenvalue
orders, producing different side-face diagonals. Sorting the two radial axes
by eigenvalue magnitude makes the polygonal interface invariant to sign
reversal. The new C++ regression and native sealed-material test pass. With 96
angular segments, the same candidate also passes analytical OCCT Boolean,
STEP/B-rep reload, and independent ten-moment checks at `1e-6`. This is a
synthetic method-validation result only; it is not an optimized solver result
or a paper-case reproduction.

## Initial solver-to-CAD integration probe (superseded below)

A nontrivial locked-material probe used a translated cube shell,
`SPINNING_TOP`, `DISCRETE_FULL_RANK_QUADRIC`, the dual-raw equation profile,
one frozen initial quadric, and the existing 60-evaluation budget.
`solve_and_seal_faster` returned `RESOURCE_LIMIT_REACHED` with no sealed
artifact, so the CAD exporter was correctly not called. This is a solver
failure for that exact method-validation fixture, not an infeasibility result.

At that point, the repository had a separate converged full-rank regression,
but no locked-material adapter for it. That gap is addressed by the explicit
adapter described in the following section. The strict CAD moment gate remains
unverified for that solver result; it is not upgraded merely because the
quadric Boolean can now be constructed. Do not change the solver budget or
substitute a candidate result to make this gate appear passed.

## Full-rank solve sealing and hyperboloid CAD probe (2026-09-27)

Added the explicit `solve_fullrank_and_seal_faster` path. It calls the existing
discrete full-rank Newton solver, records the full solver result and request
override identity (including the center-of-mass weight), re-evaluates the
solution and complement through that same discrete moment map, and only then
creates the sealed artifact. The frozen method fixture converges in 13
function evaluations and passes the independent seal checks. The older
`solve_and_seal_faster` entry remains unchanged and is not relabeled.

The pinned OCCT material protocol now also constructs two-sheet hyperboloid
lobes. The full-rank fixture returns a two-sheet hyperboloid. OCCT builds and
reloads its analytic multi-material STEP/B-rep, but the analytic result differs
from the solver's level-1 triangulated design material by
`0.11584316133993125` normalized volume/moment error (design-region raw volume
delta `-0.033688249986097596`). An exploratory level-3 attempt did not produce
a sealed converged result and therefore does not establish that refinement
closes this gap. The exporter rejects the `1e-6` gate and does not publish the
unvalidated B-rep, STEP, or identity sidecar. This is **not** a passed
solver-to-CAD validation; the discrete-to-analytic geometry discrepancy
remains a blocker.

Verification on 2026-09-27: `spin-faster-material-bridge-python` passed **1/1**
(119.57 s), including the actual full-rank solve/seal and fail-closed CAD
regression. The complete `test_occt_step_roundtrip` module passed **21/21**
(121.57 s), including the existing ellipsoid/cylinder sealed candidates and
standalone one-/two-sheet hyperboloid analytical oracles.

## Converged prescribed-waterline plane A to CAD (2026-09-27)

The dedicated plane solver's converged A result now uses this same sealed
CAD path. The seal retains the canonical request snapshot, its request digest,
the solved plane coefficients, and the target frame. Before CAD export, the
bridge reconstructs and hashes that request, confirms its outer/design meshes
match the sealed problem, re-evaluates both halfspaces, and requires their
vertices/faces to match the retained and carved sealed regions exactly. It
then transforms the physical-coordinate plane into the artifact-normalized CAD
frame before invoking OCCT. The fixture uses a low-density shell represented
as six closed solids, avoiding the existing multi-shell fixed-entity limitation.

The focused test
`test_converged_buoyant_plane_a_is_sealed_from_actual_complement_geometry`
passed after first reproducing the missing solver-operation dispatch, then
exposing and fixing the physical-to-normalized plane transform. It now checks
the actual converged plane solver result, sealed artifact binding, exact plane
material re-evaluation, analytic OCCT halfspace, reloaded STEP/B-rep, and
independent per-material moments at the existing `1e-6` gate. This remains a
synthetic prescribed-waterline method validation, not a free-waterline,
physical, author-asset, or paper-case reproduction.

## Sign-canonical saddle Boolean consistency (2026-09-27)

An oblique hyperbolic-paraboloid candidate exposed a sign-polarity issue: the
forward-polarity sealed CAD export passed the strict material-moment gate,
while the reversed polarity missed it by `1.1338246e-6`. Refining the native
quadric grid from 64 to 128 barely changed the error. Direct comparison showed
that native `q` and `-q` outputs partition the design domain to roundoff, while
the corresponding old OCCT Boolean outputs had a mass deficit near `7.08e-8`.

The pinned OCCT bridge now builds one lexicographically sign-canonical saddle
carrier for both coefficient polarities and selects Common/Cut according to
the original polynomial's sign relative to that canonical form. This retains
the requested `q>=0` region while avoiding two differently parameterized
coincident saddle surfaces in the complement pair. The regression checks
positive and reversed polarities, resolution 64 and 128, B-rep/STEP publication,
the unchanged `1e-6` per-material gate, and that the two 64-level CAD material
moment vectors sum to the analytic design-box moments within `1e-9`.

The focused oblique saddle regression passed **1/1** and the complete Faster
material bridge CTest passed **1/1** after rebuilding the pinned OCCT probe.
This closes only the tested saddle polarity consistency issue; it does not
establish full OCCT/CGAL entity equivalence, Faster-wide quadric coverage, or
paper-case reproduction.

## Elliptic-paraboloid CAD discretization audit (2026-09-27)

The existing coarse-paraboloid regression was rerun against the current build.
It correctly fails closed rather than publishing analytic STEP/B-rep whose
moments do not match the sealed native triangulation. A fixed unit-cube design,
fixed shell, coefficient vector, material densities, and all solver settings
were held constant while only `paraboloid_grid_resolution` changed:

| Native paraboloid grid | Maximum normalized material-moment error | Design volume delta |
|---:|---:|---:|
| 32 | `1.2755738841e-2` | `-9.4212338962e-4` |
| 64 | `3.3774490626e-3` | `-2.4710690590e-4` |
| 128 | `8.2780459786e-4` | `-6.0410757150e-5` |

The error falls by about four on each resolution doubling, consistent with a
second-order faceted approximation of an analytic curved surface. This
distinguishes the current failure from the saddle's complement-Boolean defect:
the strict `1e-6` mismatch is resolution-dependent and the carrier is not
shown to be wrong. A naive second-order extrapolation would require grid
resolution on the order of 3.7k to reach `1e-6`; this is only a cost warning,
not a measured resolution or accepted workaround. Do not relax the material
gate or treat a rejected coarse result as CAD support. The elliptic-paraboloid
solver-to-analytic-CAD path remains blocked until the project can either
produce the required geometric accuracy within declared resource limits or
define and independently validate a representation that preserves the sealed
material without claiming analytic equivalence.

Verification: the full `spin-faster-material-bridge-python` CTest passed
**1/1** (270.71 s); the focused coarse-paraboloid fail-closed test passed
**1/1** (10.85 s); a same-input 32/64/128 resolution diagnostic reproduced
the errors above. These checks validate rejection and the observed convergence
trend only, not analytic CAD export for the paraboloid.
