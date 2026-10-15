# Faster elliptic-cylinder stratum motion — 2026-09-24

## Result

Added a local source-vertex differential for a simple-spectrum elliptic
cylinder, expressed in the cylinder's seven-dimensional coefficient tangent
chart. The API returns a `10 x 7` coefficient basis and one `3 x 7` velocity
matrix per native source vertex. A reduced-coordinate polyhedral moment
transport now consumes these velocities and returns a `10 x 7` raw-moment
Jacobian. The existing full-rank `10 x 10` interface still delegates to the
same generalized transport implementation. Neither API is connected to the
paper's full-rank Equation-14 optimizer as a restricted-stratum solve.

The test retracts each tangent direction back to the exact rank-two,
elliptic-cylinder stratum, holds the native triangulation topology and
truncation radius fixed, and compares all source-vertex velocities against
central finite differences at four step sizes. Relative errors were
`8.09e-8`, `7.28e-9`, `8.09e-10`, and `7.22e-11` for steps
`1e-4`, `3e-5`, `1e-5`, and `3e-6` respectively. The test also verifies
deterministic hashing and fail-closed behavior for a repeated-eigenvalue
circular cylinder and for a paraboloid, which is not yet supported by this
chart. Along the same seven retracted directions, the raw ten-moment
derivative of the closed, triangulated finite source boundary matched
independent signed-tetrahedron finite differences with relative errors
`8.98e-8`, `8.08e-9`, and `8.97e-10` at steps `1e-4`, `3e-5`, and `1e-5`.

A separate exact-containment fixture puts a design cube strictly inside a
larger elliptic-cylinder solid, with the artificial cylinder end caps outside
the design. The unchanged cube boundary is design-born, so its quadric moment
sensitivity must be zero; the reduced `10 x 7` transport returns zero within
`1e-12`. These are source-motion and fixed-provenance transport checks, not a
clipped paper material case or an end-to-end Boolean solve.

## Defect found and corrected

The stratum-tangent differential used half-valued off-diagonal entries for
the derivative of the quadratic matrix with respect to cross coefficients.
That disagreed with `QuadricCoefficients10::quadraticMatrix()`, which stores a
cross coefficient directly in both symmetric off-diagonal entries. A
fixed-cylinder retraction test exposed the mismatch: perturbations left the
cylinder stratum and the native eigenvector signs produced discontinuous
vertex correspondences. The mapping now uses unit off-diagonal derivatives.

The finite paraboloid/cylinder triangulators now canonicalize each simple
eigenvector sign by making its largest-magnitude component nonnegative. This
removes the sign-flip discontinuity on a fixed simple-spectrum branch. It does
not resolve repeated eigenspaces; those remain explicitly unsupported by the
cylinder differential.

## Verification

Build directory: `build/strict-algorithms-exact-current-20260922`.

Rebuilt the affected targets and passed these eleven focused CTests:

- `spin-faster-canonical-integration`
- `spin-faster-canonical-modes`
- `spin-faster-quadric-triangulator`
- `spin-faster-inverse-gauss-map`
- `spin-faster-source-motion`
- `spin-faster-discrete-moments`
- `spin-faster-plane-family`
- `spin-faster-quadric-moments`
- `spin-faster-stratum-tangent`
- `spin-faster-polyhedral-derivative`
- `spin-faster-paper-figures-miniature`

All eleven passed on this working tree. The miniature paper-figure CTest ran
three native Python tests, including the fixed-seed rejected-arrangement event
regression. This is a focused regression result, not a full-suite result and
not a paper-case reproduction. The historical full-rank interface and
reduced-coordinate path both passed their regressions.

## Remaining boundary

- Paraboloid source motion remains unsupported.
- Repeated-eigenvalue cylinders remain unsupported because the principal
  basis is non-unique.
- A clipped cylinder/design Boolean material has not passed end-to-end
  differentiation. A deterministic reconstructed Heart replay localized one
  `Input mesh is not PWN!` failure to its iteration-1 `MINUS` trial (evaluation
  2): design `98V/192F`, source `280V/552F`, and rejected Boolean output
  `235V/466F` with 466 birth faces. The coefficients are regenerated from a
  fixed case-derived seed in the regression test. In the pinned libigl
  implementation, this `false` return
  comes from winding propagation rejecting the stitched arrangement; it is
  not an empty-region certificate. The solver records the typed failure and
  rejects the trial; the same fixed-seed solve subsequently converges and
  passes its existing final verification. This is not proof that the rejected
  arrangement is valid or that clipped-cylinder Booleans work. The exact
  containment fixture does not exercise that branch.
- A standalone replay of that coefficient vector, after applying the Heart
  case's target-frame transform, reproduces the same `280V/552F` quadric
  source size. The source boundary itself has no boundary edges, nonmanifold
  edges or nonmanifold vertices, and no detected geometric self-intersections,
  but it has **15 inconsistently oriented edges**, two face components, and
  fails PWN. The first witness is source edge `(83,79)`; source faces `105` and
  `306` both traverse it in the same direction, and both are classified as
  interface faces (not artificial caps). After exact libigl
  self-intersection remeshing/stitching, the first arrangement witness is edge
  `(55,64)`, also with signed incidence `+2`; both incident faces descend from
  quadric faces (`306` and `307`). This narrows the failure to an invalidly
  oriented quadric triangulation/arrangement input, rather than establishing
  the vendored remesher's documented edge-on-face bug as the cause. No
  automatic orientation repair has been applied.
- A controlled source-tessellation diagnostic on the same fixed coefficient
  vector gave `15`, `12`, `9`, and `0` inconsistent source edges at icosphere
  subdivisions `3`, `4`, `5`, and `6`, respectively. The sub-3/4/5 source
  meshes also failed arrangement PWN; a BFS orientation experiment could make
  them topologically coherent only by reversing `5/4/3` interface facets
  against the local implicit-gradient direction, so it is **not** a repair.
  At subdivision `6`, the source boundary passed the strict/PWN audit and the
  stitched arrangement passed PWN. The libigl `MINUS` call returned
  `2396V/4788F`; all reported edge/vertex topology and winding checks passed,
  but the ordinary audit counted two very small-area faces and its independent
  self-intersection subcheck was incomplete. The exact represented-solid
  importer accepted this libigl result, and its exact self-relation returned
  equality with zero symmetric-difference volume. This establishes neither an
  exact operand-minus relation nor the strict artifact gate, so the result is
  still only a diagnostic candidate.
- The pinned Cherchi probe on the identical subdivision-6 design/source OBJ
  operands returned `2396V/4788F` twice with identical output SHA-256
  `54531cf748d482970c825cc0f1ecc2a8bc19b685858b5c43ff206b47336c61e5`.
  Its mass/first-moment values were close to libigl's, but exact material
  comparison rejected the Cherchi mesh for five self-intersection pairs; the
  cross-backend exact-set comparison therefore remains `UNKNOWN`, not a
  pass. The exported operands and libigl result are under
  `build/heart-boolean-refinement-sub6-20260924/cherchi/`; the replay log is
  `build/heart-boolean-refinement-sub6-20260924/diagnostic.log` (SHA-256
  `9aea3f823c0fd009835c349ba9f920ee58fe5f75d94603a61a0e63d31137a50c`).
- The native full-rank discrete evaluator was then run on the exact same Heart
  request and coefficients at subdivisions `3` and `6`. At `3` it failed
  closed with the same typed arrangement-PWN error. At `6` it rejected the
  Boolean material boundary because its conservative audit found two tiny
  faces and could not complete its own self-intersection subcheck. In a
  separate Python/CGAL exact-region diagnostic, the libigl output imported and
  its self-relation completed as equal with zero symmetric-difference volume;
  the Cherchi output's exact self-relation rejected five intersection pairs.
  Thus subdivision `6` is useful diagnostic evidence, but does not pass the
  production full-rank evaluator's stricter sealing gate or close the
  Cherchi/CGAL cross-check.
- Follow-up exact operand checks changed the interpretation: self-relation of
  the libigl output is not a Boolean-correctness test. The exact-boundary audit
  of the subdivision-6 `MINUS` mesh itself completed with zero exact
  degenerate faces, zero intersecting pairs, valid closure/manifoldness and a
  positive exact component volume. But exact Nef relations against the actual
  exported operands did **not** prove the expected set relations: the candidate
  difference was not certified as a subset of the design, the candidate
  intersection was not certified as a subset of either operand, and the pair
  did not form an exact inclusive-set partition of the design
  (`complete_partition=false`, `missing_whole_content=true`,
  `excess_union_content=true`). The earlier interpretation of the
  `0.8319354948209017` symmetric-difference volume between the candidate
  difference and design as an error was incorrect: exact component-volume
  integration shows that it equals the intended removed intersection volume.
  The candidate difference lies inside the design up to an exact represented
  sliver of `1.4081342796481462e-17`; its overlap with the source quadric is
  `1.758249934452556e-17` by exact-volume arithmetic. Double volume checks give
  `vol(design)=3.3917980843488778`, `vol(MINUS)=2.559862589527979`, and
  `vol(INTERSECT)=0.8319354948209002`, whose sum matches the design to floating
  point roundoff. Thus the Boolean is volume-correct far below the project's
  `1e-6` geometry/moment tolerance, while the exact inclusive-set checker
  detects a nonzero binary64 sliver. All four regions imported individually
  as valid. This points to a represented-coordinate precision mismatch
  (likely exact Boolean intersections converted back to binary64; not yet
  proven), rather than the earlier suspected large Boolean error. The quadric
  operand has two disjoint positive root components (`parents=[-1,-1]`,
  `depths=[0,0]`, signs `[1,1]`), so nested-shell orientation is not the
  explanation. Source inspection of the vendored libigl `mesh_boolean.cpp`
  identifies a likely precision mechanism: the two-input wrapper builds its
  combined vertex matrix using the first operand's scalar type, and the core
  emits to the requested output scalar type. This diagnostic instantiates
  both input and output as `double`, so exact intersection constructions are
  rounded back to binary64 before they become the represented material mesh.
  The observed `1.4e-17`–`1.8e-17` exact-volume slivers are consistent with that
  conversion. That control has now been run with the same binary64 input
  operands and Boolean operation, changing only the output vertex scalar from
  `double` to libigl's exact Epeck scalar. The double-output run has `3,093`
  non-coplanar face-vertex incidences among `14,364`, with maximum perpendicular
  residual `1.9572685690175516e-16`; the exact-output run has `0` such
  incidences and zero measured residual. Both runs report the same `2396V/4788F`
  result topology and the same birth-face count. This confirms that returning
  constructed intersections through a binary64 output matrix is sufficient to
  explain the observed source-plane residual on this fixture. It does **not**
  establish that an exact-coordinate output passes the project's full material
  containment, self-intersection, serialization, or sealing checks, nor does it
  change the production evaluator, whose mesh and independent integration path
  currently consume binary64 coordinates. Converting the Epeck output back to
  binary64 and running the existing strict boundary audit still finds the same
  two degenerate faces and cannot complete its self-intersection subcheck; so
  exact output coordinates alone do not make the current production audit pass.
  The diagnostic utility was rebuilt;
  `spin-faster-quadric-triangulator` passed `1/1`. This focused test is a source
  triangulator regression, not a full Boolean solve or paper case. The repeated
  MINUS/INTERSECT outputs retain hashes `2ecb94e388bc1c53153c0baa3923127b496456f19ee02af446228d59f5b7c385` and
  `473dcfd19f568908420e3690ea02376fe142ca2832a62eaa68c8c288fb32bbdb`.
  A first
  diagnostic invocation accidentally transformed an already target-frame
  design twice and is excluded; the recorded run below uses the same target
  frame for both operands. The diagnostic utility now emits both MINUS and
  INTERSECT meshes for such operand-level checks. Its output hashes are:
  design `94b70085432fe36ca7997763282d235affc4998bbcde7ad9426a9310a837edaa`,
  quadric `5a5626db8eebabbc7d1bc2dcf2f24568934d35fde85ce96a54511815c5219987`,
  MINUS `2ecb94e388bc1c53153c0baa3923127b496456f19ee02af446228d59f5b7c385`,
  INTERSECT `473dcfd19f568908420e3690ea02376fe142ca2832a62eaa68c8c288fb32bbdb`.
  The earlier MINUS export hash differs, but exact Nef comparison between the
  two independently exported candidate meshes completed with `equal=true` and
  exact symmetric-difference volume `0`; therefore the byte difference is
  triangulation/order/serialization-level, not a different represented solid.
  No production audit gate was relaxed, and no candidate material is promoted
  to a sealed artifact. The diagnostic utility was rebuilt and the focused
  `spin-faster-quadric-triangulator` CTest passed `1/1`; this does not convert
  the Boolean candidate into a passing solver result.
- `spin-faster-quadric-triangulator` now includes the fixed Heart negative
  control: it verifies the pinned level-3 source still has its 15-edge winding
  witness and separately verifies that level 6 yields a strictly audited
  source boundary. The focused CTest passed `1/1`; this is a source-mesh
  triangulation regression, not a full Boolean solve or paper case.
- The evaluator now reports this as
  `BOOLEAN_ARRANGEMENT_PWN_FAILURE`, including operation and input/output mesh
  sizes in the solver event. It does not repair or accept the failed output;
  this is better failure evidence, not a geometry-kernel fix.
- A focused integration audit confirms that `FullRankMomentEvaluator` forces
  the discrete derivative even when the source-compatible settings default is
  legacy continuous. Five existing canonical-reference tests passed:
  adapter, derivative chain, solver, off-target yo-yo seed, and replay (`5/5`).
  The derivative-chain relative errors were `5.2047264006447444e-8`,
  `4.6861598988372465e-9`, and `5.1429535093897235e-10` at step sizes
  `1e-4`, `3e-5`, and `1e-5`, respectively. This verifies the fixed-topology
  full-rank Boolean-to-moment chain on its registered method fixtures; the
  off-target yo-yo record is `KNOWN_SOLUTION_VALIDATION`, not an unknown-answer
  paper experiment.
- This does not yet produce a successful full-rank solve on the Heart Boolean
  candidate: the candidate is rejected by the strict binary64 geometry audit
  before Newton can use the derivative. Exact Epeck coordinates remove the
  source-plane residual, but converting them to the current binary64 mesh still
  leaves two degenerate faces and an incomplete self-intersection subcheck.
  The evaluated mapping remains blocked at this geometry boundary; no audit
  threshold was relaxed.
- No general restricted-stratum Newton optimizer, paper Eq. 14 rank-deficient
  extension, Cherchi/CGAL agreement, author-provided cases, or full application
  matrix has been verified.

## Follow-up: clipped elliptic-cylinder Boolean differential

The earlier clipped-material gap above is now closed at the **differential
method-validation** level (not at the optimizer or paper-case level). The
source-motion regression now evaluates a fixed analytic cylinder against the
registered cube design domain, obtains the actual CGAL Boolean boundary and
birth-face provenance, and compares the reduced `10 x 7` Jacobian against
central finite differences of the evaluator's independently recomputed
Boolean material moments. The fixture has 100 constructed Boolean
intersection vertices. Relative errors were `3.2867048569e-8`,
`2.9612999206e-9`, and `3.3045183691e-10` for steps `1e-4`, `3e-5`, and
`1e-5`; topology remained fixed at all tested perturbations. The same focused
`spin-faster-source-motion` CTest passed.

This result verifies only the supported simple-spectrum elliptic-cylinder
stratum, on this fixed-topology synthetic Boolean fixture. It does not connect
the reduced chart to Equation 14 or the full-rank Newton solver; the evaluator's
normal solver path also remains unchanged. Paraboloid source motion, repeated
eigenspaces, cross-stratum transitions, and paper application cases remain
unsupported or unverified.

Project claim remains `PARTIAL_BLOCKED`.
