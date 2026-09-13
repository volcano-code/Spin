# Task5z controller protocol (fixed before region-output inspection)

Scope: actual finite-polyhedral material regions and outer/carved clearance.
This is KNOWN_SOLUTION_VALIDATION, never original-paper optimization, a seal,
full wall/web thickness, a feasibility certificate or physical validation.

Reuse the immutable Task5y fixed four-input protocol without changing geometry,
state, density, objective, tolerance or pose. Source outer cube [.3125,.6875]^3,
design [.425,.575]^3, source maximumDepth4, domain4..4, cage2, fine3, all26BBW
handles, density1.17. Full affine scale1.02 and translation(.011,-.007,.019),
angle.23 in top cases, zero objective weights and known stationary initialization.
Relaxed beta .25/.75/.5 give empty/full/tie; the separate static .4/.6 control
has its fixed axis derived from the relaxed state and fails after binary rounding.
This known-solution choice is inherited unchanged, not made to fit Task5z output.

Every actual binary input enters the new native gate at minimumOuterWall=.05.
The full case also enters with .13 (too large) and Nef complexity cap1. These
failures must remain explicit; no geometry repair or tolerance changes allowed.
Native module tests separately cover two disjoint cavities and exactly
representable threshold equality versus nextafter-above. The controller does
not relabel its four cases as those distinct native fixtures.

The independent Python path imports no native evaluator/checker. It reconstructs
full BBW interpolation and actual selected oriented boundary chains, integrates
signed rational tetrahedra, checks the fixed-frame constraints/inertia, verifies
input and selected-mesh byte identities, and recomputes all outer-to-carved
triangle-feature distances as fractions of the actual binary64 coordinates.
It verifies the reported witness pair realizes the minimum, not just the scalar.
The convex fixture's represented primitive containment is checked against all
outer triangle halfspaces. This limited independent fixture checker is not a
general replacement for production Nef set algebra or arbitrary-shape topology.

The producer's component/partition evidence is independently compared with this
fixed connected cube / one-cavity protocol; primitive-union oriented cancellation
and signed volume subtraction must agree. Empty cavity requires a real explicit
empty operand, no fabricated distance, contact, witness or threshold result.
The product must retain the original actual binary object; it never reoptimizes.

Verification order: compile actual exporter; preserve preview and any failures;
freeze effective native/compiler/dependency/checker inputs; two fresh processes
under24GiB/60min guards; byte-compare scientific JSONL excluding resource times;
independent checks and actual-output corruption negatives; verify frozen inputs
again. Follow with related registered tests and a fresh scoped code review.
Raw logs and failure attempts are retained. Main and foreign Faster work remain
untouched. No package-wide completion or paper count can be granted by this file.
