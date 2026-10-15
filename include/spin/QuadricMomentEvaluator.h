#pragma once

#include "physics/RawMassMoments10.h"
#include "spin/QuadricCoefficients.h"
#include "spin/QuadricDiscreteMomentDifferential.h"
#include "spin/QuadricTriangulator.h"
#include "spin/SpinDesignRequest.h"
#include "spin/PolyhedralMomentDifferential.h"

#include <Eigen/Core>

#include <string>
#include <array>
#include <vector>
#include <functional>
#include <memory>

namespace manufacturing { struct BoundaryAudit; }
namespace verification { class ExactMaterialRegion; }

namespace spin
{


enum class QuadricMomentDerivativeMethod
{
    LegacyContinuousImplicitShape,
    DiscreteFullRankQuadric,
    /** Fixed-stratum derivative for supported rank-deficient quadrics. This
     * is returned separately and is never presented as an ambient 10-D
     * Jacobian for the Eq. 14 Newton solver. */
    DiscreteReducedStratum,
    /** Source-compatible legacy spelling; now covers supported ellipsoids
     * and hyperboloids. */
    DiscreteEllipsoid = DiscreteFullRankQuadric
};

const char* quadricMomentDerivativeMethodToken(QuadricMomentDerivativeMethod method) noexcept;

/** Content-bound identity for the Boolean implementation used by the native
 * quadric/design-domain evaluator.  This identifies libigl's CGAL mesh
 * Boolean call and the selected libigl source checkout/tree; it does not
 * identify or claim equivalence with the separate Cherchi reference backend.
 */
std::string quadricBooleanBackendIdentity();

/** Classify which sign of q is inside a closed, outward-oriented quadric
 * source using its non-cap interface normals. A volume centroid is not a
 * reliable interior witness for non-convex solids. Throws when interface
 * orientation does not give one unambiguous sign. */
bool quadricSolidInteriorIsPositiveRegion(
    const QuadricCoefficients10& quadric,
    const QuadricTriangulation& triangulation);

/** Exact strict margin for the INACTIVE centroid rule and uniform oriented
 * interface fallback on the supplied indexed source and coefficients.
 * Establishes an open joint mesh/coefficient neighborhood of this selector
 * branch in real arithmetic. It does not certify the coefficient-to-source
 * map, its finite topology/radius branch, Boolean motion or any Jacobian.
 * Callers must separately require a regular matching source differential.
 * Equality at the availability threshold, mixed/ambiguous interface signs
 * and a selected side inconsistent with the exact fallback return false.
 */
bool quadricSelectorHasStrictInactiveCentroidMargin(
    const QuadricCoefficients10&, const QuadricTriangulation&,
    bool selectedPositiveRegion);

enum class QuadricEmptyRegionPolicy
{
    RejectEmpty,
    CertifyStrictlySeparatedV1
};

/** Policy for artificial caps on unbounded quadrics.  The default remains
 * fail-closed; the opt-in policy grows the effective truncation radius until
 * the exact convex-hull separation certificate succeeds or the registered
 * radius-growth budget is exhausted. */
enum class QuadricCapRadiusPolicy
{
    RejectUncertified,
    GrowToCertifiedV1
};

enum class QuadricMaterialRegion
{
    Nonempty,
    CertifiedEmpty
};

struct QuadricMomentEvaluatorSettings
{
    /** Optional diagnostic-only stage callback. Not a scientific setting,
     * not hashed, and never grants a derivative or convergence status. */
    std::function<void(const std::string&)> diagnosticStageObserver;
    QuadricTriangulationSettings triangulation;
    double minimumInterfaceGradient = 1e-12;
    QuadricMomentDerivativeMethod derivativeMethod =
        QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape;
    QuadricEmptyRegionPolicy emptyRegionPolicy = QuadricEmptyRegionPolicy::RejectEmpty;
    QuadricCapRadiusPolicy capRadiusPolicy = QuadricCapRadiusPolicy::RejectUncertified;
    /** Explicit exact indexed-simplex self-intersection predicates on ORIGINAL
     * binary64 coordinates. Legacy area/volume/topology gates remain unchanged.
     * Default false preserves historical floating shared-corner heuristics. */
    bool exactBinary64SelfIntersection = false;
    /** Opt-in for canonical full-rank source and input boundary admission. Requires
     * the complete exact indexed-simplex audit of original binary64 source
     * coordinates, including positive component winding and nonzero areas.
     * Input outer/design boundaries use the same exact predicates. Final material
     * and differential conditioning gates do not change. This is not a seal.
     * False preserves all historical admission and identity behavior. */
    bool exactBinary64SourceBoundary = false;
    /** Explicit solver-material proof on ORIGINAL binary64 triangles. Replaces
     * only the heuristic geometric boundary classification with a complete
     * signed-shell exact material proof (including correctly inward cavities).
     * No repair, manufacturing qualification or
     * relaxation of motion/feature conditioning, mass or equation tolerances.
     * Default false retains historical value/identity behavior. */
    bool exactBinary64MaterialBoundary = false;
};

manufacturing::BoundaryAudit auditQuadricMaterialBoundary(
    const TriangleMesh&, bool exactBinary64SelfIntersection);
manufacturing::BoundaryAudit auditQuadricMaterialBoundary(
    const TriangleMesh&, bool exactBinary64SelfIntersection, bool exactBinary64MaterialBoundary);

/** Effective request identity including evaluator controls.
 * The historical control tuple (including an eight-step cap-growth budget)
 * preserves its request hash byte for byte; later behavior-affecting defaults
 * receive a distinct identity. */
std::string computeQuadricEvaluatorRequestHash(
    const SpinDesignRequest& request,
    const QuadricMomentEvaluatorSettings& settings);

void validateQuadricMomentEvaluatorSettings(const QuadricMomentEvaluatorSettings& settings);


struct QuadricMomentEvaluation
{
    TriangleMesh materialBoundary;
    TriangleMesh materialBoundaryTarget;
    QuadricTriangulation quadricTriangulation;
    physics::RawMassMoments10 variableMomentsTarget;
    physics::RawMassMoments10 totalMomentsTarget;
    Eigen::Matrix<double, 10, 10> momentJacobian =
        Eigen::Matrix<double, 10, 10>::Zero();
    int interfaceTriangleCount = 0;
    /** True only when `momentJacobian` is the ambient 10 x 10 differential.
     * A successful reduced-stratum result is reported exclusively through
     * `stratumDifferential.defined`. */
    bool derivativeDefined = false;
    bool quadricSolidIsPositiveRegion = false;
    double residualEvaluationSeconds = 0.0;
    double derivativeEvaluationSeconds = 0.0;
    std::string geometryHash;
    std::string topologyHash;
    QuadricMomentDerivativeMethod derivativeMethod =
        QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape;
    std::string derivativeDiagnostic;
    /** Populated only for DiscreteReducedStratum. The reduced Jacobian and
     * coefficient basis keep their true stratum dimension; consumers must not
     * treat this result as a full-rank Eq. 14 Jacobian. */
    QuadricStratumDiscreteMomentDifferential stratumDifferential;
    /** v1 geometry/topology hashes above retain their historical bytes.
     * This v2 scientific identity includes method/applicability and source
     * provenance; it is not a build/run attestation. */
    std::string derivativeScientificHash;
    std::string provenanceTopologyHash;
    /** Actual value-side Boolean backend identity, included in scientific
     * evaluation hashes. Empty is never a valid completed evaluation. */
    std::string booleanBackendIdentity;
    std::vector<std::string> vertexProvenance;
    /** Boolean output face -> [design faces, quadric-source faces]. Retained so
     * reduced stratum differentials can audit and replay the same operation. */
    Eigen::VectorXi materialBirthFaces;
    /** SAME operation's rational target coordinates before binary64 assignment.
     * Populated only by the explicit exact source+material full-rank profile.
     * This is motion provenance, not an independent material acceptance seal. */
    std::vector<std::array<std::string,3>> materialExactVerticesTarget;
    /** Rational-source profile: actual input operand and density-weighted
     * variable moments before the final numerical consumer conversion. */
    std::vector<std::array<std::string,3>> designExactVerticesTarget;
    std::vector<std::string> variableMomentsExactTarget;
    int constructedIntersectionVertexCount = 0;
    double maximumConstraintResidual = 0;
    double minimumConstraintSingularRatio = 0;
    /** Empty variable region is not a final fixed-shell fabrication mesh. */
    QuadricMaterialRegion materialRegion = QuadricMaterialRegion::Nonempty;
    std::string emptyRegionCertificate;
};

/** Scoped, thread-local diagnostic tap on the actual production construction.
 * Read-only snapshots only; never changes admission or any scientific identity.
 * The observer must not recursively evaluate geometry. Restoration is automatic
 * even when the evaluator rejects a boundary. Existing settings ABI is unchanged.
 */
class ScopedQuadricGeometryObserver
{
public:
    using Observer = std::function<void(const std::string&, const QuadricMomentEvaluation&)>;
    explicit ScopedQuadricGeometryObserver(Observer observer);
    ~ScopedQuadricGeometryObserver();
    ScopedQuadricGeometryObserver(const ScopedQuadricGeometryObserver&) = delete;
    ScopedQuadricGeometryObserver& operator=(const ScopedQuadricGeometryObserver&) = delete;
private:
    Observer previous_;
};


/**
 * Evaluate Omega intersect {beta>0} with exact CGAL mesh Boolean and closed
 * surface moment integration. The default derivative is legacy continuous
 * implicit-shape quadrature. The opt-in discrete mode differentiates actual
 * supported full-rank quadric and Boolean vertices, then polynomial moments.
 * `DiscreteReducedStratum` instead returns a separate lower-dimensional
 * differential for supported planes, elliptic/hyperbolic/parabolic cylinders,
 * and paraboloids;
 * it is deliberately not represented as the ambient Jacobian consumed by
 * the Equation-14 solver. Unsupported discrete constructions preserve
 * geometry and fail explicitly.
 */
QuadricMomentEvaluation evaluateQuadricMaterialMoments(
    const SpinDesignRequest& request,
    const QuadricCoefficients10& quadricTarget,
    const QuadricMomentEvaluatorSettings& settings =
        QuadricMomentEvaluatorSettings()
);

/** Authority-driven same-operation value path. Exact material coordinates
 * and density-weighted rational moments are retained. Until authority source
 * motion is independently validated, derivativeDefined is false; this entry
 * must not be treated as a Newton-ready evaluation. */
QuadricMomentEvaluation evaluateQuadricMaterialMoments(
    const SpinDesignRequest&, const QuadricCoefficientAuthority10&,
    const QuadricMomentEvaluatorSettings&);

/** Explicit same-operation authority derivative attachment. This retains the
 * represented rational Jacobian and source enclosures. Local derivative
 * applicability is not a neighborhood/topology, Newton convergence or final
 * material acceptance certificate. Legacy result layouts remain unchanged. */
struct QuadricAuthorityMomentAttachment {
    bool defined=false;
    std::string coefficientPolicy,sourceConstructionIdentity,sourceMotionScientificHash;
    unsigned precisionBits=0;
    QuadricCoefficientAuthority10::RationalVector coefficientAuthority;
    RationalPolyhedralMomentDifferential differential;
    std::vector<RationalSourceVertexJacobian> sourceDerivativeLower,sourceDerivativeUpper;
    std::string scientificHash;
};
QuadricMomentEvaluation evaluateQuadricMaterialMoments(
    const SpinDesignRequest&,const QuadricCoefficientAuthority10&,
    const QuadricMomentEvaluatorSettings&,QuadricAuthorityMomentAttachment&);

/** Geometry-owned fixed material is OUTER minus DESIGN, independently of
 * the caller's numeric fixedShellMoments compatibility field. The factory
 * derives that view; numeric values cannot grant fixed-material authority.
 * This is input geometry/mass admission, not a fabrication or solver seal.
 */
class QuadricAuthorityMaterialContext {
public:
    using Affine = std::array<std::array<std::string,4>,3>;
    const SpinDesignRequest& request() const {return request_;}
    const std::array<std::string,10>& fixedMomentsTarget() const {return fixedMomentsTarget_;}
    const Affine& modelToTarget() const {return modelToTarget_;}
    const Affine& targetToModel() const {return targetToModel_;}
    const std::string& identity() const {return identity_;}
    double fixedDensity() const {return fixedDensity_;}
private:
    QuadricAuthorityMaterialContext()=default;
    SpinDesignRequest request_;
    double fixedDensity_=0;
    std::array<std::string,10> fixedMomentsTarget_;
    Affine modelToTarget_,targetToModel_;
    std::shared_ptr<const verification::ExactMaterialRegion> outerRegion_,designRegion_;
    std::string identity_;
    friend QuadricAuthorityMaterialContext makeQuadricAuthorityMaterialContext(const SpinDesignRequest&,double);
};
QuadricAuthorityMaterialContext makeQuadricAuthorityMaterialContext(
    const SpinDesignRequest& geometricRequest,double fixedDensity);

/** Same actual Boolean/Jacobian plus geometry-owned exact TOTAL moments.
 * Compatibility material.totalMomentsTarget remains a legacy double view;
 * new high-precision consumers must use totalMomentsTarget below.
 */
struct QuadricAuthorityTotalMomentEvaluation {
    QuadricMomentEvaluation material;
    QuadricAuthorityMomentAttachment derivative;
    std::array<std::string,10> totalMomentsTarget;
    std::string contextIdentity,scientificHash;
};
QuadricAuthorityTotalMomentEvaluation evaluateQuadricMaterialMoments(
    const QuadricAuthorityMaterialContext&,const QuadricCoefficientAuthority10&,
    const QuadricMomentEvaluatorSettings&);


} // namespace spin
