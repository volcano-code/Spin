#pragma once

#include "spin/QuadricCoefficients.h"
#include "spin/QuadricCoefficientAuthority.h"
#include "spin/RankTwoSource.h"
#include "spin/SpinDesignRequest.h"

#include <Eigen/Core>

#include <array>
#include <string>
#include <vector>


namespace spin
{


enum class QuadricSurfaceType
{
    Ellipsoid,
    Hyperboloid,
    Paraboloid,
    EllipticCylinder,
    Plane,
    Cone,
    HyperbolicCylinder,
    ParabolicCylinder,
    IntersectingPlanes
};


/** Euclidean subtype used by the paper's Pareto-front analysis and by
 * explicit classification of degenerate inputs. Families without a supported
 * source triangulation remain fail-closed; supported cylinder strata are
 * routed through their dedicated triangulators. */
enum class QuadricSurfaceSubtype
{
    Plane,
    EllipticCylinder,
    Ellipsoid,
    OneSheetHyperboloid,
    TwoSheetHyperboloid,
    Paraboloid,
    Cone,
    HyperbolicCylinder,
    ParabolicCylinder,
    IntersectingPlanes
};


const char* quadricSurfaceTypeToken(QuadricSurfaceType type) noexcept;
const char* quadricSurfaceSubtypeToken(QuadricSurfaceSubtype type) noexcept;


QuadricSurfaceSubtype classifyQuadricSurfaceSubtype(
    const QuadricCoefficients10& quadric
);


struct QuadricTriangulationSettings
{
    int icosphereSubdivisions = 3;
    /** Optional bounded fallback for ambiguous hyperboloid cap topology.
     * -1 preserves the requested resolution exactly. When enabled, failed
     * levels are retried one level at a time through this inclusive maximum.
     * This is an implementation policy, not a paper-prescribed solver step. */
    int maximumIcosphereSubdivisions = -1;
    int paraboloidGridResolution = 24;
    int cylinderAngularResolution = 96;
    int hyperbolicCylinderGridResolution = 24;
    int parabolicCylinderGridResolution = 24;
    double truncationRadius = 0.0;
    int maximumRadiusDoublings = 16;
    /** Metadata only: does not alter triangulation or its legacy hash. */
    bool retainFixedNormalProvenance = false;
    /** Explicit paper finite/ideal triangle construction. Legacy analytic
     * quadric/radius clipping remains the default. Unsupported normal-face
     * connectivity fails closed; no automatic normal-sphere refinement. */
    bool paperProjectiveHyperboloid = false;
    /** Explicit algebraic finite-point policy. Currently admitted only for
     * full-rank source meshes; rank-deficient source charts are not silently
     * substituted. Geometry and request hashes distinguish this policy. */
    bool algebraicFiniteMap = false;
    /** Explicit rational-enclosed source policy v1. Coordinates are generated
     * before binary64 assignment; the accompanying intervals enclose analytic
     * source points with absolute width <=2^-384. Unsupported charts reject.
     * The display boundary is derived and is not the authority operand. */
    bool rationalSourceCoordinates = false;
    /** Enclosed-dyadic-1024 v1; opt-in, requires rational source coordinates.
     * Cap centres remain exact means of the same rational ring. */
    bool compactRationalSourceCoordinates = false;
    /** Opt-in bounded finite fixed-normal paraboloid patch. Does not implement
     * ideal/equator faces. Requires certified patch coverage and exterior caps;
     * otherwise rejects without grid fallback or normal refinement. */
    bool fixedNormalParaboloid = false;
    /** Independent exact rank-two arrangement profile; never silently routed
     * through a legacy triangulator. Fixed radius and normal grid required. */
    bool rankTwoArrangementSource = false;
    RankTwoSourceLimits rankTwoLimits;
};

/** Deterministic local tangent frame used to triangulate a plane patch.
 * `chart` is 0 for an x-axis reference and 1 for a y-axis reference. The chart
 * switch occurs at |normal.x|=0.9. The triangulator may construct the frame on
 * the seam, but differential consumers must require `chartStable`.
 */
struct QuadricPlaneFrame
{
    Eigen::Vector3d first = Eigen::Vector3d::Zero();
    Eigen::Vector3d second = Eigen::Vector3d::Zero();
    int chart = -1;
    bool chartStable = false;
    std::string diagnostic;
};

QuadricPlaneFrame makeQuadricPlaneFrame(const Eigen::Vector3d& unitNormal);


struct QuadricTriangulation
{
    QuadricSurfaceType type = QuadricSurfaceType::Plane;
    TriangleMesh closedBoundary;
    /** Authority coordinate payload for rationalSourceCoordinates only.
     * Same indexing as closedBoundary; bounds refer to analytic source points,
     * while exactVertices define the actual rational polyhedral operand. */
    std::vector<std::array<std::string,3>> exactVertices;
    std::vector<std::array<std::string,3>> coordinateLower;
    std::vector<std::array<std::string,3>> coordinateUpper;
    Eigen::VectorXi faceIsCap;
    int interfaceVertexCount = 0;
    int interfaceFaceCount = 0;
    int capFaceCount = 0;
    int boundaryLoopCount = 0;
    int idealNormalEdgeIntersectionCount = 0;
    int initialIcosphereSubdivisions = -1;
    int effectiveIcosphereSubdivisions = -1;
    double truncationRadius = 0.0;
    std::string hash;
    /** Native content/recipe binding for an explicitly rational source. It
     * includes coefficients and provenance, unlike the legacy geometry hash.
     * Enables same-construction derivative reuse; never a material seal. */
    std::string rationalConstructionIdentity;
    bool compactRationalSourceCoordinates = false;
    /** Row i is the fixed icosphere normal that generated vertex i.
     * Empty for unsupported types or when metadata was not requested. */
    Eigen::MatrixXd fixedNormals;
    /** For capped hyperboloids, normal-sphere edge identity for each output
     * vertex. {-1,-1} denotes an original fixed-normal vertex or cap center.
     * Metadata only; excluded from the legacy geometry hash. */
    std::vector<std::array<int, 2>> sourceNormalEdgeByVertex;
    /** For capped hyperboloids, fixed-normal vertex identity for each output
     * vertex. -1 denotes an edge crossing or cap center. Metadata only. */
    std::vector<int> sourceNormalVertexByVertex;
    bool paperProjectiveHyperboloid = false;
    bool fixedNormalParaboloid = false;
    /** Paper ray identity: {finite-normal, ideal-edge-first, ideal-edge-second,
     * ideal-root-ordinal}. All -1 for finite normal vertices and cap centres.
     * Remote planar-cap vertex: {-2, interface-ring-vertex, loop, ordinal}.
     * Finite root-interval midpoint: {-3, edge-first, edge-second, 0}; its
     * normal row is coefficient-dependent and differentiated as such.
     * Always populated in paper mode, never interpreted as legacy crossings. */
    std::vector<std::array<int,4>> sourceProjectiveRayByVertex;
};


/**
 * Deterministically triangulate beta=0 for the paper's Boolean operation.
 * Hyperboloids and cylinders are truncated beyond the design domain and
 * capped. Plane, paraboloid, and parabolic-cylinder interfaces are closed
 * with finite extruded caps beyond the domain.
 */
QuadricTriangulation triangulateQuadric(
    const QuadricCoefficients10& quadric,
    const Eigen::MatrixXd& designDomainVertices,
    const QuadricTriangulationSettings& settings =
        QuadricTriangulationSettings()
);

/** Recompute the native binding of a supplied rational construction. Empty,
 * shape-invalid, stale or mutated operands fail closed. This checks content
 * identity, not manufacturing qualification or analytic-root exactness. */
bool rationalQuadricConstructionIdentityMatches(
    const QuadricCoefficients10&, const QuadricTriangulation&,
    const Eigen::MatrixXd& designDomainVertices,
    const QuadricTriangulationSettings&);

/** Explicit MPFR1024 authority path; requires rational source coordinates
 * and the existing fixed-normal full-rank chart. Unsupported profiles reject. */
QuadricTriangulation triangulateQuadric(
    const QuadricCoefficientAuthority10&, const Eigen::MatrixXd&,
    const QuadricTriangulationSettings&);
bool rationalQuadricConstructionIdentityMatches(
    const QuadricCoefficientAuthority10&, const QuadricTriangulation&,
    const Eigen::MatrixXd&, const QuadricTriangulationSettings&);


} // namespace spin
