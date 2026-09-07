#pragma once

#include "spin/QuadricCoefficients.h"
#include "spin/SpinDesignRequest.h"

#include <Eigen/Core>

#include <string>


namespace spin
{


enum class QuadricSurfaceType
{
    Ellipsoid,
    Hyperboloid,
    Paraboloid,
    EllipticCylinder,
    Plane
};


/** Euclidean subtype used by the paper's Pareto-front analysis. */
enum class QuadricSurfaceSubtype
{
    Plane,
    EllipticCylinder,
    Ellipsoid,
    OneSheetHyperboloid,
    TwoSheetHyperboloid,
    Paraboloid
};


const char* quadricSurfaceTypeToken(QuadricSurfaceType type) noexcept;
const char* quadricSurfaceSubtypeToken(QuadricSurfaceSubtype type) noexcept;


QuadricSurfaceSubtype classifyQuadricSurfaceSubtype(
    const QuadricCoefficients10& quadric
);


struct QuadricTriangulationSettings
{
    int icosphereSubdivisions = 3;
    int paraboloidGridResolution = 24;
    int cylinderAngularResolution = 96;
    double truncationRadius = 0.0;
    int maximumRadiusDoublings = 8;
    /** Metadata only: does not alter triangulation or its legacy hash. */
    bool retainFixedNormalProvenance = false;
};


struct QuadricTriangulation
{
    QuadricSurfaceType type = QuadricSurfaceType::Plane;
    TriangleMesh closedBoundary;
    Eigen::VectorXi faceIsCap;
    int interfaceVertexCount = 0;
    int interfaceFaceCount = 0;
    int capFaceCount = 0;
    int boundaryLoopCount = 0;
    int idealNormalEdgeIntersectionCount = 0;
    double truncationRadius = 0.0;
    std::string hash;
    /** Row i is the fixed icosphere normal that generated vertex i.
     * Empty for unsupported types or when metadata was not requested. */
    Eigen::MatrixXd fixedNormals;
};


/**
 * Deterministically triangulate beta=0 for the paper's Boolean operation.
 * Hyperboloids are truncated beyond the design domain and capped.  Plane and
 * paraboloid interfaces are closed with an extruded cap beyond the domain.
 */
QuadricTriangulation triangulateQuadric(
    const QuadricCoefficients10& quadric,
    const Eigen::MatrixXd& designDomainVertices,
    const QuadricTriangulationSettings& settings =
        QuadricTriangulationSettings()
);


} // namespace spin
