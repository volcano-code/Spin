#pragma once

#include "physics/RawMassMoments10.h"
#include "spin/QuadricCoefficients.h"
#include "spin/QuadricTriangulator.h"
#include "spin/SpinDesignRequest.h"

#include <Eigen/Core>

#include <string>
#include <vector>


namespace spin
{


enum class QuadricMomentDerivativeMethod
{
    LegacyContinuousImplicitShape,
    DiscreteEllipsoid
};

const char* quadricMomentDerivativeMethodToken(QuadricMomentDerivativeMethod method) noexcept;

struct QuadricMomentEvaluatorSettings
{
    QuadricTriangulationSettings triangulation;
    double minimumInterfaceGradient = 1e-12;
    QuadricMomentDerivativeMethod derivativeMethod =
        QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape;
};


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
    bool derivativeDefined = false;
    bool quadricSolidIsPositiveRegion = false;
    double residualEvaluationSeconds = 0.0;
    double derivativeEvaluationSeconds = 0.0;
    std::string geometryHash;
    std::string topologyHash;
    QuadricMomentDerivativeMethod derivativeMethod =
        QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape;
    std::string derivativeDiagnostic;
    /** v1 geometry/topology hashes above retain their historical bytes.
     * This v2 scientific identity includes method/applicability and source
     * provenance; it is not a build/run attestation. */
    std::string derivativeScientificHash;
    std::string provenanceTopologyHash;
    std::vector<std::string> vertexProvenance;
    int constructedIntersectionVertexCount = 0;
    double maximumConstraintResidual = 0;
    double minimumConstraintSingularRatio = 0;
};


/**
 * Evaluate Omega intersect {beta>0} with exact CGAL mesh Boolean and closed
 * surface moment integration. The default derivative is legacy continuous
 * implicit-shape quadrature. The opt-in discrete mode differentiates actual
 * finite ellipsoid and supported Boolean vertices, then polynomial moments.
 * Unsupported discrete constructions preserve geometry and fail explicitly.
 */
QuadricMomentEvaluation evaluateQuadricMaterialMoments(
    const SpinDesignRequest& request,
    const QuadricCoefficients10& quadricTarget,
    const QuadricMomentEvaluatorSettings& settings =
        QuadricMomentEvaluatorSettings()
);


} // namespace spin
