#pragma once

#include "physics/RawMassMoments10.h"
#include "spin/QuadricCoefficients.h"
#include "spin/QuadricTriangulator.h"
#include "spin/SpinDesignRequest.h"

#include <Eigen/Core>

#include <string>


namespace spin
{


struct QuadricMomentEvaluatorSettings
{
    QuadricTriangulationSettings triangulation;
    double minimumInterfaceGradient = 1e-12;
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
};


/**
 * Evaluate Omega intersect {beta>0} with exact CGAL mesh Boolean and closed
 * surface moment integration.  The Jacobian uses degree-five triangle
 * quadrature of the interface shape derivative and is tangent-projected to
 * the normalized coefficient sphere.
 */
QuadricMomentEvaluation evaluateQuadricMaterialMoments(
    const SpinDesignRequest& request,
    const QuadricCoefficients10& quadricTarget,
    const QuadricMomentEvaluatorSettings& settings =
        QuadricMomentEvaluatorSettings()
);


} // namespace spin
