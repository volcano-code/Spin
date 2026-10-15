#pragma once

#include "spin/QuadricMomentEvaluator.h"
#include <optional>

namespace spin {

struct PlaneMomentEvaluatorSettings {
    double truncationRadius = 0.0;
};

enum class PlaneMaterialRegion { Unclassified, StrictFull, StrictEmpty };

struct PlaneMomentEvaluation {
    using Jacobian = Eigen::Matrix<double,10,4>;
    /** Actual triangulated variable-material boundary from the same native
     * Boolean evaluation that produced `variableMomentsTarget`, in target
     * coordinates. An empty region remains explicitly empty. */
    TriangleMesh variableMaterialBoundaryTarget;
    physics::RawMassMoments10 variableMomentsTarget;
    physics::RawMassMoments10 totalMomentsTarget;
    Eigen::Vector4d normalizedCoefficients = Eigen::Vector4d::Zero();
    // Derivative with respect to the FOUR supplied, unnormalized target-frame
    // coefficients [a,bx,by,bz], including normalization. No full10 conversion.
    // Absence is unavailable, not a zero gradient. Fixed shell is independent.
    std::optional<Jacobian> coefficientJacobian;
    std::string diagnostic;
    std::string geometryHash;
    std::string scientificHash;
    std::string capSeparationCertificate;
    int interfaceTriangleCount = 0;
    PlaneMaterialRegion materialRegion = PlaneMaterialRegion::Unclassified;
};

// Restricted family beta(x)=a+b.x, material beta>0. Uses native Boolean values;
// differentiates the native finite facets in real arithmetic. IEEE rounding,
// Boolean topology events and unitOrthogonal branch changes are excluded.
// Analytic cubic quadrature transports the moving source facets. This is not a
// full ten-parameter quadric derivative, Newton adapter or manufacture seal.
// Invalid input/geometry/cap coverage throws. Contact returns no Jacobian.
PlaneMomentEvaluation evaluatePlaneMaterialMoments(
    const SpinDesignRequest& request,
    const Eigen::Vector4d& coefficientsTarget,
    const PlaneMomentEvaluatorSettings& settings = {});

} // namespace spin
