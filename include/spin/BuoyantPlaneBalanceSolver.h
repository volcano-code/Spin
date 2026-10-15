#pragma once

#include "spin/FasterTheoremConditions.h"
#include "spin/PlaneMomentEvaluator.h"

#include <Eigen/Core>
#include <optional>
#include <string>
#include <vector>

namespace spin {

enum class BuoyantPlaneBalanceStatus {
    Converged,
    IterationLimit,
    DerivativeUnavailable,
    SingularReducedSystem,
    LineSearchExhausted,
    GeometryEvaluationFailed,
    OptimalityAuditFailed,
    InvalidRequest
};

const char* buoyantPlaneBalanceStatusToken(BuoyantPlaneBalanceStatus) noexcept;

struct BuoyantPlaneBalanceSettings {
    int maximumIterations = 60;
    int maximumLineSearchSteps = 24;
    double constraintTolerance = 1e-10;
    double backtrackingFactor = 0.5;
    double minimumStep = 1e-12;
    PlaneMomentEvaluatorSettings planeEvaluation;
};

struct BuoyantPlaneBalanceTraceRow {
    int evaluation = 0;
    int iteration = 0;
    Eigen::Vector4d planeCoefficients = Eigen::Vector4d::Zero();
    Eigen::Vector3d constraints = Eigen::Vector3d::Zero();
    Eigen::Vector3d normalizedConstraints = Eigen::Vector3d::Zero();
    double residualNorm = 0.0;
    double step = 0.0;
    bool accepted = false;
    bool evaluationSucceeded = false;
    bool derivativeAvailable = false;
    std::string diagnostic;
};

struct BuoyantPlaneBalanceResult {
    BuoyantPlaneBalanceStatus status =
        BuoyantPlaneBalanceStatus::GeometryEvaluationFailed;
    /** Unnormalized target-frame coefficients [a,bx,by,bz], with bz fixed at -1. */
    Eigen::Vector4d planeCoefficients = Eigen::Vector4d(0.0, 0.0, 0.0, -1.0);
    /** Total raw mass moments, including the immutable fixed shell. */
    physics::RawMassMoments10 totalMoments;
    physics::RawMassMoments10 variableMoments;
    PlaneMomentEvaluation evaluation;
    Eigen::Vector3d constraints = Eigen::Vector3d::Zero();
    Eigen::Vector3d normalizedConstraints = Eigen::Vector3d::Zero();
    Eigen::Vector3d constraintScales = Eigen::Vector3d::Ones();
    double targetTotalMass = 0.0;
    double centerOfMassHeight = 0.0;
    double metacenterHeight = 0.0;
    double stabilityMargin = 0.0;
    Eigen::VectorXd equalityMultipliers;
    double multiplierMu = 0.0;
    std::optional<FasterTheoremOneNumericalAudit> theoremOneNumericalAudit;
    /** All three prescribed-waterline equilibrium constraints were met, with
     * a full-rank local plane-family Jacobian. This is not a global proof. */
    bool equilibriumConstraintsVerified = false;
    bool eq7NumericalConditionsVerified = false;
    bool stableAtPrescribedWaterline = false;
    std::vector<BuoyantPlaneBalanceTraceRow> trace;
    int iterations = 0;
    int acceptedSteps = 0;
    std::string diagnostic;
};

/**
 * Solve the prescribed-waterline feasibility problem from Faster §4.3 in the
 * restricted planar family beta(x)=a+bx*x+by*y-z>0. The target total mass is
 * rho_liquid * displacedLiquidVolume; the two horizontal first moments are
 * constrained to the center of buoyancy. The returned metacenter comparison
 * is a numerical check for the supplied fixed waterline, not a dynamic-fluid
 * simulation or a continuous-domain infeasibility certificate.
 */
BuoyantPlaneBalanceResult solveBuoyantStabilityByPlane(
    const SpinDesignRequest& request,
    const Eigen::Vector3d& initialParameters,
    const BuoyantPlaneBalanceSettings& settings = {});

} // namespace spin
