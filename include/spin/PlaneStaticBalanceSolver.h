#pragma once

#include "spin/PlaneMomentEvaluator.h"

#include <Eigen/Core>
#include <optional>
#include <string>
#include <vector>

namespace spin {

enum class PlaneStaticBalanceStatus {
    Converged,
    IterationLimit,
    DerivativeUnavailable,
    SingularReducedSystem,
    LineSearchExhausted,
    GeometryEvaluationFailed
};

const char* planeStaticBalanceStatusToken(PlaneStaticBalanceStatus) noexcept;

struct PlaneStaticBalanceSettings {
    int maximumIterations = 60;
    int maximumLineSearchSteps = 24;
    double constraintTolerance = 1e-10;
    double backtrackingFactor = 0.5;
    double minimumStep = 1e-12;
    PlaneMomentEvaluatorSettings planeEvaluation;
};

struct PlaneStaticBalanceTraceRow {
    int evaluation = 0;
    int iteration = 0;
    Eigen::Vector2d lambda = Eigen::Vector2d::Zero();
    Eigen::Vector4d coefficients = Eigen::Vector4d::Zero();
    Eigen::Vector2d constraints = Eigen::Vector2d::Zero();
    Eigen::Vector2d normalizedConstraints = Eigen::Vector2d::Zero();
    double normalizedStationarityResidual = 0.0;
    double residualNorm = 0.0;
    double variableMass = 0.0;
    double step = 0.0;
    bool accepted = false;
    bool evaluationSucceeded = false;
    bool derivativeAvailable = false;
    std::string diagnostic;
};

struct PlaneStaticBalanceResult {
    struct FullDomainOptimality {
        // A separate negative-objective-gradient plane, NOT a replacement for
        // the actual accepted coefficients or trace. Zero equality multipliers
        // are a sufficient numerical local condition, not a global certificate.
        Eigen::Vector4d coefficients = Eigen::Vector4d::Zero();
        double minimumGradientPolynomial = 0.0;
    };
    std::optional<FullDomainOptimality> boundaryOptimality;
    PlaneStaticBalanceStatus status = PlaneStaticBalanceStatus::GeometryEvaluationFailed;
    Eigen::Vector2d lambda = Eigen::Vector2d::Zero();
    /** Plane coefficients in dimensionless target coordinates x/L. */
    Eigen::Vector4d coefficients = Eigen::Vector4d::Zero();
    /** Raw target-frame first moments, then divided by constraintScale. */
    Eigen::Vector2d constraints = Eigen::Vector2d::Zero();
    Eigen::Vector2d normalizedConstraints = Eigen::Vector2d::Zero();
    double constraintScale = 1.0;
    double stationarityResidual = 0.0;
    physics::RawMassMoments10 variableMoments;
    PlaneMomentEvaluation evaluation;
    std::vector<PlaneStaticBalanceTraceRow> trace;
    int iterations = 0;
    int acceptedSteps = 0;
    /** True when the paper's plane KKT condition is met within the declared
     * numerical constraint tolerance. This is not an interval-arithmetic proof. */
    bool paperOptimalityConditionsVerified = false;
    std::string diagnostic;
};

/**
 * Solve the paper's 3D static mass-maximization case in its derived vertical
 * plane family beta(x)=1+lambda.x_h>0. This is a dedicated application solver,
 * not a reduced or padded implementation of the general 10D Eq. 14 Newton
 * system. It marks the plane optimality conditions as numerically verified
 * only when both total horizontal first moments satisfy the fixed-input
 * normalized tolerance; this is not an interval-arithmetic certificate.
 */
PlaneStaticBalanceResult solveStaticBalanceByPlane(
    const SpinDesignRequest& request,
    const Eigen::Vector2d& initialLambda,
    const PlaneStaticBalanceSettings& settings = {});

/**
 * Paper static lower-center-of-mass and Pareto plane-family solver. It solves
 * the two horizontal first-moment constraints together with the objective's
 * plane KKT relation from Spin-It Faster Sec. 4.2 / Eq. (13). Coefficients are
 * represented in dimensionless target coordinates x/L and remain four-wide;
 * no artificial ten-dimensional Jacobian is constructed. Maximum-mass and
 * Pareto kappa=1 dispatch to the dedicated vertical-plane specialization.
 * The request explicitly selects legacy mass or Figure 11 uniform-density
 * geometric volume. Strict-full-domain states may use a separately recorded
 * numerical negative-gradient sign witness instead of the interface equation.
 */
PlaneStaticBalanceResult solveStaticBalanceByPlaneFamily(
    const SpinDesignRequest& request,
    const Eigen::Vector4d& initialCoefficients,
    const PlaneStaticBalanceSettings& settings = {});

} // namespace spin
