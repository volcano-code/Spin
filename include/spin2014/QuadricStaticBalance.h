#pragma once
#include "spin2014/QuadricHollowing.h"
#include "spin/PlaneStaticBalanceSolver.h"

namespace spin2014 {
/** Fresh local squared-height verification, not a material, global optimum or
 * support certificate. The ordinary plane candidate's pass flag is not read.
 * Coordinates are target input units; coefficients act on x/targetLength.
 */
struct StaticSquaredHeightVerification {
    spin::SpinMomentEvaluation evaluation;
    Eigen::Vector2d equalityMultipliers=Eigen::Vector2d::Zero();
    double multiplier=0, stationarityResidual=0;
    bool degenerateObjective=false;
    bool squaredHeightConditionsVerified=false;
    bool boundaryWitnessVerified=false;
    Eigen::Vector4d boundaryWitnessCoefficients=Eigen::Vector4d::Zero();
    double boundaryMultiplier=0, boundaryGradientMinimum=0;
    std::string diagnostic;
};
StaticSquaredHeightVerification verifySpin2014StaticPlane(
    const Eigen::Matrix<double,10,1>& totalMomentsTarget,
    const Eigen::Vector4d& dimensionlessNormalizedPlane,double targetLength,
    const QuadricObjectiveSettings&,double constraintScale=1.,double tolerance=1e-10);

struct QuadricStaticMaterialRun {
    std::string objectiveProfile, objectiveIdentity;
    QuadricObjectiveSettings objectiveSettings;
    spin::SpinDesignRequest measuredRequest;
    double fixedMaterialDensity=0, targetLength=0;
    std::string initialization;
    Eigen::Vector4d initialCoefficients=Eigen::Vector4d::Zero();
    spin::PlaneStaticBalanceResult candidate;
    StaticSquaredHeightVerification squaredHeight;
};
QuadricStaticMaterialRun solveSpin2014StaticQuadricHollowing(
    const spin::SpinDesignRequest&,const QuadricObjectiveSettings&,double fixedDensity,
    const std::optional<Eigen::Vector4d>& initialDimensionlessCoefficients=std::nullopt,
    const spin::PlaneStaticBalanceSettings& = {});
} // namespace spin2014
