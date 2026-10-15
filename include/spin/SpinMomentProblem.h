#pragma once

#include "spin/SpinDesignRequest.h"

#include <Eigen/Core>

#include <vector>


namespace spin
{


struct SpinMomentEvaluation
{
    double objective = 0.0;
    Eigen::Matrix<double, 10, 1> objectiveGradient =
        Eigen::Matrix<double, 10, 1>::Zero();
    Eigen::Matrix<double, 10, 10> objectiveHessian =
        Eigen::Matrix<double, 10, 10>::Zero();
    Eigen::VectorXd constraints;
    Eigen::MatrixXd constraintJacobian;
    std::vector<Eigen::Matrix<double, 10, 10>> constraintHessians;
};


/**
 * Paper objective and equality-constraint system in canonical raw moments.
 *
 * Input ordering is `(v,l1,l2,l3,Q11,Q22,Q33,Q23,Q13,Q12)` and must already
 * be expressed in the request's target frame.  Derivatives are analytic with
 * respect to those ten coordinates.
 */
class SpinMomentProblem
{
public:
    using Vector10 = Eigen::Matrix<double, 10, 1>;

    static SpinMomentProblem staticBalance(StaticObjective objective, double massScale = 1.0);
    static SpinMomentProblem staticBalancePareto(double massWeight, double massScale = 1.0);
    static SpinMomentProblem buoyantStability(
        double targetMass,
        const Eigen::Vector2d& horizontalCenterOfBuoyancy
    );
    static SpinMomentProblem fidgetSpinner();
    static SpinMomentProblem yoYo();
    static SpinMomentProblem spinningTop(double centerOfMassWeight);
    static SpinMomentProblem spinningTop(double centerOfMassWeight, SpinningTopObjective objective);

    SpinMode mode() const noexcept;
    StaticObjective staticObjective() const noexcept;
    double centerOfMassWeight() const noexcept;
    double staticParetoWeight() const noexcept;

    SpinMomentEvaluation evaluate(const Vector10& rawMoments) const;

private:
    SpinMomentProblem(
        SpinMode mode,
        StaticObjective staticObjective,
        double centerOfMassWeight,
        double staticParetoWeight = 0.0,
        double targetMass = 0.0,
        const Eigen::Vector2d& horizontalCenterOfBuoyancy =
            Eigen::Vector2d::Zero()
    );

    SpinMode mode_ = SpinMode::StaticBalance;
    StaticObjective staticObjective_ = StaticObjective::MaximumMass;
    double centerOfMassWeight_ = 0.0;
    double staticParetoWeight_ = 0.0;
    double staticMassScale_ = 1.0;
    double targetMass_ = 0.0;
    Eigen::Vector2d horizontalCenterOfBuoyancy_ = Eigen::Vector2d::Zero();
    SpinningTopObjective spinningTopObjective_ = SpinningTopObjective::InertiaAndLowering;
};


/** Height of the lower metacenter in target-frame coordinates. */
double buoyantMetacenterHeight(const BuoyantStabilityData& buoyancy);


/** Positive exactly when the center of mass lies below the metacenter. */
double buoyantStabilityMargin(
    const Eigen::Matrix<double, 10, 1>& rawMoments,
    const BuoyantStabilityData& buoyancy
);


} // namespace spin
