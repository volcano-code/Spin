#pragma once

#include "spin/SpinDesignRequest.h"

#include <Eigen/Core>


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

    static SpinMomentProblem staticBalance(StaticObjective objective);
    static SpinMomentProblem fidgetSpinner();
    static SpinMomentProblem yoYo();
    static SpinMomentProblem spinningTop(double centerOfMassWeight);

    SpinMode mode() const noexcept;
    StaticObjective staticObjective() const noexcept;
    double centerOfMassWeight() const noexcept;

    SpinMomentEvaluation evaluate(const Vector10& rawMoments) const;

private:
    SpinMomentProblem(
        SpinMode mode,
        StaticObjective staticObjective,
        double centerOfMassWeight
    );

    SpinMode mode_ = SpinMode::StaticBalance;
    StaticObjective staticObjective_ = StaticObjective::MaximumMass;
    double centerOfMassWeight_ = 0.0;
};


} // namespace spin
