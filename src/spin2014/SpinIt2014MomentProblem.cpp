#include "spin2014/SpinIt2014MomentProblem.h"

#include <cmath>
#include <stdexcept>


namespace spin2014
{


SpinIt2014MomentProblem::SpinIt2014MomentProblem(
    spin::SpinMode mode,
    double centerWeight
)
    : mode_(mode), centerWeight_(centerWeight)
{
    if (!std::isfinite(centerWeight) || centerWeight < 0.0)
    {
        throw std::invalid_argument(
            "Spin-It 2014 center weight must be finite and nonnegative"
        );
    }
}


SpinIt2014MomentProblem SpinIt2014MomentProblem::staticBalance(
    double centerWeight
)
{
    return SpinIt2014MomentProblem(
        spin::SpinMode::StaticBalance,
        centerWeight
    );
}


SpinIt2014MomentProblem SpinIt2014MomentProblem::yoYo()
{
    return SpinIt2014MomentProblem(spin::SpinMode::YoYo, 0.0);
}


SpinIt2014MomentProblem SpinIt2014MomentProblem::spinningTop(
    double centerWeight
)
{
    return SpinIt2014MomentProblem(
        spin::SpinMode::SpinningTop,
        centerWeight
    );
}


spin::SpinMode SpinIt2014MomentProblem::mode() const noexcept
{
    return mode_;
}


spin::SpinMomentEvaluation SpinIt2014MomentProblem::evaluate(
    const Vector10& rawMoments
) const
{
    if (mode_ == spin::SpinMode::YoYo)
    {
        return spin::SpinMomentProblem::yoYo().evaluate(rawMoments);
    }
    if (mode_ == spin::SpinMode::StaticBalance)
    {
        spin::SpinMomentEvaluation result =
            spin::SpinMomentProblem::staticBalance(
                spin::StaticObjective::LowerCenterOfMass
            ).evaluate(rawMoments);
        const double height = result.objective;
        const Eigen::Matrix<double, 10, 1> heightGradient =
            result.objectiveGradient;
        result.objective = centerWeight_ * height * height;
        result.objectiveGradient =
            2.0 * centerWeight_ * height * heightGradient;
        result.objectiveHessian = 2.0 * centerWeight_ * (
            heightGradient * heightGradient.transpose() +
            height * result.objectiveHessian
        );
        return result;
    }

    spin::SpinMomentEvaluation result =
        spin::SpinMomentProblem::spinningTop(0.0).evaluate(rawMoments);
    // Paper Eq. (3): (ell*M)^2 = s_z^2 in the contact-point frame.
    result.objective += centerWeight_ * rawMoments(3) * rawMoments(3);
    result.objectiveGradient(3) +=
        2.0 * centerWeight_ * rawMoments(3);
    result.objectiveHessian(3, 3) += 2.0 * centerWeight_;
    return result;
}


} // namespace spin2014
