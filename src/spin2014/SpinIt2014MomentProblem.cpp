#include "spin2014/SpinIt2014MomentProblem.h"

#include <cmath>
#include <stdexcept>


namespace spin2014
{


const char* spin2014ObjectiveToken(Spin2014Objective objective) noexcept
{
    switch (objective)
    {
    case Spin2014Objective::YoYo: return "YOYO";
    case Spin2014Objective::TopInertiaOnly: return "TOP_INERTIA_ONLY";
    case Spin2014Objective::TopLoweringOnly: return "TOP_LOWERING_ONLY";
    case Spin2014Objective::TopMassAndLowering:
        return "TOP_MASS_AND_LOWERING";
    case Spin2014Objective::StaticLowering: return "STATIC_LOWERING";
    }
    return "UNKNOWN";
}


spin::SpinMode spin2014ObjectiveMode(Spin2014Objective objective) noexcept
{
    switch (objective)
    {
    case Spin2014Objective::YoYo:
        return spin::SpinMode::YoYo;
    case Spin2014Objective::TopInertiaOnly:
    case Spin2014Objective::TopLoweringOnly:
    case Spin2014Objective::TopMassAndLowering:
        return spin::SpinMode::SpinningTop;
    case Spin2014Objective::StaticLowering:
        return spin::SpinMode::StaticBalance;
    }
    return spin::SpinMode::StaticBalance;
}


Spin2014Objective defaultSpin2014Objective(spin::SpinMode mode)
{
    switch (mode)
    {
    case spin::SpinMode::StaticBalance:
        return Spin2014Objective::StaticLowering;
    case spin::SpinMode::YoYo:
        return Spin2014Objective::YoYo;
    case spin::SpinMode::SpinningTop:
        return Spin2014Objective::TopMassAndLowering;
    case spin::SpinMode::FidgetSpinner:
    case spin::SpinMode::BuoyantStability:
        break;
    }
    throw std::invalid_argument("mode has no Spin-It 2014 objective");
}


SpinIt2014MomentProblem::SpinIt2014MomentProblem(
    Spin2014Objective objective,
    double inertiaWeight,
    double centerWeight
)
    : objective_(objective),
      inertiaWeight_(inertiaWeight),
      centerWeight_(centerWeight)
{
    if (
        !std::isfinite(inertiaWeight) || inertiaWeight < 0.0 ||
        !std::isfinite(centerWeight) || centerWeight < 0.0
    )
    {
        throw std::invalid_argument(
            "Spin-It 2014 weights must be finite and nonnegative"
        );
    }
}


SpinIt2014MomentProblem SpinIt2014MomentProblem::fromObjective(
    Spin2014Objective objective,
    double inertiaWeight,
    double centerWeight
)
{
    return SpinIt2014MomentProblem(objective, inertiaWeight, centerWeight);
}


SpinIt2014MomentProblem SpinIt2014MomentProblem::staticBalance(
    double centerWeight
)
{
    return fromObjective(
        Spin2014Objective::StaticLowering, 0.0, centerWeight
    );
}


SpinIt2014MomentProblem SpinIt2014MomentProblem::yoYo()
{
    return fromObjective(Spin2014Objective::YoYo, 1.0, 0.0);
}


SpinIt2014MomentProblem SpinIt2014MomentProblem::spinningTop(
    double centerWeight
)
{
    return fromObjective(
        Spin2014Objective::TopMassAndLowering, 1.0, centerWeight
    );
}


Spin2014Objective SpinIt2014MomentProblem::objective() const noexcept
{
    return objective_;
}


spin::SpinMode SpinIt2014MomentProblem::mode() const noexcept
{
    return spin2014ObjectiveMode(objective_);
}


double SpinIt2014MomentProblem::inertiaWeight() const noexcept
{
    return inertiaWeight_;
}


double SpinIt2014MomentProblem::centerWeight() const noexcept
{
    return centerWeight_;
}


spin::SpinMomentEvaluation SpinIt2014MomentProblem::evaluate(
    const Vector10& rawMoments
) const
{
    if (objective_ == Spin2014Objective::StaticLowering)
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
        objective_ == Spin2014Objective::YoYo
        ? spin::SpinMomentProblem::yoYo().evaluate(rawMoments)
        : spin::SpinMomentProblem::spinningTop(0.0).evaluate(rawMoments);
    result.objective *= inertiaWeight_;
    result.objectiveGradient *= inertiaWeight_;
    result.objectiveHessian *= inertiaWeight_;

    if (objective_ == Spin2014Objective::TopLoweringOnly)
    {
        const spin::SpinMomentEvaluation height =
            spin::SpinMomentProblem::staticBalance(
                spin::StaticObjective::LowerCenterOfMass
            ).evaluate(rawMoments);
        result.objective += centerWeight_ * height.objective * height.objective;
        result.objectiveGradient += 2.0 * centerWeight_ * height.objective *
            height.objectiveGradient;
        result.objectiveHessian += 2.0 * centerWeight_ * (
            height.objectiveGradient * height.objectiveGradient.transpose() +
            height.objective * height.objectiveHessian
        );
    }
    else if (objective_ == Spin2014Objective::TopMassAndLowering)
    {
        // Paper Eq. (3): (ell*M)^2 = s_z^2 in the contact-point frame.
        result.objective += centerWeight_ * rawMoments(3) * rawMoments(3);
        result.objectiveGradient(3) +=
            2.0 * centerWeight_ * rawMoments(3);
        result.objectiveHessian(3, 3) += 2.0 * centerWeight_;
    }
    return result;
}


} // namespace spin2014
