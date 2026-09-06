#pragma once

#include "spin/SpinMomentProblem.h"

#include <Eigen/Core>


namespace spin2014
{


enum class Spin2014Objective
{
    YoYo,
    TopInertiaOnly,
    TopLoweringOnly,
    TopMassAndLowering,
    StaticLowering
};


const char* spin2014ObjectiveToken(Spin2014Objective objective) noexcept;
spin::SpinMode spin2014ObjectiveMode(Spin2014Objective objective) noexcept;
Spin2014Objective defaultSpin2014Objective(spin::SpinMode mode);


/**
 * Spin-It 2014 moment objective with the transverse principal angle eliminated.
 *
 * The invariant transverse term is analytically equivalent to the sum of the
 * two squared principal-moment ratios.  This class does not reproduce the
 * paper's non-unique Givens-angle optimization trajectory.
 */
class SpinIt2014MomentProblem
{
public:
    using Vector10 = Eigen::Matrix<double, 10, 1>;

    static SpinIt2014MomentProblem fromObjective(
        Spin2014Objective objective,
        double inertiaWeight = 1.0,
        double centerWeight = 1.0
    );
    static SpinIt2014MomentProblem staticBalance(double centerWeight = 1.0);
    static SpinIt2014MomentProblem yoYo();
    static SpinIt2014MomentProblem spinningTop(double centerWeight);

    Spin2014Objective objective() const noexcept;
    spin::SpinMode mode() const noexcept;
    double inertiaWeight() const noexcept;
    double centerWeight() const noexcept;
    spin::SpinMomentEvaluation evaluate(const Vector10& rawMoments) const;

private:
    SpinIt2014MomentProblem(
        Spin2014Objective objective,
        double inertiaWeight,
        double centerWeight
    );

    Spin2014Objective objective_ = Spin2014Objective::StaticLowering;
    double inertiaWeight_ = 1.0;
    double centerWeight_ = 1.0;
};


} // namespace spin2014
