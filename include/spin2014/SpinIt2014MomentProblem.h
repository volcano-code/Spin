#pragma once

#include "spin/SpinMomentProblem.h"

#include <Eigen/Core>


namespace spin2014
{


/** Exact Section 3--4 objective specialization used by the 2014 method. */
class SpinIt2014MomentProblem
{
public:
    using Vector10 = Eigen::Matrix<double, 10, 1>;

    static SpinIt2014MomentProblem staticBalance(double centerWeight = 1.0);
    static SpinIt2014MomentProblem yoYo();
    static SpinIt2014MomentProblem spinningTop(double centerWeight);

    spin::SpinMode mode() const noexcept;
    spin::SpinMomentEvaluation evaluate(const Vector10& rawMoments) const;

private:
    SpinIt2014MomentProblem(spin::SpinMode mode, double centerWeight);

    spin::SpinMode mode_ = spin::SpinMode::StaticBalance;
    double centerWeight_ = 1.0;
};


} // namespace spin2014
