#include "spin/PlanarSpinMomentProblem.h"

#include <stdexcept>

namespace spin
{


PlanarSpinMomentProblem PlanarSpinMomentProblem::fidgetSpinner()
{
    return PlanarSpinMomentProblem{};
}


PlanarSpinMomentEvaluation PlanarSpinMomentProblem::evaluate(
    const PlanarMassMoments6& moments) const
{
    const Eigen::Matrix<double, 6, 1> raw = moments.vector();
    if (!raw.allFinite() || !(moments.mass > 0.0))
        throw std::invalid_argument("planar fidget-spinner moments must be finite with positive mass");

    PlanarSpinMomentEvaluation result;
    result.objective = -(moments.secondXX + moments.secondYY);
    result.objectiveGradient << 0.0, 0.0, 0.0, -1.0, -1.0, 0.0;
    result.constraints << moments.firstX, moments.firstY;
    result.constraintJacobian(0, 1) = 1.0;
    result.constraintJacobian(1, 2) = 1.0;
    return result;
}


} // namespace spin
