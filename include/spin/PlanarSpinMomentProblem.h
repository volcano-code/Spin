#pragma once

#include "spin/PlanarMassMoments.h"

#include <Eigen/Core>

namespace spin
{


struct PlanarSpinMomentEvaluation
{
    double objective = 0.0;
    Eigen::Matrix<double, 6, 1> objectiveGradient =
        Eigen::Matrix<double, 6, 1>::Zero();
    Eigen::Vector2d constraints = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 2, 6> constraintJacobian =
        Eigen::Matrix<double, 2, 6>::Zero();
};


/** Paper's two-dimensional fidget-spinner mass-moment objective and constraints. */
class PlanarSpinMomentProblem
{
public:
    static PlanarSpinMomentProblem fidgetSpinner();

    /** Input order is `(mass,Lx,Ly,Qxx,Qyy,Qxy)`, already density weighted. */
    PlanarSpinMomentEvaluation evaluate(const PlanarMassMoments6& moments) const;
};


} // namespace spin
