#include "spin/PlanarSpinMomentProblem.h"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace
{

bool close(double actual, double expected)
{
    return std::abs(actual - expected) <=
        1e-12 * std::max(1.0, std::abs(expected));
}

} // namespace

int main()
{
    spin::PlanarMassMoments6 moments;
    moments.mass = 15.0;
    moments.firstX = 15.0;
    moments.firstY = 22.5;
    moments.secondXX = 20.0;
    moments.secondYY = 45.0;
    moments.secondXY = 22.5;

    const auto evaluation = spin::PlanarSpinMomentProblem::fidgetSpinner().evaluate(moments);
    if (!close(evaluation.objective, -65.0))
    {
        std::cerr << "2D fidget objective must maximize the polar moment\n";
        return 1;
    }

    Eigen::Matrix<double, 6, 1> expectedGradient;
    expectedGradient << 0.0, 0.0, 0.0, -1.0, -1.0, 0.0;
    if (!evaluation.objectiveGradient.isApprox(expectedGradient, 0.0))
    {
        std::cerr << "2D fidget objective gradient uses the wrong moment order\n";
        return 1;
    }

    Eigen::Vector2d expectedConstraints(15.0, 22.5);
    if (!evaluation.constraints.isApprox(expectedConstraints, 0.0))
    {
        std::cerr << "2D fidget constraints must place the center of mass at the origin\n";
        return 1;
    }

    Eigen::Matrix<double, 2, 6> expectedJacobian =
        Eigen::Matrix<double, 2, 6>::Zero();
    expectedJacobian(0, 1) = 1.0;
    expectedJacobian(1, 2) = 1.0;
    if (!evaluation.constraintJacobian.isApprox(expectedJacobian, 0.0))
    {
        std::cerr << "2D fidget constraint Jacobian uses the wrong moment coordinates\n";
        return 1;
    }

    return 0;
}
