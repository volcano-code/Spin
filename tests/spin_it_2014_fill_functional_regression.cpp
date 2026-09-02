#include "spin2014/AdaptiveFillFunctional.h"
#include "spin2014/AdaptiveVoxelDomain.h"
#include "spin2014/SpinIt2014MomentProblem.h"

#include "voxel/VoxelDomain.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>


namespace
{


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=spin-it-2014-fill-functional "
              << message << '\n';
    return EXIT_FAILURE;
}


std::shared_ptr<const voxel::VoxelDomain> fixture()
{
    voxel::GridSpec grid;
    grid.dimensions = Eigen::Vector3i(4, 4, 4);
    grid.minimum = Eigen::Vector3d(-1.0, -1.5, 0.0);
    grid.maximum = Eigen::Vector3d(1.0, 1.5, 3.0);
    std::vector<unsigned char> material(64U, 1U);
    std::vector<unsigned char> carveable(64U, 1U);
    // Keep the bottom layer as the fixed fabrication shell.
    for (int y = 0; y < 4; ++y)
    {
        for (int x = 0; x < 4; ++x)
        {
            carveable[static_cast<std::size_t>(x + 4 * y)] = 0U;
        }
    }
    return voxel::VoxelDomain::create(grid, material, {}, carveable);
}


} // namespace


int main()
{
    spin2014::AdaptiveVoxelDomain domain(fixture(), 1, 2);
    const spin2014::SpinIt2014MomentProblem problem =
        spin2014::SpinIt2014MomentProblem::spinningTop(0.37);
    const spin2014::AdaptiveFillFunctional functional(domain, problem, 0.19);

    Eigen::VectorXd beta(domain.leaves().size());
    for (Eigen::Index index = 0; index < beta.size(); ++index)
    {
        beta(index) = 0.13 + 0.71 *
            static_cast<double>(index + 1) /
            static_cast<double>(beta.size() + 1);
    }
    const spin2014::AdaptiveFillEvaluation value = functional.evaluate(beta);
    if (
        value.gradient.size() != beta.size() ||
        value.constraintJacobian.cols() != beta.size() ||
        value.constraints.size() != 4 ||
        !value.finalMoments.isFinite()
    )
    {
        return fail("evaluation dimensions or moment state are invalid");
    }

    const double step = 1e-6;
    double maximumGradientError = 0.0;
    double maximumJacobianError = 0.0;
    for (Eigen::Index column = 0; column < beta.size(); ++column)
    {
        Eigen::VectorXd plus = beta;
        Eigen::VectorXd minus = beta;
        plus(column) += step;
        minus(column) -= step;
        const auto plusValue = functional.evaluate(plus);
        const auto minusValue = functional.evaluate(minus);
        const double objectiveDerivative =
            (plusValue.objective - minusValue.objective) / (2.0 * step);
        maximumGradientError = std::max(
            maximumGradientError,
            std::abs(objectiveDerivative - value.gradient(column))
        );
        const Eigen::VectorXd constraintDerivative =
            (plusValue.constraints - minusValue.constraints) / (2.0 * step);
        maximumJacobianError = std::max(
            maximumJacobianError,
            (constraintDerivative - value.constraintJacobian.col(column))
                .lpNorm<Eigen::Infinity>()
        );
    }
    if (maximumGradientError > 1e-7 || maximumJacobianError > 1e-8)
    {
        return fail("analytic fill derivatives disagree with finite differences");
    }

    const auto filled = functional.evaluate(Eigen::VectorXd::Zero(beta.size()));
    const auto hollowed = functional.evaluate(Eigen::VectorXd::Ones(beta.size()));
    if (
        !filled.finalMoments.values().isApprox(
            domain.solidMoments().values(), 1e-12
        ) ||
        !hollowed.finalMoments.values().isApprox(
            (domain.solidMoments() - domain.variableMoments()).values(), 1e-12
        )
    )
    {
        return fail("beta=0/1 does not follow the paper's void convention");
    }

    const Eigen::VectorXd constant = Eigen::VectorXd::Constant(beta.size(), 0.4);
    if (std::abs(functional.laplacianEnergy(constant)) > 1e-13)
    {
        return fail("constant fill field has nonzero Laplacian energy");
    }

    std::cout << "[PASS] scenario=spin-it-2014-fill-functional"
              << " variables=" << beta.size()
              << " gradient_error=" << maximumGradientError
              << " jacobian_error=" << maximumJacobianError << '\n';
    return EXIT_SUCCESS;
}
