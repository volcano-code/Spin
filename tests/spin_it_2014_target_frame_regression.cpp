#include "spin2014/AdaptiveFillFunctional.h"
#include "spin2014/AdaptiveVoxelSpinSolver.h"

#include "voxel/VoxelDomain.h"

#include <Eigen/Geometry>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace
{

int fail(const char* message)
{
    std::cerr << "[FAIL] scenario=spin-it-2014-target-frame " << message << '\n';
    return EXIT_FAILURE;
}

std::shared_ptr<const voxel::VoxelDomain> fixture()
{
    voxel::GridSpec grid;
    grid.dimensions = Eigen::Vector3i(4, 4, 4);
    grid.minimum = Eigen::Vector3d(-1.0, -1.5, -0.5);
    grid.maximum = Eigen::Vector3d(1.0, 1.5, 2.5);
    std::vector<unsigned char> material(64U, 1U);
    std::vector<unsigned char> carveable(64U, 1U);
    for (int z = 0; z < 4; ++z)
    {
        for (int y = 0; y < 4; ++y)
        {
            carveable[static_cast<std::size_t>(4 * (y + 4 * z))] = 0U;
        }
    }
    return voxel::VoxelDomain::create(grid, material, {}, carveable);
}

} // namespace

int main()
{
    try
    {
        spin2014::AdaptiveVoxelDomain domain(fixture(), 1, 2);
        const auto problem = spin2014::SpinIt2014MomentProblem::yoYo();
        spin::SpinTargetFrame frame;
        frame.origin = Eigen::Vector3d(0.17, -0.23, 0.31);
        frame.axes = Eigen::AngleAxisd(
            0.43, Eigen::Vector3d(1.0, 2.0, -1.0).normalized()
        ).toRotationMatrix();
        const spin2014::AdaptiveFillFunctional functional(
            domain, problem, 0.07, 1.0, 1.0, frame
        );
        Eigen::VectorXd beta = Eigen::VectorXd::Constant(
            domain.leaves().size(), 0.35
        );
        const auto value = functional.evaluate(beta);
        const auto expectedTarget = value.finalMoments.transformed(
            frame.axes.transpose(), -frame.axes.transpose() * frame.origin
        );
        const auto expected = problem.evaluate(expectedTarget.values());
        if (!value.constraints.isApprox(expected.constraints, 1e-12) ||
            std::abs(value.paperObjective - expected.objective) > 1e-12)
        {
            return fail("target-frame objective or constraints ignore transformed mass moments");
        }

        const double h = 1e-6;
        for (Eigen::Index i = 0; i < beta.size(); ++i)
        {
            Eigen::VectorXd plus = beta, minus = beta;
            plus(i) += h;
            minus(i) -= h;
            const auto a = functional.evaluate(plus);
            const auto b = functional.evaluate(minus);
            if (std::abs((a.objective - b.objective) / (2.0 * h) -
                         value.gradient(i)) > 1e-7 ||
                !((a.constraints - b.constraints) / (2.0 * h)).isApprox(
                    value.constraintJacobian.col(i), 1e-8
                ))
            {
                return fail("target-frame fill derivative disagrees with finite differences");
            }
        }

        spin2014::AdaptiveVoxelSpinRequest request;
        request.domain = fixture();
        request.initialLevel = 0;
        request.maximumLevel = 2;
        request.maximumAdaptivePasses = 1;
        const std::string identityHash =
            spin2014::computeAdaptiveVoxelSpinRequestHash(request);
        request.targetFrame = frame;
        const std::string rotatedHash =
            spin2014::computeAdaptiveVoxelSpinRequestHash(request);
        if (identityHash == rotatedHash)
        {
            return fail("target frame is absent from the adaptive request identity");
        }
        request.targetFrame.axes(0, 0) = 2.0;
        bool rejected = false;
        try
        {
            spin2014::validateAdaptiveVoxelSpinRequest(request);
        }
        catch (const std::invalid_argument&)
        {
            rejected = true;
        }
        if (!rejected)
        {
            return fail("non-orthonormal target frame was accepted");
        }
        std::cout << "[PASS] target-frame raw-moment transform, derivatives, and request identity\n";
        return EXIT_SUCCESS;
    }
    catch (const std::exception& error)
    {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
