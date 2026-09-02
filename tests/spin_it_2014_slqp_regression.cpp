#include "spin2014/AdaptiveFillFunctional.h"
#include "spin2014/AdaptiveVoxelDomain.h"
#include "spin2014/SlqpLbfgsSolver.h"
#include "spin2014/SpinIt2014MomentProblem.h"

#include "voxel/VoxelDomain.h"

#include <Eigen/Core>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>


namespace
{


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=spin-it-2014-slqp " << message << '\n';
    return EXIT_FAILURE;
}


std::shared_ptr<const voxel::VoxelDomain> fixture()
{
    voxel::GridSpec grid;
    grid.dimensions = Eigen::Vector3i(4, 4, 4);
    grid.minimum = Eigen::Vector3d(-1.0, -1.0, 0.0);
    grid.maximum = Eigen::Vector3d(1.0, 1.0, 2.0);
    std::vector<unsigned char> material(64U, 1U);
    std::vector<unsigned char> carveable(64U, 1U);
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
    spin2014::AdaptiveVoxelDomain domain(fixture(), 2, 2);
    const spin2014::AdaptiveFillFunctional functional(
        domain,
        spin2014::SpinIt2014MomentProblem::staticBalance(1.0),
        1e-3
    );
    Eigen::VectorXd initial(domain.leaves().size());
    for (Eigen::Index index = 0; index < initial.size(); ++index)
    {
        initial(index) = 0.1 + 0.7 *
            static_cast<double>((index * 17 + 5) % initial.size()) /
            static_cast<double>(initial.size());
    }
    const auto initialValue = functional.evaluate(initial);

    spin2014::SlqpLbfgsSettings settings;
    settings.maximumIterations = 160;
    settings.equalityTolerance = 1e-9;
    settings.projectedGradientTolerance = 1e-7;
    const spin2014::SlqpLbfgsResult result =
        spin2014::solveSlqpLbfgs(functional, initial, settings);
    const spin2014::SlqpLbfgsResult repeated =
        spin2014::solveSlqpLbfgs(functional, initial, settings);

    if (
        result.status != spin2014::SlqpStatus::Converged ||
        result.evaluation.constraints.lpNorm<Eigen::Infinity>() > 1e-8 ||
        result.evaluation.objective > initialValue.objective + 1e-10 ||
        (result.beta.array() < 0.0).any() ||
        (result.beta.array() > 1.0).any() ||
        result.trace.empty() ||
        result.traceHash.size() != 64U ||
        result.traceHash != repeated.traceHash ||
        !result.beta.isApprox(repeated.beta, 0.0)
    )
    {
        std::cerr << "status=" << spin2014::slqpStatusToken(result.status)
                  << " diagnostic=" << result.diagnostic
                  << " iterations=" << result.trace.size()
                  << " residual="
                  << result.evaluation.constraints.lpNorm<Eigen::Infinity>()
                  << " initial_objective=" << initialValue.objective
                  << " final_objective=" << result.evaluation.objective
                  << " repeat_hash=" << (result.traceHash == repeated.traceHash)
                  << '\n';
        return fail("active-set SLQP/L-BFGS did not converge deterministically");
    }
    bool usedLbfgs = false;
    for (const spin2014::SlqpIteration& iteration : result.trace)
    {
        usedLbfgs = usedLbfgs || iteration.lbfgsHistorySize > 0;
    }
    if (!usedLbfgs)
    {
        return fail("solver never exercised its L-BFGS approximation");
    }

    std::cout << "[PASS] scenario=spin-it-2014-slqp"
              << " iterations=" << result.trace.size()
              << " residual="
              << result.evaluation.constraints.lpNorm<Eigen::Infinity>()
              << " hash=" << result.traceHash << '\n';
    return EXIT_SUCCESS;
}
