#include "spin2014/AdaptiveVoxelSpinSolver.h"

#include "voxel/VoxelDomain.h"

#include <Eigen/Core>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>


namespace
{


std::shared_ptr<const voxel::VoxelDomain> fixture(int side, bool centeredZ)
{
    voxel::GridSpec grid;
    grid.dimensions = Eigen::Vector3i::Constant(side);
    grid.minimum = centeredZ
        ? Eigen::Vector3d(-1.0, -1.0, -1.0)
        : Eigen::Vector3d(-1.0, -1.0, 0.0);
    grid.maximum = centeredZ
        ? Eigen::Vector3d(1.0, 1.0, 1.0)
        : Eigen::Vector3d(1.0, 1.0, 2.0);
    const std::size_t count = static_cast<std::size_t>(side * side * side);
    std::vector<unsigned char> material(count, 1U);
    std::vector<unsigned char> carveable(count, 0U);
    for (int z = 1; z < side - 1; ++z)
    {
        for (int y = 1; y < side - 1; ++y)
        {
            for (int x = 1; x < side - 1; ++x)
            {
                carveable[static_cast<std::size_t>(
                    x + side * (y + side * z)
                )] = 1U;
            }
        }
    }
    return voxel::VoxelDomain::create(grid, material, {}, carveable);
}


std::shared_ptr<const voxel::VoxelDomain> asymmetricFixture(int side)
{
    voxel::GridSpec grid;
    grid.dimensions = Eigen::Vector3i::Constant(side);
    grid.minimum = Eigen::Vector3d(-1.0, -0.8, 0.0);
    grid.maximum = Eigen::Vector3d(1.0, 0.8, 1.2);
    const std::size_t count = static_cast<std::size_t>(side * side * side);
    std::vector<unsigned char> material(count, 1U);
    for (int z = side - 2; z < side; ++z)
    {
        for (int y = 0; y < side; ++y)
        {
            for (int x = 0; x < 1; ++x)
            {
                material[static_cast<std::size_t>(
                    x + side * (y + side * z)
                )] = 0U;
            }
        }
    }
    std::vector<unsigned char> carveable(count, 0U);
    for (int z = 1; z < side - 1; ++z)
    {
        for (int y = 1; y < side - 1; ++y)
        {
            for (int x = 1; x < side - 1; ++x)
            {
                const int index = x + side * (y + side * z);
                const int offsets[] = {-1, 1, -side, side, -side * side, side * side};
                bool protectedCell = material[static_cast<std::size_t>(index)] != 0U;
                for (const int offset : offsets)
                {
                    protectedCell = protectedCell &&
                        material[static_cast<std::size_t>(index + offset)] != 0U;
                }
                carveable[static_cast<std::size_t>(index)] =
                    protectedCell ? 1U : 0U;
            }
        }
    }
    return voxel::VoxelDomain::create(grid, material, {}, carveable);
}


int levelForSide(int side)
{
    int level = 0;
    while ((1 << level) < side) ++level;
    return level;
}


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=spin-it-2014-resolution-ladder "
              << message << '\n';
    return EXIT_FAILURE;
}


} // namespace


int main()
{
    for (const spin::SpinMode mode : {
             spin::SpinMode::SpinningTop,
             spin::SpinMode::YoYo
         })
    {
        for (const int side : {4, 8, 16})
        {
            spin2014::AdaptiveVoxelSpinRequest request;
            request.domain = fixture(side, mode == spin::SpinMode::YoYo);
            request.mode = mode;
            request.initialLevel = std::min(2, levelForSide(side));
            request.maximumLevel = levelForSide(side);
            request.maximumAdaptivePasses =
                request.maximumLevel - request.initialLevel + 2;
            request.slqp.maximumIterations = 200;
            request.slqp.equalityTolerance = 1e-9;
            request.slqp.projectedGradientTolerance = 1e-7;
            const spin2014::AdaptiveVoxelSpinResult result =
                spin2014::solveAdaptiveVoxelSpin(request);
            if (!result.adaptiveTrace.empty())
            {
                const auto& observed = result.adaptiveTrace.back();
                std::cout << "mode=" << spin::spinModeToken(mode)
                          << " side=" << side
                          << " status=" << spin2014::adaptiveSpinStatusToken(result.status)
                          << " eq=" << observed.finalEqualityInfinityNorm
                          << " kkt=" << observed.finalProjectedGradientInfinityNorm
                          << " line_evals=" << observed.totalLineSearchEvaluations
                          << " active_transitions=" << observed.activeSetTransitions
                          << " iterations=" << observed.slqpTrace.size() << '\n';
            }
            if (
                result.status != spin2014::AdaptiveSpinStatus::Converged ||
                result.binaryConstraintInfinityNorm > 1e-8 ||
                result.adaptiveTrace.empty() ||
                result.binaryProjectionMaxDelta < 0.0 ||
                result.binaryProjectionChangedLeaves < 0
            )
            {
                if (!result.adaptiveTrace.empty())
                {
                    const auto& pass = result.adaptiveTrace.back();
                    std::cerr << " final_eq=" << pass.finalEqualityInfinityNorm
                              << " final_kkt="
                              << pass.finalProjectedGradientInfinityNorm
                              << " line_evals=" << pass.totalLineSearchEvaluations
                              << " active_transitions=" << pass.activeSetTransitions
                              << " trace_size=" << pass.slqpTrace.size();
                    if (!pass.slqpTrace.empty())
                    {
                        const auto& last = pass.slqpTrace.back();
                        std::cerr << " last_step=" << last.acceptedStep
                                  << " last_active=" << last.activeVariables
                                  << " last_merit=" << last.meritValue;
                        for (std::size_t index = 0;
                             index < pass.slqpTrace.size();
                             index += 25)
                        {
                            const auto& sample = pass.slqpTrace[index];
                            std::cerr << " [" << index
                                      << ":kkt=" << sample.kktInfinityNorm
                                      << ",active=" << sample.activeVariables
                                      << ",step=" << sample.acceptedStep << ']';
                        }
                    }
                    std::cerr << '\n';
                }
                return fail(
                    "mode=" + std::string(spin::spinModeToken(mode)) +
                    " side=" + std::to_string(side) +
                    " status=" + spin2014::adaptiveSpinStatusToken(result.status) +
                    " diagnostic=" + result.diagnostic
                );
            }
            const spin2014::AdaptiveSpinPass& pass = result.adaptiveTrace.back();
            if (
                pass.slqpStatus != spin2014::SlqpStatus::Converged ||
                !std::isfinite(pass.finalEqualityInfinityNorm) ||
                !std::isfinite(pass.finalProjectedGradientInfinityNorm) ||
                pass.totalLineSearchEvaluations < 0
            )
            {
                return fail("resolution trace omitted SLQP convergence evidence");
            }
        }
    }
    spin2014::AdaptiveVoxelSpinRequest refinedOutput;
    refinedOutput.domain = fixture(4, false);
    refinedOutput.mode = spin::SpinMode::SpinningTop;
    refinedOutput.initialLevel = 2;
    refinedOutput.maximumLevel = 5;
    refinedOutput.maximumAdaptivePasses = 5;
    refinedOutput.binaryConstraintTolerance = 5e-3;
    const auto refinedResult = spin2014::solveAdaptiveVoxelSpin(refinedOutput);
    if (
        refinedResult.status != spin2014::AdaptiveSpinStatus::Converged ||
        refinedResult.outputGridDimensions != Eigen::Vector3i::Constant(32)
    )
    {
        return fail("coarse source was not refined onto the requested octree grid");
    }
    spin2014::AdaptiveVoxelSpinRequest asymmetric;
    asymmetric.domain = asymmetricFixture(16);
    asymmetric.mode = spin::SpinMode::SpinningTop;
    asymmetric.initialLevel = 2;
    asymmetric.maximumLevel = 4;
    asymmetric.maximumAdaptivePasses = 6;
    // A 16^3 reconstructed proxy has a finite voxel quantization floor.  This
    // tolerance is intentionally local to the proxy regression; the paper
    // reproduction gate still requires the published level-9 experiment.
    asymmetric.binaryConstraintTolerance = 5e-3;
    asymmetric.slqp.maximumIterations = 240;
    const auto asymmetricResult = spin2014::solveAdaptiveVoxelSpin(asymmetric);
    if (asymmetricResult.status != spin2014::AdaptiveSpinStatus::Converged)
    {
        if (!asymmetricResult.adaptiveTrace.empty())
        {
            std::cerr << "asymmetric passes=" << asymmetricResult.adaptiveTrace.size() << '\n';
            for (const auto& prior : asymmetricResult.adaptiveTrace)
            {
                std::cerr << " pass=" << prior.pass
                          << " status=" << spin2014::slqpStatusToken(prior.slqpStatus)
                          << " split=" << prior.splitLeafCount
                          << " merge=" << prior.mergedParentCount
                          << " fractional=" << prior.fractionalLeafCount
                          << " eq=" << prior.finalEqualityInfinityNorm << '\n';
            }
            const auto& pass = asymmetricResult.adaptiveTrace.back();
            std::cerr << "asymmetric trace eq=" << pass.finalEqualityInfinityNorm
                      << " kkt=" << pass.finalProjectedGradientInfinityNorm
                      << " iterations=" << pass.slqpTrace.size() << '\n';
            for (const auto& entry : pass.slqpTrace)
            {
                std::cerr << "  i=" << entry.iteration
                          << " eq=" << entry.equalityInfinityNorm
                          << " kkt=" << entry.kktInfinityNorm
                          << " merit=" << entry.meritValue
                          << " step=" << entry.acceptedStep
                          << " fallback=" << entry.usedFirstOrderFallback
                          << " active=" << entry.activeVariables << '\n';
            }
        }
        return fail(
            "asymmetric top did not converge: " +
            std::string(spin2014::adaptiveSpinStatusToken(asymmetricResult.status)) +
            " binary=" + std::to_string(asymmetricResult.binaryConstraintInfinityNorm) +
            " diagnostic=" + asymmetricResult.diagnostic
        );
    }
    spin2014::AdaptiveVoxelSpinRequest bounded;
    bounded.domain = fixture(4, false);
    bounded.mode = spin::SpinMode::SpinningTop;
    bounded.initialLevel = 9;
    bounded.maximumLevel = 9;
    bounded.maximumLeafCount = 1000;
    const auto limited = spin2014::solveAdaptiveVoxelSpin(bounded);
    if (
        limited.status != spin2014::AdaptiveSpinStatus::ResourceLimitReached ||
        limited.diagnostic.find("leaf") == std::string::npos
    )
    {
        return fail("level-9 resource bound was not reported explicitly");
    }
    std::cout << "[PASS] scenario=spin-it-2014-resolution-ladder 4/8/16 top/yo-yo\n";
    return EXIT_SUCCESS;
}
