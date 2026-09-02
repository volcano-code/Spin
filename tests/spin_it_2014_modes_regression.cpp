#include "spin2014/AdaptiveVoxelSpinSolver.h"

#include "voxel/VoxelDomain.h"

#include <Eigen/Core>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>


namespace
{


std::shared_ptr<const voxel::VoxelDomain> fixture(bool centeredZ)
{
    voxel::GridSpec grid;
    grid.dimensions = Eigen::Vector3i(4, 4, 4);
    grid.minimum = centeredZ
        ? Eigen::Vector3d(-1.0, -1.0, -1.0)
        : Eigen::Vector3d(-1.0, -1.0, 0.0);
    grid.maximum = centeredZ
        ? Eigen::Vector3d(1.0, 1.0, 1.0)
        : Eigen::Vector3d(1.0, 1.0, 2.0);
    std::vector<unsigned char> material(64U, 1U);
    std::vector<unsigned char> carveable(64U, 0U);
    for (int z = 1; z < 3; ++z)
    {
        for (int y = 1; y < 3; ++y)
        {
            for (int x = 1; x < 3; ++x)
            {
                carveable[static_cast<std::size_t>(x + 4 * (y + 4 * z))] = 1U;
            }
        }
    }
    return voxel::VoxelDomain::create(grid, material, {}, carveable);
}


} // namespace


int main()
{
    const spin::SpinMode modes[] = {
        spin::SpinMode::StaticBalance,
        spin::SpinMode::SpinningTop,
        spin::SpinMode::YoYo
    };
    std::set<std::string> outputHashes;
    for (const spin::SpinMode mode : modes)
    {
        spin2014::AdaptiveVoxelSpinRequest request;
        request.domain = fixture(mode == spin::SpinMode::YoYo);
        request.mode = mode;
        request.initialLevel = 1;
        request.maximumLevel = 2;
        request.maximumAdaptivePasses = 3;
        request.slqp.maximumIterations = 120;
        const auto first = spin2014::solveAdaptiveVoxelSpin(request);
        const auto repeated = spin2014::solveAdaptiveVoxelSpin(request);
        if (
            first.status != spin2014::AdaptiveSpinStatus::Converged ||
            first.binaryConstraintInfinityNorm > 1e-10 ||
            first.outputHash != repeated.outputHash ||
            first.outputHash.size() != 64U
        )
        {
            std::cerr << "[FAIL] scenario=spin-it-2014-modes mode="
                      << static_cast<int>(mode)
                      << " status=" << spin2014::adaptiveSpinStatusToken(first.status)
                      << " residual=" << first.binaryConstraintInfinityNorm
                      << '\n';
            return EXIT_FAILURE;
        }
        outputHashes.insert(first.outputHash);
    }
    if (outputHashes.size() != 3U)
    {
        std::cerr << "[FAIL] scenario=spin-it-2014-modes hashes collided\n";
        return EXIT_FAILURE;
    }
    std::cout << "[PASS] scenario=spin-it-2014-modes static/top/yo-yo\n";
    return EXIT_SUCCESS;
}
