#include "spin2014/AdaptiveVoxelSpinSolver.h"

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
    std::cerr << "[FAIL] scenario=spin-it-2014-end-to-end "
              << message << '\n';
    return EXIT_FAILURE;
}


std::shared_ptr<const voxel::VoxelDomain> fixture()
{
    voxel::GridSpec grid;
    grid.dimensions = Eigen::Vector3i(8, 8, 8);
    grid.minimum = Eigen::Vector3d(-1.0, -1.0, 0.0);
    grid.maximum = Eigen::Vector3d(1.0, 1.0, 2.0);
    std::vector<unsigned char> material(512U, 1U);
    std::vector<unsigned char> carveable(512U, 0U);
    for (int z = 1; z < 7; ++z)
    {
        for (int y = 1; y < 7; ++y)
        {
            for (int x = 1; x < 7; ++x)
            {
                carveable[static_cast<std::size_t>(x + 8 * (y + 8 * z))] = 1U;
            }
        }
    }
    // An asymmetric fixed-shell defect forces a fractional coarse solution;
    // finest-level binary cells can compensate it exactly.
    const std::size_t defectA = static_cast<std::size_t>(1 + 8 * (3 + 8 * 4));
    const std::size_t defectB = static_cast<std::size_t>(1 + 8 * (4 + 8 * 4));
    material[defectA] = 0U;
    material[defectB] = 0U;
    carveable[defectA] = 0U;
    carveable[defectB] = 0U;
    return voxel::VoxelDomain::create(grid, material, {}, carveable);
}


} // namespace


int main()
{
    spin2014::AdaptiveVoxelSpinRequest request;
    request.domain = fixture();
    request.mode = spin::SpinMode::StaticBalance;
    request.initialLevel = 1;
    request.maximumLevel = 3;
    request.initialFill = 0.5;
    request.fractionalEpsilon = 0.05;
    request.laplacianWeight = 1e-3;
    request.maximumAdaptivePasses = 4;
    request.slqp.maximumIterations = 180;
    request.slqp.equalityTolerance = 1e-9;
    request.slqp.projectedGradientTolerance = 1e-7;

    const spin2014::AdaptiveVoxelSpinResult result =
        spin2014::solveAdaptiveVoxelSpin(request);
    const spin2014::AdaptiveVoxelSpinResult repeated =
        spin2014::solveAdaptiveVoxelSpin(request);
    if (
        result.status != spin2014::AdaptiveSpinStatus::Converged ||
        result.binaryConstraintInfinityNorm > 1e-8 ||
        result.materialBoundary.outerFaces.rows() == 0 ||
        result.materialBoundary.innerFaces.rows() == 0 ||
        result.adaptiveTrace.empty() ||
        result.inputHash.size() != 64U ||
        result.outputHash.size() != 64U ||
        result.outputHash != repeated.outputHash ||
        result.binaryVoidMask != repeated.binaryVoidMask
    )
    {
        std::cerr << "status=" << spin2014::adaptiveSpinStatusToken(result.status)
                  << " diagnostic=" << result.diagnostic
                  << " residual=" << result.binaryConstraintInfinityNorm
                  << " outer_faces=" << result.materialBoundary.outerFaces.rows()
                  << " inner_faces=" << result.materialBoundary.innerFaces.rows()
                  << " passes=" << result.adaptiveTrace.size() << '\n';
        return fail("adaptive hollowing did not produce a sealed deterministic result");
    }
    for (double value : result.binaryFill)
    {
        if (value != 0.0 && value != 1.0)
        {
            return fail("printable output contains a nonbinary fill value");
        }
    }

    std::cout << "[PASS] scenario=spin-it-2014-end-to-end"
              << " passes=" << result.adaptiveTrace.size()
              << " residual=" << result.binaryConstraintInfinityNorm
              << " hash=" << result.outputHash << '\n';
    return EXIT_SUCCESS;
}
