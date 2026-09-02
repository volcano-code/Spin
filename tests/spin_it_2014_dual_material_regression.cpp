#include "spin2014/AdaptiveVoxelSpinSolver.h"
#include "spin2014/DualMaterialSpinSolver.h"

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
    std::cerr << "[FAIL] scenario=spin-it-2014-dual-material "
              << message << '\n';
    return EXIT_FAILURE;
}


std::shared_ptr<const voxel::VoxelDomain> fixture()
{
    voxel::GridSpec grid;
    grid.dimensions = Eigen::Vector3i(4, 4, 4);
    grid.minimum = Eigen::Vector3d(-1.0, -1.0, 0.0);
    grid.maximum = Eigen::Vector3d(1.0, 1.0, 2.0);
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


spin2014::AdaptiveVoxelSpinRequest request()
{
    spin2014::AdaptiveVoxelSpinRequest result;
    result.domain = fixture();
    result.mode = spin::SpinMode::StaticBalance;
    result.initialLevel = 1;
    result.maximumLevel = 2;
    result.maximumAdaptivePasses = 3;
    result.initialFill = 0.5;
    result.laplacianWeight = 1e-3;
    result.slqp.maximumIterations = 160;
    return result;
}


} // namespace


int main()
{
    const spin2014::AdaptiveVoxelSpinResult oneMaterial =
        spin2014::solveAdaptiveVoxelSpin(request());
    spin2014::DualMaterialSpinRequest equalDensity;
    equalDensity.hollowing = request();
    equalDensity.boundaryDensity = 1.0;
    equalDensity.interiorDensity = 1.0;
    const spin2014::DualMaterialSpinResult equal =
        spin2014::solveDualMaterialSpin(equalDensity);
    if (
        oneMaterial.status != spin2014::AdaptiveSpinStatus::Converged ||
        equal.hollowing.status != spin2014::AdaptiveSpinStatus::Converged ||
        equal.hollowing.outputHash != oneMaterial.outputHash ||
        equal.hollowing.binaryVoidMask != oneMaterial.binaryVoidMask ||
        !equal.hollowing.finalMoments.values().isApprox(
            oneMaterial.finalMoments.values(), 0.0
        )
    )
    {
        return fail("equal-density limit does not reduce exactly to hollowing");
    }

    spin2014::DualMaterialSpinRequest denseInterior = equalDensity;
    denseInterior.boundaryDensity = 1.17;
    denseInterior.interiorDensity = 8.1;
    const spin2014::DualMaterialSpinResult dense =
        spin2014::solveDualMaterialSpin(denseInterior);
    if (
        dense.hollowing.status != spin2014::AdaptiveSpinStatus::Converged ||
        !dense.hollowing.finalMoments.values().isApprox(
            (dense.boundaryContribution + dense.interiorContribution).values(),
            1e-11
        ) ||
        dense.materialLabels.size() != denseInterior.hollowing.domain->cellCount() ||
        dense.outputHash.size() != 64U
    )
    {
        return fail("density-weighted material contributions are inconsistent");
    }
    bool sawBoundary = false;
    bool sawInterior = false;
    for (unsigned char label : dense.materialLabels)
    {
        sawBoundary = sawBoundary || label == 1U;
        sawInterior = sawInterior || label == 2U;
        if (label > 2U)
        {
            return fail("dual-material label is outside {0,1,2}");
        }
    }
    if (!sawBoundary || !sawInterior)
    {
        return fail("dual-material output omitted a retained material class");
    }

    std::cout << "[PASS] scenario=spin-it-2014-dual-material"
              << " mass=" << dense.hollowing.finalMoments.mass()
              << " hash=" << dense.outputHash << '\n';
    return EXIT_SUCCESS;
}
