#include "api/SpinIt2014.h"

#include "voxel/VoxelDomain.h"

#include <Eigen/Core>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>


int main()
{
    voxel::GridSpec grid;
    grid.dimensions = Eigen::Vector3i(4, 4, 4);
    grid.minimum = Eigen::Vector3d(-1.0, -1.0, 0.0);
    grid.maximum = Eigen::Vector3d(1.0, 1.0, 2.0);
    std::vector<unsigned char> material(64U, 1U);
    std::vector<unsigned char> carveable(64U, 0U);
    for (int z = 1; z < 3; ++z)
        for (int y = 1; y < 3; ++y)
            for (int x = 1; x < 3; ++x)
                carveable[static_cast<std::size_t>(x + 4 * (y + 4 * z))] = 1U;

    spin2014::AdaptiveVoxelSpinRequest request;
    request.domain = voxel::VoxelDomain::create(grid, material, {}, carveable);
    request.mode = spin::SpinMode::StaticBalance;
    request.initialLevel = 1;
    request.maximumLevel = 2;
    request.maximumAdaptivePasses = 3;
    request.slqp.maximumIterations = 160;
    const api::SpinIt2014Artifact artifact = api::solveSpinIt2014(
        request,
        api::SpinIt2014Evidence::PaperFormulaReproduction
    );
    if (
        !artifact.succeeded() || !artifact.geometrySealed() ||
        artifact.outerFaces().rows() == 0 ||
        artifact.designFaces().rows() == 0 ||
        artifact.materialFaces().rows() == 0 ||
        artifact.cavityFaces().rows() == 0 ||
        artifact.diagnosticsJson().find("spin-it-2014-artifact-v1") ==
            std::string::npos ||
        artifact.diagnosticsJson().find("\"geometry_sealed\": true") ==
            std::string::npos
    )
    {
        std::cerr << "[FAIL] scenario=spin-it-2014-api status="
                  << artifact.statusToken() << " sealed="
                  << artifact.geometrySealed() << '\n';
        return EXIT_FAILURE;
    }

    // A valid optimization request may still terminate before producing
    // physical mass properties.  The API must return an auditable failed
    // artifact instead of throwing while serializing diagnostics.
    voxel::GridSpec hardGrid;
    hardGrid.dimensions = Eigen::Vector3i(8, 8, 8);
    hardGrid.minimum = Eigen::Vector3d(-1.0, -1.0, 0.0);
    hardGrid.maximum = Eigen::Vector3d(1.0, 1.0, 2.0);
    std::vector<unsigned char> hardMaterial(512U, 1U);
    std::vector<unsigned char> hardCarveable(512U, 0U);
    for (int z = 1; z < 7; ++z)
        for (int y = 1; y < 7; ++y)
            for (int x = 1; x < 7; ++x)
                hardCarveable[static_cast<std::size_t>(x + 8 * (y + 8 * z))] = 1U;
    const std::size_t defect = static_cast<std::size_t>(1 + 8 * (3 + 8 * 4));
    hardMaterial[defect] = 0U;
    hardCarveable[defect] = 0U;
    spin2014::AdaptiveVoxelSpinRequest unfinished;
    unfinished.domain = voxel::VoxelDomain::create(
        hardGrid, hardMaterial, {}, hardCarveable
    );
    unfinished.mode = spin::SpinMode::StaticBalance;
    unfinished.initialLevel = 1;
    unfinished.maximumLevel = 3;
    unfinished.maximumAdaptivePasses = 1;
    unfinished.slqp.maximumIterations = 1;
    const api::SpinIt2014Artifact failed = api::solveSpinIt2014(unfinished);
    if (
        failed.succeeded() || failed.geometrySealed() ||
        failed.diagnosticsJson().find("\"mass_properties_available\": false") ==
            std::string::npos ||
        failed.diagnosticsJson().find("\"center_of_mass_model\": null") ==
            std::string::npos
    )
    {
        std::cerr << "[FAIL] scenario=spin-it-2014-api failed artifact contract\n";
        return EXIT_FAILURE;
    }
    std::cout << "[PASS] scenario=spin-it-2014-api status="
              << artifact.statusToken() << '\n';
    return EXIT_SUCCESS;
}
