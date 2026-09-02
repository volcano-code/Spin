#include "api/BalanceModel.h"
#include "manufacturing/ManufacturingGeometry.h"
#include "voxel/MaterialOccupancy.h"
#include "voxel/VoxelDomain.h"
#include "voxel/VoxelGrid.h"
#include "voxel/VoxelTopology.h"

#include <Eigen/Core>

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>


namespace
{


std::shared_ptr<const voxel::VoxelDomain> cubeWithOptionalCavity(bool cavity)
{
    voxel::GridSpec grid;
    grid.dimensions = Eigen::Vector3i(5, 5, 5);
    grid.minimum = Eigen::Vector3d(-1.0, 2.0, 4.0);
    grid.maximum = Eigen::Vector3d(4.0, 7.0, 9.0);
    std::vector<unsigned char> material(125U, 1U);
    if (cavity)
    {
        const std::size_t center = static_cast<std::size_t>(
            2 + 5 * (2 + 5 * 2)
        );
        material[center] = 0U;
    }
    return voxel::VoxelDomain::create(grid, material);
}


int domainScenario()
{
    const auto domain = cubeWithOptionalCavity(false);
    if (
        domain->cellCount() != 125U ||
        domain->cornerCount() != 216U ||
        domain->cellIndex(1, 2, 3) != 86U ||
        domain->cellCoordinates(86U) != Eigen::Vector3i(1, 2, 3) ||
        domain->artifactId().size() != 64U ||
        domain->resourceEstimate().materialCellCount != 125U ||
        domain->resourceEstimate().fillStateBytes != 16U
    )
    {
        std::cerr << "immutable voxel-domain contract mismatch\n";
        return 1;
    }
    voxel::VoxelFillState first(domain);
    voxel::VoxelFillState second = first;
    first.setFilled(domain->cellIndex(2, 2, 2), false);
    if (
        first.filledCount() != 124U ||
        second.filledCount() != 125U ||
        first.contentHash() == second.contentHash()
    )
    {
        std::cerr << "fill snapshots are not independent compact states\n";
        return 1;
    }
    first.refillMaterial();
    if (first.filledCount() != 125U ||
        first.contentHash() != second.contentHash())
    {
        std::cerr << "material refill did not restore the baseline\n";
        return 1;
    }
    voxel::GridSpec shiftedGrid = domain->grid();
    shiftedGrid.minimum.array() += 1.0;
    shiftedGrid.maximum.array() += 1.0;
    const auto shiftedDomain = voxel::VoxelDomain::create(
        shiftedGrid, domain->materialMask()
    );
    const voxel::VoxelFillState shifted(shiftedDomain);
    if (shifted.contentHash() == second.contentHash())
    {
        std::cerr << "fill hash did not bind the immutable domain\n";
        return 1;
    }
    return 0;
}


int boundaryScenario()
{
    const auto domain = cubeWithOptionalCavity(true);
    const voxel::MaterialBoundaryMesh boundary =
        voxel::extractMaterialBoundaries(domain);
    const manufacturing::BoundaryAudit outer =
        manufacturing::auditTriangleBoundary(
            boundary.outerVertices,
            boundary.outerFaces,
            manufacturing::BoundaryOrientation::OutwardPositive
        );
    const manufacturing::BoundaryAudit inner =
        manufacturing::auditTriangleBoundary(
            boundary.innerVertices,
            boundary.innerFaces,
            manufacturing::BoundaryOrientation::InwardNegative
        );
    if (
        boundary.materialComponentCount != 1 ||
        boundary.cavityComponentCount != 1 ||
        !outer.strictBoundaryValid() ||
        !inner.strictBoundaryValid()
    )
    {
        std::cerr << "material boundary extraction failed strict audit\n";
        return 1;
    }

    VoxelGrid grid = VoxelGrid::fromMaterialOccupancy(domain, 0.0, 0);
    const VoxelTopologyReport before = analyzeVoxelTopology(grid);
    if (
        before.materialComponentCount != 1 ||
        before.carvedVoidComponentCount != 1 ||
        before.materialVoxelCount != 124 ||
        before.carvedVoidVoxelCount != 1 ||
        grid(2, 2, 2).isFilled ||
        !grid(2, 2, 2).isInner
    )
    {
        std::cerr << "material occupancy facade lost its initial cavity\n";
        return 1;
    }

    // A normal carving edit is reversible, while the source cavity remains.
    grid(1, 1, 1).isFilled = false;
    const int changed = grid.refillAll();
    if (changed != 1 || grid(2, 2, 2).isFilled ||
        !grid(1, 1, 1).isFilled)
    {
        std::cerr << "baseline refill incorrectly filled a source cavity\n";
        return 1;
    }
    return 0;
}


int invalidScenario()
{
    voxel::GridSpec grid;
    grid.dimensions = Eigen::Vector3i(2, 2, 2);
    std::vector<unsigned char> material(8U, 0U);
    try
    {
        static_cast<void>(voxel::VoxelDomain::create(grid, material));
    }
    catch (const std::invalid_argument&)
    {
        material[0] = 1U;
        std::vector<unsigned char> carveable(8U, 0U);
        carveable[1] = 1U;
        try
        {
            static_cast<void>(voxel::VoxelDomain::create(
                grid, material, material, carveable
            ));
        }
        catch (const std::invalid_argument&)
        {
            return 0;
        }
    }
    std::cerr << "invalid material domains were not rejected\n";
    return 1;
}


int modelScenario()
{
    api::MaterialOccupancyData occupancy;
    occupancy.dimensions = Eigen::Vector3i(5, 5, 5);
    occupancy.minimum = Eigen::Vector3d::Zero();
    occupancy.maximum = Eigen::Vector3d::Ones();
    occupancy.materialMaskXFast.assign(125U, 1U);
    occupancy.materialMaskXFast[static_cast<std::size_t>(
        2 + 5 * (2 + 5 * 2)
    )] = 0U;
    occupancy.protectedBoundaryLayers = 0;
    api::BalanceModelConfig config;
    config.voxelResolution = 5;
    std::unique_ptr<api::BalanceModel> model =
        api::BalanceModel::fromMaterialOccupancy(occupancy, config);
    model->setGravity(Eigen::Vector3d(0.0, 0.0, -1.0));
    model->setSupportFromMesh(0.0, 0.0);
    api::BalanceRunConfig run;
    run.optimizer.enableDeformation = false;
    run.optimizer.maximumCarvingBlocks = 1;
    const api::BalanceResult result = model->optimize(run);
    if (
        !result.succeeded() ||
        !result.manufacturingReady() ||
        result.innerFaces().rows() != 12
    )
    {
        std::cerr << "material occupancy did not produce a sealed cavity result\n";
        for (const std::string& reason : result.manufacturingFailures())
        {
            std::cerr << "  " << reason << '\n';
        }
        return 1;
    }
    return 0;
}


} // namespace


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        return 2;
    }
    const std::string scenario(argv[1]);
    if (scenario == "domain") return domainScenario();
    if (scenario == "boundary") return boundaryScenario();
    if (scenario == "invalid") return invalidScenario();
    if (scenario == "model") return modelScenario();
    return 2;
}
