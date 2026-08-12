#include "physics/MassProperties.h"
#include "voxel/VoxelGrid.h"
#include "voxel/VoxelTopology.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>


namespace
{

void makeBox(
    Eigen::MatrixXd& vertices,
    Eigen::MatrixXi& faces
)
{
    vertices.resize(8, 3);
    vertices <<
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        1.0, 1.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        1.0, 0.0, 1.0,
        1.0, 1.0, 1.0,
        0.0, 1.0, 1.0;

    faces.resize(12, 3);
    faces <<
        0, 2, 1, 0, 3, 2,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5;
}


int fail(
    const std::string& scenario,
    const std::string& message
)
{
    std::cerr
        << "[FAIL] scenario=" << scenario
        << ' ' << message << '\n';
    return EXIT_FAILURE;
}


bool validateWallGrid(
    const VoxelGrid& grid,
    double minimumWallThickness,
    int& carveableCount,
    int& protectedInsideCount
)
{
    carveableCount = 0;
    protectedInsideCount = 0;
    const Eigen::Vector3i dimensions = grid.dimensions();
    const double halfDiagonal = 0.5 * grid.voxelSize().norm();

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                const Voxel& voxel = grid(x, y, z);
                const double exactCubeSurfaceDistance = std::min(
                    std::min(
                        std::min(voxel.center(0), 1.0 - voxel.center(0)),
                        std::min(voxel.center(1), 1.0 - voxel.center(1))
                    ),
                    std::min(voxel.center(2), 1.0 - voxel.center(2))
                );
                const double exactBoxClearance =
                    exactCubeSurfaceDistance - halfDiagonal;

                if (!voxel.isInner)
                {
                    if (
                        voxel.distanceToMesh < 0.0 &&
                        exactBoxClearance + 1e-12 < minimumWallThickness
                    )
                    {
                        ++protectedInsideCount;
                    }
                    continue;
                }

                ++carveableCount;

                if (exactBoxClearance + 1e-12 < minimumWallThickness)
                {
                    return false;
                }
            }
        }
    }

    return true;
}


int runWall()
{
    const std::string scenario = "wall";
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    makeBox(vertices, faces);

    const double thickness = 0.12;
    VoxelGrid coarse(vertices, faces, thickness, 24);
    VoxelGrid fine(vertices, faces, thickness, 48);

    int coarseCount = 0;
    int fineCount = 0;
    int coarseProtectedCount = 0;
    int fineProtectedCount = 0;

    if (
        !validateWallGrid(
            coarse,
            thickness,
            coarseCount,
            coarseProtectedCount
        ) ||
        !validateWallGrid(
            fine,
            thickness,
            fineCount,
            fineProtectedCount
        ) ||
        coarseCount <= 0 ||
        fineCount <= coarseCount ||
        coarseProtectedCount <= 0 ||
        fineProtectedCount <= 0 ||
        coarse.longestAxisResolution() != 24 ||
        fine.longestAxisResolution() != 48 ||
        coarse.minimumWallThickness() != thickness ||
        fine.minimumWallThickness() != thickness
    )
    {
        return fail(
            scenario,
            "model-space t_min contract failed across resolutions"
        );
    }

    bool rejectedThickness = false;
    bool rejectedResolution = false;
    const double retainedThickness = coarse.minimumWallThickness();
    const int retainedResolution = coarse.longestAxisResolution();

    try
    {
        VoxelGrid invalid(vertices, faces, -1.0, 24);
        static_cast<void>(invalid);
    }
    catch (const std::invalid_argument&)
    {
        rejectedThickness = true;
    }

    try
    {
        coarse.build(vertices, faces, 0.25, 2);
    }
    catch (const std::invalid_argument&)
    {
        rejectedResolution = true;
    }

    if (
        !rejectedThickness ||
        !rejectedResolution ||
        coarse.minimumWallThickness() != retainedThickness ||
        coarse.longestAxisResolution() != retainedResolution
    )
    {
        return fail(scenario, "invalid manufacturing input was accepted");
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " coarse=" << coarseCount
        << " fine=" << fineCount << '\n';
    return EXIT_SUCCESS;
}


int runTopology()
{
    const std::string scenario = "topology";
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
    makeBox(vertices, faces);
    VoxelGrid grid(vertices, faces, 0.0, 24);

    const VoxelTopologyReport initial =
        analyzeVoxelTopology(grid);

    if (
        initial.materialComponentCount != 1 ||
        initial.carvedVoidComponentCount != 0 ||
        initial.materialVoxelCount <= 0
    )
    {
        return fail(scenario, "initial component oracle failed");
    }

    bool carvedFirst = false;
    bool carvedSecond = false;
    const Eigen::Vector3i dimensions = grid.dimensions();

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                Voxel& voxel = grid(x, y, z);

                if (!voxel.isInner)
                {
                    continue;
                }

                if (!carvedFirst)
                {
                    voxel.isFilled = false;
                    carvedFirst = true;
                    continue;
                }

                if (!carvedSecond)
                {
                    voxel.isFilled = false;
                    carvedSecond = true;
                }
            }
        }
    }

    const VoxelTopologyReport carved =
        analyzeVoxelTopology(grid);

    if (
        !carvedFirst ||
        !carvedSecond ||
        carved.materialComponentCount != 1 ||
        carved.carvedVoidComponentCount < 1 ||
        carved.carvedVoidVoxelCount != 2
    )
    {
        return fail(scenario, "six-neighbor component report failed");
    }

    VoxelGrid disconnected(vertices, faces, 0.0, 24);
    bool separatedFirst = false;
    bool separatedSecond = false;
    Eigen::Vector3i firstKey = Eigen::Vector3i::Constant(-1);

    for (int x = 1; x < dimensions(0) - 1; ++x)
    {
        for (int y = 1; y < dimensions(1) - 1; ++y)
        {
            for (int z = 1; z < dimensions(2) - 1; ++z)
            {
                Voxel& voxel = disconnected(x, y, z);

                if (!voxel.isInner)
                {
                    continue;
                }

                if (!separatedFirst)
                {
                    voxel.isFilled = false;
                    firstKey = Eigen::Vector3i(x, y, z);
                    separatedFirst = true;
                    continue;
                }

                if (
                    !separatedSecond &&
                    (Eigen::Vector3i(x, y, z) - firstKey)
                        .cwiseAbs()
                        .maxCoeff() >= 2
                )
                {
                    voxel.isFilled = false;
                    separatedSecond = true;
                }
            }
        }
    }

    const VoxelTopologyReport separated =
        analyzeVoxelTopology(disconnected);

    if (
        !separatedFirst ||
        !separatedSecond ||
        separated.carvedVoidVoxelCount != 2 ||
        separated.carvedVoidComponentCount != 2
    )
    {
        return fail(scenario, "disconnected void components were not distinguished");
    }

    VoxelGrid separatedMaterial(vertices, faces, 0.0, 24);
    Eigen::Vector3i keptFirst = Eigen::Vector3i::Constant(-1);
    Eigen::Vector3i keptSecond = Eigen::Vector3i::Constant(-1);

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                Voxel& voxel = separatedMaterial(x, y, z);

                if (voxel.distanceToMesh >= 0.0)
                {
                    continue;
                }

                const Eigen::Vector3i key(x, y, z);

                if ((keptFirst.array() < 0).any())
                {
                    keptFirst = key;
                    continue;
                }

                if (
                    (keptSecond.array() < 0).any() &&
                    (key - keptFirst).cwiseAbs().maxCoeff() >= 2
                )
                {
                    keptSecond = key;
                    continue;
                }

                voxel.isFilled = false;
            }
        }
    }

    if ((keptFirst.array() < 0).any() || (keptSecond.array() < 0).any())
    {
        return fail(scenario, "could not construct separated material fixture");
    }

    const VoxelTopologyReport material =
        analyzeVoxelTopology(separatedMaterial);

    if (
        material.materialVoxelCount != 2 ||
        material.materialComponentCount != 2
    )
    {
        return fail(scenario, "disconnected material components were not distinguished");
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " material_components=" << carved.materialComponentCount
        << " void_components=" << carved.carvedVoidComponentCount
        << '\n';
    return EXIT_SUCCESS;
}

} // namespace


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        return fail("unknown", "usage: voxel-constraints-regression <scenario>");
    }

    const std::string scenario = argv[1];

    if (scenario == "wall")
    {
        return runWall();
    }
    if (scenario == "topology")
    {
        return runTopology();
    }

    return fail(scenario, "unknown scenario");
}
