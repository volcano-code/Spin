#include "voxel/VoxelTopology.h"

#include "voxel/VoxelGrid.h"

#include <cstddef>
#include <limits>
#include <queue>
#include <stdexcept>
#include <vector>


namespace
{

enum class Domain
{
    Material,
    CarvedVoid
};


std::size_t linearIndex(
    int x,
    int y,
    int z,
    const Eigen::Vector3i& dimensions
)
{
    return static_cast<std::size_t>(
        (x * dimensions(1) + y) * dimensions(2) + z
    );
}


bool belongsToDomain(
    const Voxel& voxel,
    Domain domain
)
{
    if (domain == Domain::Material)
    {
        return voxel.distanceToMesh < 0.0 && voxel.isFilled;
    }

    return voxel.isInner && !voxel.isFilled;
}


int countComponents(
    const VoxelGrid& voxelGrid,
    Domain domain,
    int& voxelCount
)
{
    const Eigen::Vector3i& dimensions =
        voxelGrid.dimensions();

    std::size_t totalVoxelCount = 1;

    for (int axis = 0; axis < 3; ++axis)
    {
        if (dimensions(axis) <= 0)
        {
            throw std::invalid_argument(
                "topology requires positive grid dimensions"
            );
        }

        const std::size_t dimension =
            static_cast<std::size_t>(dimensions(axis));

        if (
            totalVoxelCount >
                std::numeric_limits<std::size_t>::max() / dimension
        )
        {
            throw std::overflow_error(
                "voxel topology grid size overflowed"
            );
        }

        totalVoxelCount *= dimension;
    }

    std::vector<unsigned char> visited(
        totalVoxelCount,
        0U
    );

    const int offsets[6][3] =
    {
        {-1, 0, 0}, {1, 0, 0},
        {0, -1, 0}, {0, 1, 0},
        {0, 0, -1}, {0, 0, 1}
    };

    int components = 0;
    voxelCount = 0;

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                if (belongsToDomain(voxelGrid(x, y, z), domain))
                {
                    ++voxelCount;
                }

                const std::size_t start =
                    linearIndex(x, y, z, dimensions);

                if (
                    visited[start] != 0U ||
                    !belongsToDomain(voxelGrid(x, y, z), domain)
                )
                {
                    continue;
                }

                ++components;
                visited[start] = 1U;

                std::queue<Eigen::Vector3i> pending;
                pending.push(Eigen::Vector3i(x, y, z));

                while (!pending.empty())
                {
                    const Eigen::Vector3i current =
                        pending.front();
                    pending.pop();

                    for (const auto& offset : offsets)
                    {
                        const Eigen::Vector3i neighbor(
                            current(0) + offset[0],
                            current(1) + offset[1],
                            current(2) + offset[2]
                        );

                        if (
                            (neighbor.array() < 0).any() ||
                            (neighbor.array() >= dimensions.array()).any()
                        )
                        {
                            continue;
                        }

                        const std::size_t index = linearIndex(
                            neighbor(0),
                            neighbor(1),
                            neighbor(2),
                            dimensions
                        );

                        if (
                            visited[index] == 0U &&
                            belongsToDomain(
                                voxelGrid(
                                    neighbor(0),
                                    neighbor(1),
                                    neighbor(2)
                                ),
                                domain
                            )
                        )
                        {
                            visited[index] = 1U;
                            pending.push(neighbor);
                        }
                    }
                }
            }
        }
    }

    return components;
}

} // namespace


VoxelTopologyReport::VoxelTopologyReport()
    :
    materialComponentCount(0),
    carvedVoidComponentCount(0),
    materialVoxelCount(0),
    carvedVoidVoxelCount(0)
{
}


VoxelTopologyReport analyzeVoxelTopology(
    const VoxelGrid& voxelGrid
)
{
    VoxelTopologyReport result;
    result.materialComponentCount = countComponents(
        voxelGrid,
        Domain::Material,
        result.materialVoxelCount
    );
    result.carvedVoidComponentCount = countComponents(
        voxelGrid,
        Domain::CarvedVoid,
        result.carvedVoidVoxelCount
    );
    return result;
}
