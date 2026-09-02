#include "voxel/MaterialOccupancy.h"

#include <array>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>


namespace voxel
{
namespace
{


struct FaceCollector
{
    std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>>
        faces;
};


bool inBounds(const Eigen::Vector3i& cell, const Eigen::Vector3i& dimensions)
{
    return
        (cell.array() >= 0).all() &&
        (cell.array() < dimensions.array()).all();
}


void appendFace(
    FaceCollector& collector,
    const VoxelDomain& domain,
    int x,
    int y,
    int z,
    int direction
)
{
    const auto corner = [&](int dx, int dy, int dz)
    {
        return static_cast<int>(domain.cornerIndex(x + dx, y + dy, z + dz));
    };

    std::array<int, 4> quad{};
    switch (direction)
    {
        case 0: quad = {corner(0,0,0), corner(0,0,1),
                        corner(0,1,0), corner(0,1,1)}; break; // -x
        case 1: quad = {corner(1,0,0), corner(1,1,0),
                        corner(1,0,1), corner(1,1,1)}; break; // +x
        case 2: quad = {corner(0,0,0), corner(1,0,0),
                        corner(0,0,1), corner(1,0,1)}; break; // -y
        case 3: quad = {corner(0,1,0), corner(0,1,1),
                        corner(1,1,0), corner(1,1,1)}; break; // +y
        case 4: quad = {corner(0,0,0), corner(0,1,0),
                        corner(1,0,0), corner(1,1,0)}; break; // -z
        case 5: quad = {corner(0,0,1), corner(1,0,1),
                        corner(0,1,1), corner(1,1,1)}; break; // +z
        default: throw std::logic_error("invalid voxel face direction");
    }
    collector.faces.emplace_back(quad[0], quad[1], quad[2]);
    collector.faces.emplace_back(quad[2], quad[1], quad[3]);
}


void compact(
    const VoxelDomain& domain,
    const FaceCollector& source,
    Eigen::MatrixXd& vertices,
    Eigen::MatrixXi& faces
)
{
    if (source.faces.empty())
    {
        vertices.resize(0, 3);
        faces.resize(0, 3);
        return;
    }
    std::vector<int> remap(domain.cornerCount(), -1);
    int vertexCount = 0;
    for (const auto& face : source.faces)
    {
        for (int column = 0; column < 3; ++column)
        {
            int& mapped = remap[static_cast<std::size_t>(face(column))];
            if (mapped < 0)
            {
                mapped = vertexCount++;
            }
        }
    }
    vertices.resize(vertexCount, 3);
    const Eigen::Vector3i cornerDimensions =
        domain.grid().dimensions + Eigen::Vector3i::Ones();
    for (std::size_t global = 0; global < remap.size(); ++global)
    {
        if (remap[global] < 0)
        {
            continue;
        }
        std::size_t index = global;
        const int x = static_cast<int>(
            index % static_cast<std::size_t>(cornerDimensions(0))
        );
        index /= static_cast<std::size_t>(cornerDimensions(0));
        const int y = static_cast<int>(
            index % static_cast<std::size_t>(cornerDimensions(1))
        );
        const int z = static_cast<int>(
            index / static_cast<std::size_t>(cornerDimensions(1))
        );
        vertices.row(remap[global]) =
            domain.cornerPosition(x, y, z).transpose();
    }
    faces.resize(static_cast<Eigen::Index>(source.faces.size()), 3);
    for (std::size_t row = 0; row < source.faces.size(); ++row)
    {
        for (int column = 0; column < 3; ++column)
        {
            faces(static_cast<Eigen::Index>(row), column) = remap[
                static_cast<std::size_t>(source.faces[row](column))
            ];
        }
    }
}


int countComponents(
    const VoxelDomain& domain,
    const std::vector<unsigned char>& membership
)
{
    static constexpr int offsets[6][3] = {
        {-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}
    };
    std::vector<unsigned char> visited(domain.cellCount(), 0U);
    int components = 0;
    for (std::size_t start = 0; start < domain.cellCount(); ++start)
    {
        if (membership[start] == 0U || visited[start] != 0U)
        {
            continue;
        }
        ++components;
        visited[start] = 1U;
        std::queue<std::size_t> pending;
        pending.push(start);
        while (!pending.empty())
        {
            const Eigen::Vector3i cell =
                domain.cellCoordinates(pending.front());
            pending.pop();
            for (const auto& offset : offsets)
            {
                const Eigen::Vector3i neighbor(
                    cell(0) + offset[0],
                    cell(1) + offset[1],
                    cell(2) + offset[2]
                );
                if (!inBounds(neighbor, domain.grid().dimensions))
                {
                    continue;
                }
                const std::size_t index = domain.cellIndex(
                    neighbor(0), neighbor(1), neighbor(2)
                );
                if (membership[index] != 0U && visited[index] == 0U)
                {
                    visited[index] = 1U;
                    pending.push(index);
                }
            }
        }
    }
    return components;
}


} // namespace


MaterialBoundaryMesh extractMaterialBoundaries(
    const std::shared_ptr<const VoxelDomain>& domain
)
{
    if (!domain)
    {
        throw std::invalid_argument("material extraction requires a domain");
    }
    const Eigen::Vector3i dimensions = domain->grid().dimensions;
    static constexpr int offsets[6][3] = {
        {-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}
    };

    // Flood the complement from the lattice boundary. Any unvisited empty
    // cell is an enclosed cavity and therefore receives an inner boundary.
    std::vector<unsigned char> exterior(domain->cellCount(), 0U);
    std::queue<std::size_t> pending;
    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                if (
                    x != 0 && x != dimensions(0) - 1 &&
                    y != 0 && y != dimensions(1) - 1 &&
                    z != 0 && z != dimensions(2) - 1
                )
                {
                    continue;
                }
                const std::size_t index = domain->cellIndex(x, y, z);
                if (!domain->material(index) && exterior[index] == 0U)
                {
                    exterior[index] = 1U;
                    pending.push(index);
                }
            }
        }
    }
    while (!pending.empty())
    {
        const Eigen::Vector3i cell =
            domain->cellCoordinates(pending.front());
        pending.pop();
        for (const auto& offset : offsets)
        {
            const Eigen::Vector3i neighbor(
                cell(0) + offset[0],
                cell(1) + offset[1],
                cell(2) + offset[2]
            );
            if (!inBounds(neighbor, dimensions))
            {
                continue;
            }
            const std::size_t index = domain->cellIndex(
                neighbor(0), neighbor(1), neighbor(2)
            );
            if (!domain->material(index) && exterior[index] == 0U)
            {
                exterior[index] = 1U;
                pending.push(index);
            }
        }
    }

    FaceCollector outer;
    FaceCollector inner;
    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                if (!domain->material(domain->cellIndex(x, y, z)))
                {
                    continue;
                }
                for (int direction = 0; direction < 6; ++direction)
                {
                    const Eigen::Vector3i neighbor(
                        x + offsets[direction][0],
                        y + offsets[direction][1],
                        z + offsets[direction][2]
                    );
                    if (!inBounds(neighbor, dimensions))
                    {
                        appendFace(outer, *domain, x, y, z, direction);
                        continue;
                    }
                    const std::size_t neighborIndex = domain->cellIndex(
                        neighbor(0), neighbor(1), neighbor(2)
                    );
                    if (domain->material(neighborIndex))
                    {
                        continue;
                    }
                    appendFace(
                        exterior[neighborIndex] != 0U ? outer : inner,
                        *domain,
                        x,
                        y,
                        z,
                        direction
                    );
                }
            }
        }
    }

    MaterialBoundaryMesh result;
    compact(*domain, outer, result.outerVertices, result.outerFaces);
    compact(*domain, inner, result.innerVertices, result.innerFaces);
    result.materialComponentCount = countComponents(
        *domain,
        domain->materialMask()
    );
    std::vector<unsigned char> cavities(domain->cellCount(), 0U);
    for (std::size_t index = 0; index < cavities.size(); ++index)
    {
        cavities[index] =
            !domain->material(index) && exterior[index] == 0U ? 1U : 0U;
    }
    result.cavityComponentCount = countComponents(*domain, cavities);
    return result;
}


} // namespace voxel
