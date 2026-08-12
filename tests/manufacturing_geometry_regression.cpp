#include "manufacturing/ManufacturingGeometry.h"

#include "deformation/GridDeformationDomain.h"
#include "extraction/VoidMeshExtractor.h"
#include "optimization/AlternatingBalanceOptimizer.h"
#include "voxel/VoxelGrid.h"

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <unistd.h>


namespace
{


using manufacturing::BoundaryAudit;
using manufacturing::BoundaryOrientation;
using manufacturing::ManufacturingGeometry;


struct TriangleMesh
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
};


int fail(const std::string& scenario, const std::string& reason)
{
    std::cerr
        << "[FAIL] scenario=" << scenario
        << " reason=" << reason << '\n';
    return EXIT_FAILURE;
}


int pass(const std::string& scenario)
{
    std::cout << "[PASS] scenario=" << scenario << '\n';
    return EXIT_SUCCESS;
}


TriangleMesh makeBox(
    const Eigen::Vector3d& minimum,
    const Eigen::Vector3d& maximum,
    bool outward = true
)
{
    TriangleMesh result;
    result.vertices.resize(8, 3);
    result.vertices <<
        minimum(0), minimum(1), minimum(2),
        maximum(0), minimum(1), minimum(2),
        maximum(0), maximum(1), minimum(2),
        minimum(0), maximum(1), minimum(2),
        minimum(0), minimum(1), maximum(2),
        maximum(0), minimum(1), maximum(2),
        maximum(0), maximum(1), maximum(2),
        minimum(0), maximum(1), maximum(2);

    result.faces.resize(12, 3);
    result.faces <<
        0, 2, 1, 0, 3, 2,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5;

    if (!outward)
    {
        result.faces.col(1).swap(result.faces.col(2));
    }
    return result;
}


TriangleMesh join(const TriangleMesh& first, const TriangleMesh& second)
{
    TriangleMesh result;
    result.vertices.resize(
        first.vertices.rows() + second.vertices.rows(),
        3
    );
    result.vertices.topRows(first.vertices.rows()) = first.vertices;
    result.vertices.bottomRows(second.vertices.rows()) = second.vertices;

    result.faces.resize(first.faces.rows() + second.faces.rows(), 3);
    result.faces.topRows(first.faces.rows()) = first.faces;
    result.faces.bottomRows(second.faces.rows()) =
        (second.faces.array() + first.vertices.rows()).matrix();
    return result;
}


BoundaryAudit auditOutward(const TriangleMesh& mesh)
{
    return manufacturing::auditTriangleBoundary(
        mesh.vertices,
        mesh.faces,
        BoundaryOrientation::OutwardPositive
    );
}


bool hasFailure(
    const manufacturing::ManufacturingReadinessReport& report,
    const std::string& expected
)
{
    return std::find(
        report.failures.begin(),
        report.failures.end(),
        expected
    ) != report.failures.end();
}


deformation::DeformedGeometry identityGeometry(
    const TriangleMesh& outer,
    const VoxelGrid& grid
)
{
    deformation::DeformedGeometry result;
    result.outerVertices = outer.vertices;
    result.outerFaces = outer.faces;
    result.gridCorners = grid.corners();

    const Eigen::Vector3i dimensions = grid.dimensions();
    result.voxelCenters.resize(dimensions.prod(), 3);
    Eigen::Index row = 0;
    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                result.voxelCenters.row(row++) =
                    grid(x, y, z).center.transpose();
            }
        }
    }

    VoidMeshExtractor::extract(
        grid,
        result.innerVertices,
        result.innerFaces
    );
    return result;
}


struct GridCell
{
    int x = -1;
    int y = -1;
    int z = -1;

    bool valid() const noexcept
    {
        return x >= 0 && y >= 0 && z >= 0;
    }
};


GridCell nearestEligibleCell(
    const VoxelGrid& grid,
    const Eigen::Vector3d& target = Eigen::Vector3d::Constant(0.5)
)
{
    GridCell result;
    double best = std::numeric_limits<double>::infinity();
    const Eigen::Vector3i dimensions = grid.dimensions();
    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                const Voxel& voxel = grid(x, y, z);
                if (!voxel.isInner)
                {
                    continue;
                }
                const double distance = (voxel.center - target).squaredNorm();
                if (distance < best)
                {
                    best = distance;
                    result = GridCell{x, y, z};
                }
            }
        }
    }
    return result;
}


int firstActiveCell(const VoxelGrid& grid)
{
    const deformation::GridDeformationDomain domain =
        deformation::GridDeformationDomain::fromVoxelGrid(grid);
    for (int cell = 0; cell < domain.cellCount(); ++cell)
    {
        if (domain.isActiveCell(cell))
        {
            return cell;
        }
    }
    return -1;
}


int runPrepareEmpty()
{
    const std::string scenario = "prepare-empty";
    const TriangleMesh outer = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );
    const VoxelGrid grid(outer.vertices, outer.faces, 0.1, 8);
    const ManufacturingGeometry prepared =
        manufacturing::prepareManufacturingGeometry(
            identityGeometry(outer, grid),
            grid
        );

    if (
        prepared.report.ready ||
        prepared.report.sourceBound ||
        !hasFailure(prepared.report, "r7_source_result_not_bound") ||
        !prepared.report.synchronizedInnerGeometry ||
        !prepared.report.deformedCellsValid ||
        !prepared.report.massPropertiesConsistent ||
        !prepared.report.singleMaterialComponent ||
        !prepared.report.innerBoundary.empty ||
        prepared.report.wallClearance.applicable ||
        !prepared.report.wallClearance.completed ||
        !prepared.report.containment.completed ||
        !prepared.report.containment.allCavitiesInsideOuter ||
        std::abs(prepared.report.signedMaterialVolume - 1.0) > 1.0e-10 ||
        (prepared.report.centerOfMass - Eigen::Vector3d::Constant(0.5)).norm() >
            1.0e-10
    )
    {
        return fail(scenario, "identity solid did not pass the final gate");
    }
    return pass(scenario);
}


int runPrepareCavity()
{
    const std::string scenario = "prepare-cavity";
    const TriangleMesh outer = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );
    VoxelGrid grid(outer.vertices, outer.faces, 0.1, 8);
    const GridCell carved = nearestEligibleCell(grid);
    if (!carved.valid())
    {
        return fail(scenario, "fixture has no eligible interior voxel");
    }
    grid(carved.x, carved.y, carved.z).isFilled = false;

    const ManufacturingGeometry prepared =
        manufacturing::prepareManufacturingGeometry(
            identityGeometry(outer, grid),
            grid
        );
    const double expectedVolume = 1.0 - grid.voxelSize().prod();
    if (
        prepared.report.ready ||
        prepared.report.sourceBound ||
        !hasFailure(prepared.report, "r7_source_result_not_bound") ||
        prepared.report.innerBoundary.empty ||
        prepared.report.innerBoundary.faceComponentCount != 1 ||
        !prepared.report.wallClearance.applicable ||
        !prepared.report.wallClearance.intersectionFree ||
        !prepared.report.wallClearance.meetsRequestedThickness ||
        !prepared.report.containment.allCavitiesInsideOuter ||
        prepared.report.voxelTopology.carvedVoidComponentCount != 1 ||
        prepared.report.voxelTopology.carvedVoidVoxelCount != 1 ||
        std::abs(prepared.report.signedMaterialVolume - expectedVolume) >
            1.0e-10
    )
    {
        return fail(scenario, "identity single cavity did not pass the final gate");
    }
    return pass(scenario);
}


int runStaleInner()
{
    const std::string scenario = "stale-inner";
    const TriangleMesh outer = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );
    VoxelGrid grid(outer.vertices, outer.faces, 0.0, 8);
    const GridCell carved = nearestEligibleCell(grid);
    if (!carved.valid())
    {
        return fail(scenario, "fixture has no eligible interior voxel");
    }
    grid(carved.x, carved.y, carved.z).isFilled = false;
    deformation::DeformedGeometry geometry = identityGeometry(outer, grid);
    geometry.innerFaces.resize(0, 3);

    const ManufacturingGeometry prepared =
        manufacturing::prepareManufacturingGeometry(geometry, grid);
    if (
        prepared.report.ready ||
        prepared.report.synchronizedInnerGeometry ||
        !hasFailure(prepared.report, "inner_geometry_not_synchronized")
    )
    {
        return fail(scenario, "stale cavity connectivity was not rejected");
    }
    return pass(scenario);
}


int runCellMutation(bool collapse)
{
    const std::string scenario = collapse ? "cell-collapse" : "cell-swap";
    const TriangleMesh outer = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );
    const VoxelGrid grid(outer.vertices, outer.faces, 0.0, 8);
    deformation::DeformedGeometry geometry = identityGeometry(outer, grid);
    const deformation::GridDeformationDomain domain =
        deformation::GridDeformationDomain::fromVoxelGrid(grid);
    const int cell = firstActiveCell(grid);
    if (cell < 0)
    {
        return fail(scenario, "fixture has no active deformation cell");
    }

    const std::array<int, 8> corners = domain.cellGlobalCorners(cell);
    if (collapse)
    {
        geometry.gridCorners.row(corners[6]) =
            geometry.gridCorners.row(corners[7]);
    }
    else
    {
        geometry.gridCorners.row(corners[6]).swap(
            geometry.gridCorners.row(corners[7])
        );
    }
    geometry.innerVertices = geometry.gridCorners;

    const ManufacturingGeometry prepared =
        manufacturing::prepareManufacturingGeometry(geometry, grid);
    if (
        prepared.report.ready ||
        prepared.report.deformedCellsValid ||
        prepared.report.deformedCells.invalidCellCount == 0 ||
        !hasFailure(prepared.report, "deformed_cells_invalid")
    )
    {
        return fail(scenario, "invalid local cell deformation escaped the gate");
    }
    return pass(scenario);
}


GridCell cellWithEligibleSixNeighbors(const VoxelGrid& grid)
{
    const Eigen::Vector3i dimensions = grid.dimensions();
    const int offsets[6][3] =
    {
        {-1, 0, 0}, {1, 0, 0},
        {0, -1, 0}, {0, 1, 0},
        {0, 0, -1}, {0, 0, 1}
    };
    for (int x = 1; x < dimensions(0) - 1; ++x)
    {
        for (int y = 1; y < dimensions(1) - 1; ++y)
        {
            for (int z = 1; z < dimensions(2) - 1; ++z)
            {
                if (!grid(x, y, z).isInner)
                {
                    continue;
                }
                bool eligible = true;
                for (const auto& offset : offsets)
                {
                    eligible = eligible && grid(
                        x + offset[0],
                        y + offset[1],
                        z + offset[2]
                    ).isInner;
                }
                if (eligible)
                {
                    return GridCell{x, y, z};
                }
            }
        }
    }
    return GridCell{};
}


int runMaterialComponents()
{
    const std::string scenario = "material-components";
    const TriangleMesh outer = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );
    VoxelGrid grid(outer.vertices, outer.faces, 0.0, 12);
    const GridCell island = cellWithEligibleSixNeighbors(grid);
    if (!island.valid())
    {
        return fail(scenario, "fixture lacks an eligible six-neighbor star");
    }

    const int offsets[6][3] =
    {
        {-1, 0, 0}, {1, 0, 0},
        {0, -1, 0}, {0, 1, 0},
        {0, 0, -1}, {0, 0, 1}
    };
    for (const auto& offset : offsets)
    {
        grid(
            island.x + offset[0],
            island.y + offset[1],
            island.z + offset[2]
        ).isFilled = false;
    }

    const ManufacturingGeometry prepared =
        manufacturing::prepareManufacturingGeometry(
            identityGeometry(outer, grid),
            grid
        );
    if (
        prepared.report.ready ||
        prepared.report.singleMaterialComponent ||
        prepared.report.voxelTopology.materialComponentCount != 2 ||
        !hasFailure(prepared.report, "material_not_single_component")
    )
    {
        return fail(scenario, "floating material island escaped connectivity gate");
    }
    return pass(scenario);
}


int runContainment()
{
    const std::string scenario = "containment";
    const TriangleMesh outer = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );
    VoxelGrid grid(outer.vertices, outer.faces, 0.0, 8);
    const GridCell carved = nearestEligibleCell(grid);
    if (!carved.valid())
    {
        return fail(scenario, "fixture has no eligible interior voxel");
    }
    grid(carved.x, carved.y, carved.z).isFilled = false;
    deformation::DeformedGeometry geometry = identityGeometry(outer, grid);
    geometry.gridCorners.rowwise() += Eigen::RowVector3d(2.0, 0.0, 0.0);
    geometry.innerVertices = geometry.gridCorners;

    const ManufacturingGeometry prepared =
        manufacturing::prepareManufacturingGeometry(geometry, grid);
    if (
        prepared.report.ready ||
        !prepared.report.containment.completed ||
        prepared.report.containment.allCavitiesInsideOuter ||
        prepared.report.containment.failedComponent < 0 ||
        !hasFailure(
            prepared.report,
            "cavity_not_contained_by_outer_boundary"
        )
    )
    {
        return fail(scenario, "external cavity escaped containment gate");
    }
    return pass(scenario);
}


ManufacturingGeometry prepareScaledCavity(
    const TriangleMesh& outer,
    const VoxelGrid& grid,
    const GridCell& carved,
    double targetClearance
)
{
    deformation::DeformedGeometry geometry = identityGeometry(outer, grid);
    const Eigen::Vector3d center = grid(
        carved.x,
        carved.y,
        carved.z
    ).center;
    const double nearestOuterPlane = std::min(
        {
            center(0), 1.0 - center(0),
            center(1), 1.0 - center(1),
            center(2), 1.0 - center(2)
        }
    );
    const double restHalfExtent = 0.5 * grid.voxelSize().minCoeff();
    const double factor =
        (nearestOuterPlane - targetClearance) / restHalfExtent;

    for (Eigen::Index row = 0; row < geometry.gridCorners.rows(); ++row)
    {
        geometry.gridCorners.row(row) =
            (center + factor *
                (geometry.gridCorners.row(row).transpose() - center)
            ).transpose();
    }
    geometry.innerVertices = geometry.gridCorners;
    return manufacturing::prepareManufacturingGeometry(geometry, grid);
}


int runWallThreshold()
{
    const std::string scenario = "wall-threshold";
    const TriangleMesh outer = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );
    constexpr double requested = 0.2;
    VoxelGrid grid(outer.vertices, outer.faces, requested, 8);
    const GridCell carved = nearestEligibleCell(grid);
    if (!carved.valid())
    {
        return fail(scenario, "fixture has no requested-wall eligible voxel");
    }
    grid(carved.x, carved.y, carved.z).isFilled = false;

    const ManufacturingGeometry above = prepareScaledCavity(
        outer, grid, carved, requested + 0.01
    );
    const ManufacturingGeometry boundary = prepareScaledCavity(
        outer, grid, carved, requested
    );
    const ManufacturingGeometry below = prepareScaledCavity(
        outer, grid, carved, requested - 0.01
    );
    if (
        above.report.ready ||
        above.report.sourceBound ||
        !above.report.wallClearance.meetsRequestedThickness ||
        boundary.report.ready ||
        boundary.report.sourceBound ||
        !boundary.report.wallClearance.meetsRequestedThickness ||
        below.report.ready ||
        below.report.wallClearance.meetsRequestedThickness ||
        !hasFailure(
            below.report,
            "deformed_wall_below_requested_thickness"
        ) ||
        std::abs(
            boundary.report.wallClearance.measuredMinimumClearance - requested
        ) > 1.0e-10
    )
    {
        std::cerr
            << " measured_above="
            << above.report.wallClearance.measuredMinimumClearance
            << " measured_boundary="
            << boundary.report.wallClearance.measuredMinimumClearance
            << " measured_below="
            << below.report.wallClearance.measuredMinimumClearance;
        return fail(scenario, "deformed wall threshold contract failed");
    }
    return pass(scenario);
}


class ScopedTemporaryDirectory
{
public:
    explicit ScopedTemporaryDirectory(const std::string& scenario)
    {
        const std::filesystem::path base =
            std::filesystem::temp_directory_path();
        for (int attempt = 0; attempt < 100; ++attempt)
        {
            path_ = base / (
                "make_it_stand_r8_" + scenario + "_" +
                std::to_string(static_cast<long long>(::getpid())) + "_" +
                std::to_string(attempt)
            );
            std::error_code error;
            if (std::filesystem::create_directory(path_, error))
            {
                return;
            }
            if (error)
            {
                throw std::runtime_error(
                    "failed to create scoped R8 test directory"
                );
            }
        }
        throw std::runtime_error("could not allocate scoped R8 test directory");
    }

    ~ScopedTemporaryDirectory()
    {
        std::error_code ignored;
        std::filesystem::remove_all(path_, ignored);
    }

    const std::filesystem::path& path() const noexcept
    {
        return path_;
    }

private:
    std::filesystem::path path_;
};


std::vector<unsigned char> readBytes(const std::filesystem::path& path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input)
    {
        throw std::runtime_error("failed to read R8 test output");
    }
    return std::vector<unsigned char>(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>()
    );
}


std::string readText(const std::filesystem::path& path)
{
    const std::vector<unsigned char> bytes = readBytes(path);
    return std::string(bytes.begin(), bytes.end());
}


void writeSentinel(const std::filesystem::path& path)
{
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    output << "preexisting-sentinel\n";
    if (!output)
    {
        throw std::runtime_error("failed to create R8 sentinel");
    }
}


std::size_t countOccurrences(
    const std::string& text,
    const std::string& token
)
{
    std::size_t result = 0;
    std::size_t position = 0;
    while ((position = text.find(token, position)) != std::string::npos)
    {
        ++result;
        position += token.size();
    }
    return result;
}


std::uint32_t readLittleEndian32(
    const std::vector<unsigned char>& bytes,
    std::size_t offset
)
{
    if (offset + 4 > bytes.size())
    {
        throw std::runtime_error("truncated binary STL fixture");
    }
    return
        static_cast<std::uint32_t>(bytes[offset]) |
        (static_cast<std::uint32_t>(bytes[offset + 1]) << 8U) |
        (static_cast<std::uint32_t>(bytes[offset + 2]) << 16U) |
        (static_cast<std::uint32_t>(bytes[offset + 3]) << 24U);
}


float readLittleEndianFloat(
    const std::vector<unsigned char>& bytes,
    std::size_t offset
)
{
    const std::uint32_t bits = readLittleEndian32(bytes, offset);
    float result = 0.0F;
    static_assert(sizeof(result) == sizeof(bits), "STL float must be 32-bit");
    std::memcpy(&result, &bits, sizeof(result));
    return result;
}


TriangleMesh readBinaryStlMesh(const std::vector<unsigned char>& bytes)
{
    if (bytes.size() < 84)
    {
        throw std::runtime_error("binary STL header is truncated");
    }
    const std::uint32_t faceCount = readLittleEndian32(bytes, 80);
    if (bytes.size() != 84U + 50U * faceCount)
    {
        throw std::runtime_error("binary STL byte count is inconsistent");
    }

    using VertexBits = std::array<std::uint32_t, 3>;
    std::map<VertexBits, int> unique;
    std::vector<Eigen::Vector3d> vertices;
    TriangleMesh result;
    result.faces.resize(static_cast<Eigen::Index>(faceCount), 3);

    for (std::uint32_t face = 0; face < faceCount; ++face)
    {
        const std::size_t record = 84U + 50U * face;
        for (int corner = 0; corner < 3; ++corner)
        {
            const std::size_t vertexOffset = record + 12U + 12U * corner;
            VertexBits bits{};
            Eigen::Vector3d point;
            for (int axis = 0; axis < 3; ++axis)
            {
                bits[static_cast<std::size_t>(axis)] = readLittleEndian32(
                    bytes,
                    vertexOffset + 4U * axis
                );
                point(axis) = static_cast<double>(readLittleEndianFloat(
                    bytes,
                    vertexOffset + 4U * axis
                ));
            }

            const auto inserted = unique.emplace(
                bits,
                static_cast<int>(vertices.size())
            );
            if (inserted.second)
            {
                vertices.push_back(point);
            }
            result.faces(face, corner) = inserted.first->second;
        }
    }

    result.vertices.resize(static_cast<Eigen::Index>(vertices.size()), 3);
    for (std::size_t index = 0; index < vertices.size(); ++index)
    {
        result.vertices.row(static_cast<Eigen::Index>(index)) =
            vertices[index].transpose();
    }
    return result;
}


std::uint64_t testOccupancyHash(const VoxelGrid& grid)
{
    std::uint64_t hash = UINT64_C(1469598103934665603);
    std::uint64_t index = 0U;
    const Eigen::Vector3i dimensions = grid.dimensions();
    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                hash ^= index++;
                hash *= UINT64_C(1099511628211);
                hash ^= grid(x, y, z).isFilled ? UINT64_C(1) : UINT64_C(0);
                hash *= UINT64_C(1099511628211);
            }
        }
    }
    return hash;
}


optimization::AlternatingBalanceOptimizationResult sourceResultFixture()
{
    const TriangleMesh outer = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );
    const VoxelGrid grid(outer.vertices, outer.faces, 0.1, 8);
    optimization::AlternatingBalanceOptimizationResult result;
    result.voxelFill = grid;
    result.geometry = identityGeometry(outer, grid);
    result.massProperties = computeSurfaceMassProperties(
        outer.vertices,
        outer.faces
    );
    if (!result.massProperties.tryCenterOfMass(result.centerOfMass, 0.0))
    {
        throw std::runtime_error("source fixture has invalid analytic mass");
    }
    result.geometryEpoch = 3;
    result.targets.allReached = true;
    result.targets.generation = 7;
    optimization::AlternatingCarvingRecord record;
    record.geometryEpoch = result.geometryEpoch;
    record.alphaHash = testOccupancyHash(result.voxelFill);
    result.carvingBlocks.push_back(record);
    return result;
}


ManufacturingGeometry preparedSolidFixture()
{
    return manufacturing::prepareManufacturingGeometry(sourceResultFixture());
}


int runExportRoundTrip()
{
    const std::string scenario = "export-roundtrip";
    const ManufacturingGeometry prepared = preparedSolidFixture();
    if (!prepared.report.ready)
    {
        return fail(scenario, "positive export fixture is not ready");
    }

    const ScopedTemporaryDirectory temporary(scenario);
    const std::filesystem::path obj = temporary.path() / "mesh.obj";
    const std::filesystem::path ascii = temporary.path() / "mesh_ascii.stl";
    const std::filesystem::path binary = temporary.path() / "mesh_binary.stl";
    const std::filesystem::path binarySecond =
        temporary.path() / "mesh_binary_second.stl";
    const std::filesystem::path json = temporary.path() / "readiness.json";

    manufacturing::exportReadyManufacturingMesh(
        prepared, obj.string(), manufacturing::ManufacturingMeshFormat::Obj
    );
    manufacturing::exportReadyManufacturingMesh(
        prepared,
        ascii.string(),
        manufacturing::ManufacturingMeshFormat::AsciiStl
    );
    manufacturing::exportReadyManufacturingMesh(
        prepared,
        binary.string(),
        manufacturing::ManufacturingMeshFormat::BinaryStl
    );
    manufacturing::exportReadyManufacturingMesh(
        prepared,
        binarySecond.string(),
        manufacturing::ManufacturingMeshFormat::BinaryStl
    );
    manufacturing::writeManufacturingReadinessReport(
        prepared.report,
        json.string()
    );

    const std::string objText = readText(obj);
    const std::string asciiText = readText(ascii);
    const std::vector<unsigned char> binaryBytes = readBytes(binary);
    const TriangleMesh binaryMesh = readBinaryStlMesh(binaryBytes);
    const BoundaryAudit binaryAudit = auditOutward(binaryMesh);
    const std::string jsonText = readText(json);

    if (
        countOccurrences(objText, "\nv ") !=
            static_cast<std::size_t>(prepared.boundary.vertices.rows()) ||
        countOccurrences(objText, "\nf ") !=
            static_cast<std::size_t>(prepared.boundary.faces.rows()) ||
        asciiText.rfind("solid make_it_stand_r8\n", 0) != 0 ||
        countOccurrences(asciiText, "facet normal") !=
            static_cast<std::size_t>(prepared.boundary.faces.rows()) ||
        asciiText.find(temporary.path().string()) != std::string::npos ||
        asciiText.find(ascii.filename().string()) != std::string::npos ||
        readLittleEndian32(binaryBytes, 80) !=
            static_cast<std::uint32_t>(prepared.boundary.faces.rows()) ||
        binaryBytes.size() !=
            84U + 50U * static_cast<std::size_t>(prepared.boundary.faces.rows()) ||
        readBytes(binarySecond) != binaryBytes ||
        !binaryAudit.strictBoundaryValid() ||
        binaryAudit.faceComponentCount != 1 ||
        std::abs(binaryAudit.signedVolume - 1.0) > 1.0e-6 ||
        jsonText != prepared.report.toJson() ||
        jsonText.find(temporary.path().string()) != std::string::npos
    )
    {
        return fail(scenario, "export format or round-trip contract failed");
    }
    return pass(scenario);
}


int runExportTransaction()
{
    const std::string scenario = "export-transaction";
    const ManufacturingGeometry prepared = preparedSolidFixture();
    const ScopedTemporaryDirectory temporary(scenario);
    const std::filesystem::path target = temporary.path() / "preserve.obj";
    writeSentinel(target);
    const std::vector<unsigned char> sentinel = readBytes(target);

    ManufacturingGeometry rejected = prepared;
    rejected.report.ready = false;
    bool rejectedNotReady = false;
    try
    {
        manufacturing::exportReadyManufacturingMesh(
            rejected,
            target.string(),
            manufacturing::ManufacturingMeshFormat::Obj
        );
    }
    catch (const std::logic_error&)
    {
        rejectedNotReady = true;
    }

    bool rejectedFormat = false;
    try
    {
        manufacturing::exportReadyManufacturingMesh(
            prepared,
            target.string(),
            static_cast<manufacturing::ManufacturingMeshFormat>(999)
        );
    }
    catch (const std::invalid_argument&)
    {
        rejectedFormat = true;
    }

    const std::filesystem::path occupied = temporary.path() / "occupied";
    std::filesystem::create_directory(occupied);
    bool rejectedDirectoryTarget = false;
    try
    {
        manufacturing::writeManufacturingReadinessReport(
            prepared.report,
            occupied.string()
        );
    }
    catch (const std::runtime_error&)
    {
        rejectedDirectoryTarget = true;
    }

    bool temporaryLeak = false;
    for (const auto& entry : std::filesystem::directory_iterator(temporary.path()))
    {
        if (entry.path().filename().string().find(".r8tmp.") !=
            std::string::npos)
        {
            temporaryLeak = true;
        }
    }

    if (
        !rejectedNotReady ||
        !rejectedFormat ||
        !rejectedDirectoryTarget ||
        readBytes(target) != sentinel ||
        !std::filesystem::is_directory(occupied) ||
        temporaryLeak
    )
    {
        return fail(scenario, "failed export modified state or leaked a temp file");
    }
    return pass(scenario);
}


int runSourceBound()
{
    const std::string scenario = "source-bound";
    const ManufacturingGeometry prepared = preparedSolidFixture();
    if (
        !prepared.report.ready ||
        !prepared.report.sourceBound ||
        !prepared.report.targetsAllReached ||
        !prepared.report.sourceTraceHashConsistent ||
        !prepared.report.sourceMassPropertiesConsistent ||
        prepared.report.geometryEpoch != 3 ||
        prepared.report.targetGeneration != 7 ||
        !prepared.report.failures.empty()
    )
    {
        return fail(scenario, "valid committed R7 result did not bind and seal");
    }
    return pass(scenario);
}


int runExportSeal()
{
    const std::string scenario = "export-seal";
    const ScopedTemporaryDirectory temporary(scenario);
    const std::filesystem::path target = temporary.path() / "preserve.obj";
    writeSentinel(target);
    const std::vector<unsigned char> sentinel = readBytes(target);

    ManufacturingGeometry forged;
    forged.report.ready = true;
    forged.report.sourceBound = true;
    bool forgedRejected = false;
    try
    {
        manufacturing::exportReadyManufacturingMesh(
            forged,
            target.string(),
            manufacturing::ManufacturingMeshFormat::Obj
        );
    }
    catch (const std::logic_error&)
    {
        forgedRejected = true;
    }

    ManufacturingGeometry mutatedBoundary = preparedSolidFixture();
    mutatedBoundary.boundary.vertices(0, 0) += 0.125;
    bool boundaryRejected = false;
    try
    {
        manufacturing::exportReadyManufacturingMesh(
            mutatedBoundary,
            target.string(),
            manufacturing::ManufacturingMeshFormat::Obj
        );
    }
    catch (const std::logic_error&)
    {
        boundaryRejected = true;
    }

    ManufacturingGeometry mutatedReport = preparedSolidFixture();
    mutatedReport.report.centerOfMass(0) += 0.125;
    bool reportRejected = false;
    try
    {
        manufacturing::exportReadyManufacturingMesh(
            mutatedReport,
            target.string(),
            manufacturing::ManufacturingMeshFormat::Obj
        );
    }
    catch (const std::logic_error&)
    {
        reportRejected = true;
    }

    if (
        !forgedRejected ||
        !boundaryRejected ||
        !reportRejected ||
        readBytes(target) != sentinel
    )
    {
        return fail(scenario, "forged or mutated export escaped integrity seal");
    }
    return pass(scenario);
}


int runSourceRejection()
{
    const std::string scenario = "source-rejection";

    auto expect = [&scenario](
        optimization::AlternatingBalanceOptimizationResult source,
        const std::string& expected
    ) -> bool
    {
        const ManufacturingGeometry prepared =
            manufacturing::prepareManufacturingGeometry(source);
        if (
            prepared.report.ready ||
            !prepared.report.sourceBound ||
            !hasFailure(prepared.report, expected)
        )
        {
            std::cerr
                << "[FAIL] scenario=" << scenario
                << " expected=" << expected << '\n';
            return false;
        }
        return true;
    };

    optimization::AlternatingBalanceOptimizationResult invalidEpoch =
        sourceResultFixture();
    invalidEpoch.geometryEpoch = -1;
    invalidEpoch.carvingBlocks.back().geometryEpoch = -1;

    optimization::AlternatingBalanceOptimizationResult wrongHash =
        sourceResultFixture();
    ++wrongHash.carvingBlocks.back().alphaHash;

    optimization::AlternatingBalanceOptimizationResult wrongTraceEpoch =
        sourceResultFixture();
    ++wrongTraceEpoch.carvingBlocks.back().geometryEpoch;

    optimization::AlternatingBalanceOptimizationResult nonfinite =
        sourceResultFixture();
    nonfinite.parameters = Eigen::VectorXd::Constant(
        1,
        std::numeric_limits<double>::quiet_NaN()
    );

    optimization::AlternatingBalanceOptimizationResult unsafe =
        sourceResultFixture();
    unsafe.targets.allReached = false;

    optimization::AlternatingBalanceOptimizationResult wrongMass =
        sourceResultFixture();
    wrongMass.massProperties.mass += 0.125;

    optimization::AlternatingBalanceOptimizationResult wrongMoment =
        sourceResultFixture();
    wrongMoment.massProperties.firstMoment(1) += 0.125;

    optimization::AlternatingBalanceOptimizationResult wrongCenter =
        sourceResultFixture();
    wrongCenter.centerOfMass(2) += 0.125;

    if (
        !expect(invalidEpoch, "source_geometry_epoch_invalid") ||
        !expect(wrongHash, "source_occupancy_trace_hash_mismatch") ||
        !expect(wrongTraceEpoch, "source_geometry_epoch_trace_mismatch") ||
        !expect(nonfinite, "source_parameters_nonfinite") ||
        !expect(unsafe, "source_targets_not_all_reached") ||
        !expect(wrongMass, "source_mass_properties_disagree") ||
        !expect(wrongMoment, "source_mass_properties_disagree") ||
        !expect(wrongCenter, "source_mass_properties_disagree")
    )
    {
        return EXIT_FAILURE;
    }
    return pass(scenario);
}


int runRawSubstitution()
{
    const std::string scenario = "raw-substitution";
    optimization::AlternatingBalanceOptimizationResult source =
        sourceResultFixture();

    TriangleMesh substitute = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d(0.9, 1.0, 1.0)
    );
    source.geometry.outerVertices = substitute.vertices;
    source.geometry.outerFaces = substitute.faces;

    const ManufacturingGeometry raw =
        manufacturing::prepareManufacturingGeometry(
            source.geometry,
            source.voxelFill
        );
    if (
        raw.report.ready ||
        raw.report.sourceBound ||
        !hasFailure(raw.report, "r7_source_result_not_bound")
    )
    {
        return fail(scenario, "raw substituted outer was incorrectly source-bound");
    }

    const ScopedTemporaryDirectory temporary(scenario);
    const std::filesystem::path target = temporary.path() / "preserve.obj";
    writeSentinel(target);
    const std::vector<unsigned char> sentinel = readBytes(target);
    bool rejected = false;
    try
    {
        manufacturing::exportReadyManufacturingMesh(
            raw,
            target.string(),
            manufacturing::ManufacturingMeshFormat::Obj
        );
    }
    catch (const std::logic_error&)
    {
        rejected = true;
    }

    if (!rejected || readBytes(target) != sentinel)
    {
        return fail(scenario, "raw substituted outer escaped export guard");
    }
    return pass(scenario);
}


int runScaleTranslation()
{
    const std::string scenario = "scale-translation";
    const TriangleMesh reference = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );
    const BoundaryAudit baseline = auditOutward(reference);
    if (!baseline.strictBoundaryValid())
    {
        return fail(scenario, "reference cube is invalid");
    }

    const double scales[3] = {1.0e-4, 1.0, 1.0e4};
    const Eigen::Vector3d translations[3] =
    {
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d(1.0e6, -2.0e6, 3.0e6),
        Eigen::Vector3d(-7.0e8, 4.0e8, 2.0e8)
    };

    for (int index = 0; index < 3; ++index)
    {
        TriangleMesh transformed = reference;
        transformed.vertices *= scales[index];
        transformed.vertices.rowwise() += translations[index].transpose();
        const BoundaryAudit audit = auditOutward(transformed);
        const double expectedVolume =
            scales[index] * scales[index] * scales[index];
        const double relativeError = std::abs(
            audit.signedVolume - expectedVolume
        ) / std::max(expectedVolume, std::numeric_limits<double>::min());
        if (
            !audit.strictBoundaryValid() ||
            audit.faceComponentCount != 1 ||
            relativeError > 1.0e-8
        )
        {
            std::cerr
                << " transformed_index=" << index
                << " measured_volume=" << audit.signedVolume
                << " expected_volume=" << expectedVolume;
            return fail(
                scenario,
                "boundary verdict changed under scale/translation"
            );
        }
    }
    return pass(scenario);
}


int runCube()
{
    const std::string scenario = "cube";
    const TriangleMesh outward = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );
    const BoundaryAudit outerAudit = auditOutward(outward);

    const TriangleMesh inward = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones(),
        false
    );
    const BoundaryAudit innerAudit = manufacturing::auditTriangleBoundary(
        inward.vertices,
        inward.faces,
        BoundaryOrientation::InwardNegative
    );

    if (
        !outerAudit.strictBoundaryValid() ||
        outerAudit.faceComponentCount != 1 ||
        std::abs(outerAudit.signedVolume - 1.0) > 1.0e-12 ||
        !innerAudit.strictBoundaryValid() ||
        innerAudit.faceComponentCount != 1 ||
        std::abs(innerAudit.signedVolume + 1.0) > 1.0e-12
    )
    {
        return fail(scenario, "oriented analytic cube contract failed");
    }
    return pass(scenario);
}


int runMissing()
{
    const std::string scenario = "missing";
    TriangleMesh mesh = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );
    mesh.faces.conservativeResize(mesh.faces.rows() - 1, 3);
    const BoundaryAudit audit = auditOutward(mesh);
    if (audit.strictBoundaryValid() || audit.closed || audit.boundaryEdgeCount == 0)
    {
        return fail(scenario, "open boundary was not rejected");
    }
    return pass(scenario);
}


int runReversed()
{
    const std::string scenario = "reversed";
    const TriangleMesh mesh = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones(),
        false
    );
    const BoundaryAudit audit = auditOutward(mesh);
    if (
        audit.strictBoundaryValid() ||
        !audit.consistentlyOriented ||
        audit.componentOrientationsValid ||
        !(audit.signedVolume < 0.0)
    )
    {
        return fail(scenario, "globally reversed shell was not diagnosed");
    }
    return pass(scenario);
}


int runDuplicate()
{
    const std::string scenario = "duplicate";
    TriangleMesh mesh = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );
    const Eigen::Index oldRows = mesh.faces.rows();
    mesh.faces.conservativeResize(oldRows + 1, 3);
    mesh.faces.row(oldRows) = mesh.faces.row(0);
    const BoundaryAudit audit = auditOutward(mesh);
    if (
        audit.strictBoundaryValid() ||
        audit.nonManifoldEdgeCount == 0 ||
        audit.edgeManifold
    )
    {
        return fail(scenario, "duplicate facet was not rejected");
    }
    return pass(scenario);
}


int runDegenerate()
{
    const std::string scenario = "degenerate";
    TriangleMesh mesh = makeBox(
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Ones()
    );
    mesh.faces(0, 2) = mesh.faces(0, 1);
    const BoundaryAudit audit = auditOutward(mesh);
    if (
        audit.strictBoundaryValid() ||
        audit.nondegenerateFaces ||
        audit.degenerateFaceCount == 0
    )
    {
        return fail(scenario, "zero-area facet was not rejected");
    }
    return pass(scenario);
}


int runDoubleOuter()
{
    const std::string scenario = "double-outer";
    const TriangleMesh mesh = join(
        makeBox(Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones()),
        makeBox(
            Eigen::Vector3d(2.0, 0.0, 0.0),
            Eigen::Vector3d(3.0, 1.0, 1.0)
        )
    );
    const BoundaryAudit audit = auditOutward(mesh);
    const bool acceptedAsSingleOuter =
        audit.strictBoundaryValid() && audit.faceComponentCount == 1;
    if (
        acceptedAsSingleOuter ||
        audit.faceComponentCount != 2 ||
        !audit.selfIntersections.intersectionFree
    )
    {
        return fail(scenario, "multiple outer components were not diagnosed");
    }
    return pass(scenario);
}


int runSelfIntersection()
{
    const std::string scenario = "self-intersection";
    const TriangleMesh mesh = join(
        makeBox(Eigen::Vector3d::Zero(), Eigen::Vector3d::Ones()),
        makeBox(
            Eigen::Vector3d(0.5, 0.5, 0.5),
            Eigen::Vector3d(1.5, 1.5, 1.5)
        )
    );
    const BoundaryAudit audit = auditOutward(mesh);
    if (
        audit.strictBoundaryValid() ||
        !audit.selfIntersections.completed ||
        audit.selfIntersections.intersectionFree ||
        audit.selfIntersections.intersectingPairCount == 0
    )
    {
        return fail(scenario, "crossing closed shells escaped intersection audit");
    }
    return pass(scenario);
}


} // namespace


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cerr
            << "usage: manufacturing-geometry-regression "
            << "<cube|missing|reversed|duplicate|degenerate|double-outer|"
            << "self-intersection|prepare-empty|prepare-cavity|stale-inner|"
            << "cell-collapse|cell-swap|material-components|containment|"
            << "wall-threshold|source-bound|export-roundtrip|"
            << "export-transaction|export-seal|source-rejection|"
            << "raw-substitution|scale-translation>\n";
        return EXIT_FAILURE;
    }

    const std::string scenario = argv[1];
    if (scenario == "cube") return runCube();
    if (scenario == "missing") return runMissing();
    if (scenario == "reversed") return runReversed();
    if (scenario == "duplicate") return runDuplicate();
    if (scenario == "degenerate") return runDegenerate();
    if (scenario == "double-outer") return runDoubleOuter();
    if (scenario == "self-intersection") return runSelfIntersection();
    if (scenario == "prepare-empty") return runPrepareEmpty();
    if (scenario == "prepare-cavity") return runPrepareCavity();
    if (scenario == "stale-inner") return runStaleInner();
    if (scenario == "cell-collapse") return runCellMutation(true);
    if (scenario == "cell-swap") return runCellMutation(false);
    if (scenario == "material-components") return runMaterialComponents();
    if (scenario == "containment") return runContainment();
    if (scenario == "wall-threshold") return runWallThreshold();
    if (scenario == "source-bound") return runSourceBound();
    if (scenario == "export-roundtrip") return runExportRoundTrip();
    if (scenario == "export-transaction") return runExportTransaction();
    if (scenario == "export-seal") return runExportSeal();
    if (scenario == "source-rejection") return runSourceRejection();
    if (scenario == "raw-substitution") return runRawSubstitution();
    if (scenario == "scale-translation") return runScaleTranslation();

    return fail(scenario, "unknown scenario");
}
