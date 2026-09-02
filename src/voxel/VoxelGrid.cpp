#include "voxel/VoxelGrid.h"

#include "grid_util.h"
#include "voxel/MaterialOccupancy.h"
#include "voxel/VoxelDomain.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

#include <igl/signed_distance.h>


namespace
{

int getBinIndex(
    double minValue,
    double maxValue,
    double binSize,
    double currentLocation
)
{
    assert(currentLocation >= minValue);
    assert(currentLocation <= maxValue);

    return static_cast<int>(
        std::floor(
            (currentLocation - minValue) / binSize
        )
    );
}


void getNeighborIndices(
    double minValue,
    double maxValue,
    double binSize,
    double currentBorder,
    double epsilon,
    std::vector<int>& output
)
{
    assert(output.size() >= 2);

    assert(
        currentBorder >= minValue - epsilon &&
        currentBorder <= maxValue + epsilon
    );

    const int borderIndex = static_cast<int>(
        std::round(
            (currentBorder - minValue) / binSize
        )
    );

    assert(
        std::abs(
            currentBorder -
            (minValue + borderIndex * binSize)
        ) < epsilon
    );

    const int maxIndex = static_cast<int>(
        std::round(
            (maxValue - minValue) / binSize
        )
    );

    if (borderIndex == maxIndex)
    {
        output[0] = borderIndex - 1;
        output[1] = -1;
    }
    else if (borderIndex == 0)
    {
        output[0] = 0;
        output[1] = -1;
    }
    else
    {
        output[0] = borderIndex - 1;
        output[1] = borderIndex;
    }
}

} // namespace


VoxelGrid VoxelGrid::fromMaterialOccupancy(
    const std::shared_ptr<const voxel::VoxelDomain>& domain,
    double minimumWallThickness,
    int protectedBoundaryLayers
)
{
    if (!domain)
    {
        throw std::invalid_argument("material occupancy requires a domain");
    }
    if (
        !std::isfinite(minimumWallThickness) ||
        minimumWallThickness < 0.0 ||
        protectedBoundaryLayers < 0
    )
    {
        throw std::invalid_argument(
            "material occupancy wall settings are invalid"
        );
    }
    const voxel::MaterialBoundaryMesh boundaries =
        voxel::extractMaterialBoundaries(domain);
    if (boundaries.materialComponentCount != 1)
    {
        throw std::invalid_argument(
            "material occupancy must contain one six-connected component"
        );
    }

    const Eigen::Vector3i dimensions = domain->grid().dimensions;
    const Eigen::Vector3i cornerDimensions =
        dimensions + Eigen::Vector3i::Ones();
    const std::size_t cellCount = domain->cellCount();
    Eigen::MatrixXd corners(
        static_cast<Eigen::Index>(domain->cornerCount()),
        3
    );
    for (int x = 0; x < cornerDimensions(0); ++x)
    {
        for (int y = 0; y < cornerDimensions(1); ++y)
        {
            for (int z = 0; z < cornerDimensions(2); ++z)
            {
                corners.row(static_cast<Eigen::Index>(
                    domain->cornerIndex(x, y, z)
                )) = domain->cornerPosition(x, y, z).transpose();
            }
        }
    }

    static constexpr int offsets[6][3] = {
        {-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}
    };
    const auto inBounds = [&dimensions](int x, int y, int z)
    {
        return
            x >= 0 && x < dimensions(0) &&
            y >= 0 && y < dimensions(1) &&
            z >= 0 && z < dimensions(2);
    };

    // Identify the exterior complement so enclosed complement cells can be
    // represented as pre-existing cavities in the legacy facade.
    std::vector<unsigned char> exterior(cellCount, 0U);
    std::queue<Eigen::Vector3i> pending;
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
                    pending.emplace(x, y, z);
                }
            }
        }
    }
    while (!pending.empty())
    {
        const Eigen::Vector3i cell = pending.front();
        pending.pop();
        for (const auto& offset : offsets)
        {
            const int nx = cell(0) + offset[0];
            const int ny = cell(1) + offset[1];
            const int nz = cell(2) + offset[2];
            if (!inBounds(nx, ny, nz))
            {
                continue;
            }
            const std::size_t index = domain->cellIndex(nx, ny, nz);
            if (!domain->material(index) && exterior[index] == 0U)
            {
                exterior[index] = 1U;
                pending.emplace(nx, ny, nz);
            }
        }
    }

    // Six-neighbour material depth protects both the outer surface and every
    // pre-existing cavity. A cell can be carved only beyond both the explicit
    // layer count and requested metric wall thickness.
    std::vector<int> depth(cellCount, -1);
    std::queue<Eigen::Vector3i> depthQueue;
    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                const std::size_t index = domain->cellIndex(x, y, z);
                if (!domain->material(index))
                {
                    continue;
                }
                bool boundary = false;
                for (const auto& offset : offsets)
                {
                    const int nx = x + offset[0];
                    const int ny = y + offset[1];
                    const int nz = z + offset[2];
                    if (
                        !inBounds(nx, ny, nz) ||
                        !domain->material(domain->cellIndex(nx, ny, nz))
                    )
                    {
                        boundary = true;
                        break;
                    }
                }
                if (boundary)
                {
                    depth[index] = 0;
                    depthQueue.emplace(x, y, z);
                }
            }
        }
    }
    while (!depthQueue.empty())
    {
        const Eigen::Vector3i cell = depthQueue.front();
        depthQueue.pop();
        const int current = depth[domain->cellIndex(
            cell(0), cell(1), cell(2)
        )];
        for (const auto& offset : offsets)
        {
            const int nx = cell(0) + offset[0];
            const int ny = cell(1) + offset[1];
            const int nz = cell(2) + offset[2];
            if (!inBounds(nx, ny, nz))
            {
                continue;
            }
            const std::size_t index = domain->cellIndex(nx, ny, nz);
            if (domain->material(index) && depth[index] < 0)
            {
                depth[index] = current + 1;
                depthQueue.emplace(nx, ny, nz);
            }
        }
    }

    List3d<Voxel> voxels(
        static_cast<std::size_t>(dimensions(0)),
        static_cast<std::size_t>(dimensions(1)),
        static_cast<std::size_t>(dimensions(2))
    );
    std::vector<unsigned char> baseline(cellCount, 1U);
    const Eigen::Vector3d spacing = domain->spacing();
    const double halfDiagonal = 0.5 * spacing.norm();
    const double metricLayer = spacing.minCoeff();
    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                Voxel& cell = voxels(
                    static_cast<std::size_t>(x),
                    static_cast<std::size_t>(y),
                    static_cast<std::size_t>(z)
                );
                const std::size_t canonical = domain->cellIndex(x, y, z);
                const std::size_t legacy = static_cast<std::size_t>(
                    (x * dimensions(1) + y) * dimensions(2) + z
                );
                const bool material = domain->material(canonical);
                const bool cavity = !material && exterior[canonical] == 0U;
                cell.center = domain->cellCenter(x, y, z);
                cell.xIdx = x;
                cell.yIdx = y;
                cell.zIdx = z;
                cell.distanceToMesh = (material || cavity)
                    ? -halfDiagonal
                    : 2.0 * halfDiagonal;
                cell.isFilled = !cavity;
                baseline[legacy] = cell.isFilled ? 1U : 0U;
                const double conservativeClearance =
                    depth[canonical] < 0
                        ? 0.0
                        : static_cast<double>(depth[canonical]) * metricLayer;
                cell.isInner = cavity || (
                    material &&
                    domain->carveable(canonical) &&
                    depth[canonical] > protectedBoundaryLayers &&
                    conservativeClearance >= minimumWallThickness &&
                    x > 0 && x < dimensions(0) - 1 &&
                    y > 0 && y < dimensions(1) - 1 &&
                    z > 0 && z < dimensions(2) - 1
                );
                for (int dx = 0; dx <= 1; ++dx)
                {
                    for (int dy = 0; dy <= 1; ++dy)
                    {
                        for (int dz = 0; dz <= 1; ++dz)
                        {
                            cell.cornerIndices.insert(static_cast<int>(
                                domain->cornerIndex(x + dx, y + dy, z + dz)
                            ));
                        }
                    }
                }
            }
        }
    }

    VoxelGrid result;
    result.voxels_ = std::move(voxels);
    result.corners_ = std::move(corners);
    result.dimensions_ = dimensions;
    result.voxelSize_ = spacing;
    result.minimumWallThickness_ = minimumWallThickness;
    result.longestAxisResolution_ = dimensions.maxCoeff();
    result.fixedUnitCubeDomain_ =
        domain->grid().minimum.isApprox(Eigen::Vector3d::Zero(), 0.0) &&
        domain->grid().maximum.isApprox(Eigen::Vector3d::Ones(), 0.0) &&
        dimensions(0) == dimensions(1) && dimensions(1) == dimensions(2);
    result.baselineFilled_ = std::move(baseline);
    return result;
}


VoxelGrid VoxelGrid::fromReleasedBetaOccupancy(
    int resolution,
    const std::vector<int>& cellIdsZFast,
    int hullDepth
)
{
    if (resolution < 3)
    {
        throw std::invalid_argument(
            "released-beta voxel resolution must be at least three"
        );
    }
    if (hullDepth < 0)
    {
        throw std::invalid_argument(
            "released-beta hull depth must be nonnegative"
        );
    }

    const long long side64 = static_cast<long long>(resolution);
    const long long maximumIndex =
        static_cast<long long>(std::numeric_limits<int>::max());
    if (
        side64 > maximumIndex / side64 ||
        side64 * side64 > maximumIndex / side64
    )
    {
        throw std::overflow_error(
            "released-beta voxel resolution exceeds integer indexing"
        );
    }
    const long long expectedCount64 = side64 * side64 * side64;
    if (
        expectedCount64 >
            static_cast<long long>(std::numeric_limits<int>::max()) ||
        cellIdsZFast.size() != static_cast<std::size_t>(expectedCount64)
    )
    {
        throw std::invalid_argument(
            "released-beta .vox cache must contain exactly N^3 entries"
        );
    }

    std::size_t occupiedCount = 0U;
    for (const int id : cellIdsZFast)
    {
        if (id < -1)
        {
            throw std::invalid_argument(
                "released-beta .vox entries must be -1 or nonnegative IDs"
            );
        }
        occupiedCount += id >= 0 ? 1U : 0U;
    }
    if (occupiedCount == 0U)
    {
        throw std::invalid_argument(
            "released-beta .vox cache contains no occupied cells"
        );
    }
    if (
        occupiedCount >
        static_cast<std::size_t>(std::numeric_limits<int>::max())
    )
    {
        throw std::overflow_error(
            "released-beta occupied-cell count exceeds integer indexing"
        );
    }

    std::vector<unsigned char> seen(occupiedCount, 0U);
    for (const int id : cellIdsZFast)
    {
        if (id < 0)
        {
            continue;
        }
        if (
            static_cast<std::size_t>(id) >= occupiedCount ||
            seen[static_cast<std::size_t>(id)] != 0U
        )
        {
            throw std::invalid_argument(
                "released-beta occupied IDs must be unique and contiguous"
            );
        }
        seen[static_cast<std::size_t>(id)] = 1U;
    }

    Eigen::MatrixXd centers;
    Eigen::MatrixXd corners;
    Eigen::Vector3i dimensions;
    Eigen::Vector3d voxelSize;
    createUnitCubeVoxelGrid(
        centers,
        corners,
        dimensions,
        voxelSize,
        resolution
    );

    List3d<Voxel> voxels(
        static_cast<std::size_t>(resolution),
        static_cast<std::size_t>(resolution),
        static_cast<std::size_t>(resolution)
    );
    const int cornerResolution = resolution + 1;
    const double halfDiagonal = 0.5 * voxelSize.norm();
    const double occupiedDistance = -halfDiagonal;
    const double inactiveDistance = 2.0 * halfDiagonal;

    const auto betaCellIndex = [resolution](int x, int y, int z)
    {
        return static_cast<std::size_t>(
            (x * resolution + y) * resolution + z
        );
    };
    const auto currentCornerIndex = [cornerResolution](int x, int y, int z)
    {
        return
            x + cornerResolution *
            (y + cornerResolution * z);
    };

    for (int x = 0; x < resolution; ++x)
    {
        for (int y = 0; y < resolution; ++y)
        {
            for (int z = 0; z < resolution; ++z)
            {
                Voxel& voxel = voxels(
                    static_cast<std::size_t>(x),
                    static_cast<std::size_t>(y),
                    static_cast<std::size_t>(z)
                );
                voxel.xIdx = x;
                voxel.yIdx = y;
                voxel.zIdx = z;
                voxel.center = Eigen::Vector3d(
                    (static_cast<double>(x) + 0.5) * voxelSize(0),
                    (static_cast<double>(y) + 0.5) * voxelSize(1),
                    (static_cast<double>(z) + 0.5) * voxelSize(2)
                );
                const bool occupied =
                    cellIdsZFast[betaCellIndex(x, y, z)] >= 0;
                voxel.distanceToMesh = occupied
                    ? occupiedDistance
                    : inactiveDistance;
                voxel.isFilled = true;
                voxel.isInner = false;
                for (int dx = 0; dx <= 1; ++dx)
                {
                    for (int dy = 0; dy <= 1; ++dy)
                    {
                        for (int dz = 0; dz <= 1; ++dz)
                        {
                            voxel.cornerIndices.insert(currentCornerIndex(
                                x + dx,
                                y + dy,
                                z + dz
                            ));
                        }
                    }
                }
            }
        }
    }

    /*
     * Match BoxGrid::initStructure(): a 26-neighbour exposure seeds depth
     * zero, then six-neighbour BFS assigns the interior depth. Only cells at
     * depth strictly greater than hullDepth may be carved.
     */
    std::vector<int> depth(
        static_cast<std::size_t>(expectedCount64),
        -1
    );
    std::queue<Eigen::Vector3i> frontier;
    const auto occupiedAt = [&](int x, int y, int z)
    {
        return
            x >= 0 && x < resolution &&
            y >= 0 && y < resolution &&
            z >= 0 && z < resolution &&
            cellIdsZFast[betaCellIndex(x, y, z)] >= 0;
    };
    for (int x = 0; x < resolution; ++x)
    {
        for (int y = 0; y < resolution; ++y)
        {
            for (int z = 0; z < resolution; ++z)
            {
                if (!occupiedAt(x, y, z))
                {
                    continue;
                }
                bool hull = false;
                for (int dx = -1; dx <= 1 && !hull; ++dx)
                {
                    for (int dy = -1; dy <= 1 && !hull; ++dy)
                    {
                        for (int dz = -1; dz <= 1; ++dz)
                        {
                            if (dx == 0 && dy == 0 && dz == 0)
                            {
                                continue;
                            }
                            if (!occupiedAt(x + dx, y + dy, z + dz))
                            {
                                hull = true;
                                break;
                            }
                        }
                    }
                }
                if (hull)
                {
                    depth[betaCellIndex(x, y, z)] = 0;
                    frontier.emplace(x, y, z);
                }
            }
        }
    }

    constexpr int neighborOffsets[6][3] =
    {
        {-1, 0, 0}, {1, 0, 0},
        {0, -1, 0}, {0, 1, 0},
        {0, 0, -1}, {0, 0, 1}
    };
    while (!frontier.empty())
    {
        const Eigen::Vector3i cell = frontier.front();
        frontier.pop();
        const int cellDepth = depth[betaCellIndex(
            cell(0), cell(1), cell(2)
        )];
        for (const auto& offset : neighborOffsets)
        {
            const int nx = cell(0) + offset[0];
            const int ny = cell(1) + offset[1];
            const int nz = cell(2) + offset[2];
            if (!occupiedAt(nx, ny, nz))
            {
                continue;
            }
            int& neighborDepth = depth[betaCellIndex(nx, ny, nz)];
            if (neighborDepth < 0)
            {
                neighborDepth = cellDepth + 1;
                frontier.emplace(nx, ny, nz);
            }
        }
    }

    for (int x = 0; x < resolution; ++x)
    {
        for (int y = 0; y < resolution; ++y)
        {
            for (int z = 0; z < resolution; ++z)
            {
                if (!occupiedAt(x, y, z))
                {
                    continue;
                }
                const int cellDepth = depth[betaCellIndex(x, y, z)];
                if (cellDepth < 0)
                {
                    throw std::invalid_argument(
                        "released-beta occupied component has no exposed hull"
                    );
                }
                voxels(
                    static_cast<std::size_t>(x),
                    static_cast<std::size_t>(y),
                    static_cast<std::size_t>(z)
                ).isInner = cellDepth > hullDepth;
            }
        }
    }

    VoxelGrid result;
    result.voxels_ = std::move(voxels);
    result.corners_ = std::move(corners);
    result.dimensions_ = dimensions;
    result.voxelSize_ = voxelSize;
    result.minimumWallThickness_ = 0.0;
    result.longestAxisResolution_ = resolution;
    result.fixedUnitCubeDomain_ = true;
    result.baselineFilled_.assign(
        static_cast<std::size_t>(expectedCount64),
        1U
    );
    return result;
}


VoxelGrid::VoxelGrid()
    :
    voxels_(),
    corners_(0, 3),
    dimensions_(Eigen::Vector3i::Zero()),
    voxelSize_(Eigen::Vector3d::Zero()),
    minimumWallThickness_(0.0),
    longestAxisResolution_(50),
    fixedUnitCubeDomain_(false)
{
}


VoxelGrid::VoxelGrid(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    double minimumWallThickness,
    int longestAxisResolution,
    bool fixedUnitCubeDomain
)
    :
    voxels_(),
    corners_(0, 3),
    dimensions_(Eigen::Vector3i::Zero()),
    voxelSize_(Eigen::Vector3d::Zero()),
    minimumWallThickness_(0.0),
    longestAxisResolution_(50),
    fixedUnitCubeDomain_(false)
{
    build(
        V,
        F,
        minimumWallThickness,
        longestAxisResolution,
        fixedUnitCubeDomain
    );
}


void VoxelGrid::build(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    double minimumWallThickness,
    int longestAxisResolution,
    bool fixedUnitCubeDomain
)
{
    if (
        !std::isfinite(minimumWallThickness) ||
        minimumWallThickness < 0.0
    )
    {
        throw std::invalid_argument(
            "minimum wall thickness must be finite and nonnegative"
        );
    }

    if (longestAxisResolution < 3)
    {
        throw std::invalid_argument(
            "voxel resolution must be at least three"
        );
    }

    Eigen::MatrixXd centers;
    Eigen::MatrixXd newCorners;
    Eigen::Vector3i newDimensions;
    Eigen::Vector3d newVoxelSize;

    if (fixedUnitCubeDomain)
    {
        const double coordinateTolerance =
            64.0 * std::numeric_limits<double>::epsilon();
        if (
            V.rows() == 0 ||
            V.cols() != 3 ||
            !V.allFinite() ||
            (V.array() < -coordinateTolerance).any() ||
            (V.array() > 1.0 + coordinateTolerance).any()
        )
        {
            throw std::invalid_argument(
                "fixed unit-cube voxel grid requires mesh vertices in [0,1]^3"
            );
        }
        createUnitCubeVoxelGrid(
            centers,
            newCorners,
            newDimensions,
            newVoxelSize,
            longestAxisResolution
        );
    }
    else
    {
        createVoxelGrid(
            V,
            centers,
            newCorners,
            newDimensions,
            newVoxelSize,
            longestAxisResolution
        );
    }

    const double minX = newCorners.col(0).minCoeff();
    const double maxX = newCorners.col(0).maxCoeff();

    const double minY = newCorners.col(1).minCoeff();
    const double maxY = newCorners.col(1).maxCoeff();

    const double minZ = newCorners.col(2).minCoeff();
    const double maxZ = newCorners.col(2).maxCoeff();

    Eigen::VectorXd distances;
    Eigen::VectorXd closestFaces;
    Eigen::MatrixXd closestPoints;
    Eigen::MatrixXd closestNormals;

    igl::signed_distance(
        centers,
        V,
        F,
        igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER,
        std::numeric_limits<double>::min(),
        std::numeric_limits<double>::max(),
        distances,
        closestFaces,
        closestPoints,
        closestNormals
    );

    if (
        distances.size() != centers.rows() ||
        !distances.allFinite()
    )
    {
        throw std::runtime_error(
            "voxel signed-distance evaluation returned invalid values"
        );
    }

    List3d<Voxel> newVoxels(
        static_cast<std::size_t>(newDimensions(0)),
        static_cast<std::size_t>(newDimensions(1)),
        static_cast<std::size_t>(newDimensions(2))
    );

    for (int rowIndex = 0; rowIndex < centers.rows(); ++rowIndex)
    {
        const Eigen::RowVector3d currentCenter =
            centers.row(rowIndex);

        const int xIndex = getBinIndex(
            minX,
            maxX,
            newVoxelSize(0),
            currentCenter(0)
        );

        const int yIndex = getBinIndex(
            minY,
            maxY,
            newVoxelSize(1),
            currentCenter(1)
        );

        const int zIndex = getBinIndex(
            minZ,
            maxZ,
            newVoxelSize(2),
            currentCenter(2)
        );

        Voxel& voxel =
            newVoxels(xIndex, yIndex, zIndex);

        voxel.center =
            currentCenter.transpose();

        voxel.xIdx = xIndex;
        voxel.yIdx = yIndex;
        voxel.zIdx = zIndex;

        voxel.distanceToMesh =
            distances(rowIndex);

        /*
         * The signed-distance function is 1-Lipschitz. For an axis-aligned
         * voxel with half diagonal h, -sdf(center)-h is therefore a
         * conservative lower bound on every point's clearance from M_O.
         * Keeping only boxes whose lower bound is at least t_min protects the
         * requested model-space shell during fixed-geometry carving.
         */
        const double conservativeClearance =
            -distances(rowIndex) - 0.5 * newVoxelSize.norm();

        if (
            distances(rowIndex) < 0.0 &&
            conservativeClearance >= minimumWallThickness
        )
        {
            // Border voxels are excluded because carving uses six neighbors.
            if (
                xIndex > 0 &&
                xIndex < newDimensions(0) - 1 &&
                yIndex > 0 &&
                yIndex < newDimensions(1) - 1 &&
                zIndex > 0 &&
                zIndex < newDimensions(2) - 1
            )
            {
                voxel.isInner = true;
            }
        }
    }

    const double epsilon =
        newVoxelSize(0) / 10.0;

    std::vector<int> xNeighbors(2, -1);
    std::vector<int> yNeighbors(2, -1);
    std::vector<int> zNeighbors(2, -1);

    for (int cornerIndex = 0; cornerIndex < newCorners.rows(); ++cornerIndex)
    {
        const Eigen::RowVector3d currentCorner =
            newCorners.row(cornerIndex);

        getNeighborIndices(
            minX,
            maxX,
            newVoxelSize(0),
            currentCorner(0),
            epsilon,
            xNeighbors
        );

        getNeighborIndices(
            minY,
            maxY,
            newVoxelSize(1),
            currentCorner(1),
            epsilon,
            yNeighbors
        );

        getNeighborIndices(
            minZ,
            maxZ,
            newVoxelSize(2),
            currentCorner(2),
            epsilon,
            zNeighbors
        );

        for (int i = 0; i < 2; ++i)
        {
            for (int j = 0; j < 2; ++j)
            {
                for (int k = 0; k < 2; ++k)
                {
                    if (
                        xNeighbors[i] != -1 &&
                        yNeighbors[j] != -1 &&
                        zNeighbors[k] != -1
                    )
                    {
                        newVoxels(
                            xNeighbors[i],
                            yNeighbors[j],
                            zNeighbors[k]
                        ).cornerIndices.insert(cornerIndex);
                    }
                }
            }
        }
    }

    voxels_ = std::move(newVoxels);
    corners_ = std::move(newCorners);
    dimensions_ = newDimensions;
    voxelSize_ = newVoxelSize;
    minimumWallThickness_ = minimumWallThickness;
    longestAxisResolution_ = longestAxisResolution;
    fixedUnitCubeDomain_ = fixedUnitCubeDomain;
    baselineFilled_.assign(
        static_cast<std::size_t>(
            newDimensions(0) * newDimensions(1) * newDimensions(2)
        ),
        1U
    );
}


Voxel& VoxelGrid::operator()(
    int x,
    int y,
    int z
)
{
    return voxels_(
        static_cast<std::size_t>(x),
        static_cast<std::size_t>(y),
        static_cast<std::size_t>(z)
    );
}


const Voxel& VoxelGrid::operator()(
    int x,
    int y,
    int z
) const
{
    return voxels_(
        static_cast<std::size_t>(x),
        static_cast<std::size_t>(y),
        static_cast<std::size_t>(z)
    );
}


const Eigen::MatrixXd& VoxelGrid::corners() const
{
    return corners_;
}


const Eigen::Vector3i& VoxelGrid::dimensions() const
{
    return dimensions_;
}


const Eigen::Vector3d& VoxelGrid::voxelSize() const
{
    return voxelSize_;
}


double VoxelGrid::minimumWallThickness() const
{
    return minimumWallThickness_;
}


int VoxelGrid::longestAxisResolution() const
{
    return longestAxisResolution_;
}


bool VoxelGrid::usesFixedUnitCubeDomain() const
{
    return fixedUnitCubeDomain_;
}


int VoxelGrid::refillAll()
{
    int refilled = 0;

    for (int x = 0; x < dimensions_(0); ++x)
    {
        for (int y = 0; y < dimensions_(1); ++y)
        {
            for (int z = 0; z < dimensions_(2); ++z)
            {
                Voxel& voxel = (*this)(x, y, z);
                const std::size_t index = static_cast<std::size_t>(
                    (x * dimensions_(1) + y) * dimensions_(2) + z
                );
                const bool baseline = baselineFilled_.empty() ||
                    baselineFilled_[index] != 0U;
                if (voxel.isFilled != baseline)
                {
                    voxel.isFilled = baseline;
                    ++refilled;
                }
            }
        }
    }

    return refilled;
}
