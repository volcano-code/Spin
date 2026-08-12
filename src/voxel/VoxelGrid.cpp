#include "voxel/VoxelGrid.h"

#include "grid_util.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
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


VoxelGrid::VoxelGrid()
    :
    voxels_(),
    corners_(0, 3),
    dimensions_(Eigen::Vector3i::Zero()),
    voxelSize_(Eigen::Vector3d::Zero()),
    minimumWallThickness_(0.0),
    longestAxisResolution_(50)
{
}


VoxelGrid::VoxelGrid(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    double minimumWallThickness,
    int longestAxisResolution
)
    :
    voxels_(),
    corners_(0, 3),
    dimensions_(Eigen::Vector3i::Zero()),
    voxelSize_(Eigen::Vector3d::Zero()),
    minimumWallThickness_(0.0),
    longestAxisResolution_(50)
{
    build(V, F, minimumWallThickness, longestAxisResolution);
}


void VoxelGrid::build(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    double minimumWallThickness,
    int longestAxisResolution
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

    createVoxelGrid(
        V,
        centers,
        newCorners,
        newDimensions,
        newVoxelSize,
        longestAxisResolution
    );

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
