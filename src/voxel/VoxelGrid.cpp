#include "voxel/VoxelGrid.h"

#include "grid_util.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
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
    voxelSize_(Eigen::Vector3d::Zero())
{
}


VoxelGrid::VoxelGrid(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F
)
    :
    voxels_(),
    corners_(0, 3),
    dimensions_(Eigen::Vector3i::Zero()),
    voxelSize_(Eigen::Vector3d::Zero())
{
    build(V, F);
}


void VoxelGrid::build(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F
)
{
    Eigen::MatrixXd centers;

    createVoxelGrid(
        V,
        centers,
        corners_,
        dimensions_,
        voxelSize_
    );

    const double minX = corners_.col(0).minCoeff();
    const double maxX = corners_.col(0).maxCoeff();

    const double minY = corners_.col(1).minCoeff();
    const double maxY = corners_.col(1).maxCoeff();

    const double minZ = corners_.col(2).minCoeff();
    const double maxZ = corners_.col(2).maxCoeff();

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

    // Recreate the container so rebuilding resets all Voxel fields.
    voxels_ = List3d<Voxel>(
        static_cast<std::size_t>(dimensions_(0)),
        static_cast<std::size_t>(dimensions_(1)),
        static_cast<std::size_t>(dimensions_(2))
    );

    for (int rowIndex = 0; rowIndex < centers.rows(); ++rowIndex)
    {
        const Eigen::RowVector3d currentCenter =
            centers.row(rowIndex);

        const int xIndex = getBinIndex(
            minX,
            maxX,
            voxelSize_(0),
            currentCenter(0)
        );

        const int yIndex = getBinIndex(
            minY,
            maxY,
            voxelSize_(1),
            currentCenter(1)
        );

        const int zIndex = getBinIndex(
            minZ,
            maxZ,
            voxelSize_(2),
            currentCenter(2)
        );

        Voxel& voxel =
            voxels_(xIndex, yIndex, zIndex);

        voxel.center =
            currentCenter.transpose();

        voxel.xIdx = xIndex;
        voxel.yIdx = yIndex;
        voxel.zIdx = zIndex;

        voxel.distanceToMesh =
            distances(rowIndex);

        // Retain a one-voxel shell, matching the original implementation.
        if (
            distances(rowIndex) < 0.0 &&
            std::abs(distances(rowIndex)) >= voxelSize_(0)
        )
        {
            // Border voxels are excluded because carving uses six neighbors.
            if (
                xIndex > 0 &&
                xIndex < dimensions_(0) - 1 &&
                yIndex > 0 &&
                yIndex < dimensions_(1) - 1 &&
                zIndex > 0 &&
                zIndex < dimensions_(2) - 1
            )
            {
                voxel.isInner = true;
            }
        }
    }

    const double epsilon =
        voxelSize_(0) / 10.0;

    std::vector<int> xNeighbors(2, -1);
    std::vector<int> yNeighbors(2, -1);
    std::vector<int> zNeighbors(2, -1);

    for (int cornerIndex = 0; cornerIndex < corners_.rows(); ++cornerIndex)
    {
        const Eigen::RowVector3d currentCorner =
            corners_.row(cornerIndex);

        getNeighborIndices(
            minX,
            maxX,
            voxelSize_(0),
            currentCorner(0),
            epsilon,
            xNeighbors
        );

        getNeighborIndices(
            minY,
            maxY,
            voxelSize_(1),
            currentCorner(1),
            epsilon,
            yNeighbors
        );

        getNeighborIndices(
            minZ,
            maxZ,
            voxelSize_(2),
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
                        voxels_(
                            xNeighbors[i],
                            yNeighbors[j],
                            zNeighbors[k]
                        ).cornerIndices.insert(cornerIndex);
                    }
                }
            }
        }
    }
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
