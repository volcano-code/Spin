#include "extraction/VoidMeshExtractor.h"

#include "voxel/VoxelGrid.h"
#include "voxel/VoxelSurface.h"


void VoidMeshExtractor::extract(
    const VoxelGrid& voxelGrid,
    Eigen::MatrixXd& vertices,
    Eigen::MatrixXi& faces
)
{
    // Preserve the original representation: retain all shared grid corners.
    vertices = voxelGrid.corners();

    extractFaces(voxelGrid, faces);
}


void VoidMeshExtractor::extractFaces(
    const VoxelGrid& voxelGrid,
    Eigen::MatrixXi& faces
)
{
    const Eigen::MatrixXd& corners =
        voxelGrid.corners();

    const Eigen::Vector3i& dimensions =
        voxelGrid.dimensions();

    faces.resize(0, 3);

    for (int x = 1; x < dimensions(0) - 1; ++x)
    {
        for (int y = 1; y < dimensions(1) - 1; ++y)
        {
            for (int z = 1; z < dimensions(2) - 1; ++z)
            {
                const Voxel& currentVoxel =
                    voxelGrid(x, y, z);

                if (!voxel_surface::isVoid(currentVoxel))
                {
                    continue;
                }

                Eigen::MatrixXi currentFaces(0, 3);

                voxel_surface::appendCommonSide(
                    voxelGrid(x - 1, y, z),
                    currentVoxel,
                    corners,
                    currentFaces
                );

                voxel_surface::appendCommonSide(
                    voxelGrid(x + 1, y, z),
                    currentVoxel,
                    corners,
                    currentFaces
                );

                voxel_surface::appendCommonSide(
                    voxelGrid(x, y - 1, z),
                    currentVoxel,
                    corners,
                    currentFaces
                );

                voxel_surface::appendCommonSide(
                    voxelGrid(x, y + 1, z),
                    currentVoxel,
                    corners,
                    currentFaces
                );

                voxel_surface::appendCommonSide(
                    voxelGrid(x, y, z - 1),
                    currentVoxel,
                    corners,
                    currentFaces
                );

                voxel_surface::appendCommonSide(
                    voxelGrid(x, y, z + 1),
                    currentVoxel,
                    corners,
                    currentFaces
                );

                const Eigen::Index newFaceCount =
                    currentFaces.rows();

                if (newFaceCount == 0)
                {
                    continue;
                }

                const Eigen::Index oldFaceCount =
                    faces.rows();

                faces.conservativeResize(
                    oldFaceCount + newFaceCount,
                    3
                );

                faces.block(
                    oldFaceCount,
                    0,
                    newFaceCount,
                    3
                ) = currentFaces;
            }
        }
    }
}
