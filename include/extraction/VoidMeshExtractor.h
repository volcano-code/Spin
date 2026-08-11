#pragma once

#include <Eigen/Core>


class VoxelGrid;


/** Convert the carved voxels of a VoxelGrid into a cavity triangle mesh. */
class VoidMeshExtractor
{
public:
    static void extract(
        const VoxelGrid& voxelGrid,
        Eigen::MatrixXd& vertices,
        Eigen::MatrixXi& faces
    );
};
