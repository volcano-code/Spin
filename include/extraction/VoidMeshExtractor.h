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

    /**
     * Extract only the fixed corner-index connectivity of the cavity.
     *
     * Face orientation is determined in the VoxelGrid's rest geometry, so
     * callers may replace the returned full corner array by a synchronously
     * deformed one without recomputing topology.
     */
    static void extractFaces(
        const VoxelGrid& voxelGrid,
        Eigen::MatrixXi& faces
    );
};
