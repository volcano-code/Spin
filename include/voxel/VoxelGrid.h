#pragma once

#include <Eigen/Core>

#include "list3d.h"
#include "voxel/Voxel.h"


/**
 * Regular volumetric representation of an input triangle mesh.
 *
 * This class owns the voxel states, the shared corner positions, the grid
 * dimensions, and the voxel size in model coordinates. It deliberately does not contain
 * any carving policy.
 */
class VoxelGrid
{
public:
    VoxelGrid();

    VoxelGrid(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        double minimumWallThickness = 0.0,
        int longestAxisResolution = 50
    );

    /** Build or rebuild the voxel representation from a triangle mesh. */
    void build(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        double minimumWallThickness = 0.0,
        int longestAxisResolution = 50
    );

    Voxel& operator()(
        int x,
        int y,
        int z
    );

    const Voxel& operator()(
        int x,
        int y,
        int z
    ) const;

    const Eigen::MatrixXd& corners() const;

    const Eigen::Vector3i& dimensions() const;

    const Eigen::Vector3d& voxelSize() const;

    /** Requested model-space shell thickness for the fixed outer surface. */
    double minimumWallThickness() const;

    int longestAxisResolution() const;

private:
    List3d<Voxel> voxels_;

    Eigen::MatrixXd corners_;

    Eigen::Vector3i dimensions_;

    Eigen::Vector3d voxelSize_;

    double minimumWallThickness_;

    int longestAxisResolution_;
};
