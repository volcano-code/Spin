#pragma once

#include <Eigen/Core>

#include "voxel/Voxel.h"


namespace voxel_surface
{

/** Return true when a voxel represents carved internal empty space. */
bool isVoid(const Voxel& voxel);


/**
 * Triangulate the common face of two neighboring voxels.
 *
 * When onlyNecessary is true, two triangles are appended only when exactly
 * one of the two voxels is carved void. This produces the cavity boundary
 * rather than internal faces between two void voxels.
 */
void appendCommonSide(
    const Voxel& voxel1,
    const Voxel& voxel2,
    const Eigen::MatrixXd& vertices,
    Eigen::MatrixXi& outFaces,
    bool onlyNecessary = true
);

} // namespace voxel_surface
