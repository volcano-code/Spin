#pragma once

#include "voxel/VoxelDomain.h"

#include <Eigen/Core>

#include <memory>


namespace voxel
{


/** Outer and cavity boundaries reconstructed deterministically from cells. */
struct MaterialBoundaryMesh
{
    Eigen::MatrixXd outerVertices;
    Eigen::MatrixXi outerFaces;
    Eigen::MatrixXd innerVertices;
    Eigen::MatrixXi innerFaces;
    int materialComponentCount = 0;
    int cavityComponentCount = 0;
};


/**
 * Extract the exact axis-aligned boundary represented by a cell-centered
 * material mask. Outer faces point out of material; cavity faces consequently
 * carry the inward signed-volume orientation required by the manufacturing
 * audit.
 */
MaterialBoundaryMesh extractMaterialBoundaries(
    const std::shared_ptr<const VoxelDomain>& domain
);


} // namespace voxel
