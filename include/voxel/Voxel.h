#pragma once

#include <Eigen/Core>

#include <set>


class Voxel
{
public:
    Eigen::Vector3d center;

    // Signed distance from the voxel center to the closest mesh surface.
    double distanceToMesh;

    // Position of this voxel inside the regular grid.
    int xIdx;
    int yIdx;
    int zIdx;

    // Indices into the shared voxel-grid corner matrix.
    std::set<int> cornerIndices;

    // true: material remains; false: the voxel has been carved away.
    bool isFilled;

    // Whether the voxel is a valid interior carving candidate.
    bool isInner;

    Voxel()
        :
        center(Eigen::Vector3d::Zero()),
        distanceToMesh(0.0),
        xIdx(-1),
        yIdx(-1),
        zIdx(-1),
        isFilled(true),
        isInner(false)
    {
    }
};
