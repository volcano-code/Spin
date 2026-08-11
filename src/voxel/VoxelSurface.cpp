#include "voxel/VoxelSurface.h"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <vector>

#include <igl/per_face_normals.h>


namespace
{

bool areFacingSameDirection(
    const Eigen::RowVector3d& first,
    const Eigen::RowVector3d& second
)
{
    return first.dot(second) > 0.0;
}

} // namespace


namespace voxel_surface
{

bool isVoid(const Voxel& voxel)
{
    return voxel.isInner && !voxel.isFilled;
}


void appendCommonSide(
    const Voxel& voxel1,
    const Voxel& voxel2,
    const Eigen::MatrixXd& vertices,
    Eigen::MatrixXi& outFaces,
    bool onlyNecessary
)
{
    const Voxel* voidVoxel = nullptr;

    if (onlyNecessary)
    {
        const bool voxel1IsVoid = isVoid(voxel1);
        const bool voxel2IsVoid = isVoid(voxel2);

        // Both void or both non-void means there is no cavity boundary here.
        if (voxel1IsVoid == voxel2IsVoid)
        {
            return;
        }

        voidVoxel = voxel1IsVoid
            ? &voxel1
            : &voxel2;
    }
    else
    {
        // Preserve the behavior of the original helper.
        voidVoxel = &voxel1;
    }

    std::vector<int> commonCorners;

    std::set_intersection(
        voxel1.cornerIndices.begin(),
        voxel1.cornerIndices.end(),
        voxel2.cornerIndices.begin(),
        voxel2.cornerIndices.end(),
        std::back_inserter(commonCorners)
    );

    assert(commonCorners.size() == 4);

    Eigen::MatrixXi sideFaces(2, 3);

    sideFaces <<
        Eigen::RowVector3i(
            commonCorners[0],
            commonCorners[1],
            commonCorners[2]
        ),
        Eigen::RowVector3i(
            commonCorners[2],
            commonCorners[1],
            commonCorners[3]
        );

    Eigen::MatrixXd normals;

    igl::per_face_normals(
        vertices,
        sideFaces,
        normals
    );

    assert(
        areFacingSameDirection(
            normals.row(0),
            normals.row(1)
        )
    );

    const Eigen::RowVector3d vectorToSide =
        vertices.row(sideFaces(0, 0)) -
        voidVoxel->center.transpose();

    // The cavity surface is oriented inward relative to remaining material.
    if (
        areFacingSameDirection(
            normals.row(0),
            vectorToSide
        )
    )
    {
        sideFaces <<
            Eigen::RowVector3i(
                commonCorners[2],
                commonCorners[1],
                commonCorners[0]
            ),
            Eigen::RowVector3i(
                commonCorners[3],
                commonCorners[1],
                commonCorners[2]
            );
    }

    const Eigen::Index oldFaceCount =
        outFaces.rows();

    outFaces.conservativeResize(
        oldFaceCount + 2,
        3
    );

    outFaces.block(
        oldFaceCount,
        0,
        2,
        3
    ) = sideFaces;
}

} // namespace voxel_surface
