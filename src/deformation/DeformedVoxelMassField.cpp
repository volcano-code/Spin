#include "deformation/DeformedVoxelMassField.h"

#include "deformation/GridDeformationDomain.h"

#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>


namespace deformation
{
namespace
{


/*
 * Corner bit order is x + 2*y + 4*z. The face diagonals match the sorted-
 * global-corner convention used by VoxelSurface/VoidMeshExtractor, while
 * every triangle below points outward from an orientation-preserving cell.
 */
constexpr int OUTWARD_TRIANGLES[12][3] =
{
    {2, 1, 0}, {3, 1, 2}, // z = 0, -z
    {4, 5, 6}, {6, 5, 7}, // z = 1, +z
    {0, 1, 4}, {4, 1, 5}, // y = 0, -y
    {6, 3, 2}, {7, 3, 6}, // y = 1, +y
    {4, 2, 0}, {6, 2, 4}, // x = 0, -x
    {1, 3, 5}, {5, 3, 7}  // x = 1, +x
};


MassMoments integrateCell(
    const std::array<Eigen::Vector3d, 8>& corners,
    double density
)
{
    Eigen::Vector3d reference = Eigen::Vector3d::Zero();
    for (const Eigen::Vector3d& corner : corners)
    {
        reference += corner;
    }
    reference /= 8.0;

    MassMoments result;

    for (const auto& triangle : OUTWARD_TRIANGLES)
    {
        const Eigen::Vector3d first =
            corners[static_cast<std::size_t>(triangle[0])] - reference;
        const Eigen::Vector3d second =
            corners[static_cast<std::size_t>(triangle[1])] - reference;
        const Eigen::Vector3d third =
            corners[static_cast<std::size_t>(triangle[2])] - reference;

        const double signedMass =
            density * first.dot(second.cross(third)) / 6.0;
        const Eigen::Vector3d tetrahedronCenter =
            reference + (first + second + third) / 4.0;

        const Eigen::Vector3d tetrahedronVertices[4] =
        {
            reference,
            corners[static_cast<std::size_t>(triangle[0])],
            corners[static_cast<std::size_t>(triangle[1])],
            corners[static_cast<std::size_t>(triangle[2])]
        };
        Eigen::Vector3d vertexSum = Eigen::Vector3d::Zero();
        Eigen::Matrix3d selfOuter = Eigen::Matrix3d::Zero();
        for (const Eigen::Vector3d& vertex : tetrahedronVertices)
        {
            vertexSum += vertex;
            selfOuter += vertex * vertex.transpose();
        }

        result.mass += signedMass;
        result.firstMoment += signedMass * tetrahedronCenter;
        result.secondMoment += signedMass * (
            vertexSum * vertexSum.transpose() + selfOuter
        ) / 20.0;
    }

    return result;
}


void requireCellIndex(
    int linearCell,
    Eigen::Index cellCount
)
{
    if (linearCell < 0 || linearCell >= cellCount)
    {
        throw std::out_of_range(
            "linear cell index is outside the deformed mass field"
        );
    }
}


} // namespace


DeformedVoxelMassField::DeformedVoxelMassField(
    const GridDeformationDomain& domain,
    const Eigen::MatrixXd& fullDeformedCorners,
    double density,
    double minimumPositiveCellMass
)
    :
    dimensions_(domain.dimensions()),
    activeMask_(domain.activeCellMask()),
    masses_(Eigen::VectorXd::Zero(domain.cellCount())),
    firstMoments_(Eigen::MatrixXd::Zero(domain.cellCount(), 3)),
    secondMoments_(Eigen::MatrixXd::Zero(domain.cellCount(), 9)),
    centroids_(Eigen::MatrixXd::Zero(domain.cellCount(), 3)),
    cornerAverages_(Eigen::MatrixXd::Zero(domain.cellCount(), 3)),
    activeCellCount_(domain.activeCellCount()),
    density_(density),
    minimumPositiveCellMass_(minimumPositiveCellMass)
{
    if (
        fullDeformedCorners.rows() != domain.fullRestCorners().rows() ||
        fullDeformedCorners.cols() != 3 ||
        !fullDeformedCorners.allFinite()
    )
    {
        throw std::invalid_argument(
            "deformed mass field requires one finite 3D position per full grid corner"
        );
    }

    if (
        !std::isfinite(density_) ||
        density_ <= 0.0 ||
        !std::isfinite(minimumPositiveCellMass_) ||
        minimumPositiveCellMass_ < 0.0
    )
    {
        throw std::invalid_argument(
            "deformed mass-field density must be positive and its mass floor nonnegative"
        );
    }

    const double restCellMass =
        density_ * domain.voxelSize().prod();
    const double numericalMassFloor =
        1024.0 * std::numeric_limits<double>::epsilon() *
        std::abs(restCellMass);
    const double requiredMass = std::max(
        minimumPositiveCellMass_,
        numericalMassFloor
    );

    for (Eigen::Index cell = 0; cell < domain.cellCount(); ++cell)
    {
        if (activeMask_(cell) == 0)
        {
            continue;
        }

        const std::array<int, 8> globalCorners =
            domain.cellGlobalCorners(static_cast<int>(cell));
        std::array<Eigen::Vector3d, 8> corners{};

        for (int bit = 0; bit < 8; ++bit)
        {
            const int global = globalCorners[static_cast<std::size_t>(bit)];
            if (global < 0 || global >= fullDeformedCorners.rows())
            {
                throw std::logic_error(
                    "deformation cell references an invalid full-grid corner"
                );
            }
            corners[static_cast<std::size_t>(bit)] =
                fullDeformedCorners.row(global).transpose();
            cornerAverages_.row(cell) +=
                corners[static_cast<std::size_t>(bit)].transpose() / 8.0;
        }

        const MassMoments moments = integrateCell(corners, density_);
        const MassProperties properties(moments.mass, moments.firstMoment);
        Eigen::Vector3d center = Eigen::Vector3d::Zero();

        if (
            !properties.isFinite() ||
            !(properties.mass > requiredMass) ||
            !properties.tryCenterOfMass(center, requiredMass)
        )
        {
            std::ostringstream message;
            message
                << "active deformed voxel " << cell
                << " is inverted, collapsed, or has non-finite mass properties";
            throw InfeasibleDeformedVoxel(message.str());
        }

        masses_(cell) = properties.mass;
        firstMoments_.row(cell) = properties.firstMoment.transpose();
        for (int row = 0; row < 3; ++row)
        {
            for (int column = 0; column < 3; ++column)
            {
                secondMoments_(cell, 3 * row + column) =
                    moments.secondMoment(row, column);
            }
        }
        centroids_.row(cell) = center.transpose();
    }
}


MassMoments DeformedVoxelMassField::massMoments(int linearCell) const
{
    requireCellIndex(linearCell, cellCount());
    if (activeMask_(linearCell) == 0)
    {
        throw std::invalid_argument(
            "inactive grid cell has no deformed mass moments"
        );
    }

    MassMoments result;
    result.mass = masses_(linearCell);
    result.firstMoment = firstMoments_.row(linearCell).transpose();
    for (int row = 0; row < 3; ++row)
    {
        for (int column = 0; column < 3; ++column)
        {
            result.secondMoment(row, column) =
                secondMoments_(linearCell, 3 * row + column);
        }
    }
    return result;
}


Eigen::Index DeformedVoxelMassField::cellCount() const noexcept
{
    return masses_.size();
}


Eigen::Index DeformedVoxelMassField::activeCellCount() const noexcept
{
    return activeCellCount_;
}


const Eigen::Vector3i&
DeformedVoxelMassField::dimensions() const noexcept
{
    return dimensions_;
}


bool DeformedVoxelMassField::contains(int linearCell) const
{
    requireCellIndex(linearCell, cellCount());
    return activeMask_(linearCell) != 0;
}


MassProperties DeformedVoxelMassField::massProperties(
    int linearCell
) const
{
    requireCellIndex(linearCell, cellCount());
    if (activeMask_(linearCell) == 0)
    {
        throw std::invalid_argument(
            "inactive grid cell has no deformed mass properties"
        );
    }

    return MassProperties(
        masses_(linearCell),
        firstMoments_.row(linearCell).transpose()
    );
}


Eigen::Vector3d DeformedVoxelMassField::centroid(
    int linearCell
) const
{
    requireCellIndex(linearCell, cellCount());
    if (activeMask_(linearCell) == 0)
    {
        throw std::invalid_argument(
            "inactive grid cell has no deformed centroid"
        );
    }
    return centroids_.row(linearCell).transpose();
}


Eigen::Vector3d DeformedVoxelMassField::cornerAverage(
    int linearCell
) const
{
    requireCellIndex(linearCell, cellCount());
    if (activeMask_(linearCell) == 0)
    {
        throw std::invalid_argument(
            "inactive grid cell has no deformed corner average"
        );
    }
    return cornerAverages_.row(linearCell).transpose();
}


double DeformedVoxelMassField::density() const noexcept
{
    return density_;
}


double DeformedVoxelMassField::minimumPositiveCellMass() const noexcept
{
    return minimumPositiveCellMass_;
}


} // namespace deformation
