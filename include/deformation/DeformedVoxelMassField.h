#pragma once

#include "physics/MassProperties.h"

#include <Eigen/Core>

#include <stdexcept>


namespace deformation
{


class GridDeformationDomain;


/** A finite deformed cell whose oriented volume is non-positive or collapsed. */
class InfeasibleDeformedVoxel : public std::domain_error
{
public:
    using std::domain_error::domain_error;
};


/**
 * Exact mass sidecar for the active cells of a fixed deformation lattice.
 *
 * VoxelGrid remains the immutable topology/fill container. This class binds
 * the domain's fixed eight-corner connectivity to one synchronously deformed
 * full-corner array and integrates every active cell as a closed, outward-
 * oriented 12-triangle hexahedron. The resulting centroid is the exact
 * uniform-volume centroid of that piecewise-linear boundary, not the average
 * of its eight corners.
 */
class DeformedVoxelMassField
{
public:
    DeformedVoxelMassField(
        const GridDeformationDomain& domain,
        const Eigen::MatrixXd& fullDeformedCorners,
        double density = 1.0,
        double minimumPositiveCellMass = 0.0
    );

    Eigen::Index cellCount() const noexcept;
    Eigen::Index activeCellCount() const noexcept;
    const Eigen::Vector3i& dimensions() const noexcept;

    /** True only for cells represented by the deformation object domain. */
    bool contains(int linearCell) const;

    /** Exact positive mass and first moment of one active deformed cell. */
    MassProperties massProperties(int linearCell) const;

    /** Exact uniform-volume centroid, equal to firstMoment / mass. */
    Eigen::Vector3d centroid(int linearCell) const;

    double density() const noexcept;
    double minimumPositiveCellMass() const noexcept;

private:
    Eigen::Vector3i dimensions_ = Eigen::Vector3i::Zero();
    Eigen::VectorXi activeMask_;
    Eigen::VectorXd masses_;
    Eigen::MatrixXd firstMoments_;
    Eigen::MatrixXd centroids_;
    Eigen::Index activeCellCount_ = 0;
    double density_ = 1.0;
    double minimumPositiveCellMass_ = 0.0;
};


} // namespace deformation
