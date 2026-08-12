#pragma once

#include <Eigen/Core>


/** Return a unit gravity direction without overflow or underflow. */
Eigen::Vector3d normalizedGravityDirection(
    const Eigen::Vector3d& gravity
);


/** Project a displacement onto the plane perpendicular to gravity. */
Eigen::Vector3d projectPerpendicularToGravity(
    const Eigen::Vector3d& displacement,
    const Eigen::Vector3d& gravity
);


/** Paper Eq. (2) projected center-of-mass distance (before squaring). */
double projectedCenterOfMassDistance(
    const Eigen::Vector3d& centerOfMass,
    const Eigen::Vector3d& targetCenterOfMass,
    const Eigen::Vector3d& gravity
);


/** Paper Eq. (2): 0.5 * squared projected center-of-mass distance. */
double projectedCenterOfMassEnergy(
    const Eigen::Vector3d& centerOfMass,
    const Eigen::Vector3d& targetCenterOfMass,
    const Eigen::Vector3d& gravity
);
