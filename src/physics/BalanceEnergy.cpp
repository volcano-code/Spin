#include "physics/BalanceEnergy.h"

#include <cmath>
#include <stdexcept>


Eigen::Vector3d normalizedGravityDirection(
    const Eigen::Vector3d& gravity
)
{
    if (!gravity.allFinite())
    {
        throw std::invalid_argument(
            "gravity must be finite and non-zero"
        );
    }

    /*
     * Scaling by the largest component keeps the subsequent norm in
     * [1,sqrt(3)]. This accepts every finite non-zero direction without
     * underflowing tiny inputs or overflowing near DBL_MAX.
     */
    const double scale = gravity.cwiseAbs().maxCoeff();

    if (!(scale > 0.0))
    {
        throw std::invalid_argument(
            "gravity must be finite and non-zero"
        );
    }

    const Eigen::Vector3d scaledGravity = gravity / scale;
    const double scaledNorm = scaledGravity.norm();

    if (!std::isfinite(scaledNorm) || !(scaledNorm > 0.0))
    {
        throw std::invalid_argument(
            "gravity must be finite and non-zero"
        );
    }

    return scaledGravity / scaledNorm;
}


Eigen::Vector3d projectPerpendicularToGravity(
    const Eigen::Vector3d& displacement,
    const Eigen::Vector3d& gravity
)
{
    if (!displacement.allFinite() || !gravity.allFinite())
    {
        throw std::invalid_argument(
            "projected center-of-mass inputs must be finite"
        );
    }

    const Eigen::Vector3d normalizedGravity =
        normalizedGravityDirection(gravity);

    return displacement -
        displacement.dot(normalizedGravity) * normalizedGravity;
}


double projectedCenterOfMassDistance(
    const Eigen::Vector3d& centerOfMass,
    const Eigen::Vector3d& targetCenterOfMass,
    const Eigen::Vector3d& gravity
)
{
    return projectPerpendicularToGravity(
        centerOfMass - targetCenterOfMass,
        gravity
    ).norm();
}


double projectedCenterOfMassEnergy(
    const Eigen::Vector3d& centerOfMass,
    const Eigen::Vector3d& targetCenterOfMass,
    const Eigen::Vector3d& gravity
)
{
    const Eigen::Vector3d projectedDifference =
        projectPerpendicularToGravity(
            centerOfMass - targetCenterOfMass,
            gravity
        );

    return 0.5 * projectedDifference.squaredNorm();
}
