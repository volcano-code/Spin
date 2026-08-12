#include "physics/BalanceObjective.h"

#include "physics/BalanceEnergy.h"

#include <cmath>
#include <limits>
#include <stdexcept>


namespace
{

void validateAggregateInputs(
    const Eigen::Vector3d& centerOfMass,
    const BalanceObjectives& objectives
)
{
    if (!centerOfMass.allFinite())
    {
        throw std::invalid_argument(
            "center of mass must be finite"
        );
    }

    if (objectives.empty())
    {
        throw std::invalid_argument(
            "at least one balance objective is required"
        );
    }
}


double aggregateSquaredProjectedDistance(
    const Eigen::Vector3d& centerOfMass,
    const BalanceObjectives& objectives
)
{
    validateAggregateInputs(centerOfMass, objectives);

    double squaredDistanceSum = 0.0;

    for (const BalanceObjective& objective : objectives)
    {
        const Eigen::Vector3d residual =
            projectPerpendicularToGravity(
                centerOfMass - objective.target(),
                objective.gravity()
            );

        const double distance = residual.norm();

        if (!std::isfinite(distance))
        {
            throw std::overflow_error(
                "projected balance objective overflowed"
            );
        }

        const double squaredDistance = distance * distance;

        if (
            !std::isfinite(squaredDistance) ||
            squaredDistanceSum >
                std::numeric_limits<double>::max() -
                squaredDistance
        )
        {
            throw std::overflow_error(
                "aggregate balance objective overflowed"
            );
        }

        squaredDistanceSum += squaredDistance;
    }

    return squaredDistanceSum;
}

} // namespace


BalanceObjective::BalanceObjective(
    const Eigen::Vector3d& targetPoint,
    const Eigen::Vector3d& gravityDirection
)
    :
    target_(targetPoint),
    gravity_(Eigen::Vector3d::Zero())
{
    if (!target_.allFinite())
    {
        throw std::invalid_argument(
            "balance target must be finite"
        );
    }

    gravity_ = normalizedGravityDirection(
        gravityDirection
    );
}


const Eigen::Vector3d& BalanceObjective::target() const
{
    return target_;
}


const Eigen::Vector3d& BalanceObjective::gravity() const
{
    return gravity_;
}


BalanceObjective makeSuspensionObjective(
    const Eigen::Vector3d& suspensionPoint,
    const Eigen::Vector3d& gravityDirection
)
{
    return BalanceObjective(
        suspensionPoint,
        gravityDirection
    );
}


double aggregateProjectedCenterOfMassDistance(
    const Eigen::Vector3d& centerOfMass,
    const BalanceObjectives& objectives
)
{
    return std::sqrt(
        aggregateSquaredProjectedDistance(
            centerOfMass,
            objectives
        )
    );
}


double aggregateProjectedCenterOfMassEnergy(
    const Eigen::Vector3d& centerOfMass,
    const BalanceObjectives& objectives
)
{
    return 0.5 * aggregateSquaredProjectedDistance(
        centerOfMass,
        objectives
    );
}
