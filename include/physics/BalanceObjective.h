#pragma once

#include <Eigen/Core>
#include <Eigen/StdVector>

#include <vector>


/** One projected center-of-mass objective with a normalized gravity. */
class BalanceObjective
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    BalanceObjective(
        const Eigen::Vector3d& targetPoint,
        const Eigen::Vector3d& gravityDirection
    );

    const Eigen::Vector3d& target() const;
    const Eigen::Vector3d& gravity() const;

private:
    Eigen::Vector3d target_;
    Eigen::Vector3d gravity_;
};


using BalanceObjectives = std::vector<
    BalanceObjective,
    Eigen::aligned_allocator<BalanceObjective>
>;


/** Construct the projected CoM objective for a selected suspension point. */
BalanceObjective makeSuspensionObjective(
    const Eigen::Vector3d& suspensionPoint,
    const Eigen::Vector3d& gravityDirection
);


/** Square root of the summed squared projected distances. */
double aggregateProjectedCenterOfMassDistance(
    const Eigen::Vector3d& centerOfMass,
    const BalanceObjectives& objectives
);


/** Paper Eq. (5) CoM term: one half the summed squared distances. */
double aggregateProjectedCenterOfMassEnergy(
    const Eigen::Vector3d& centerOfMass,
    const BalanceObjectives& objectives
);
