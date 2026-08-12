#pragma once

#include <Eigen/Core>

#include "physics/BalanceObjective.h"
#include "physics/MassProperties.h"

#include <cstddef>
#include <vector>


/** Stable grid identity used to break exact score ties deterministically. */
struct VoxelKey
{
    int x;
    int y;
    int z;
};


/** A removable voxel together with its exact density-one mass properties. */
struct PlaneCarvingCandidate
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    VoxelKey key;
    Eigen::Vector3d center;
    MassProperties massProperties;
    double score;

    PlaneCarvingCandidate();
};


using PlaneCarvingCandidates = std::vector<
    PlaneCarvingCandidate,
    Eigen::aligned_allocator<PlaneCarvingCandidate>
>;


/**
 * Paper Section 4 cutting-plane state for one target.
 *
 * The raw projected residual is the plane normal in printed Equation (8).
 * A zero residual means that this target is already satisfied and contributes
 * neither a plane nor a candidate score.
 */
struct CuttingPlane
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Eigen::Vector3d point;
    Eigen::Vector3d normal;
    Eigen::Vector3d gravity;
    bool active;

    CuttingPlane();
};


using CuttingPlanes = std::vector<
    CuttingPlane,
    Eigen::aligned_allocator<CuttingPlane>
>;


/** Result of scanning prefix 0 and every valid ordered deletion prefix. */
struct PlaneCarvingPrefixResult
{
    std::size_t evaluatedPrefixCount;
    std::size_t bestPrefixLength;
    double initialEnergy;
    double bestEnergy;
    MassProperties bestMassProperties;
    Eigen::Vector3d bestCenterOfMass;

    PlaneCarvingPrefixResult();
};


/** Construct one paper Equation (8) plane per objective. */
CuttingPlanes makeCuttingPlanes(
    const Eigen::Vector3d& currentCenterOfMass,
    const BalanceObjectives& objectives
);


/** Sum only strictly positive per-plane Equation (8) values. */
double aggregatePositivePlaneScore(
    const Eigen::Vector3d& point,
    const CuttingPlanes& planes
);


/** Sort score descending, then grid key lexicographically ascending. */
void sortPlaneCarvingCandidates(
    PlaneCarvingCandidates& candidates
);


/**
 * Simulate the complete ordered queue without mutating production state.
 * Prefix 0 participates and equal energies retain the shorter prefix.
 */
PlaneCarvingPrefixResult selectBestCarvingPrefix(
    const MassProperties& initialMassProperties,
    const BalanceObjectives& objectives,
    const PlaneCarvingCandidates& orderedCandidates
);
