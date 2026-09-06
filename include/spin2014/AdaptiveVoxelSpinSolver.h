#pragma once

#include "physics/RawMassMoments10.h"
#include "spin/SpinDesignRequest.h"
#include "spin2014/SlqpLbfgsSolver.h"
#include "spin2014/SpinIt2014MomentProblem.h"
#include "voxel/MaterialOccupancy.h"
#include "voxel/VoxelDomain.h"

#include <Eigen/Core>

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>


namespace spin2014
{


enum class AdaptiveSpinStatus
{
    Converged,
    RelaxationDidNotConverge,
    AdaptivePassLimit,
    BinaryValidationFailed,
    ResourceLimitReached,
    InvalidRequest
};


const char* adaptiveSpinStatusToken(AdaptiveSpinStatus status) noexcept;


struct AdaptiveVoxelSpinRequest
{
    std::shared_ptr<const voxel::VoxelDomain> domain;
    spin::SpinMode mode = spin::SpinMode::SpinningTop;
    std::optional<Spin2014Objective> objective;
    int initialLevel = 3;
    int maximumLevel = 7;
    int maximumAdaptivePasses = 8;
    int maximumLeafCount = 1000000;
    std::uint64_t maximumOutputCellCount = 150000000U;
    double initialFill = 0.5;
    double fractionalEpsilon = 0.05;
    double laplacianWeight = 1e-3;
    double centerOfMassWeight = 1.0;
    double inertiaWeight = 1.0;
    double boundaryDensity = 1.0;
    double interiorDensity = 1.0;
    double binaryConstraintTolerance = 1e-6;
    SlqpLbfgsSettings slqp;
};


struct AdaptiveSpinPass
{
    int pass = 0;
    int leafCountBefore = 0;
    int fractionalLeafCount = 0;
    int splitLeafCount = 0;
    int mergedParentCount = 0;
    SlqpStatus slqpStatus = SlqpStatus::InvalidProblem;
    double finalEqualityInfinityNorm = 0.0;
    double finalProjectedGradientInfinityNorm = 0.0;
    int totalLineSearchEvaluations = 0;
    int activeSetTransitions = 0;
    std::string partitionHash;
    std::string slqpTraceHash;
    std::vector<SlqpIteration> slqpTrace;
};


struct AdaptiveVoxelSpinResult
{
    AdaptiveSpinStatus status = AdaptiveSpinStatus::InvalidRequest;
    std::string diagnostic;
    std::string inputHash;
    std::string adaptiveTraceHash;
    std::string geometryHash;
    std::string outputHash;
    Eigen::VectorXd relaxedFill;
    std::vector<double> binaryFill;
    std::vector<unsigned char> binaryVoidMask;
    physics::RawMassMoments10 finalMoments;
    Eigen::VectorXd binaryConstraints;
    double binaryConstraintInfinityNorm = 0.0;
    int binaryProjectionChangedLeaves = 0;
    double binaryProjectionMaxDelta = 0.0;
    Eigen::Vector3i outputGridDimensions = Eigen::Vector3i::Zero();
    voxel::MaterialBoundaryMesh sourceBoundary;
    voxel::MaterialBoundaryMesh designBoundary;
    voxel::MaterialBoundaryMesh materialBoundary;
    std::vector<AdaptiveSpinPass> adaptiveTrace;
};


void validateAdaptiveVoxelSpinRequest(const AdaptiveVoxelSpinRequest& request);
std::string computeAdaptiveVoxelSpinRequestHash(
    const AdaptiveVoxelSpinRequest& request
);
AdaptiveVoxelSpinResult solveAdaptiveVoxelSpin(
    const AdaptiveVoxelSpinRequest& request
);


} // namespace spin2014
