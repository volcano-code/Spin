#pragma once

#include "physics/RawMassMoments10.h"
#include "spin/SpinDesignRequest.h"
#include "spin2014/SlqpLbfgsSolver.h"
#include "voxel/MaterialOccupancy.h"
#include "voxel/VoxelDomain.h"

#include <Eigen/Core>

#include <memory>
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
    InvalidRequest
};


const char* adaptiveSpinStatusToken(AdaptiveSpinStatus status) noexcept;


struct AdaptiveVoxelSpinRequest
{
    std::shared_ptr<const voxel::VoxelDomain> domain;
    spin::SpinMode mode = spin::SpinMode::SpinningTop;
    int initialLevel = 3;
    int maximumLevel = 7;
    int maximumAdaptivePasses = 8;
    double initialFill = 0.5;
    double fractionalEpsilon = 0.05;
    double laplacianWeight = 1e-3;
    double centerOfMassWeight = 1.0;
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
    std::string partitionHash;
    std::string slqpTraceHash;
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
