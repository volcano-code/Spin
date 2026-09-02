#pragma once

#include "deformation/BoundedBiharmonicWeights.h"
#include "deformation/HandleSet.h"
#include "voxel/VoxelDomain.h"

#include <Eigen/Core>

#include <memory>
#include <string>
#include <vector>


namespace spin2014
{


/** Inputs for the deterministic reconstructed cage/BBW binding stage. */
struct AutomaticCageBindingRequest
{
    std::shared_ptr<const voxel::VoxelDomain> domain;
    Eigen::MatrixXd outerVertices;
    std::vector<Eigen::MatrixXd> voidCellVertices;
    int maximumHandleCount = 16;
    deformation::WeightSolveOptions weightOptions;
};


/**
 * A sparse deformation cage sampled from the material-grid boundary and its
 * bounded-biharmonic weights.  This is a deterministic reconstruction of the
 * paper's cage-binding stage; author-selected handles are not public assets.
 */
struct AutomaticCageBinding
{
    deformation::HandleVector handles;
    Eigen::MatrixXd nodeWeights;
    Eigen::MatrixXd outerWeights;
    std::vector<Eigen::MatrixXd> voidCellWeights;
    deformation::WeightSolveDiagnostics weightDiagnostics;
    int boundaryCandidateCount = 0;
    std::string bindingHash;
};


AutomaticCageBinding buildAutomaticCageBinding(
    const AutomaticCageBindingRequest& request
);


} // namespace spin2014
