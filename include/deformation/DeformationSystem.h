#pragma once

#include <Eigen/Core>

#include "deformation/BoundedBiharmonicWeights.h"
#include "deformation/GridDeformationDomain.h"
#include "deformation/HandleSet.h"
#include "voxel/VoxelGrid.h"


namespace deformation
{


/** A coherent, immutable-result snapshot of one Equation (9) evaluation. */
struct DeformedGeometry
{
    Eigen::MatrixXd outerVertices;
    Eigen::MatrixXi outerFaces;

    /** Full regular-grid corner array in the VoxelGrid's global indexing. */
    Eigen::MatrixXd gridCorners;

    /** One deformed eight-corner barycenter per VoxelGrid cell. */
    Eigen::MatrixXd voxelCenters;

    /** Full corner array plus the current carved-cavity connectivity. */
    Eigen::MatrixXd innerVertices;
    Eigen::MatrixXi innerFaces;
};


/**
 * Precomputed Section 5 deformation representation.
 *
 * The object stores immutable rest geometry and weights. evaluate() returns a
 * value snapshot and never mutates VoxelGrid, so the R4 axis-aligned carving
 * state cannot become partially or silently stale during R5.
 */
class DeformationSystem
{
public:
    DeformationSystem(
        const VoxelGrid& voxelGrid,
        const Eigen::MatrixXd& outerRestVertices,
        const Eigen::MatrixXi& outerFaces,
        const HandleSet& handles,
        const WeightSolveOptions& options = WeightSolveOptions()
    );

    /** Build a system around already verified node weights (test/import API). */
    DeformationSystem(
        const VoxelGrid& voxelGrid,
        const Eigen::MatrixXd& outerRestVertices,
        const Eigen::MatrixXi& outerFaces,
        const HandleSet& handles,
        const Eigen::MatrixXd& nodeWeights,
        bool verifiedWeightImport
    );

    const GridDeformationDomain& domain() const;
    const HandleSet& bindingHandles() const;
    const Eigen::MatrixXd& nodeWeights() const;
    const Eigen::MatrixXd& fullCornerWeights() const;
    const Eigen::MatrixXd& outerWeights() const;
    const Eigen::MatrixXd& outerRestVertices() const;
    const Eigen::MatrixXi& outerFaces() const;
    const WeightSolveDiagnostics& weightDiagnostics() const;

    /**
     * Validate that a fill snapshot still belongs to the immutable rest grid.
     * This is the allocation-free contract check used by fixed-alpha setup.
     */
    void validateVoxelFill(const VoxelGrid& currentVoxelFill) const;

    /**
     * Evaluate M_O, M_I and the voxel lattice from the same immutable rest
     * state. currentVoxelFill contributes only the cavity faces/fill topology.
     */
    DeformedGeometry evaluate(
        const HandleSet& currentHandles,
        const VoxelGrid& currentVoxelFill
    ) const;

private:
    void initializeVerifiedWeights(
        const Eigen::MatrixXd& nodeWeights
    );

    GridDeformationDomain domain_;
    VoxelGrid restVoxelGrid_;
    HandleSet bindingHandles_;
    Eigen::MatrixXd outerRestVertices_;
    Eigen::MatrixXi outerFaces_;
    Eigen::MatrixXd nodeWeights_;
    Eigen::MatrixXd fullCornerWeights_;
    Eigen::MatrixXd outerWeights_;
    WeightSolveDiagnostics weightDiagnostics_;
};


} // namespace deformation
