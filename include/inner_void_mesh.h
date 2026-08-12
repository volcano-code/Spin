#pragma once

#include <Eigen/Core>

#include "carving/CarvingSolver.h"
#include "voxel/VoxelGrid.h"


/**
 * Compatibility facade used by the existing viewer.
 *
 * VoxelGrid owns the volumetric data, CarvingSolver owns optimization state,
 * and VoidMeshExtractor (called in the implementation file) creates the
 * cavity triangle mesh. Keeping this facade preserves the inherited viewer
 * API while the solver is replaced in audited stages.
 */
class InnerVoidMesh
{
private:
    // Declaration order is important: CarvingSolver stores a reference to it.
    VoxelGrid voxelGrid;

    CarvingSolver carvingSolver;

public:
    /** Construct the facade with a non-empty balance-objective collection. */
    InnerVoidMesh(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const BalanceObjectives& objectives,
        double minimumWallThickness = 0.0
    );

    /** Compatibility constructor for the existing single-target viewer. */
    InnerVoidMesh(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXd& planeV,
        const Eigen::MatrixXi& planeF,
        const Eigen::Vector3d& gravity,
        const Eigen::Vector3d& balancePoint
    );

    void convertToMesh(
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F
    );

    void carveInnerMeshStep(
        Eigen::MatrixXd& carvePlaneV,
        Eigen::MatrixXi& carvePlaneF,
        Eigen::Vector3d& com
    );

    double getCoMDistance();

    const BalanceObjectives& balanceObjectives() const;

    bool isOptimized();

    void carveInnerMesh();
};
