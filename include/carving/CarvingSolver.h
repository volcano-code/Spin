#pragma once

#include <Eigen/Core>

#include "carving/CarvingDiagnostics.h"
#include "physics/BalanceObjective.h"
#include "physics/MassProperties.h"

#include <string>
#include <vector>


class VoxelGrid;
class Voxel;


/**
 * Implements the mass-redistribution carving procedure.
 *
 * The solver owns the optimization state (mass, first moment, derived center
 * of mass, iteration count and carved-voxel count), while the referenced
 * VoxelGrid owns all voxel data. The grid must outlive the solver and must not
 * be mutated externally while the solver is active. The solver deliberately
 * does not perform visualization or cavity-mesh extraction.
 */
class CarvingSolver
{
public:
    /**
     * Construct a solver around one or more projected CoM objectives.
     *
     * Aggregate objective queries support any non-empty collection. Actual
     * multi-target carving is deliberately deferred to R4.
     */
    CarvingSolver(
        VoxelGrid& voxelGrid,
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const BalanceObjectives& objectives
    );

    /** Compatibility constructor for the existing single-target viewer. */
    CarvingSolver(
        VoxelGrid& voxelGrid,
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXd& planeV,
        const Eigen::MatrixXi& planeF,
        const Eigen::Vector3d& gravity,
        const Eigen::Vector3d& targetCenterOfMass
    );

    CarvingSolver(const CarvingSolver&) = delete;
    CarvingSolver& operator=(const CarvingSolver&) = delete;

    /**
     * Perform one carving iteration.
     *
     * Callers must check !isOptimized() first. The viewer and carveAll()
     * enforce this precondition; direct low-level callers own it.
     */
    void carveStep(
        Eigen::MatrixXd& carvePlaneV,
        Eigen::MatrixXi& carvePlaneF,
        Eigen::Vector3d& centerOfMass
    );

    /** Continue carving until an active stopping condition is reached. */
    void carveAll();

    /** Root-sum-square projected distance; Eq. 2 distance for one target. */
    double getCoMDistance() const;

    /** Summed projected CoM energy: paper Eq. 2 for K=1, Eq. 5 for K>1. */
    double getCoMEnergy() const;

    /** Compatibility alias for getCoMDistance(). */
    double getGeometricProjectedDistance() const;

    /** Return the one-target stopping condition; K>1 carving is deferred. */
    CarvingStopReason stopReason() const;

    /** Return true for a one-target stop; K>1 carving is deferred. */
    bool isOptimized() const;

    const Eigen::Vector3d& currentCenterOfMass() const;

    /** Compatibility accessor for the first (and historically only) target. */
    const Eigen::Vector3d& targetCenterOfMass() const;

    /** The non-empty objective collection used by every energy query. */
    const BalanceObjectives& balanceObjectives() const;

    const MassProperties& initialMassProperties() const;
    const MassProperties& currentMassProperties() const;

    /** Current remaining signed mass, not the initial mass. */
    double objectMass() const;
    int iterationCount() const;
    int carvedVoxelCount() const;

    /** Enable or disable one-line diagnostic output after every step. */
    void setDiagnosticsEnabled(bool enabled);
    bool diagnosticsEnabled() const;

    /** All step records collected since solver construction. */
    const std::vector<CarvingStepDiagnostics>& diagnosticsHistory() const;

    /** Null when no carving step has been executed yet. */
    const CarvingStepDiagnostics* lastDiagnostics() const;

    /** Write the collected diagnostics to a CSV file. */
    bool writeDiagnosticsCsv(
        const std::string& filePath
    ) const;

    int initialInsideVoxelCount() const;
    int initialCarveableVoxelCount() const;

private:
    struct RecomputedCarvingState
    {
        MassProperties massProperties;
        Eigen::Vector3d centerOfMass;
        int filledInsideVoxelCount;

        RecomputedCarvingState();
    };

    void ensureDiagnosticsInitialized();

    RecomputedCarvingState
    recomputeCarvingState() const;

    void subtractVoxelMassProperties(
        const Voxel& voxel
    );

    void refreshCenterOfMass();

    double carvedGridRatio() const;
    double carvedCarveableRatio() const;

    VoxelGrid& voxelGrid_;

    BalanceObjectives objectives_;
    Eigen::Vector3d currCOM_;

    MassProperties initialMassProperties_;
    MassProperties currentMassProperties_;

    int carvedVoxelCount_;
    int iterations_;

    bool candidateSearchPerformed_;
    int lastCandidateVoxelCount_;

    bool diagnosticsEnabled_;
    bool diagnosticsInitialized_;

    int initialInsideVoxelCount_;
    int initialCarveableVoxelCount_;

    std::vector<CarvingStepDiagnostics> diagnosticsHistory_;
};
