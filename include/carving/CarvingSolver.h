#pragma once

#include <Eigen/Core>

#include "carving/CarvingDiagnostics.h"
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

    /** Paper Eq. (2) projected distance before the 0.5*square operation. */
    double getCoMDistance() const;

    /** Paper Eq. (2): 0.5 times the squared projected CoM distance. */
    double getCoMEnergy() const;

    /** Compatibility alias for getCoMDistance(). */
    double getGeometricProjectedDistance() const;

    /** Return the currently active stopping condition. */
    CarvingStopReason stopReason() const;

    /** Return true when any active stopping condition is satisfied. */
    bool isOptimized() const;

    const Eigen::Vector3d& currentCenterOfMass() const;
    const Eigen::Vector3d& targetCenterOfMass() const;

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

    Eigen::Vector3d gravity_;
    Eigen::Vector3d targetCOM_;
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
