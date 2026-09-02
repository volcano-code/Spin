#pragma once

#include <Eigen/Core>

#include "carving/CarvingDiagnostics.h"
#include "carving/PlaneCarvingPolicy.h"
#include "physics/BalanceObjective.h"
#include "physics/MassProperties.h"

#include <functional>
#include <string>
#include <vector>


class VoxelGrid;


namespace deformation
{
class DeformedVoxelMassField;
}


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
     * Every objective contributes one Section 4 cutting plane and one
     * Equation (5) term.
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

    /**
     * Construct a carving solver for one synchronously deformed lattice.
     *
     * The outer mass properties and per-cell sidecar must describe the same
     * deformed geometry. VoxelGrid contributes only immutable topology and
     * transactional fill flags; its rest centers and axis-aligned voxel size
     * are deliberately not used for carving mass or ordering in this mode.
     * The mass field is referenced, not copied, and must outlive the solver.
     */
    CarvingSolver(
        VoxelGrid& voxelGrid,
        const MassProperties& deformedOuterMassProperties,
        const deformation::DeformedVoxelMassField& deformedVoxelMassField,
        const BalanceObjectives& objectives
    );

    CarvingSolver(
        VoxelGrid&,
        const MassProperties&,
        deformation::DeformedVoxelMassField&&,
        const BalanceObjectives&
    ) = delete;

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

    /** Return the active paper-progress or local iteration-guard condition. */
    CarvingStopReason stopReason() const;

    /** Return true once any active stopping condition has been reached. */
    bool isOptimized() const;

    const Eigen::Vector3d& currentCenterOfMass() const;

    /** Compatibility accessor for the first (and historically only) target. */
    const Eigen::Vector3d& targetCenterOfMass() const;

    /** The non-empty objective collection used by every energy query. */
    const BalanceObjectives& balanceObjectives() const;

    /**
     * Replace the targets used by the next carving step.
     *
     * This is a scheduling operation between transactional carveStep() calls:
     * it changes no voxel flag, mass property, center of mass, counter or
     * terminal state. Existing diagnostics remain historical records under
     * the targets active when their steps were committed; they are not
     * reweighted retroactively. A terminal solver is deliberately not
     * reopened. Reentrant updates during carveStep() are rejected. The current
     * energy/distance accessors immediately use the new targets. Empty input
     * is rejected with the strong exception guarantee.
     */
    void updateBalanceObjectives(
        const BalanceObjectives& objectives
    );

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

    /** Local operationalization of the paper's unspecified "too few" rule. */
    void setMinimumRemovalsToContinue(int count);
    int minimumRemovalsToContinue() const;

    /**
     * Reproduce the released beta's candidate ordering for deformed cells.
     * Plane scores then use the mean of eight deformed corners, while
     * committed mass and first moments remain exactly integrated.
     */
    void setUseDeformedCornerAverageCandidatePositions(bool enabled);
    bool usesDeformedCornerAverageCandidatePositions() const noexcept;

private:
    friend struct CarvingSolverTestAccess;

    void setPreCommitHookForTesting(
        const std::function<void()>& hook
    );
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

    void refreshCenterOfMass();

    void validateInitialGridState() const;

    void initializeDeformedVoxelMassField(
        const deformation::DeformedVoxelMassField& deformedVoxelMassField
    );

    int linearCellIndex(
        int x,
        int y,
        int z
    ) const;

    MassProperties voxelMassProperties(
        int x,
        int y,
        int z
    ) const;

    Eigen::Vector3d voxelMassCentroid(
        int x,
        int y,
        int z
    ) const;

    void buildDebugPlaneMesh(
        const CuttingPlanes& planes,
        Eigen::MatrixXd& vertices,
        Eigen::MatrixXi& faces
    ) const;

    double carvedGridRatio() const;
    double carvedCarveableRatio() const;

    VoxelGrid& voxelGrid_;

    BalanceObjectives objectives_;
    Eigen::Vector3d currCOM_;

    MassProperties initialMassProperties_;
    MassProperties currentMassProperties_;

    /** Null selects the legacy axis-aligned rest-cell mass oracle. */
    const deformation::DeformedVoxelMassField* deformedVoxelMassField_;
    bool useDeformedCornerAverageCandidatePositions_;

    int carvedVoxelCount_;
    int iterations_;

    int minimumRemovalsToContinue_;
    CarvingStopReason terminalReason_;

    bool diagnosticsEnabled_;
    bool diagnosticsInitialized_;

    int initialInsideVoxelCount_;
    int initialCarveableVoxelCount_;

    std::vector<CarvingStepDiagnostics> diagnosticsHistory_;
    std::function<void()> preCommitHookForTesting_;
    bool carveStepInProgress_;
};
