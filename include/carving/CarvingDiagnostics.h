#pragma once

#include <Eigen/Core>

#include <iosfwd>


/**
 * Stopping conditions exposed by the carving solver.
 *
 * R1 adds an explicit no-candidate terminal state so a completed candidate
 * search cannot leave the production loop spinning through no-op iterations.
 * The legacy target, empty-ratio and iteration thresholds remain unchanged.
 */
enum class CarvingStopReason
{
    None,
    TargetToleranceReached,
    MaximumEmptyRatioReached,
    NoCandidatesAvailable,
    MaximumIterationsReached
};


/**
 * Diagnostic record for one call to CarvingSolver::carveStep().
 *
 * R2 records both the incremental mass-property state and a full recompute
 * from the immutable outer surface minus every carved voxel. Their
 * discrepancies are expected to remain near floating-point roundoff.
 */
struct CarvingStepDiagnostics
{
    int iteration;

    Eigen::Vector3d centerOfMassBefore;
    Eigen::Vector3d centerOfMassAfter;
    Eigen::Vector3d recomputedCenterOfMassAfter;

    Eigen::Vector3d firstMomentBefore;
    Eigen::Vector3d firstMomentAfter;
    Eigen::Vector3d recomputedFirstMomentAfter;

    double projectedDistanceBefore;
    double projectedDistanceAfter;

    double projectedEnergyBefore;
    double projectedEnergyAfter;

    double centerOfMassDiscrepancyAfter;
    double firstMomentDiscrepancyAfter;

    /** Density-one signed volume proxy, not a calibrated SI mass. */
    double initialMass;
    double remainingMassBefore;
    double remainingMassAfter;
    double recomputedRemainingMassAfter;
    double massDiscrepancyAfter;

    int candidateVoxelCount;
    int carvedVoxelCountThisStep;
    int totalCarvedVoxelCount;

    int initialInsideVoxelCount;
    int initialCarveableVoxelCount;
    int remainingInsideVoxelCount;

    double carvedGridRatio;
    double carvedCarveableRatio;

    CarvingStopReason stopReason;

    CarvingStepDiagnostics();
};


const char* carvingStopReasonName(
    CarvingStopReason reason
);


void printCarvingStepDiagnostics(
    std::ostream& output,
    const CarvingStepDiagnostics& diagnostics
);


void writeCarvingDiagnosticsCsvHeader(
    std::ostream& output
);


void writeCarvingDiagnosticsCsvRow(
    std::ostream& output,
    const CarvingStepDiagnostics& diagnostics
);
