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
 * Two projected-error values are intentionally recorded:
 *
 * - legacyProjectedError*: the value used by the current/original solver;
 * - geometricProjectedError*: a direct ground-plane projection distance.
 *
 * The second value is diagnostic only in Phase 1.6 and does not yet affect
 * stopping behavior.
 */
struct CarvingStepDiagnostics
{
    int iteration;

    Eigen::Vector3d centerOfMassBefore;
    Eigen::Vector3d centerOfMassAfter;
    Eigen::Vector3d approximateVoxelCenterOfMassAfter;

    double legacyProjectedErrorBefore;
    double legacyProjectedErrorAfter;

    double geometricProjectedErrorBefore;
    double geometricProjectedErrorAfter;

    double centerOfMassDiscrepancyAfter;

    double exactInitialMass;
    double approximateRemainingVoxelMass;

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
