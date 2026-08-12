#pragma once

#include <Eigen/Core>

#include <iosfwd>


/**
 * Stopping conditions exposed by the carving solver.
 *
 * The solver distinguishes an empty queue from a non-improving/too-small
 * best prefix.
 * The iteration terminal is an explicit local run guard, not paper convergence.
 */
enum class CarvingStopReason
{
    None,
    /** Retained for source compatibility; the R4 solver no longer emits it. */
    TargetToleranceReached,
    NoCandidatesAvailable,
    TooFewVoxelsRemoved,
    MaximumIterationsReached,
    /** Retained for source compatibility; the R4 solver no longer emits it. */
    MaximumEmptyRatioReached
};


/**
 * Diagnostic record for one call to CarvingSolver::carveStep().
 *
 * Diagnostics record both the incremental mass-property state and a full recompute
 * from the immutable outer surface minus every carved voxel. Their
 * discrepancies are expected to remain near floating-point roundoff. They also include
 * complete-prefix selection and topology-report fields.
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
    int evaluatedPrefixCount;
    int bestPrefixLength;
    int carvedVoxelCountThisStep;
    int totalCarvedVoxelCount;

    int initialInsideVoxelCount;
    int initialCarveableVoxelCount;
    int remainingInsideVoxelCount;

    double carvedGridRatio;
    double carvedCarveableRatio;

    double bestPrefixEnergy;

    int materialComponentCountAfter;
    int carvedVoidComponentCountAfter;

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
