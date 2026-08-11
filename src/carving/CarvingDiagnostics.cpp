#include "carving/CarvingDiagnostics.h"

#include <iomanip>
#include <ostream>


CarvingStepDiagnostics::CarvingStepDiagnostics()
    :
    iteration(0),
    centerOfMassBefore(Eigen::Vector3d::Zero()),
    centerOfMassAfter(Eigen::Vector3d::Zero()),
    approximateVoxelCenterOfMassAfter(Eigen::Vector3d::Zero()),
    legacyProjectedErrorBefore(0.0),
    legacyProjectedErrorAfter(0.0),
    geometricProjectedErrorBefore(0.0),
    geometricProjectedErrorAfter(0.0),
    centerOfMassDiscrepancyAfter(0.0),
    exactInitialMass(0.0),
    approximateRemainingVoxelMass(0.0),
    candidateVoxelCount(0),
    carvedVoxelCountThisStep(0),
    totalCarvedVoxelCount(0),
    initialInsideVoxelCount(0),
    initialCarveableVoxelCount(0),
    remainingInsideVoxelCount(0),
    carvedGridRatio(0.0),
    carvedCarveableRatio(0.0),
    stopReason(CarvingStopReason::None)
{
}


const char* carvingStopReasonName(
    CarvingStopReason reason
)
{
    switch (reason)
    {
        case CarvingStopReason::None:
            return "none";

        case CarvingStopReason::TargetToleranceReached:
            return "target_tolerance";

        case CarvingStopReason::MaximumEmptyRatioReached:
            return "maximum_empty_ratio";

        case CarvingStopReason::NoCandidatesAvailable:
            return "no_candidates";

        case CarvingStopReason::MaximumIterationsReached:
            return "maximum_iterations";
    }

    return "unknown";
}


void printCarvingStepDiagnostics(
    std::ostream& output,
    const CarvingStepDiagnostics& diagnostics
)
{
    output
        << std::setprecision(10)
        << "[carving] iteration=" << diagnostics.iteration
        << " candidates=" << diagnostics.candidateVoxelCount
        << " carved_step=" << diagnostics.carvedVoxelCountThisStep
        << " carved_total=" << diagnostics.totalCarvedVoxelCount
        << " legacy_error="
        << diagnostics.legacyProjectedErrorBefore
        << "->"
        << diagnostics.legacyProjectedErrorAfter
        << " geometric_error="
        << diagnostics.geometricProjectedErrorBefore
        << "->"
        << diagnostics.geometricProjectedErrorAfter
        << " com_discrepancy="
        << diagnostics.centerOfMassDiscrepancyAfter
        << " carved_grid_ratio="
        << diagnostics.carvedGridRatio
        << " carved_carveable_ratio="
        << diagnostics.carvedCarveableRatio
        << " stop="
        << carvingStopReasonName(diagnostics.stopReason)
        << '\n';
}


void writeCarvingDiagnosticsCsvHeader(
    std::ostream& output
)
{
    output
        << "iteration,"
        << "com_before_x,com_before_y,com_before_z,"
        << "com_after_x,com_after_y,com_after_z,"
        << "voxel_com_after_x,voxel_com_after_y,voxel_com_after_z,"
        << "legacy_error_before,legacy_error_after,"
        << "geometric_error_before,geometric_error_after,"
        << "com_discrepancy_after,"
        << "exact_initial_mass,approximate_remaining_voxel_mass,"
        << "candidate_voxel_count,carved_voxel_count_step,"
        << "total_carved_voxel_count,"
        << "initial_inside_voxel_count,initial_carveable_voxel_count,"
        << "remaining_inside_voxel_count,"
        << "carved_grid_ratio,carved_carveable_ratio,stop_reason\n";
}


void writeCarvingDiagnosticsCsvRow(
    std::ostream& output,
    const CarvingStepDiagnostics& diagnostics
)
{
    output << std::setprecision(17);

    output
        << diagnostics.iteration << ','
        << diagnostics.centerOfMassBefore(0) << ','
        << diagnostics.centerOfMassBefore(1) << ','
        << diagnostics.centerOfMassBefore(2) << ','
        << diagnostics.centerOfMassAfter(0) << ','
        << diagnostics.centerOfMassAfter(1) << ','
        << diagnostics.centerOfMassAfter(2) << ','
        << diagnostics.approximateVoxelCenterOfMassAfter(0) << ','
        << diagnostics.approximateVoxelCenterOfMassAfter(1) << ','
        << diagnostics.approximateVoxelCenterOfMassAfter(2) << ','
        << diagnostics.legacyProjectedErrorBefore << ','
        << diagnostics.legacyProjectedErrorAfter << ','
        << diagnostics.geometricProjectedErrorBefore << ','
        << diagnostics.geometricProjectedErrorAfter << ','
        << diagnostics.centerOfMassDiscrepancyAfter << ','
        << diagnostics.exactInitialMass << ','
        << diagnostics.approximateRemainingVoxelMass << ','
        << diagnostics.candidateVoxelCount << ','
        << diagnostics.carvedVoxelCountThisStep << ','
        << diagnostics.totalCarvedVoxelCount << ','
        << diagnostics.initialInsideVoxelCount << ','
        << diagnostics.initialCarveableVoxelCount << ','
        << diagnostics.remainingInsideVoxelCount << ','
        << diagnostics.carvedGridRatio << ','
        << diagnostics.carvedCarveableRatio << ','
        << carvingStopReasonName(diagnostics.stopReason)
        << '\n';
}
