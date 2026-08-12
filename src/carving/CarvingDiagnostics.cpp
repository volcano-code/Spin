#include "carving/CarvingDiagnostics.h"

#include <iomanip>
#include <ostream>


CarvingStepDiagnostics::CarvingStepDiagnostics()
    :
    iteration(0),
    centerOfMassBefore(Eigen::Vector3d::Zero()),
    centerOfMassAfter(Eigen::Vector3d::Zero()),
    recomputedCenterOfMassAfter(Eigen::Vector3d::Zero()),
    firstMomentBefore(Eigen::Vector3d::Zero()),
    firstMomentAfter(Eigen::Vector3d::Zero()),
    recomputedFirstMomentAfter(Eigen::Vector3d::Zero()),
    projectedDistanceBefore(0.0),
    projectedDistanceAfter(0.0),
    projectedEnergyBefore(0.0),
    projectedEnergyAfter(0.0),
    centerOfMassDiscrepancyAfter(0.0),
    firstMomentDiscrepancyAfter(0.0),
    initialMass(0.0),
    remainingMassBefore(0.0),
    remainingMassAfter(0.0),
    recomputedRemainingMassAfter(0.0),
    massDiscrepancyAfter(0.0),
    candidateVoxelCount(0),
    evaluatedPrefixCount(0),
    bestPrefixLength(0),
    carvedVoxelCountThisStep(0),
    totalCarvedVoxelCount(0),
    initialInsideVoxelCount(0),
    initialCarveableVoxelCount(0),
    remainingInsideVoxelCount(0),
    carvedGridRatio(0.0),
    carvedCarveableRatio(0.0),
    bestPrefixEnergy(0.0),
    materialComponentCountAfter(0),
    carvedVoidComponentCountAfter(0),
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

        case CarvingStopReason::NoCandidatesAvailable:
            return "no_candidates";

        case CarvingStopReason::TooFewVoxelsRemoved:
            return "too_few_voxels_removed";

        case CarvingStopReason::MaximumIterationsReached:
            return "maximum_iterations";

        case CarvingStopReason::MaximumEmptyRatioReached:
            return "maximum_empty_ratio";
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
        << " prefixes=" << diagnostics.evaluatedPrefixCount
        << " best_prefix=" << diagnostics.bestPrefixLength
        << " carved_step=" << diagnostics.carvedVoxelCountThisStep
        << " carved_total=" << diagnostics.totalCarvedVoxelCount
        << " distance="
        << diagnostics.projectedDistanceBefore
        << "->"
        << diagnostics.projectedDistanceAfter
        << " energy="
        << diagnostics.projectedEnergyBefore
        << "->"
        << diagnostics.projectedEnergyAfter
        << " best_prefix_energy="
        << diagnostics.bestPrefixEnergy
        << " material_components="
        << diagnostics.materialComponentCountAfter
        << " void_components="
        << diagnostics.carvedVoidComponentCountAfter
        << " mass="
        << diagnostics.remainingMassBefore
        << "->"
        << diagnostics.remainingMassAfter
        << " mass_discrepancy="
        << diagnostics.massDiscrepancyAfter
        << " first_moment_discrepancy="
        << diagnostics.firstMomentDiscrepancyAfter
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
        << "recomputed_com_after_x,recomputed_com_after_y,"
        << "recomputed_com_after_z,"
        << "first_moment_before_x,first_moment_before_y,"
        << "first_moment_before_z,"
        << "first_moment_after_x,first_moment_after_y,"
        << "first_moment_after_z,"
        << "recomputed_first_moment_after_x,"
        << "recomputed_first_moment_after_y,"
        << "recomputed_first_moment_after_z,"
        << "projected_distance_before,projected_distance_after,"
        << "projected_energy_before,projected_energy_after,"
        << "com_discrepancy_after,first_moment_discrepancy_after,"
        << "initial_mass,remaining_mass_before,remaining_mass_after,"
        << "recomputed_remaining_mass_after,mass_discrepancy_after,"
        << "candidate_voxel_count,evaluated_prefix_count,"
        << "best_prefix_length,best_prefix_energy,"
        << "material_component_count_after,"
        << "carved_void_component_count_after,"
        << "carved_voxel_count_step,"
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
        << diagnostics.recomputedCenterOfMassAfter(0) << ','
        << diagnostics.recomputedCenterOfMassAfter(1) << ','
        << diagnostics.recomputedCenterOfMassAfter(2) << ','
        << diagnostics.firstMomentBefore(0) << ','
        << diagnostics.firstMomentBefore(1) << ','
        << diagnostics.firstMomentBefore(2) << ','
        << diagnostics.firstMomentAfter(0) << ','
        << diagnostics.firstMomentAfter(1) << ','
        << diagnostics.firstMomentAfter(2) << ','
        << diagnostics.recomputedFirstMomentAfter(0) << ','
        << diagnostics.recomputedFirstMomentAfter(1) << ','
        << diagnostics.recomputedFirstMomentAfter(2) << ','
        << diagnostics.projectedDistanceBefore << ','
        << diagnostics.projectedDistanceAfter << ','
        << diagnostics.projectedEnergyBefore << ','
        << diagnostics.projectedEnergyAfter << ','
        << diagnostics.centerOfMassDiscrepancyAfter << ','
        << diagnostics.firstMomentDiscrepancyAfter << ','
        << diagnostics.initialMass << ','
        << diagnostics.remainingMassBefore << ','
        << diagnostics.remainingMassAfter << ','
        << diagnostics.recomputedRemainingMassAfter << ','
        << diagnostics.massDiscrepancyAfter << ','
        << diagnostics.candidateVoxelCount << ','
        << diagnostics.evaluatedPrefixCount << ','
        << diagnostics.bestPrefixLength << ','
        << diagnostics.bestPrefixEnergy << ','
        << diagnostics.materialComponentCountAfter << ','
        << diagnostics.carvedVoidComponentCountAfter << ','
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
