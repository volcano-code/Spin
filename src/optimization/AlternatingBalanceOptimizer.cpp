#include "optimization/AlternatingBalanceOptimizer.h"

#include "carving/CarvingSolver.h"
#include "deformation/DeformedVoxelMassField.h"
#include "deformation/DeformationObjective.h"
#include "manufacturing/ManufacturingGeometry.h"
#include "physics/BalanceObjective.h"
#include "physics/MassProperties.h"
#include "voxel/VoxelTopology.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>


namespace optimization
{
namespace
{


struct MaterialSnapshot
{
    VoxelGrid voxelFill;
    deformation::DeformedGeometry geometry;
    MassProperties massProperties;
    Eigen::Vector3d centerOfMass = Eigen::Vector3d::Zero();
};


deformation::HandleSet makeHandles(
    const deformation::DeformationSystem& system,
    const Eigen::VectorXd& parameters
);


std::uint64_t alphaHash(const VoxelGrid& grid);


void validateConfig(
    const AlternatingBalanceOptimizerConfig& config
)
{
    const deformation::PaperLineSearchConfig& line = config.lineSearch;

    if (
        config.maximumCarvingBlocks <= 0 ||
        config.maximumDeformationIterationsPerBlock <= 0 ||
        config.minimumCarvingRemovalsToContinue < 0 ||
        !std::isfinite(config.gradientTolerance) ||
        config.gradientTolerance < 0.0 ||
        !std::isfinite(config.shapeScaleLambda) ||
        config.shapeScaleLambda < 0.0 ||
        !std::isfinite(config.density) ||
        config.density <= 0.0 ||
        !std::isfinite(config.minimumPositiveMass) ||
        config.minimumPositiveMass < 0.0 ||
        !std::isfinite(config.minimumPositiveCellMass) ||
        config.minimumPositiveCellMass < 0.0 ||
        !std::isfinite(line.initialStep) ||
        !std::isfinite(line.reductionRatio) ||
        !std::isfinite(line.minimumStep) ||
        !std::isfinite(line.requiredRelativeDecrease) ||
        !std::isfinite(line.muDecrement) ||
        line.initialStep <= 0.0 ||
        line.minimumStep <= 0.0 ||
        line.minimumStep > line.initialStep ||
        line.reductionRatio <= 0.0 ||
        line.reductionRatio >= 1.0 ||
        line.requiredRelativeDecrease < 0.0 ||
        line.requiredRelativeDecrease >= 1.0 ||
        line.muDecrement <= 0.0 ||
        line.muDecrement > 1.0 ||
        (!config.enableCarving && !config.enableDeformation) ||
        (
            config.releasedBetaOuterLoop &&
            (
                !config.enableCarving ||
                !config.enableDeformation ||
                !config.freezeTargetsAfterInitialRefresh ||
                !line.releasedBetaSemantics
            )
        )
    )
    {
        throw std::invalid_argument(
            "alternating optimizer configuration is invalid"
        );
    }
}


bool allStandingConstraintsFeasible(
    const stability::BalanceTargetSnapshot& targets,
    bool& hasStandingConstraints
)
{
    hasStandingConstraints = false;
    bool allFeasible = true;
    for (const stability::BalanceConstraintSnapshot& constraint :
         targets.constraints)
    {
        if (constraint.type != stability::BalanceConstraintType::Standing)
        {
            continue;
        }
        hasStandingConstraints = true;
        allFeasible =
            allFeasible &&
            constraint.standingState.has_value() &&
            *constraint.standingState !=
                stability::StandingState::OutsideSupport;
    }
    return hasStandingConstraints && allFeasible;
}


void appendEvent(
    AlternatingBalanceOptimizationResult& result,
    AlternatingEventKind kind,
    int block,
    int iteration,
    int geometryEpoch,
    bool geometryChanged,
    const VoxelGrid& voxelFill,
    const deformation::DeformedGeometry& geometry,
    const Eigen::VectorXd& parameters,
    const MassProperties& massProperties,
    const Eigen::Vector3d& centerOfMass,
    const stability::BalanceTargetSnapshot& targets
)
{
    AlternatingEventRecord record;
    record.ordinal = static_cast<int>(result.events.size());
    record.kind = kind;
    record.block = block;
    record.iteration = iteration;
    record.geometryEpoch = geometryEpoch;
    record.geometryChanged = geometryChanged;
    record.allStandingFeasible = allStandingConstraintsFeasible(
        targets,
        record.hasStandingConstraints
    );
    record.allTargetsReached = targets.allReached;
    record.targetGeneration = targets.generation;
    record.targetSnapshotHash =
        stability::balanceTargetSnapshotContentHash(targets);
    record.alphaHash = alphaHash(voxelFill);
    record.geometryHash = deformation::deformedGeometryContentHash(geometry);
    record.parameters = parameters;
    record.projectedTargetEnergy = aggregateProjectedCenterOfMassEnergy(
        centerOfMass,
        targets.objectives
    );
    record.massProperties = massProperties;
    record.centerOfMass = centerOfMass;

    if (
        result.firstStandingFeasibleEvent < 0 &&
        record.allStandingFeasible
    )
    {
        result.firstStandingFeasibleEvent = record.ordinal;
    }
    if (result.firstSafeEvent < 0 && record.allTargetsReached)
    {
        result.firstSafeEvent = record.ordinal;
    }
    result.events.push_back(std::move(record));
}


void finalizeMetrics(
    const deformation::DeformationSystem& system,
    AlternatingBalanceOptimizationResult& result,
    double effectiveShapeScaleLambda,
    double density,
    double minimumPositiveMass
)
{
    const deformation::DeformationObjective objective(
        system,
        result.voxelFill,
        result.targets.objectives,
        effectiveShapeScaleLambda,
        density,
        minimumPositiveMass
    );
    const deformation::DeformationEnergyEvaluation evaluation =
        objective.evaluate(result.parameters, result.mu);
    result.finalProjectedTargetEnergy = evaluation.centerOfMassEnergy;
    result.finalUnscaledShapeEnergy = evaluation.unscaledShapeEnergy;
    result.finalShapeEnergy = evaluation.shapeEnergy;
    result.finalTotalEnergy = evaluation.totalEnergy;

    result.finalMassMoments = computeSurfaceMassMoments(
        result.geometry.outerVertices,
        result.geometry.outerFaces,
        density
    );
    const deformation::DeformedVoxelMassField finalCellMoments(
        system.domain(),
        result.geometry.gridCorners,
        density,
        0.0
    );
    const Eigen::Vector3i finalDimensions = result.voxelFill.dimensions();
    for (int x = 0; x < finalDimensions(0); ++x)
    {
        for (int y = 0; y < finalDimensions(1); ++y)
        {
            for (int z = 0; z < finalDimensions(2); ++z)
            {
                const int linear = system.domain().cellLinearIndex(x, y, z);
                if (
                    finalCellMoments.contains(linear) &&
                    result.voxelFill(x, y, z).isInner &&
                    !result.voxelFill(x, y, z).isFilled
                )
                {
                    result.finalMassMoments -=
                        finalCellMoments.massMoments(linear);
                }
            }
        }
    }
    const double momentTolerance =
        1.0e-9 * std::max(1.0, std::abs(result.massProperties.mass));
    if (
        std::abs(
            result.finalMassMoments.mass - result.massProperties.mass
        ) > momentTolerance ||
        (
            result.finalMassMoments.firstMoment -
            result.massProperties.firstMoment
        ).norm() > momentTolerance
    )
    {
        throw std::logic_error(
            "occupancy-bound mass moments disagree with final material mass"
        );
    }
    result.finalInertiaTensorAboutCenter =
        result.finalMassMoments.inertiaAboutCenter(minimumPositiveMass);
    result.finalInertiaTensorAboutCenter = 0.5 * (
        result.finalInertiaTensorAboutCenter +
        result.finalInertiaTensorAboutCenter.transpose()
    );

    const MassProperties finalOuter = computeSurfaceMassProperties(
        result.geometry.outerVertices,
        result.geometry.outerFaces,
        1.0
    );
    const MassProperties finalInner = computeSurfaceMassProperties(
        result.geometry.innerVertices,
        result.geometry.innerFaces,
        1.0
    );
    result.finalOuterVolume = finalOuter.mass;
    const double volumeTolerance =
        4096.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::abs(result.finalOuterVolume));
    if (result.geometry.innerFaces.rows() == 0)
    {
        result.finalCavityVolume = 0.0;
    }
    else if (finalInner.mass < -volumeTolerance)
    {
        result.finalCavityVolume = -finalInner.mass;
    }
    else
    {
        // A non-empty cavity with non-negative signed volume violates the
        // inward-boundary convention.  Preserve that invalidity for the
        // diagnostics/experiment layer instead of silently converting it to
        // a zero-volume cavity.
        result.finalCavityVolume =
            std::numeric_limits<double>::quiet_NaN();
    }
    result.finalCavityVolumeFraction =
        result.finalOuterVolume > 0.0 &&
        std::isfinite(result.finalCavityVolume) &&
        result.finalCavityVolume <= result.finalOuterVolume + volumeTolerance
            ? result.finalCavityVolume / result.finalOuterVolume
            : std::numeric_limits<double>::quiet_NaN();
    result.finalMaterialMassChangeFraction =
        result.initialMassProperties.mass > 0.0
            ? (result.initialMassProperties.mass -
               result.massProperties.mass) /
                result.initialMassProperties.mass
            : 0.0;

    const Eigen::MatrixXd displacement =
        result.geometry.outerVertices - system.outerRestVertices();
    const Eigen::Vector3d extent =
        system.outerRestVertices().colwise().maxCoeff().transpose() -
        system.outerRestVertices().colwise().minCoeff().transpose();
    const double normalization = extent.maxCoeff();
    if (!(normalization > 0.0) || !std::isfinite(normalization))
    {
        throw std::logic_error(
            "outer rest mesh has no finite positive normalization extent"
        );
    }

    double squaredDisplacement = 0.0;
    double maximumDisplacement = 0.0;
    for (Eigen::Index row = 0; row < displacement.rows(); ++row)
    {
        const double magnitude = displacement.row(row).norm();
        squaredDisplacement += magnitude * magnitude;
        maximumDisplacement = std::max(maximumDisplacement, magnitude);
    }
    result.normalizedRmsVertexDisplacement = std::sqrt(
        squaredDisplacement /
        static_cast<double>(std::max<Eigen::Index>(1, displacement.rows()))
    ) / normalization;
    result.normalizedMaximumVertexDisplacement =
        maximumDisplacement / normalization;

    const deformation::HandleSet finalHandles =
        makeHandles(system, result.parameters);
    for (const deformation::Handle& handle : finalHandles.handles())
    {
        if (handle.locked)
        {
            continue;
        }
        result.maximumUserHandleTranslation = std::max(
            result.maximumUserHandleTranslation,
            handle.translation.norm()
        );
        result.maximumUserHandleScaleDeviation = std::max(
            result.maximumUserHandleScaleDeviation,
            std::abs(handle.uniformScale - 1.0)
        );
    }
}


void validateInitialParameters(
    const deformation::DeformationSystem& system,
    const Eigen::VectorXd& parameters,
    double initialMu
)
{
    if (
        parameters.size() !=
            system.bindingHandles().optimizableParameterCount() ||
        !parameters.allFinite() ||
        !std::isfinite(initialMu) ||
        initialMu < 0.0 ||
        initialMu > 1.0
    )
    {
        throw std::invalid_argument(
            "alternating optimizer initial state is invalid"
        );
    }

    deformation::HandleSet handles = system.bindingHandles();
    handles.setOptimizationParameters(parameters);
}


int alphaHammingDistance(
    const VoxelGrid& first,
    const VoxelGrid& second
)
{
    if (first.dimensions() != second.dimensions())
    {
        throw std::invalid_argument(
            "alpha comparison requires identical grid dimensions"
        );
    }

    int distance = 0;
    const Eigen::Vector3i dimensions = first.dimensions();
    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                distance +=
                    first(x, y, z).isFilled !=
                    second(x, y, z).isFilled ? 1 : 0;
            }
        }
    }
    return distance;
}


std::uint64_t alphaHash(const VoxelGrid& grid)
{
    std::uint64_t hash = UINT64_C(1469598103934665603);
    std::uint64_t index = 0;
    const Eigen::Vector3i dimensions = grid.dimensions();

    for (int x = 0; x < dimensions(0); ++x)
    {
        for (int y = 0; y < dimensions(1); ++y)
        {
            for (int z = 0; z < dimensions(2); ++z)
            {
                hash ^= index++;
                hash *= UINT64_C(1099511628211);
                hash ^= grid(x, y, z).isFilled ? UINT64_C(1) : UINT64_C(0);
                hash *= UINT64_C(1099511628211);
            }
        }
    }
    return hash;
}


Eigen::Vector3d requireCenterOfMass(
    const MassProperties& properties,
    double minimumPositiveMass,
    const char* message
)
{
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    if (
        !properties.isFinite() ||
        !(properties.mass > minimumPositiveMass) ||
        !properties.tryCenterOfMass(center, minimumPositiveMass)
    )
    {
        throw deformation::InfeasibleDeformationState(message);
    }
    return center;
}


deformation::HandleSet makeHandles(
    const deformation::DeformationSystem& system,
    const Eigen::VectorXd& parameters
)
{
    deformation::HandleSet handles = system.bindingHandles();
    try
    {
        handles.setOptimizationParameters(parameters);
    }
    catch (const std::invalid_argument& error)
    {
        throw deformation::InfeasibleDeformationState(error.what());
    }
    return handles;
}


MaterialSnapshot evaluateFixedAlphaState(
    const deformation::DeformationSystem& system,
    const VoxelGrid& fixedVoxelFill,
    const Eigen::VectorXd& parameters,
    double density,
    double minimumPositiveMass,
    double minimumPositiveCellMass
)
{
    MaterialSnapshot snapshot;
    snapshot.voxelFill = fixedVoxelFill;
    const deformation::HandleSet handles = makeHandles(system, parameters);
    snapshot.geometry = system.evaluate(handles, snapshot.voxelFill);

    const manufacturing::BoundaryAudit outerAudit =
        manufacturing::auditTriangleBoundary(
            snapshot.geometry.outerVertices,
            snapshot.geometry.outerFaces,
            manufacturing::BoundaryOrientation::OutwardPositive,
            false,
            0.0
        );
    if (!outerAudit.strictBoundaryValid())
    {
        throw deformation::InfeasibleDeformationState(
            "deformation produced an invalid outer boundary"
        );
    }

    try
    {
        static_cast<void>(deformation::DeformedVoxelMassField(
            system.domain(),
            snapshot.geometry.gridCorners,
            density,
            minimumPositiveCellMass
        ));
    }
    catch (const deformation::InfeasibleDeformedVoxel& error)
    {
        throw deformation::InfeasibleDeformationState(error.what());
    }

    snapshot.massProperties =
        computeSurfaceMassProperties(
            snapshot.geometry.outerVertices,
            snapshot.geometry.outerFaces,
            density
        ) +
        computeSurfaceMassProperties(
            snapshot.geometry.innerVertices,
            snapshot.geometry.innerFaces,
            density
        );
    snapshot.centerOfMass = requireCenterOfMass(
        snapshot.massProperties,
        minimumPositiveMass,
        "deformed material has non-positive or non-finite mass"
    );
    return snapshot;
}


MaterialSnapshot evaluateSolidState(
    const deformation::DeformationSystem& system,
    const VoxelGrid& restGrid,
    const Eigen::VectorXd& parameters,
    double density,
    double minimumPositiveMass,
    double minimumPositiveCellMass
)
{
    VoxelGrid solid = restGrid;
    solid.refillAll();
    return evaluateFixedAlphaState(
        system,
        solid,
        parameters,
        density,
        minimumPositiveMass,
        minimumPositiveCellMass
    );
}


struct CarvingBlockOutcome
{
    MaterialSnapshot state;
    AlternatingCarvingRecord record;
    int carvingIterations = 0;
    int carvedVoxelCount = 0;
    CarvingStopReason stopReason = CarvingStopReason::None;
    bool stable = false;
    bool maximumIterationsReached = false;
};


struct ReleasedBetaBalanceOutcome
{
    MaterialSnapshot state;
    int carvingIterations = 0;
    int carvedVoxelCount = 0;
    CarvingStopReason stopReason = CarvingStopReason::None;
};


CarvingBlockOutcome runCarvingBlock(
    const deformation::DeformationSystem& system,
    const VoxelGrid& restGrid,
    const VoxelGrid& alphaBeforeReset,
    const Eigen::VectorXd& parameters,
    stability::BalanceTargetSchedule& schedule,
    stability::BalanceTargetSnapshot& targets,
    const AlternatingBalanceOptimizerConfig& config,
    int block,
    int& geometryEpoch,
    std::vector<AlternatingCarvingIterationRecord>& iterationTrace,
    AlternatingBalanceOptimizationResult& result
)
{
    CarvingBlockOutcome outcome;
    outcome.record.block = block;
    outcome.record.geometryEpochBeforeReset = geometryEpoch;
    outcome.record.alphaHashBeforeReset = alphaHash(alphaBeforeReset);

    VoxelGrid resetProbe = alphaBeforeReset;
    outcome.record.refilledVoxelCount = resetProbe.refillAll();
    outcome.state = evaluateSolidState(
        system,
        restGrid,
        parameters,
        config.density,
        config.minimumPositiveMass,
        config.minimumPositiveCellMass
    );
    if (outcome.record.refilledVoxelCount > 0)
    {
        ++geometryEpoch;
    }
    outcome.record.geometryEpochAfterReset = geometryEpoch;
    outcome.record.alphaHashAfterReset = alphaHash(outcome.state.voxelFill);

    outcome.record.resetTargetGenerationBefore = targets.generation;
    targets = schedule.refresh(outcome.state.centerOfMass);
    outcome.record.resetTargetGenerationAfter = targets.generation;
    outcome.record.resetObjectiveContextChanged = targets.contextChanged;
    outcome.record.resetAllTargetsReached = targets.allReached;
    outcome.record.targetGenerationBefore = targets.generation;
    outcome.record.resetTargetEnergy =
        aggregateProjectedCenterOfMassEnergy(
            outcome.state.centerOfMass,
            targets.objectives
        );

    if (outcome.record.refilledVoxelCount > 0)
    {
        appendEvent(
            result,
            AlternatingEventKind::SolidRefill,
            block,
            -1,
            geometryEpoch,
            true,
            outcome.state.voxelFill,
            outcome.state.geometry,
            parameters,
            outcome.state.massProperties,
            outcome.state.centerOfMass,
            targets
        );
    }

    if (targets.allReached)
    {
        outcome.stable = true;
        outcome.record.geometryEpoch = geometryEpoch;
        outcome.record.alphaHammingDistance = alphaHammingDistance(
            alphaBeforeReset,
            outcome.state.voxelFill
        );
        outcome.record.targetGenerationAfter = targets.generation;
        outcome.record.allTargetsReached = true;
        outcome.record.finalTargetEnergy = outcome.record.resetTargetEnergy;
        outcome.record.alphaHash = alphaHash(outcome.state.voxelFill);
        outcome.record.geometryHash =
            deformation::deformedGeometryContentHash(outcome.state.geometry);
        outcome.record.massProperties = outcome.state.massProperties;
        outcome.record.centerOfMass = outcome.state.centerOfMass;
        return outcome;
    }

    try
    {
        const deformation::DeformedVoxelMassField field(
            system.domain(),
            outcome.state.geometry.gridCorners,
            config.density,
            config.minimumPositiveCellMass
        );
        CarvingSolver solver(
            outcome.state.voxelFill,
            outcome.state.massProperties,
            field,
            targets.objectives
        );
        if (config.releasedBetaOuterLoop)
        {
            solver.setUseDeformedCornerAverageCandidatePositions(true);
        }
        const deformation::HandleSet handles = makeHandles(
            system,
            parameters
        );
        solver.setMinimumRemovalsToContinue(
            config.minimumCarvingRemovalsToContinue
        );

        while (!solver.isOptimized())
        {
            const VoxelGrid alphaBeforeStep = outcome.state.voxelFill;
            const std::uint64_t alphaBeforeHash = alphaHash(alphaBeforeStep);
            const std::size_t generationBefore = targets.generation;
            const double energyBefore = solver.getCoMEnergy();
            const int carvedBefore = solver.carvedVoxelCount();
            const int epochBefore = geometryEpoch;

            Eigen::MatrixXd temporaryPlaneVertices;
            Eigen::MatrixXi temporaryPlaneFaces;
            Eigen::Vector3d temporaryCenterOfMass;
            solver.carveStep(
                temporaryPlaneVertices,
                temporaryPlaneFaces,
                temporaryCenterOfMass
            );

            const int carvedThisStep =
                solver.carvedVoxelCount() - carvedBefore;
            if (carvedThisStep > 0)
            {
                ++geometryEpoch;
            }

            /*
             * A committed event is a complete state snapshot.  Carving owns
             * alpha while the block's deformation parameters stay frozen, so
             * rebuild the extracted inner boundary before hashing the event.
             * Otherwise alpha/mass/CoM would describe the post-carve state
             * while geometryHash still described the block-start solid.
             */
            outcome.state.geometry = system.evaluate(
                handles,
                outcome.state.voxelFill
            );

            targets = schedule.refresh(solver.currentCenterOfMass());

            AlternatingCarvingIterationRecord iteration;
            iteration.block = block;
            iteration.iteration = solver.iterationCount() - 1;
            iteration.geometryEpochBefore = epochBefore;
            iteration.geometryEpochAfter = geometryEpoch;
            iteration.carvedVoxelCount = carvedThisStep;
            iteration.totalCarvedVoxelCount = solver.carvedVoxelCount();
            iteration.stopReason = solver.stopReason();
            iteration.frozenTargetEnergyBefore = energyBefore;
            iteration.frozenTargetEnergyAfter =
                solver.diagnosticsHistory().back().projectedEnergyAfter;
            iteration.refreshedTargetEnergyAfter =
                aggregateProjectedCenterOfMassEnergy(
                    solver.currentCenterOfMass(),
                    targets.objectives
                );
            iteration.targetGenerationBefore = generationBefore;
            iteration.targetGenerationAfter = targets.generation;
            iteration.objectiveContextChanged = targets.contextChanged;
            iteration.allTargetsReached = targets.allReached;
            iteration.alphaHashBefore = alphaBeforeHash;
            iteration.alphaHashAfter = alphaHash(outcome.state.voxelFill);
            iteration.massProperties = solver.currentMassProperties();
            iteration.centerOfMass = solver.currentCenterOfMass();

            const int observedHamming = alphaHammingDistance(
                alphaBeforeStep,
                outcome.state.voxelFill
            );
            const double energyTolerance =
                1024.0 * std::numeric_limits<double>::epsilon() *
                std::max({
                    1.0,
                    std::abs(iteration.frozenTargetEnergyBefore),
                    std::abs(iteration.frozenTargetEnergyAfter)
                });
            if (
                observedHamming != carvedThisStep ||
                iteration.frozenTargetEnergyAfter >
                    iteration.frozenTargetEnergyBefore + energyTolerance
            )
            {
                throw std::logic_error(
                    "committed carving iteration violated its frozen context"
                );
            }
            iterationTrace.push_back(iteration);
            appendEvent(
                result,
                AlternatingEventKind::CarvingIteration,
                block,
                iteration.iteration,
                geometryEpoch,
                carvedThisStep > 0,
                outcome.state.voxelFill,
                outcome.state.geometry,
                parameters,
                iteration.massProperties,
                iteration.centerOfMass,
                targets
            );

            if (targets.allReached)
            {
                outcome.stable = true;
                break;
            }

            /*
             * The released beta's balanceByPlaneCarving() performs one
             * candidate sort and one best-prefix commit.  It is called once
             * when plane-carving mode is selected and once after every
             * deformation trial.  Do not turn that call into the paper
             * profile's repeated carving block.
             */
            if (config.releasedBetaOuterLoop)
            {
                break;
            }

            if (solver.isOptimized())
            {
                break;
            }

            solver.updateBalanceObjectives(targets.objectives);
        }

        outcome.carvingIterations = solver.iterationCount();
        outcome.carvedVoxelCount = solver.carvedVoxelCount();
        outcome.stopReason = solver.stopReason();
        outcome.maximumIterationsReached =
            !outcome.stable &&
            outcome.stopReason == CarvingStopReason::MaximumIterationsReached;
        outcome.state.geometry = system.evaluate(
            handles,
            outcome.state.voxelFill
        );
        outcome.state.massProperties = solver.currentMassProperties();
        outcome.state.centerOfMass = solver.currentCenterOfMass();

        const MassProperties surfaceOracle =
            computeSurfaceMassProperties(
                outcome.state.geometry.outerVertices,
                outcome.state.geometry.outerFaces,
                config.density
            ) +
            computeSurfaceMassProperties(
                outcome.state.geometry.innerVertices,
                outcome.state.geometry.innerFaces,
                config.density
            );
        const double scale = std::max({
            1.0,
            surfaceOracle.mass,
            outcome.state.massProperties.mass,
            surfaceOracle.firstMoment.norm(),
            outcome.state.massProperties.firstMoment.norm()
        });
        const double tolerance =
            4096.0 * std::numeric_limits<double>::epsilon() * scale;
        if (
            std::abs(
                surfaceOracle.mass - outcome.state.massProperties.mass
            ) > tolerance ||
            (
                surfaceOracle.firstMoment -
                outcome.state.massProperties.firstMoment
            ).norm() > tolerance
        )
        {
            throw std::logic_error(
                "deformed carving mass disagrees with the cavity surface"
            );
        }
        outcome.record.geometryEpoch = geometryEpoch;
        outcome.record.alphaHammingDistance = alphaHammingDistance(
            alphaBeforeReset,
            outcome.state.voxelFill
        );
        outcome.record.carvingIterations = outcome.carvingIterations;
        outcome.record.carvedVoxelCount = outcome.carvedVoxelCount;
        outcome.record.stopReason = outcome.stopReason;
        outcome.record.targetGenerationAfter = targets.generation;
        outcome.record.objectiveContextChanged =
            targets.generation != outcome.record.targetGenerationBefore;
        outcome.record.allTargetsReached = targets.allReached;
        outcome.record.finalTargetEnergy =
            aggregateProjectedCenterOfMassEnergy(
                outcome.state.centerOfMass,
                targets.objectives
            );
        outcome.record.alphaHash = alphaHash(outcome.state.voxelFill);
        outcome.record.geometryHash =
            deformation::deformedGeometryContentHash(outcome.state.geometry);
        outcome.record.massProperties = outcome.state.massProperties;
        outcome.record.centerOfMass = outcome.state.centerOfMass;
        return outcome;
    }
    catch (const deformation::InfeasibleDeformedVoxel& error)
    {
        throw deformation::InfeasibleDeformationState(error.what());
    }
}


ReleasedBetaBalanceOutcome runReleasedBetaBalance(
    const deformation::DeformationSystem& system,
    const VoxelGrid& restGrid,
    const Eigen::VectorXd& parameters,
    const BalanceObjectives& frozenObjectives,
    const AlternatingBalanceOptimizerConfig& config
)
{
    ReleasedBetaBalanceOutcome outcome;
    outcome.state = evaluateSolidState(
        system,
        restGrid,
        parameters,
        config.density,
        config.minimumPositiveMass,
        config.minimumPositiveCellMass
    );

    try
    {
        const deformation::DeformedVoxelMassField field(
            system.domain(),
            outcome.state.geometry.gridCorners,
            config.density,
            config.minimumPositiveCellMass
        );
        CarvingSolver solver(
            outcome.state.voxelFill,
            outcome.state.massProperties,
            field,
            frozenObjectives
        );
        solver.setUseDeformedCornerAverageCandidatePositions(true);
        solver.setMinimumRemovalsToContinue(
            config.minimumCarvingRemovalsToContinue
        );
        if (!solver.isOptimized())
        {
            Eigen::MatrixXd temporaryPlaneVertices;
            Eigen::MatrixXi temporaryPlaneFaces;
            Eigen::Vector3d temporaryCenterOfMass;
            solver.carveStep(
                temporaryPlaneVertices,
                temporaryPlaneFaces,
                temporaryCenterOfMass
            );
        }
        outcome.carvingIterations = solver.iterationCount();
        outcome.carvedVoxelCount = solver.carvedVoxelCount();
        outcome.stopReason = solver.stopReason();
        if (outcome.stopReason == CarvingStopReason::MaximumIterationsReached)
        {
            throw deformation::InfeasibleDeformationState(
                "released-beta balance hit the local carving iteration guard"
            );
        }

        const deformation::HandleSet handles = makeHandles(system, parameters);
        outcome.state.geometry = system.evaluate(
            handles,
            outcome.state.voxelFill
        );
        outcome.state.massProperties = solver.currentMassProperties();
        outcome.state.centerOfMass = solver.currentCenterOfMass();

        const MassProperties surfaceOracle =
            computeSurfaceMassProperties(
                outcome.state.geometry.outerVertices,
                outcome.state.geometry.outerFaces,
                config.density
            ) +
            computeSurfaceMassProperties(
                outcome.state.geometry.innerVertices,
                outcome.state.geometry.innerFaces,
                config.density
            );
        const double scale = std::max({
            1.0,
            surfaceOracle.mass,
            outcome.state.massProperties.mass,
            surfaceOracle.firstMoment.norm(),
            outcome.state.massProperties.firstMoment.norm()
        });
        const double tolerance =
            4096.0 * std::numeric_limits<double>::epsilon() * scale;
        if (
            std::abs(surfaceOracle.mass - outcome.state.massProperties.mass) >
                tolerance ||
            (surfaceOracle.firstMoment -
             outcome.state.massProperties.firstMoment).norm() > tolerance
        )
        {
            throw std::logic_error(
                "released-beta recarving mass disagrees with cavity surface"
            );
        }
        return outcome;
    }
    catch (const deformation::InfeasibleDeformedVoxel& error)
    {
        throw deformation::InfeasibleDeformationState(error.what());
    }
}


void copyStateToResult(
    const MaterialSnapshot& state,
    AlternatingBalanceOptimizationResult& result
)
{
    result.voxelFill = state.voxelFill;
    result.geometry = state.geometry;
    result.massProperties = state.massProperties;
    result.centerOfMass = state.centerOfMass;
    result.geometryHash = deformation::deformedGeometryContentHash(
        state.geometry
    );
}


void completeResult(
    const deformation::DeformationSystem& system,
    const MaterialSnapshot& state,
    const AlternatingBalanceOptimizerConfig& config,
    AlternatingBalanceOptimizationResult& result
)
{
    copyStateToResult(state, result);
    result.mu = result.lineSearchState.mu;
    finalizeMetrics(
        system,
        result,
        result.effectiveShapeScaleLambda,
        config.density,
        config.minimumPositiveMass
    );
}


} // namespace


AlternatingBalanceOptimizer::AlternatingBalanceOptimizer(
    const deformation::DeformationSystem& deformationSystem,
    const VoxelGrid& restVoxelGrid,
    const stability::BalanceTargetSchedule& targetSchedule,
    const AlternatingBalanceOptimizerConfig& config
)
    :
    deformationSystem_(deformationSystem),
    restVoxelGrid_(restVoxelGrid),
    targetSpecifications_(targetSchedule.specifications()),
    config_(config)
{
    validateConfig(config_);
    deformationSystem_.validateVoxelFill(restVoxelGrid_);
    restVoxelGrid_.refillAll();
}


AlternatingBalanceOptimizationResult
AlternatingBalanceOptimizer::optimize(
    const Eigen::VectorXd& initialParameters,
    double initialMu
) const
{
    validateInitialParameters(
        deformationSystem_,
        initialParameters,
        initialMu
    );

    AlternatingBalanceOptimizationResult result;
    result.carvingEnabled = config_.enableCarving;
    result.deformationEnabled = config_.enableDeformation;
    result.laplacianEnabled = config_.enableLaplacian;
    result.releasedBetaOuterLoopEnabled = config_.releasedBetaOuterLoop;
    result.parameters = initialParameters;
    result.mu = initialMu;
    result.effectiveShapeScaleLambda = config_.enableLaplacian
        ? config_.shapeScaleLambda
        : 0.0;
    result.gridDimensions = restVoxelGrid_.dimensions();
    result.voxelSize = restVoxelGrid_.voxelSize();
    result.activeDeformationNodeCount =
        static_cast<int>(deformationSystem_.domain().nodeCount());
    std::vector<int> fixedNodes;
    for (const deformation::Handle& handle :
         deformationSystem_.bindingHandles().handles())
    {
        if (handle.locked)
        {
            fixedNodes.insert(
                fixedNodes.end(),
                handle.constrainedNodes.begin(),
                handle.constrainedNodes.end()
            );
        }
        else
        {
            ++result.userHandleCount;
        }
    }
    std::sort(fixedNodes.begin(), fixedNodes.end());
    fixedNodes.erase(
        std::unique(fixedNodes.begin(), fixedNodes.end()),
        fixedNodes.end()
    );
    result.fixedDeformationNodeCount = static_cast<int>(fixedNodes.size());
    result.lineSearchState.step = config_.lineSearch.initialStep;
    result.lineSearchState.mu = initialMu;
    result.usesHypotheticalFlatteningSupport = std::any_of(
        targetSpecifications_.begin(),
        targetSpecifications_.end(),
        [](const stability::BalanceTargetSpecification& specification)
        {
            const auto* standing =
                std::get_if<stability::StandingTargetSpecification>(
                    &specification
                );
            return standing != nullptr &&
                standing->support.contactThreshold() > 0.0;
        }
    );

    stability::BalanceTargetSchedule schedule(
        targetSpecifications_,
        config_.freezeTargetsAfterInitialRefresh
    );
    MaterialSnapshot current = evaluateSolidState(
        deformationSystem_,
        restVoxelGrid_,
        result.parameters,
        config_.density,
        config_.minimumPositiveMass,
        config_.minimumPositiveCellMass
    );

    const VoxelTopologyReport initialTopology = analyzeVoxelTopology(
        current.voxelFill
    );
    result.initialMaterialVoxelCount = initialTopology.materialVoxelCount;
    result.initialCarveableVoxelCount = 0;
    for (int x = 0; x < result.gridDimensions(0); ++x)
    {
        for (int y = 0; y < result.gridDimensions(1); ++y)
        {
            for (int z = 0; z < result.gridDimensions(2); ++z)
            {
                result.initialCarveableVoxelCount +=
                    current.voxelFill(x, y, z).isInner ? 1 : 0;
            }
        }
    }
    result.initialMassProperties = current.massProperties;
    result.initialOuterVolume = computeSurfaceMassProperties(
        current.geometry.outerVertices,
        current.geometry.outerFaces,
        1.0
    ).mass;
    result.initialCenterOfMass = current.centerOfMass;
    result.initialTargets = schedule.refresh(current.centerOfMass);
    result.targets = result.initialTargets;
    appendEvent(
        result,
        AlternatingEventKind::Initial,
        -1,
        -1,
        result.geometryEpoch,
        false,
        current.voxelFill,
        current.geometry,
        result.parameters,
        current.massProperties,
        current.centerOfMass,
        result.targets
    );

    if (result.targets.allReached && !config_.enableCarving)
    {
        result.status = config_.enableCarving
            ? AlternatingOptimizationStatus::StableAfterCarving
            : AlternatingOptimizationStatus::StableWithoutCarving;
        completeResult(deformationSystem_, current, config_, result);
        return result;
    }

    for (int block = 0;
         block < config_.maximumCarvingBlocks;
         ++block)
    {
        CarvingBlockOutcome carved;
        if (config_.enableCarving)
        {
            carved = runCarvingBlock(
                deformationSystem_,
                restVoxelGrid_,
                current.voxelFill,
                result.parameters,
                schedule,
                result.targets,
                config_,
                block,
                result.geometryEpoch,
                result.carvingIterations,
                result
            );
        }
        else
        {
            carved.state = current;
            result.targets = schedule.refresh(current.centerOfMass);
        }

        current = std::move(carved.state);
        if (config_.enableCarving)
        {
            result.carvingBlocks.push_back(carved.record);
        }

        if (carved.stable)
        {
            result.status = AlternatingOptimizationStatus::StableAfterCarving;
            completeResult(deformationSystem_, current, config_, result);
            return result;
        }

        if (carved.maximumIterationsReached)
        {
            result.status =
                AlternatingOptimizationStatus::MaximumCarvingIterations;
            completeResult(deformationSystem_, current, config_, result);
            return result;
        }

        if (!config_.enableDeformation)
        {
            result.status =
                AlternatingOptimizationStatus::DeformationDisabledUnstable;
            completeResult(deformationSystem_, current, config_, result);
            return result;
        }

        if (result.parameters.size() == 0)
        {
            result.status =
                AlternatingOptimizationStatus::NoDeformationDegreesOfFreedom;
            completeResult(deformationSystem_, current, config_, result);
            return result;
        }

        if (config_.releasedBetaOuterLoop)
        {
            for (int iteration = 0;
                 iteration < config_.maximumDeformationIterationsPerBlock;
                 ++iteration)
            {
                const VoxelGrid alphaBefore = current.voxelFill;
                const Eigen::VectorXd parametersBefore = result.parameters;
                const int geometryEpochBefore = result.geometryEpoch;
                const std::size_t targetGenerationBefore =
                    result.targets.generation;
                const std::uint64_t alphaBeforeHash = alphaHash(alphaBefore);
                const std::uint64_t geometryBeforeHash =
                    deformation::deformedGeometryContentHash(current.geometry);

                deformation::DeformationObjective objective(
                    deformationSystem_,
                    current.voxelFill,
                    result.targets.objectives,
                    result.effectiveShapeScaleLambda,
                    config_.density,
                    config_.minimumPositiveMass
                );
                const deformation::DeformationEnergyEvaluation evaluation =
                    objective.evaluate(
                        result.parameters,
                        result.lineSearchState.mu
                    );
                const double gradientNorm = evaluation.gradient.norm();
                if (!std::isfinite(gradientNorm))
                {
                    throw std::overflow_error(
                        "released-beta deformation gradient norm overflowed"
                    );
                }
                if (gradientNorm <= config_.gradientTolerance)
                {
                    result.status = AlternatingOptimizationStatus::
                        DeformationConvergedUnstable;
                    completeResult(deformationSystem_, current, config_, result);
                    return result;
                }

                std::optional<ReleasedBetaBalanceOutcome> trialBalance;
                bool trialBalanceAttempted = false;
                int energyEvaluationCalls = 0;
                const deformation::PaperLineSearchResult trial =
                    deformation::paperGradientDescentTrial(
                        result.parameters,
                        evaluation.gradient,
                        result.lineSearchState,
                        config_.lineSearch,
                        [&](const Eigen::VectorXd& parameters, double mu)
                        {
                            ++energyEvaluationCalls;
                            if (energyEvaluationCalls == 1)
                            {
                                if (
                                    !parameters.isApprox(parametersBefore, 0.0) ||
                                    mu != result.lineSearchState.mu
                                )
                                {
                                    throw std::logic_error(
                                        "released-beta energy callback did not "
                                        "begin from the committed state"
                                    );
                                }
                                return evaluation.totalEnergy;
                            }
                            if (energyEvaluationCalls != 2)
                            {
                                throw std::logic_error(
                                    "released-beta line search requested an "
                                    "unexpected objective reweight"
                                );
                            }

                            trialBalanceAttempted = true;
                            trialBalance = runReleasedBetaBalance(
                                deformationSystem_,
                                restVoxelGrid_,
                                parameters,
                                result.targets.objectives,
                                config_
                            );
                            const deformation::DeformationObjective trialObjective(
                                deformationSystem_,
                                trialBalance->state.voxelFill,
                                result.targets.objectives,
                                result.effectiveShapeScaleLambda,
                                config_.density,
                                config_.minimumPositiveMass
                            );
                            return trialObjective.evaluate(parameters, mu).totalEnergy;
                        }
                    );

                AlternatingDeformationRecord record;
                record.block = block;
                record.iteration = iteration;
                record.geometryEpochBefore = geometryEpochBefore;
                record.decision = trial.decision;
                record.accepted = trial.accepted;
                record.rolledBack = trial.rolledBack;
                record.energyBefore = trial.energyBefore;
                record.trialEnergy = trial.trialEnergy;
                record.gradientNorm = gradientNorm;
                record.attemptedStep = trial.attemptedStep;
                record.muBefore = trial.muBefore;
                record.objectiveReweighted = trial.objectiveReweighted;
                record.stateAfter = trial.nextState;
                record.targetGenerationBefore = targetGenerationBefore;
                record.alphaHashBefore = alphaBeforeHash;
                record.parametersBefore = parametersBefore;
                record.geometryHashBefore = geometryBeforeHash;
                record.releasedBetaOuterLoop = true;
                record.trialBalanceCalls = trialBalanceAttempted ? 1 : 0;
                if (trialBalance.has_value())
                {
                    record.trialCarvingIterations =
                        trialBalance->carvingIterations;
                    record.trialCarvedVoxelCount =
                        trialBalance->carvedVoxelCount;
                    record.trialAlphaHash =
                        alphaHash(trialBalance->state.voxelFill);
                }

                if (trial.accepted)
                {
                    if (!trialBalance.has_value())
                    {
                        throw std::logic_error(
                            "accepted released-beta trial has no balanced state"
                        );
                    }
                    result.parameters = trial.parameters;
                    current = std::move(trialBalance->state);
                    ++result.geometryEpoch;
                }
                else
                {
                    ReleasedBetaBalanceOutcome rollback =
                        runReleasedBetaBalance(
                            deformationSystem_,
                            restVoxelGrid_,
                            parametersBefore,
                            result.targets.objectives,
                            config_
                        );
                    record.rollbackBalanceCalls = 1;
                    record.rollbackCarvingIterations =
                        rollback.carvingIterations;
                    record.rollbackCarvedVoxelCount = rollback.carvedVoxelCount;
                    record.rollbackAlphaHash =
                        alphaHash(rollback.state.voxelFill);
                    current = std::move(rollback.state);

                    if (
                        !result.parameters.isApprox(parametersBefore, 0.0) ||
                        record.rollbackAlphaHash != alphaBeforeHash ||
                        deformation::deformedGeometryContentHash(
                            current.geometry
                        ) != geometryBeforeHash
                    )
                    {
                        throw std::logic_error(
                            "released-beta rollback did not restore and "
                            "rebalance the committed geometry"
                        );
                    }
                }

                record.refillCalls =
                    record.trialBalanceCalls + record.rollbackBalanceCalls;
                record.carvingCalls = record.refillCalls;
                result.lineSearchState = trial.nextState;
                result.mu = result.lineSearchState.mu;
                result.targets = schedule.refresh(current.centerOfMass);

                if (trial.accepted)
                {
                    appendEvent(
                        result,
                        AlternatingEventKind::AcceptedDeformation,
                        block,
                        iteration,
                        result.geometryEpoch,
                        true,
                        current.voxelFill,
                        current.geometry,
                        result.parameters,
                        current.massProperties,
                        current.centerOfMass,
                        result.targets
                    );
                }

                const deformation::DeformationObjective retainedObjective(
                    deformationSystem_,
                    current.voxelFill,
                    result.targets.objectives,
                    result.effectiveShapeScaleLambda,
                    config_.density,
                    config_.minimumPositiveMass
                );
                record.energyAfter = retainedObjective.evaluate(
                    result.parameters,
                    trial.muBefore
                ).totalEnergy;
                record.geometryEpochAfter = result.geometryEpoch;
                record.targetGenerationAfter = result.targets.generation;
                record.objectiveContextChanged = result.targets.contextChanged;
                record.allTargetsReached = result.targets.allReached;
                record.refreshedTargetEnergyAfter =
                    aggregateProjectedCenterOfMassEnergy(
                        current.centerOfMass,
                        result.targets.objectives
                    );
                record.alphaHashAfter = alphaHash(current.voxelFill);
                record.alphaHammingDistance = alphaHammingDistance(
                    alphaBefore,
                    current.voxelFill
                );
                record.parametersAfter = result.parameters;
                record.geometryHashAfter =
                    deformation::deformedGeometryContentHash(current.geometry);
                record.massPropertiesAfter = current.massProperties;
                record.centerOfMassAfter = current.centerOfMass;
                result.deformationIterations.push_back(record);

                if (result.targets.allReached)
                {
                    result.status =
                        AlternatingOptimizationStatus::StableAfterDeformation;
                    completeResult(deformationSystem_, current, config_, result);
                    return result;
                }
                if (
                    trial.decision ==
                        deformation::PaperTrialDecision::MuExhausted
                )
                {
                    result.status =
                        AlternatingOptimizationStatus::MuExhaustedUnstable;
                    completeResult(deformationSystem_, current, config_, result);
                    return result;
                }
            }

            result.status =
                AlternatingOptimizationStatus::MaximumDeformationIterations;
            completeResult(deformationSystem_, current, config_, result);
            return result;
        }

        const std::uint64_t blockAlphaHash = alphaHash(current.voxelFill);
        const VoxelGrid fixedAlpha = current.voxelFill;
        const std::size_t deformationTraceBegin =
            result.deformationIterations.size();
        bool gradientConverged = false;

        for (int iteration = 0;
             iteration < config_.maximumDeformationIterationsPerBlock;
             ++iteration)
        {
            deformation::DeformationObjective objective(
                deformationSystem_,
                fixedAlpha,
                result.targets.objectives,
                result.effectiveShapeScaleLambda,
                config_.density,
                config_.minimumPositiveMass
            );
            const deformation::DeformationEnergyEvaluation evaluation =
                objective.evaluate(
                    result.parameters,
                    result.lineSearchState.mu
                );
            const double gradientNorm = evaluation.gradient.norm();
            if (!std::isfinite(gradientNorm))
            {
                throw std::overflow_error(
                    "alternating deformation gradient norm overflowed"
                );
            }

            if (gradientNorm <= config_.gradientTolerance)
            {
                gradientConverged = true;
                break;
            }

            const Eigen::VectorXd parametersBefore = result.parameters;
            const int geometryEpochBefore = result.geometryEpoch;
            const std::size_t targetGenerationBefore =
                result.targets.generation;
            const std::uint64_t alphaBeforeHash =
                alphaHash(current.voxelFill);
            const std::uint64_t geometryBeforeHash =
                deformation::deformedGeometryContentHash(current.geometry);

            const deformation::PaperLineSearchResult trial =
                deformation::paperGradientDescentTrial(
                    result.parameters,
                    evaluation.gradient,
                    result.lineSearchState,
                    config_.lineSearch,
                    [&](const Eigen::VectorXd& parameters, double mu)
                    {
                        const deformation::DeformationEnergyEvaluation value =
                            objective.evaluate(parameters, mu);
                        if (!parameters.isApprox(result.parameters, 0.0))
                        {
                            static_cast<void>(evaluateFixedAlphaState(
                                deformationSystem_,
                                fixedAlpha,
                                parameters,
                                config_.density,
                                config_.minimumPositiveMass,
                                config_.minimumPositiveCellMass
                            ));
                        }
                        return value.totalEnergy;
                    }
                );

            AlternatingDeformationRecord record;
            record.block = block;
            record.iteration = iteration;
            record.geometryEpochBefore = geometryEpochBefore;
            record.decision = trial.decision;
            record.accepted = trial.accepted;
            record.rolledBack = trial.rolledBack;
            record.energyBefore = trial.energyBefore;
            record.trialEnergy = trial.trialEnergy;
            record.energyAfter = trial.energyAfter;
            record.gradientNorm = gradientNorm;
            record.attemptedStep = trial.attemptedStep;
            record.muBefore = trial.muBefore;
            record.objectiveReweighted = trial.objectiveReweighted;
            record.stateAfter = trial.nextState;
            record.targetGenerationBefore = targetGenerationBefore;
            record.alphaHashBefore = alphaBeforeHash;
            record.parametersBefore = parametersBefore;
            record.geometryHashBefore = geometryBeforeHash;

            if (trial.accepted)
            {
                MaterialSnapshot accepted = evaluateFixedAlphaState(
                    deformationSystem_,
                    fixedAlpha,
                    trial.parameters,
                    config_.density,
                    config_.minimumPositiveMass,
                    config_.minimumPositiveCellMass
                );
                result.parameters = trial.parameters;
                current = std::move(accepted);
                ++result.geometryEpoch;
                result.targets = schedule.refresh(current.centerOfMass);
                appendEvent(
                    result,
                    AlternatingEventKind::AcceptedDeformation,
                    block,
                    iteration,
                    result.geometryEpoch,
                    true,
                    current.voxelFill,
                    current.geometry,
                    result.parameters,
                    current.massProperties,
                    current.centerOfMass,
                    result.targets
                );
            }
            else
            {
                if (
                    !trial.parameters.isApprox(parametersBefore, 0.0) ||
                    !result.parameters.isApprox(parametersBefore, 0.0)
                )
                {
                    throw std::logic_error(
                        "rejected deformation trial changed live parameters"
                    );
                }
            }

            result.lineSearchState = trial.nextState;
            result.mu = result.lineSearchState.mu;
            record.geometryEpochAfter = result.geometryEpoch;
            record.targetGenerationAfter = result.targets.generation;
            record.objectiveContextChanged =
                trial.accepted && result.targets.contextChanged;
            record.allTargetsReached = result.targets.allReached;
            record.refreshedTargetEnergyAfter =
                aggregateProjectedCenterOfMassEnergy(
                    current.centerOfMass,
                    result.targets.objectives
                );
            record.alphaHashAfter = alphaHash(current.voxelFill);
            record.alphaHammingDistance = alphaHammingDistance(
                fixedAlpha,
                current.voxelFill
            );
            record.parametersAfter = result.parameters;
            record.geometryHashAfter =
                deformation::deformedGeometryContentHash(current.geometry);
            record.massPropertiesAfter = current.massProperties;
            record.centerOfMassAfter = current.centerOfMass;

            if (
                record.alphaHashBefore != blockAlphaHash ||
                record.alphaHashAfter != blockAlphaHash ||
                record.alphaHammingDistance != 0
            )
            {
                throw std::logic_error(
                    "fixed-alpha deformation block changed voxel fill"
                );
            }

            result.deformationIterations.push_back(record);

            if (result.targets.allReached)
            {
                result.status =
                    config_.enableCarving
                        ? AlternatingOptimizationStatus::StableAfterDeformation
                        : AlternatingOptimizationStatus::StableWithoutCarving;
                completeResult(deformationSystem_, current, config_, result);
                return result;
            }

            if (
                trial.decision ==
                    deformation::PaperTrialDecision::MuExhausted
            )
            {
                result.status =
                    AlternatingOptimizationStatus::MuExhaustedUnstable;
                completeResult(deformationSystem_, current, config_, result);
                return result;
            }
        }

        if (!gradientConverged)
        {
            result.status =
                AlternatingOptimizationStatus::MaximumDeformationIterations;
            completeResult(deformationSystem_, current, config_, result);
            return result;
        }

        if (!config_.enableCarving)
        {
            result.status =
                AlternatingOptimizationStatus::DeformationConvergedUnstable;
            completeResult(deformationSystem_, current, config_, result);
            return result;
        }

        const bool acceptedDeformationThisBlock = std::any_of(
            result.deformationIterations.begin() +
                static_cast<std::ptrdiff_t>(deformationTraceBegin),
            result.deformationIterations.end(),
            [](const AlternatingDeformationRecord& record)
            {
                return record.accepted;
            }
        );
        if (
            result.carvingBlocks.back().alphaHammingDistance == 0 &&
            !acceptedDeformationThisBlock
        )
        {
            result.status = AlternatingOptimizationStatus::AlternatingStalled;
            completeResult(deformationSystem_, current, config_, result);
            return result;
        }
    }

    result.status = AlternatingOptimizationStatus::MaximumCarvingBlocks;
    completeResult(deformationSystem_, current, config_, result);
    return result;
}


} // namespace optimization
