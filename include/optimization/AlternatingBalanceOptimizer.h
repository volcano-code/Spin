#pragma once

#include "carving/CarvingDiagnostics.h"
#include "deformation/DeformationOptimizer.h"
#include "deformation/DeformationSystem.h"
#include "physics/MassProperties.h"
#include "stability/BalanceTargetSchedule.h"
#include "voxel/VoxelGrid.h"

#include <Eigen/Core>

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>


namespace optimization
{


struct AlternatingBalanceOptimizerConfig
{
    deformation::PaperLineSearchConfig lineSearch;
    int maximumCarvingBlocks = 20;
    int maximumDeformationIterationsPerBlock = 200;
    double gradientTolerance = 1.0e-8;
    double shapeScaleLambda = 1.0;
    double density = 1.0;
    double minimumPositiveMass = 0.0;
    double minimumPositiveCellMass = 0.0;
    int minimumCarvingRemovalsToContinue = 1;

    /** R10 ablation switches; the paper profile keeps all three enabled. */
    bool enableCarving = true;
    bool enableDeformation = true;
    bool enableLaplacian = true;
};


enum class AlternatingOptimizationStatus
{
    StableAfterCarving = 0,
    StableAfterDeformation = 1,
    NoDeformationDegreesOfFreedom = 2,
    MuExhaustedUnstable = 3,
    AlternatingStalled = 4,
    MaximumCarvingBlocks = 5,
    MaximumCarvingIterations = 6,
    MaximumDeformationIterations = 7,

    /** R10 additions are appended so the pre-R10 enum ordinals stay stable. */
    StableWithoutCarving = 8,
    DeformationDisabledUnstable = 9,
    DeformationConvergedUnstable = 10
};


/** Ordered committed-state events used by R10 convergence accounting. */
enum class AlternatingEventKind
{
    Initial,
    SolidRefill,
    CarvingIteration,
    AcceptedDeformation
};


struct AlternatingEventRecord
{
    int ordinal = 0;
    AlternatingEventKind kind = AlternatingEventKind::Initial;
    int block = -1;
    int iteration = -1;
    int geometryEpoch = 0;
    bool geometryChanged = false;
    bool hasStandingConstraints = false;
    bool allStandingFeasible = false;
    bool allTargetsReached = false;
    double projectedTargetEnergy = 0.0;
    std::size_t targetGeneration = 0;
    std::uint64_t targetSnapshotHash = 0;
    std::uint64_t alphaHash = 0;
    std::uint64_t geometryHash = 0;
    Eigen::VectorXd parameters;
    MassProperties massProperties;
    Eigen::Vector3d centerOfMass = Eigen::Vector3d::Zero();
};


/** One committed Section 4 iteration with a frozen target context. */
struct AlternatingCarvingIterationRecord
{
    int block = 0;
    int iteration = 0;
    int geometryEpochBefore = 0;
    int geometryEpochAfter = 0;
    int carvedVoxelCount = 0;
    int totalCarvedVoxelCount = 0;
    CarvingStopReason stopReason = CarvingStopReason::None;
    double frozenTargetEnergyBefore = 0.0;
    double frozenTargetEnergyAfter = 0.0;
    double refreshedTargetEnergyAfter = 0.0;
    std::size_t targetGenerationBefore = 0;
    std::size_t targetGenerationAfter = 0;
    bool objectiveContextChanged = false;
    bool allTargetsReached = false;
    std::uint64_t alphaHashBefore = 0;
    std::uint64_t alphaHashAfter = 0;
    MassProperties massProperties;
    Eigen::Vector3d centerOfMass = Eigen::Vector3d::Zero();
};


/** One paper Section 3.2 carving block on a frozen deformed lattice. */
struct AlternatingCarvingRecord
{
    int block = 0;
    int geometryEpochBeforeReset = 0;
    int geometryEpochAfterReset = 0;
    int geometryEpoch = 0;
    int refilledVoxelCount = 0;
    int alphaHammingDistance = 0;
    int carvingIterations = 0;
    int carvedVoxelCount = 0;
    CarvingStopReason stopReason = CarvingStopReason::None;
    double resetTargetEnergy = 0.0;
    double finalTargetEnergy = 0.0;
    std::size_t resetTargetGenerationBefore = 0;
    std::size_t resetTargetGenerationAfter = 0;
    bool resetObjectiveContextChanged = false;
    bool resetAllTargetsReached = false;
    std::size_t targetGenerationBefore = 0;
    std::size_t targetGenerationAfter = 0;
    /** True when a post-reset committed carving refresh changed targets. */
    bool objectiveContextChanged = false;
    bool allTargetsReached = false;
    std::uint64_t alphaHashBeforeReset = 0;
    std::uint64_t alphaHashAfterReset = 0;
    std::uint64_t alphaHash = 0;
    std::uint64_t geometryHash = 0;
    MassProperties massProperties;
    Eigen::Vector3d centerOfMass = Eigen::Vector3d::Zero();
};


/** One fixed-alpha R6 line-search trial inside an R7 deformation block. */
struct AlternatingDeformationRecord
{
    int block = 0;
    int iteration = 0;
    int geometryEpochBefore = 0;
    int geometryEpochAfter = 0;
    deformation::PaperTrialDecision decision =
        deformation::PaperTrialDecision::StepReduced;
    bool accepted = false;
    bool rolledBack = true;
    double energyBefore = 0.0;
    double trialEnergy = 0.0;
    double energyAfter = 0.0;
    /** Energy of the committed/retained CoM after refreshing its targets. */
    double refreshedTargetEnergyAfter = 0.0;
    double gradientNorm = 0.0;
    double attemptedStep = 0.0;
    double muBefore = 0.0;
    bool objectiveReweighted = false;
    deformation::PaperLineSearchState stateAfter;
    std::size_t targetGenerationBefore = 0;
    std::size_t targetGenerationAfter = 0;
    bool objectiveContextChanged = false;
    bool allTargetsReached = false;

    /** Paper profile discriminator: both remain zero inside this block. */
    int refillCalls = 0;
    int carvingCalls = 0;
    int alphaHammingDistance = 0;
    std::uint64_t alphaHashBefore = 0;
    std::uint64_t alphaHashAfter = 0;
    Eigen::VectorXd parametersBefore;
    Eigen::VectorXd parametersAfter;

    /** Content hashes bind this trace record to its geometry snapshots. */
    std::uint64_t geometryHashBefore = 0;
    std::uint64_t geometryHashAfter = 0;

    MassProperties massPropertiesAfter;
    Eigen::Vector3d centerOfMassAfter = Eigen::Vector3d::Zero();
};


struct AlternatingBalanceOptimizationResult
{
    /** Immutable execution profile copied from the validated optimizer config. */
    bool carvingEnabled = true;
    bool deformationEnabled = true;
    bool laplacianEnabled = true;
    AlternatingOptimizationStatus status =
        AlternatingOptimizationStatus::MaximumCarvingBlocks;
    Eigen::VectorXd parameters;
    double mu = 0.0;
    deformation::PaperLineSearchState lineSearchState;
    VoxelGrid voxelFill;
    deformation::DeformedGeometry geometry;
    MassProperties massProperties;
    Eigen::Vector3d centerOfMass = Eigen::Vector3d::Zero();

    /** R10 initial/final quantities use local paper-profile semantics. */
    Eigen::Vector3i gridDimensions = Eigen::Vector3i::Zero();
    Eigen::Vector3d voxelSize = Eigen::Vector3d::Zero();
    int initialMaterialVoxelCount = 0;
    int initialCarveableVoxelCount = 0;
    int activeDeformationNodeCount = 0;
    int fixedDeformationNodeCount = 0;
    int userHandleCount = 0;
    MassProperties initialMassProperties;
    Eigen::Vector3d initialCenterOfMass = Eigen::Vector3d::Zero();
    stability::BalanceTargetSnapshot initialTargets{
        0,
        false,
        false,
        BalanceObjectives(),
        {}
    };
    double effectiveShapeScaleLambda = 0.0;
    double finalProjectedTargetEnergy = 0.0;
    double finalUnscaledShapeEnergy = 0.0;
    double finalShapeEnergy = 0.0;
    double finalTotalEnergy = 0.0;
    double initialOuterVolume = 0.0;
    double finalOuterVolume = 0.0;
    double finalCavityVolume = 0.0;
    double finalCavityVolumeFraction = 0.0;
    double finalMaterialMassChangeFraction = 0.0;
    double normalizedRmsVertexDisplacement = 0.0;
    double normalizedMaximumVertexDisplacement = 0.0;
    int firstStandingFeasibleEvent = -1;
    int firstSafeEvent = -1;
    double maximumUserHandleTranslation = 0.0;
    double maximumUserHandleScaleDeviation = 0.0;
    /** Content hash of the returned outer/inner/lattice geometry snapshot. */
    std::uint64_t geometryHash = 0;
    /** True when a standing target used a positive, unapplied flatten band. */
    bool usesHypotheticalFlatteningSupport = false;
    stability::BalanceTargetSnapshot targets{
        0,
        false,
        false,
        BalanceObjectives(),
        {}
    };
    int geometryEpoch = 0;
    std::vector<AlternatingCarvingRecord> carvingBlocks;
    std::vector<AlternatingCarvingIterationRecord> carvingIterations;
    std::vector<AlternatingDeformationRecord> deformationIterations;
    std::vector<AlternatingEventRecord> events;
};


/**
 * Paper-governed block-coordinate controller for Sections 3.2--5.3.
 *
 * Every carving block begins from an all-filled copy and uses the currently
 * deformed hexahedral mass field.  A standing target is refreshed after that
 * reset and after every committed carving iteration.  The following
 * deformation block freezes alpha until its fixed-alpha gradient descent
 * converges.  Fixed supports and suspension attachments are never rebuilt.
 */
class AlternatingBalanceOptimizer
{
public:
    AlternatingBalanceOptimizer(
        const deformation::DeformationSystem& deformationSystem,
        const VoxelGrid& restVoxelGrid,
        const stability::BalanceTargetSchedule& targetSchedule,
        const AlternatingBalanceOptimizerConfig& config =
            AlternatingBalanceOptimizerConfig()
    );

    AlternatingBalanceOptimizationResult optimize(
        const Eigen::VectorXd& initialParameters,
        double initialMu
    ) const;

private:
    const deformation::DeformationSystem& deformationSystem_;
    VoxelGrid restVoxelGrid_;
    stability::BalanceTargetSpecifications targetSpecifications_;
    AlternatingBalanceOptimizerConfig config_;
};


} // namespace optimization
