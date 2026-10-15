#pragma once
#include "spin2014/SourceBinaryMoments.h"
#include "spin2014/SlqpLbfgsSolver.h"
#include "spin/SpinDesignRequest.h"
#include <optional>
#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
enum class SourceAdaptiveSpinStatus {BinaryMomentValidated,BinaryValidationFailed,InfeasibleCertified,RelaxationUnknown,AdaptivePassLimit,ResourceLimitReached,StageBoundaryPaused};
const char* sourceAdaptiveSpinStatusToken(SourceAdaptiveSpinStatus) noexcept;
/** A certificate scoped to the current finite source partition.  It proves
 * that an affine equality row cannot reach zero for any beta in [0,1]^n;
 * it deliberately says nothing about a refined or continuous source domain.
 */
struct SourceAffineConstraintCertificate {
    int constraintIndex=-1;
    double lower=0.,upper=0.,tolerance=0.;
    std::string scope="DISCRETE_SOURCE_PARTITION_AFFINE_CONSTRAINT_BOX_V1";
};
enum class SourceAdaptiveTermination { UnchangedPartitionKktV1, PaperPostTransitionV1 };
const char* sourceAdaptiveTerminationToken(SourceAdaptiveTermination) noexcept;
/** PAPER_GIVENS + AUDITED_KKT_RESTORATION_V1. All numerical defaults and
 * post-merge resolve policy are IMPLEMENTATION_DEFAULT, not author trajectory.
 * Single density, native z axis/origin, no automatic geometry rotation.
 */
struct SourceAdaptiveSpinRequest {
    std::shared_ptr<const SourceSurfaceOctree> source;
    Spin2014Objective objective=Spin2014Objective::StaticLowering;
    int initialLevel=2,maximumOptimizationLevel=5,maximumAdaptivePasses=16;
    std::size_t maximumLeafCount=1000000,maximumOutputCuboids=1000000;
    double initialFill=.5,initialPhi=0,fractionalEpsilon=.05;
    double laplacianWeight=1e-3,inertiaWeight=1,centerWeight=1,density=1,binaryConstraintTolerance=1e-6;
    // Optional source-surface dual material densities. When both are absent,
    // legacy `density` is used for both regions; specifying only one is invalid.
    std::optional<double> boundaryDensity,interiorDensity;
    double resolvedBoundaryDensity() const noexcept {return boundaryDensity.value_or(density);}
    double resolvedInteriorDensity() const noexcept {return interiorDensity.value_or(density);}
    /** Objective/constraints frame; source geometry and emitted cuboids remain
     * in model coordinates. */
    spin::SpinTargetFrame targetFrame;
    SlqpLbfgsSettings slqp=[] {SlqpLbfgsSettings s;s.algorithm=SlqpAlgorithm::AuditedKktRestorationV1;return s;}();
    // Preserve legacy callers. Paper Sec. 5.3 is an explicit opt-in, not a
    // claim that the refined, inherited state has been optimized again.
    SourceAdaptiveTermination termination=SourceAdaptiveTermination::UnchangedPartitionKktV1;
    // Exclusive directory, one complete iteration artifact per attempted solve.
    // This is a cumulative source-run budget, excluding resident pass metadata.
    std::string iterationStoreDirectory;
    std::uint64_t iterationStoreMaximumBytes=1073741824;
};
struct SourceAdaptivePass {
    int index=0;
    std::vector<SourceAdaptiveLeaf> before,after;
    std::string beforeHash,afterHash;
    Eigen::VectorXd inputBeta,outputBeta;
    double inputPhi=0;
    SlqpLbfgsResult solve;
    SourcePartitionTransition split,merge;
    std::vector<std::size_t> levelCounts,fractionalLevelCounts;
    std::size_t fractionalCount=0;
    std::vector<std::size_t> outputLevelCounts,outputFractionalLevelCounts;
    std::size_t outputFractionalCount=0;
    bool paperResolutionComplete=false;
};
struct SourceAdaptiveSpinResult {
    SourceAdaptiveSpinStatus status=SourceAdaptiveSpinStatus::RelaxationUnknown;
    std::string diagnostic,inputHash,traceHash,outputHash,finalPartitionHash;
    int geometryDepth=0,maximumOptimizationLevel=0;
    double finalPhi=0;
    Eigen::VectorXd relaxedFill,binaryFill,roundingDelta;
    std::vector<SourceAdaptiveLeaf> finalPartition;
    std::vector<SourceAdaptivePass> passes;
    std::optional<physics::RawMassMoments10> binaryMoments;
    std::optional<SourceBinaryCuboids> binaryCuboids;
    std::optional<SourceBinaryMomentCheck> binaryCheck;
    std::optional<SourceAffineConstraintCertificate> infeasibilityCertificate;
    SourceAdaptiveTermination termination=SourceAdaptiveTermination::UnchangedPartitionKktV1;
    // Last solve's KKT describes before/solve.beta, never the rounded object.
    bool terminalRelaxationKktApplies=false;
};
void validateSourceAdaptiveSpinRequest(const SourceAdaptiveSpinRequest&);
std::string computeSourceAdaptiveSpinRequestHash(const SourceAdaptiveSpinRequest&);
SourceAdaptiveSpinResult solveSourceAdaptiveSpin(const SourceAdaptiveSpinRequest&);
}
#endif
