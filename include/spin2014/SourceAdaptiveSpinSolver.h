#pragma once
#include "spin2014/SourceBinaryMoments.h"
#include "spin2014/SlqpLbfgsSolver.h"
#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
enum class SourceAdaptiveSpinStatus {BinaryMomentValidated,BinaryValidationFailed,RelaxationUnknown,AdaptivePassLimit,ResourceLimitReached};
const char* sourceAdaptiveSpinStatusToken(SourceAdaptiveSpinStatus) noexcept;
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
    SlqpLbfgsSettings slqp=[] {SlqpLbfgsSettings s;s.algorithm=SlqpAlgorithm::AuditedKktRestorationV1;return s;}();
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
};
void validateSourceAdaptiveSpinRequest(const SourceAdaptiveSpinRequest&);
std::string computeSourceAdaptiveSpinRequestHash(const SourceAdaptiveSpinRequest&);
SourceAdaptiveSpinResult solveSourceAdaptiveSpin(const SourceAdaptiveSpinRequest&);
}
#endif
