#pragma once
#include "spin2014/SourceCageGeometry.h"

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
struct SourceCageAdaptiveRebindBudget {
    std::size_t maximumOwnedBytes=256*1024*1024;
};
struct SourceCageAdaptiveRebindSettings {
    double fractionalEpsilon=.05; // IMPLEMENTATION_DEFAULT, not a paper parameter.
    SourceCageGeometrySettings geometry;
    SourceCageAdaptiveRebindBudget budget;
};
enum class SourceCageAdaptiveRebindOutcome {Published,RejectedGeometry,ResourceLimit,Unknown};
enum class SourceCageAdaptiveRebindStage {
    Preflight,InputGeometry,InputMaterial,Split,Merge,Rebuild,OutputGeometry,OutputMaterial
};
enum class SourceCageAdaptiveRebindStageOutcome {NotAttempted,Attempted,Passed,RejectedGeometry,ResourceLimit,Unknown};
struct SourceCageAdaptiveRebindRecord {
    SourceCageAdaptiveRebindStageOutcome outcome=SourceCageAdaptiveRebindStageOutcome::NotAttempted;
    std::array<char,160> reason{};
    std::array<char,65> reasonSha256{};
};
struct SourceCageAdaptiveMaterialFacts {
    bool available=false;
    physics::RawMassMoments10 worldMoments,spinMoments;
    std::string outputSha256,worldOutputSha256; // Full derivative evaluation.
};
struct SourceCageAdaptiveRebindFacts {
    SourceCageAdaptiveRebindOutcome outcome=SourceCageAdaptiveRebindOutcome::Unknown;
    SourceCageAdaptiveRebindStage lastStage=SourceCageAdaptiveRebindStage::Preflight;
    std::array<SourceCageAdaptiveRebindRecord,8> stages{};
    SourceCageAdaptiveRebindSettings settings;
    // Bounds transaction evidence/state/evaluation payload, separately from
    // upstream domain/rebuild/checker budgets; never a native RSS guarantee.
    std::size_t plannedOwnedBytes=0;
    std::vector<SourcePartitionKey> oldLeaves,newLeaves;
    Eigen::VectorXd oldBeta,newBeta;
    SourcePartitionTransition split,merge; // removedBefore/After are REST moments.
    std::shared_ptr<const SourceCageGeometryAudit> inputGeometry,outputGeometry;
    SourceCageAdaptiveMaterialFacts beforeMaterial,afterMaterial;
    Eigen::Matrix<double,10,1> deformedMomentDifference=Eigen::Matrix<double,10,1>::Zero();
    bool deformedMomentDifferenceAvailable=false;
    std::string oldStateSha256,newStateSha256,oldFunctionalSha256,newFunctionalSha256;
    std::string oldPartitionSha256,newPartitionSha256,oldBindingSha256,newBindingSha256;
    std::string configSha256,scientificSha256;
};
class SourceCageAdaptiveRebindAudit;
/** Only an entirely successful actual transaction can construct this immutable
 * pair. No numerical history, KKT, solver feasibility or geometry conclusion is
 * transferred into a subsequent solveSourceCageGuardedSlqp call.
 */
class SourceCageAdaptiveCheckpoint {
public:
    const SourceCagePaperFunctional& functional() const {return functional_;}
    const Eigen::VectorXd& state() const {return state_;}
private:
    SourceCageAdaptiveCheckpoint(const SourceCagePaperFunctional& f,Eigen::VectorXd&& x)
        :functional_(f),state_(std::move(x)) {}
    const SourceCagePaperFunctional functional_;
    const Eigen::VectorXd state_;
    friend SourceCageAdaptiveRebindAudit rebindSourceCageAdaptive(
        const SourceCagePaperFunctional&,const Eigen::VectorXd&,const SourceCageAdaptiveRebindSettings&);
};
/** Shared const evidence and checkpoint storage does not move its scalar facts
 * away from previously obtained references. A moved-from accessor throws.
 * A failed audit always has a null checkpoint; copied facts cannot authorize one.
 */
class SourceCageAdaptiveRebindAudit {
public:
    const SourceCageAdaptiveRebindFacts& facts() const;
    std::shared_ptr<const SourceCageAdaptiveCheckpoint> checkpoint() const;
private:
    struct Data;
    explicit SourceCageAdaptiveRebindAudit(std::shared_ptr<const Data> data):data_(std::move(data)) {}
    std::shared_ptr<const Data> data_;
    friend SourceCageAdaptiveRebindAudit rebindSourceCageAdaptive(
        const SourceCagePaperFunctional&,const Eigen::VectorXd&,const SourceCageAdaptiveRebindSettings&);
};
// Always checks actual input and rebuilt output, including unchanged partitions.
// Split then merge run only on a local copy. All upstream configs/caps, density,
// frame, source and full BBW are retained. Invalid API inputs may throw; resource
// failures return nonpublishing evidence unless diagnostic allocation also fails.
// Use the 24 GiB / 60 minute external process guards.
SourceCageAdaptiveRebindAudit rebindSourceCageAdaptive(const SourceCagePaperFunctional&,
    const Eigen::VectorXd&,const SourceCageAdaptiveRebindSettings& = {});
}
#endif
