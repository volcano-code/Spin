#pragma once
#include "spin2014/SourceCageSlqp.h"
#include "spin2014/SourceCageGeometry.h"
#include <stdexcept>
#include <utility>

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
struct SourceCageGuardedSlqpBudget {
    std::size_t maximumGuardCalls=1024;
    std::size_t maximumRetainedStageRecords=65536;
    std::size_t maximumGuardedBytes=std::size_t(2)*1024*1024*1024;
};
struct SourceCageGuardedSlqpSettings {
    SourceCageSlqpSettings source;
    SourceCageGeometrySettings geometry;
    SourceCageGuardedSlqpBudget budget;
};
struct SourceCageGuardedSlqpResources {
    SourceCageSlqpResources adapter;
    std::size_t stagesPerGuard=0,stageCapacity=0,attemptCapacity=0;
    std::size_t retainedEvidenceBytes=0,checkerScratchBytes=0,requiredBytes=0;
};
enum class SourceCageGuardedFinalStatus {Unavailable,ChecksPassed,RejectedGeometry,Unknown,ResourceLimit};
struct SourceCageSlqpGuardAttempt {
    // Index into source.records and source.callbackStates, including final check.
    std::size_t callbackRecord=0,firstStage=0,stageCount=0;
    bool checkerInvoked=false,checkCompleted=false,sourceValuesInvoked=false;
    SourceCageGeometryOutcome outcome=SourceCageGeometryOutcome::Unknown;
    std::array<char,65> inputSha256{},sourceOutputSha256{},scientificSha256{};
    std::size_t plannedCheckerBytes=0,pairsVisited=0,exactPairRelations=0;
};
struct SourceCageGuardedSlqpFacts {
    // This legacy subresult retains NOT_EVALUATED. Geometry authorization lives
    // only in this privately constructed audit's finalStatus and check history.
    SourceCageSlqpResult source;
    SourceCageGuardedSlqpSettings settings;
    SourceCageGuardedSlqpResources resources;
    std::vector<SourceCageSlqpGuardAttempt> attempts;
    std::vector<SourceCageGeometryRecord> stages;
    std::size_t checkerInvocationCount=0,geometryValueEvaluationCount=0,finalRecheckCount=0;
    // False if a checker threw before publishing its actual stage tape. Never
    // describe the sum of known counts as a total in that case.
    bool geometryValueEvaluationCountComplete=true;
    SourceCageGuardedFinalStatus finalStatus=SourceCageGuardedFinalStatus::Unavailable;
    SourceCageSlqpFeasibility feasibility=SourceCageSlqpFeasibility::Unknown;
    std::string configSha256,scientificSha256;
};
class SourceCageGuardedSlqpAudit {
public:
    SourceCageGuardedSlqpAudit(SourceCageGuardedSlqpAudit&& other)
        :facts_(std::move(other.facts_)),available_(other.available_) {
        other.available_=false;
        // Also invalidate references obtained before the move, not only future
        // accessor calls. Numerical data transferred away is no longer present.
        other.facts_.finalStatus=SourceCageGuardedFinalStatus::Unavailable;
        other.facts_.source.finalEvaluationAvailable=false;
        other.facts_.source.numericalResultAvailable=false;
    }
    SourceCageGuardedSlqpAudit& operator=(SourceCageGuardedSlqpAudit&&)=delete;
    SourceCageGuardedSlqpAudit(const SourceCageGuardedSlqpAudit&)=delete;
    SourceCageGuardedSlqpAudit& operator=(const SourceCageGuardedSlqpAudit&)=delete;
    const SourceCageGuardedSlqpFacts& facts() const {
        if(!available_)throw std::logic_error("moved-from guarded audit has no evidence");
        return facts_;
    }
private:
    SourceCageGuardedSlqpAudit()=default;
    SourceCageGuardedSlqpFacts facts_;
    bool available_=true;
    friend SourceCageGuardedSlqpAudit solveSourceCageGuardedSlqp(
        const SourceCagePaperFunctional&,const Eigen::VectorXd&,const SourceCageGuardedSlqpSettings&);
};
SourceCageGuardedSlqpResources planSourceCageGuardedSlqpResources(
    const SourceCagePaperFunctional&,const SourceCageGuardedSlqpSettings& = {});
// The actual immutable source and actual checker are the only admission path.
// Numerical convergence, represented geometry and final material recomputation
// remain separate facts. None grants solver feasibility, a seal, or physics.
SourceCageGuardedSlqpAudit solveSourceCageGuardedSlqp(
    const SourceCagePaperFunctional&,const Eigen::VectorXd&,const SourceCageGuardedSlqpSettings& = {});
}
#endif
