#pragma once
#include "spin2014/SourceCageGuardedSlqp.h"
#include "spin2014/SourceCageAdaptiveRebind.h"

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
struct SourceCageAdaptiveSolverBudget {
    // Cumulative conservative reservations, never refunded. NOT object-graph
    // bytes or RSS: rebuilt immutable layers retain separate upstream caps.
    std::size_t maximumKnownBufferBytes=std::size_t(8)*1024*1024*1024;
    std::size_t maximumRetainedSourcePairs=16;
};
struct SourceCageAdaptiveSolverSettings {
    std::size_t maximumPasses=4;
    SourceCageGuardedSlqpSettings solve;
    SourceCageAdaptiveRebindSettings rebind;
    SourceCageAdaptiveSolverBudget budget;
};
enum class SourceCageAdaptiveSolverOutcome {
    RelaxedAdaptiveComplete,MaximumPasses,SolveNotConverged,SolveGeometryRejected,
    SolveUnknown,RebindGeometryRejected,RebindUnknown,IdentityMismatch,ResourceLimit
};
struct SourceCageAdaptiveSolverResources {
    std::size_t cumulativeKnownBufferBytes=0,guardedReservedBytes=0,rebindReservedBytes=0;
    // Admitted construction slots versus successfully retained immutable pairs.
    // At most maximumPasses+1 and maximumRetainedSourcePairs slots are admitted.
    std::size_t admittedSourcePairs=0,retainedSourcePairs=0;
};
struct SourceCageAdaptiveSolverPass {
    std::size_t passNumber=0,startingStateValues=0,startingBetaValues=0,startingHandleValues=0;
    bool hasAngle=false,freshHistory=false; // Derived from actual trace, never a grant.
    // A pre-admission resource failure can leave source null and initialState
    // empty. Dimensions/hashes still describe the proposed pass; no solve ran.
    std::shared_ptr<const SourceCagePaperFunctional> source;
    Eigen::VectorXd initialState;
    std::shared_ptr<const SourceCageGuardedSlqpAudit> solve;
    std::shared_ptr<const SourceCageAdaptiveRebindAudit> rebind;
    std::string oldPartitionSha256,newPartitionSha256,oldStateSha256,solvedStateSha256,newStateSha256;
    std::string oldFunctionalSha256,newFunctionalSha256,guardedScientificSha256,rebindScientificSha256;
    std::size_t cumulativeKnownBufferBytes=0;
};
struct SourceCageAdaptiveSolverFacts {
    SourceCageAdaptiveSolverSettings settings;
    SourceCageAdaptiveSolverOutcome outcome=SourceCageAdaptiveSolverOutcome::SolveUnknown;
    SourceCageAdaptiveSolverResources resources;
    std::vector<SourceCageAdaptiveSolverPass> passes;
    std::size_t remainingFractionalLeaves=0,refinableFractionalLeaves=0;
    bool fractionalCountsAvailable=false;
    std::string initialFunctionalSha256,initialStateSha256,configSha256,scientificSha256;
    std::array<char,256> terminalReason{};
    std::string terminalReasonSha256;
};
class SourceCageAdaptiveSolverAudit;
/** An owned immutable full-state endpoint whose own fresh guarded solve has
 * numerically converged and freshly passed represented geometry and material.
 * This does not establish feasibility, binary material, spin stability, a seal
 * or original-paper completion. Copied facts cannot construct an endpoint.
 */
class SourceCageAdaptiveSolvedEndpoint {
public:
    const SourceCagePaperFunctional& functional() const {return *source_;}
    const Eigen::VectorXd& state() const {return solve_->facts().source.numerical.finalState;}
    const SourceCageGuardedSlqpAudit& solve() const {return *solve_;}
private:
    SourceCageAdaptiveSolvedEndpoint(std::shared_ptr<const SourceCagePaperFunctional> source,
        std::shared_ptr<const SourceCageGuardedSlqpAudit> solve):source_(std::move(source)),solve_(std::move(solve)) {}
    const std::shared_ptr<const SourceCagePaperFunctional> source_;
    const std::shared_ptr<const SourceCageGuardedSlqpAudit> solve_;
    friend SourceCageAdaptiveSolverAudit solveSourceCageAdaptive(const SourceCagePaperFunctional&,
        const Eigen::VectorXd&,const SourceCageAdaptiveSolverSettings&);
};
/** Shared immutable ownership keeps references valid after move while an owner
 * remains; moved-from accessors throw. Pending is a changed, successfully rebound
 * checkpoint not yet successfully solved, and stays separate on limits/failure.
 * RelaxedAdaptiveComplete permits maximum-depth fractional leaves. It is never
 * binary validation, a feasible witness, geometry seal, physics or paper grant.
 */
class SourceCageAdaptiveSolverAudit {
public:
    const SourceCageAdaptiveSolverFacts& facts() const;
    std::shared_ptr<const SourceCageAdaptiveSolvedEndpoint> lastSolved() const;
    std::shared_ptr<const SourceCageAdaptiveCheckpoint> pendingCheckpoint() const;
private:
    struct Data;
    explicit SourceCageAdaptiveSolverAudit(std::shared_ptr<const Data> data):data_(std::move(data)) {}
    std::shared_ptr<const Data> data_;
    friend SourceCageAdaptiveSolverAudit solveSourceCageAdaptive(const SourceCagePaperFunctional&,
        const Eigen::VectorXd&,const SourceCageAdaptiveSolverSettings&);
};
SourceCageAdaptiveSolverAudit solveSourceCageAdaptive(const SourceCagePaperFunctional&,
    const Eigen::VectorXd&,const SourceCageAdaptiveSolverSettings& = {});
}
#endif
