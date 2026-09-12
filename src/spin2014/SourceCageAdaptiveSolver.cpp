#include "spin2014/SourceCageAdaptiveSolver.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <cmath>
#include <cstring>

namespace spin2014 {
namespace {
using Outcome=SourceCageAdaptiveSolverOutcome;
[[noreturn]] void resource(const char* reason) {throw std::length_error(reason);}
std::size_t add(std::size_t a,std::size_t b) {if(a>SIZE_MAX-b)resource("aggregate addition overflow");return a+b;}
std::size_t mul(std::size_t a,std::size_t b) {if(b&&a>SIZE_MAX/b)resource("aggregate multiplication overflow");return a*b;}
class Hash {
    verification::Sha256 h_;
public:
    void integer(std::uint64_t n) {unsigned char b[8];for(int i=0;i<8;++i)b[7-i]=n>>(8*i);h_.update(b,8);}
    void string(const std::string& s) {integer(s.size());h_.update(reinterpret_cast<const unsigned char*>(s.data()),s.size());}
    void real(double d) {std::uint64_t n;std::memcpy(&n,&d,8);integer(n);}
    void values(const Eigen::VectorXd& x) {integer(x.size());for(Eigen::Index i=0;i<x.size();++i)real(x[i]);}
    std::string finish() {return h_.finish();}
};
// Same identity domains as the actual rebind interface, so adjacent pass links
// can be compared without treating configuration hash alone as full identity.
std::string stateHash(const Eigen::VectorXd& x) {Hash h;h.string("source-adaptive-full-state-v1");h.values(x);return h.finish();}
std::string functionalHash(const SourceCagePaperFunctional& f) {
    Hash h;h.string("source-adaptive-functional-v1");
    for(const auto* s:{&f.configSha256(),&f.upstreamSha256(),&f.source().configSha256(),
        &f.source().upstreamSha256(),&f.source().boundary().configSha256(),&f.source().boundary().outputSha256(),
        &f.source().boundary().binding().configSha256(),&f.source().boundary().binding().outputSha256(),
        &f.source().boundary().binding().domain().partitionHash()})h.string(*s);
    return h.finish();
}
void geometryHash(Hash& h,const SourceCageGeometrySettings& g) {
    h.real(g.minimumOuterDomainClearance);
    for(auto n:{g.maximumPoints,g.maximumTriangles,g.maximumPrimitives,g.maximumPairs,g.maximumRecords,
        g.maximumCheckerBytes,g.maximumNefComplexity,g.maximumConnectedShells,g.maximumExactStringBytes})h.integer(n);
}
std::string configHash(const SourceCageAdaptiveSolverSettings& s) {
    Hash h;h.string("actual-source-adaptive-solver-v1 cumulative-known-buffer-reservations-v1 no-feasibility-binary-physics-seal-paper-grant");
    h.integer(s.maximumPasses);h.integer(s.budget.maximumKnownBufferBytes);h.integer(s.budget.maximumRetainedSourcePairs);
    const auto& n=s.solve.source.solver;
    for(auto i:{int(n.algorithm),n.maximumIterations,n.maximumLineSearchSteps,n.lbfgsMemory,n.audit.lpIterationLimit})h.integer(i);
    for(auto d:{n.equalityTolerance,n.projectedGradientTolerance,n.activeSetTolerance,n.meritPenalty,n.lineSearchReduction,
        n.armijo,n.audit.restorationTrustBound,n.audit.lpTimeLimit,n.audit.lpFeasibilityTolerance})h.real(d);
    const auto& b=s.solve.source.budget;
    for(auto i:{b.maximumStateValues,b.maximumSourceEvaluations,b.maximumTraceRows,b.maximumCallbackRecords,b.maximumAdapterBytes,
        s.solve.budget.maximumGuardCalls,s.solve.budget.maximumRetainedStageRecords,s.solve.budget.maximumGuardedBytes,
        s.rebind.budget.maximumOwnedBytes})h.integer(i);
    geometryHash(h,s.solve.geometry);geometryHash(h,s.rebind.geometry);h.real(s.rebind.fractionalEpsilon);return h.finish();
}
void terminal(SourceCageAdaptiveSolverFacts& f,Outcome outcome,const char* reason) {
    f.outcome=outcome;const auto n=std::min(f.terminalReason.size()-1,std::strlen(reason));
    std::memcpy(f.terminalReason.data(),reason,n);f.terminalReason[n]=0;
    Hash h;h.string(reason);f.terminalReasonSha256=h.finish();
}
void reserve(SourceCageAdaptiveSolverFacts& f,std::size_t bytes) {
    const auto total=add(f.resources.cumulativeKnownBufferBytes,bytes);
    if(total>f.settings.budget.maximumKnownBufferBytes)resource("aggregate known-buffer reservation cap before work");
    f.resources.cumulativeKnownBufferBytes=total;
    if(!f.passes.empty())f.passes.back().cumulativeKnownBufferBytes=total;
}
void fractional(SourceCageAdaptiveSolverFacts& f,const SourceCagePaperFunctional& source,const Eigen::VectorXd& x) {
    f.remainingFractionalLeaves=0;f.refinableFractionalLeaves=0;
    const auto& d=source.source().boundary().binding().domain();const auto eps=f.settings.rebind.fractionalEpsilon;
    for(std::size_t i=0;i<d.leaves().size();++i)if(x[i]>=eps&&x[i]<=1.-eps) {
        ++f.remainingFractionalLeaves;if(d.leaves()[i].level<d.maximumOptimizationLevel())++f.refinableFractionalLeaves;
    }
    f.fractionalCountsAvailable=true;
}
bool sameState(const Eigen::VectorXd& a,const Eigen::VectorXd& b) {
    return a.size()==b.size()&&std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0;
}
bool solved(const SourceCageGuardedSlqpFacts& g) {
    const auto& s=g.source;
    return s.outcome==SourceCageSlqpOutcome::Completed&&s.numericalResultAvailable&&s.numerical.status==SlqpStatus::Converged&&
        g.finalStatus==SourceCageGuardedFinalStatus::ChecksPassed&&g.finalRecheckCount==1&&s.finalEvaluationAvailable&&
        s.recordCount>0&&s.records[s.recordCount-1].purpose==SourceCageSlqpEvaluationPurpose::FinalAcceptedState&&
        s.records[s.recordCount-1].status==SourceCageSlqpEvaluationStatus::Evaluated&&
        !g.attempts.empty()&&g.attempts.back().checkCompleted&&g.attempts.back().outcome==SourceCageGeometryOutcome::ChecksPassed;
}
void finish(SourceCageAdaptiveSolverFacts& f,bool hasSolved,bool hasPending) {
    Hash h;h.string("actual-source-adaptive-solver-scientific-v1 relaxed-only");
    for(const auto* s:{&f.configSha256,&f.initialFunctionalSha256,&f.initialStateSha256,&f.terminalReasonSha256})h.string(*s);
    h.integer(unsigned(f.outcome));h.string(f.terminalReason.data());h.integer(hasSolved);h.integer(hasPending);
    const auto& r=f.resources;
    for(auto n:{r.cumulativeKnownBufferBytes,r.guardedReservedBytes,r.rebindReservedBytes,r.admittedSourcePairs,r.retainedSourcePairs,
        f.remainingFractionalLeaves,f.refinableFractionalLeaves})h.integer(n);
    h.integer(f.fractionalCountsAvailable);h.integer(f.passes.size());
    for(const auto& p:f.passes) {
        for(auto n:{p.passNumber,p.startingStateValues,p.startingBetaValues,p.startingHandleValues,p.cumulativeKnownBufferBytes})h.integer(n);
        h.integer(p.hasAngle);h.integer(p.freshHistory);h.integer(bool(p.solve));h.integer(bool(p.rebind));
        for(const auto* s:{&p.oldPartitionSha256,&p.newPartitionSha256,&p.oldStateSha256,&p.solvedStateSha256,&p.newStateSha256,
            &p.oldFunctionalSha256,&p.newFunctionalSha256,&p.guardedScientificSha256,&p.rebindScientificSha256})h.string(*s);
        h.values(p.initialState);
    }
    f.scientificSha256=h.finish();
}
}
struct SourceCageAdaptiveSolverAudit::Data {
    SourceCageAdaptiveSolverFacts facts;
    std::shared_ptr<const SourceCageAdaptiveSolvedEndpoint> solved;
    std::shared_ptr<const SourceCageAdaptiveCheckpoint> pending;
};
const SourceCageAdaptiveSolverFacts& SourceCageAdaptiveSolverAudit::facts() const {
    if(!data_)throw std::logic_error("moved-from adaptive solver audit has no evidence");return data_->facts;
}
std::shared_ptr<const SourceCageAdaptiveSolvedEndpoint> SourceCageAdaptiveSolverAudit::lastSolved() const {(void)facts();return data_->solved;}
std::shared_ptr<const SourceCageAdaptiveCheckpoint> SourceCageAdaptiveSolverAudit::pendingCheckpoint() const {(void)facts();return data_->pending;}

SourceCageAdaptiveSolverAudit solveSourceCageAdaptive(const SourceCagePaperFunctional& original,
    const Eigen::VectorXd& initial,const SourceCageAdaptiveSolverSettings& settings) {
    if(!settings.maximumPasses||initial.size()!=Eigen::Index(original.stateCount())||!initial.allFinite())
        throw std::invalid_argument("adaptive solver requires positive passes and finite full-state dimensions");
    if((initial.head(original.betaCount()).array()<0.).any()||(initial.head(original.betaCount()).array()>1.).any())
        throw std::invalid_argument("adaptive solver beta must be in [0,1]");
    for(std::size_t h=0;h<original.handleCount();++h)if(initial[original.handleOffset()+4*h]<=0.)
        throw std::invalid_argument("adaptive solver handle scales must be positive");
    if(!std::isfinite(settings.rebind.fractionalEpsilon)||settings.rebind.fractionalEpsilon<=0.||settings.rebind.fractionalEpsilon>=.5)
        throw std::invalid_argument("adaptive solver epsilon must be in (0,.5)");
    if(!settings.rebind.budget.maximumOwnedBytes||settings.rebind.budget.maximumOwnedBytes>std::size_t(PTRDIFF_MAX))
        throw std::invalid_argument("invalid rebind owned byte cap");
    const SourceCageGeometryChecker validateRebind(original,settings.rebind.geometry);
    auto data=std::make_shared<SourceCageAdaptiveSolverAudit::Data>();auto& out=data->facts;
    out.settings=settings;out.configSha256=configHash(settings);
    out.initialFunctionalSha256=functionalHash(original);out.initialStateSha256=stateHash(initial);
    auto run=[&] {
        const auto pairBound=add(settings.maximumPasses,1);
        if(!settings.budget.maximumRetainedSourcePairs)resource("zero retained source pair cap");
        // Reserve every pass's fixed tape entry up front, including vector growth,
        // shared controls, hashes and fixed failure diagnostics. No maxPasses loop
        // or capacity allocation is admitted before checked arithmetic/cap checks.
        reserve(out,add(sizeof(SourceCageAdaptiveSolverFacts)+8192,
            mul(settings.maximumPasses,mul(2,sizeof(SourceCageAdaptiveSolverPass)+4096))));
        out.passes.reserve(settings.maximumPasses);
        std::shared_ptr<const SourceCagePaperFunctional> current;
        const Eigen::VectorXd* state=&initial;
        for(std::size_t i=0;i<settings.maximumPasses;++i) {
            out.passes.emplace_back();auto& p=out.passes.back();p.passNumber=i+1;
            p.cumulativeKnownBufferBytes=out.resources.cumulativeKnownBufferBytes;
            const auto& f=current?*current:original;
            p.startingStateValues=f.stateCount();p.startingBetaValues=f.betaCount();
            p.startingHandleValues=mul(f.handleCount(),4);p.hasAngle=f.hasAngle();
            p.oldPartitionSha256=f.source().boundary().binding().domain().partitionHash();
            p.oldFunctionalSha256=functionalHash(f);p.oldStateSha256=stateHash(*state);
            // Initial state snapshot plus immutable source/solved wrappers. Guarded
            // output owns the solved vector; endpoint access aliases that vector.
            reserve(out,add(mul(f.stateCount(),sizeof(double)),4096));
            if(!current) {
                current=std::make_shared<const SourceCagePaperFunctional>(original);
                out.resources.admittedSourcePairs=out.resources.retainedSourcePairs=1;
            }
            p.source=current;p.initialState=*state;
            const auto planned=planSourceCageGuardedSlqpResources(f,settings.solve);
            reserve(out,planned.requiredBytes);
            out.resources.guardedReservedBytes=add(out.resources.guardedReservedBytes,planned.requiredBytes);
            p.solve=std::make_shared<const SourceCageGuardedSlqpAudit>(solveSourceCageGuardedSlqp(f,*state,settings.solve));
            const auto& g=p.solve->facts();const auto& raw=g.source;
            p.guardedScientificSha256=g.scientificSha256;
            p.freshHistory=raw.numericalResultAvailable&&!raw.numerical.trace.empty()&&
                raw.numerical.trace.front().historyBefore==0&&raw.numerical.trace.front().lbfgsHistorySize==0;
            if(!solved(g)) {
                if(raw.outcome==SourceCageSlqpOutcome::ResourceLimit||g.finalStatus==SourceCageGuardedFinalStatus::ResourceLimit)
                    terminal(out,Outcome::ResourceLimit,"guarded solve resource limit; no adaptation");
                else if(g.finalStatus==SourceCageGuardedFinalStatus::RejectedGeometry||
                    (raw.outcome==SourceCageSlqpOutcome::InitialEvaluationFailure&&!g.attempts.empty()&&
                     g.attempts.front().checkCompleted&&g.attempts.front().outcome==SourceCageGeometryOutcome::RejectedGeometry))
                    terminal(out,Outcome::SolveGeometryRejected,"actual guarded solve geometry rejected; no adaptation");
                else if(raw.numericalResultAvailable&&raw.numerical.status!=SlqpStatus::Converged)
                    terminal(out,Outcome::SolveNotConverged,"raw numerical termination is not Converged; no adaptation or infeasibility conclusion");
                else terminal(out,Outcome::SolveUnknown,"guarded solve lacks fresh completed final geometry/material evidence");
                return;
            }
            p.solvedStateSha256=stateHash(raw.numerical.finalState);
            data->solved=std::shared_ptr<const SourceCageAdaptiveSolvedEndpoint>(new SourceCageAdaptiveSolvedEndpoint(current,p.solve));
            data->pending.reset();
            fractional(out,f,raw.numerical.finalState);
            // Each actual transaction can retain another distinct immutable graph,
            // including an unchanged rebuild. Reserve its slot before any rebuild.
            if(out.resources.admittedSourcePairs>=std::min(pairBound,settings.budget.maximumRetainedSourcePairs))
                resource("retained source/checkpoint pair cap before rebind");
            reserve(out,add(settings.rebind.budget.maximumOwnedBytes,4096));
            out.resources.rebindReservedBytes=add(out.resources.rebindReservedBytes,settings.rebind.budget.maximumOwnedBytes);
            ++out.resources.admittedSourcePairs;
            p.rebind=std::make_shared<const SourceCageAdaptiveRebindAudit>(rebindSourceCageAdaptive(f,raw.numerical.finalState,settings.rebind));
            const auto& r=p.rebind->facts();p.rebindScientificSha256=r.scientificSha256;
            p.newPartitionSha256=r.newPartitionSha256;p.newStateSha256=r.newStateSha256;p.newFunctionalSha256=r.newFunctionalSha256;
            const auto checkpoint=p.rebind->checkpoint();
            if(!checkpoint||r.outcome!=SourceCageAdaptiveRebindOutcome::Published) {
                if(r.outcome==SourceCageAdaptiveRebindOutcome::ResourceLimit)terminal(out,Outcome::ResourceLimit,"actual rebind resource limit; no candidate checkpoint published");
                else if(r.outcome==SourceCageAdaptiveRebindOutcome::RejectedGeometry)terminal(out,Outcome::RebindGeometryRejected,"actual rebind geometry rejected; no candidate checkpoint published");
                else terminal(out,Outcome::RebindUnknown,"actual rebind unknown; no candidate checkpoint published");
                return;
            }
            ++out.resources.retainedSourcePairs;
            if(r.oldPartitionSha256!=p.oldPartitionSha256||r.oldFunctionalSha256!=p.oldFunctionalSha256||r.oldStateSha256!=p.solvedStateSha256||
                r.newPartitionSha256!=checkpoint->functional().source().boundary().binding().domain().partitionHash()||
                r.newFunctionalSha256!=functionalHash(checkpoint->functional())||r.newStateSha256!=stateHash(checkpoint->state())) {
                terminal(out,Outcome::IdentityMismatch,"actual rebind input/output identity mismatch");return;
            }
            if(p.newPartitionSha256==p.oldPartitionSha256) {
                if(p.newFunctionalSha256!=p.oldFunctionalSha256||p.newStateSha256!=p.solvedStateSha256||
                    !sameState(checkpoint->state(),raw.numerical.finalState)||out.refinableFractionalLeaves) {
                    terminal(out,Outcome::IdentityMismatch,"unchanged partition lacks identity or has refinable fractional leaves");return;
                }
                terminal(out,Outcome::RelaxedAdaptiveComplete,"unchanged converged checked partition; maximum-depth fractional leaves may remain; relaxed completion only");return;
            }
            data->pending=checkpoint;
            // Aliasing ownership keeps the checkpoint's immutable functional alive
            // without cloning any graph or transferring old numerical history.
            current=std::shared_ptr<const SourceCagePaperFunctional>(checkpoint,&checkpoint->functional());
            state=&checkpoint->state();
        }
        terminal(out,Outcome::MaximumPasses,"maximum passes reached with changed rebound checkpoint not yet solved");
    };
    try {run();}
    catch(const std::bad_alloc&) {terminal(out,Outcome::ResourceLimit,"allocation failure; prior evidence and endpoints retained");}
    catch(const std::length_error& e) {terminal(out,Outcome::ResourceLimit,e.what());}
    // Invalid public settings may throw. Diagnostic allocation can also propagate;
    // no incomplete audit or candidate endpoint escapes in those exceptional paths.
    finish(out,bool(data->solved),bool(data->pending));
    return SourceCageAdaptiveSolverAudit(std::move(data));
}
}
