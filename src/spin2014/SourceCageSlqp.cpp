#include "spin2014/SourceCageSlqp.h"
#include "spin2014/SourceCageGuardedSlqp.h"
#include "spin2014/SlqpAuditTrace.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>

namespace spin2014 {
namespace {
[[noreturn]] void resource(const char* reason) {
    throw std::length_error(std::string("RESOURCE_LIMIT: source cage SLQP ")+reason);
}
std::size_t add(std::size_t a,std::size_t b) {
    if(a>std::numeric_limits<std::size_t>::max()-b)resource("size addition overflow");return a+b;
}
std::size_t mul(std::size_t a,std::size_t b) {
    if(b&&a>std::numeric_limits<std::size_t>::max()/b)resource("size multiplication overflow");return a*b;
}
std::size_t sum(std::initializer_list<std::size_t> list) {
    std::size_t n=0;for(auto v:list)n=add(n,v);return n;
}
class Digest {
    verification::Sha256 hash_;
public:
    void integer(std::uint64_t x) {
        unsigned char b[8];for(int i=0;i<8;++i)b[7-i]=static_cast<unsigned char>(x>>(8*i));hash_.update(b,8);
    }
    void string(const char* x) {
        const auto n=std::strlen(x);integer(n);hash_.update(reinterpret_cast<const unsigned char*>(x),n);
    }
    void string(const std::string& x) {integer(x.size());hash_.update(reinterpret_cast<const unsigned char*>(x.data()),x.size());}
    void real(double x) {std::uint64_t bits;static_assert(sizeof(x)==8,"IEEE double required");std::memcpy(&bits,&x,8);integer(bits);}
    template<class D> void matrix(const Eigen::MatrixBase<D>& m) {
        integer(m.rows());integer(m.cols());for(Eigen::Index i=0;i<m.rows();++i)for(Eigen::Index j=0;j<m.cols();++j)real(m(i,j));
    }
    std::string finish() {return hash_.finish();}
};
template<std::size_t N> void text(std::array<char,N>& out,const char* in) {
    std::size_t i=0;for(;i+1<N&&in[i];++i)out[i]=in[i];out[i]='\0';
}
void failure(SourceCageSlqpEvaluationRecord& event,SourceCageSlqpEvaluationStatus status,const char* message) {
    event.status=status;text(event.reason,message);Digest hash;hash.string(message);text(event.reasonSha256,hash.finish().c_str());
}
void markBookkeepingResource(SourceCageSlqpResult& out) noexcept {
    out.outcome=SourceCageSlqpOutcome::ResourceLimit;
    if(!out.recordCount)return;
    auto& event=out.records[out.recordCount-1];
    event.status=SourceCageSlqpEvaluationStatus::ResourceLimit;
    text(event.reason,"RESOURCE_LIMIT: callback provenance allocation failure");
    event.reasonSha256.fill('\0'); // Recovery hash is deferred until memory is available.
}
void validateSettings(const SourceCageSlqpSettings& settings) {
    if(!settings.solver.iterationStorePath.empty())
        throw std::invalid_argument("disk iteration storage is unsupported for cage post-solve phi mutation");
    const auto& s=settings.solver;
    if(s.directionPolicy!=SlqpDirectionPolicy::ProjectedGradientV1)
        throw std::invalid_argument("metric tangent rescue is not yet admitted for source cage deformation");
    if(s.algorithm!=SlqpAlgorithm::AuditedKktRestorationV1)
        throw std::invalid_argument("actual source cage SLQP requires AuditedKktRestorationV1");
    validateSlqpAuditSettings(s.audit);
    if(s.maximumIterations<=0||s.maximumLineSearchSteps<=0||s.lbfgsMemory<=0||
       !std::isfinite(s.equalityTolerance)||s.equalityTolerance<=0||
       !std::isfinite(s.projectedGradientTolerance)||s.projectedGradientTolerance<=0||
       !std::isfinite(s.activeSetTolerance)||s.activeSetTolerance<0||s.activeSetTolerance>=.5||
       !std::isfinite(s.meritPenalty)||s.meritPenalty<=0||
       !std::isfinite(s.lineSearchReduction)||s.lineSearchReduction<=0||s.lineSearchReduction>=1||
       !std::isfinite(s.armijo)||s.armijo<=0||s.armijo>=1)
        throw std::invalid_argument("invalid actual source cage SLQP settings");
    const auto& b=settings.budget;
    if(!b.maximumStateValues||!b.maximumSourceEvaluations||!b.maximumTraceRows||!b.maximumCallbackRecords||!b.maximumAdapterBytes)
        resource("zero cap");
    // The unchanged generic trace serializer indexes every K x N coefficient
    // with int; K can be six. Cap the public state budget conservatively too.
    if(b.maximumStateValues>std::size_t(INT_MAX/6)||b.maximumSourceEvaluations>std::size_t(INT_MAX)||
       b.maximumTraceRows>std::size_t(INT_MAX)||b.maximumCallbackRecords>std::size_t(INT_MAX)||
       b.maximumAdapterBytes>std::size_t(std::numeric_limits<std::ptrdiff_t>::max()))resource("unrepresentable cap");
}
std::string configHash(const SourceCagePaperFunctional& f,const Eigen::VectorXd& x,const SourceCageSlqpSettings& settings) {
    Digest h;h.string("source-cage-slqp-v1 IMPLEMENTATION_DEFAULT fixed-source full-state first-order-only");
    h.string("ACCOUNTING_V1; geometry=NOT_EVALUATED; feasibility=UNKNOWN; no physical release; HIGHS threads=1 parallel=off seed=0 simplex presolve=on");
    h.string(f.configSha256());h.string(f.upstreamSha256());h.string(f.source().upstreamSha256());
    h.string(f.source().boundary().binding().bbw().outputSha256());
    h.integer(f.betaCount());h.integer(f.handleCount());h.integer(f.hasAngle());h.matrix(x);
    const auto& s=settings.solver;
    h.integer(static_cast<unsigned>(s.algorithm));h.integer(s.maximumIterations);h.integer(s.maximumLineSearchSteps);h.integer(s.lbfgsMemory);
    h.real(s.equalityTolerance);h.real(s.projectedGradientTolerance);h.real(s.activeSetTolerance);h.real(s.meritPenalty);
    h.real(s.lineSearchReduction);h.real(s.armijo);h.real(s.audit.restorationTrustBound);h.integer(s.audit.lpIterationLimit);
    h.real(s.audit.lpTimeLimit);h.real(s.audit.lpFeasibilityTolerance);
    const auto& b=settings.budget;
    h.integer(b.maximumStateValues);h.integer(b.maximumSourceEvaluations);h.integer(b.maximumTraceRows);
    h.integer(b.maximumCallbackRecords);h.integer(b.maximumAdapterBytes);return h.finish();
}
void finishIdentity(SourceCageSlqpResult& out) {
    Digest tape;tape.string("source-cage-slqp-call-tape-v1");
    tape.integer(out.recordCount);tape.integer(out.callbackCount);tape.integer(out.sourceEvaluationCount);
    for(std::size_t i=0;i<out.recordCount;++i) {
        auto& r=out.records[i];
        if(r.status==SourceCageSlqpEvaluationStatus::Unused)resource("incomplete attempted callback record");
        if(std::strlen(r.stateSha256.data())!=64) {
            Digest state;state.matrix(out.callbackStates.col(i));text(r.stateSha256,state.finish().c_str());
        }
        if(r.status!=SourceCageSlqpEvaluationStatus::Evaluated&&std::strlen(r.reasonSha256.data())!=64) {
            if(!r.reason[0])resource("missing callback failure reason");
            Digest reason;reason.string(r.reason.data());text(r.reasonSha256,reason.finish().c_str());
        }
        if(std::strlen(r.stateSha256.data())!=64||
           (r.status==SourceCageSlqpEvaluationStatus::Evaluated?std::strlen(r.outputSha256.data())!=64:
            std::strlen(r.reasonSha256.data())!=64))resource("incomplete callback digest after recovery");
        tape.integer(r.sequence);tape.integer(static_cast<unsigned>(r.purpose));
        tape.integer(static_cast<unsigned>(r.status));tape.integer(r.sourceInvoked);
        tape.matrix(out.callbackStates.col(i)); // Bind the retained raw state, including recovery paths.
        tape.string(r.stateSha256.data());tape.string(r.outputSha256.data());tape.string(r.reasonSha256.data());
    }
    out.callbackTapeSha256=tape.finish();
    Digest h;h.string("source-cage-slqp-result-v1");h.string(out.configSha256);h.string(out.callbackTapeSha256);
    h.integer(static_cast<unsigned>(out.outcome));h.integer(static_cast<unsigned>(out.geometryStatus));h.integer(static_cast<unsigned>(out.feasibility));
    h.integer(out.numericalResultAvailable);
    h.string(out.numerical.traceHash);h.matrix(out.numerical.finalState);h.matrix(out.numerical.beta);
    h.real(out.numerical.initialPhi);h.real(out.numerical.finalPhi);h.integer(static_cast<unsigned>(out.numerical.formulation));
    // Generic serializer binds membershipHash; additionally bind raw membership
    // and its diagnostic counts instead of trusting the supplied hash alone.
    for(const auto& t:out.numerical.trace) {
        h.integer(t.activeMembership.size());for(auto v:t.activeMembership)h.integer(v);
        h.integer(t.activeVariables);h.integer(t.lbfgsHistorySize);
    }
    h.integer(out.finalEvaluationAvailable);
    if(out.finalEvaluationAvailable) {h.string(out.finalEvaluation.outputSha256);h.string(out.finalEvaluation.world.outputSha256);}
    out.scientificSha256=h.finish();
}
// File-private gate: no user callback, evaluation, or boolean can supply a pass.
struct GuardContext {
    const SourceCageGeometryChecker checker;
    SourceCageGuardedSlqpFacts& facts;
    void begin(const SourceCageSlqpEvaluationRecord& event) {
        if(facts.attempts.size()>=facts.resources.attemptCapacity)resource("guard attempt capacity exceeded");
        auto& a=facts.attempts.emplace_back();a.callbackRecord=event.sequence;a.firstStage=facts.stages.size();
    }
    void bookkeepingResource() noexcept {
        if(!facts.attempts.empty()&&!facts.attempts.back().checkCompleted)
            facts.attempts.back().outcome=SourceCageGeometryOutcome::ResourceLimit;
    }
    void check(const Eigen::VectorXd& x,SourceCageSlqpResult& out,SourceCageSlqpEvaluationRecord& event) {
        auto& a=facts.attempts.back();
        if(out.outcome==SourceCageSlqpOutcome::ResourceLimit||
           facts.checkerInvocationCount>=facts.settings.budget.maximumGuardCalls) {
            a.outcome=SourceCageGeometryOutcome::ResourceLimit;
            out.outcome=SourceCageSlqpOutcome::ResourceLimit;
            failure(event,SourceCageSlqpEvaluationStatus::ResourceLimit,"RESOURCE_LIMIT: guard call cap or prior resource latch");
            resource("guard call cap or prior resource latch");
        }
        a.checkerInvoked=true;++facts.checkerInvocationCount;
        if(event.purpose==SourceCageSlqpEvaluationPurpose::FinalAcceptedState)++facts.finalRecheckCount;
        bool returned=false;
        try {
            const auto audit=checker.check(x);returned=true;const auto& g=audit.facts();
            for(const auto& stage:g.records)if(stage.stage==SourceCageGeometryStage::SourceEvaluation) {
                a.sourceValuesInvoked=true;++facts.geometryValueEvaluationCount;break;
            }
            a.outcome=g.outcome;a.plannedCheckerBytes=g.plannedCheckerBytes;
            a.pairsVisited=g.pairsVisited;a.exactPairRelations=g.exactPairRelations;
            a.exactPairCacheHits=g.exactPairCacheHits;
            a.exactPairCacheMisses=g.exactPairCacheMisses;
            a.exactPairBooleanEvaluations=g.exactPairBooleanEvaluations;
            facts.exactPairRelations=add(facts.exactPairRelations,g.exactPairRelations);
            facts.exactPairCacheHits=add(facts.exactPairCacheHits,g.exactPairCacheHits);
            facts.exactPairCacheMisses=add(facts.exactPairCacheMisses,g.exactPairCacheMisses);
            facts.exactPairBooleanEvaluations=add(facts.exactPairBooleanEvaluations,g.exactPairBooleanEvaluations);
            text(a.inputSha256,g.inputSha256.c_str());text(a.sourceOutputSha256,g.sourceOutputSha256.c_str());
            text(a.scientificSha256,g.scientificSha256.c_str());
            if(g.records.size()>facts.resources.stageCapacity-facts.stages.size())resource("retained geometry stages exceeded");
            for(const auto& stage:g.records) {
                if(stage.squaredDistanceExact.size()>facts.settings.geometry.maximumExactStringBytes)resource("retained exact string exceeded");
                facts.stages.push_back(stage);++a.stageCount;
            }
            a.checkCompleted=true;
            if(g.outcome==SourceCageGeometryOutcome::ResourceLimit) {
                out.outcome=SourceCageSlqpOutcome::ResourceLimit;
                failure(event,SourceCageSlqpEvaluationStatus::ResourceLimit,g.records.back().reason.data());
                resource("actual geometry resource failure");
            }
            if(g.outcome!=SourceCageGeometryOutcome::ChecksPassed) {
                failure(event,SourceCageSlqpEvaluationStatus::InvalidNumericalTrial,g.records.back().reason.data());
                throw std::runtime_error("actual geometry trial was rejected or incomplete");
            }
        }catch(const std::bad_alloc&) {
            if(!returned)facts.geometryValueEvaluationCountComplete=false;
            a.checkCompleted=false;a.outcome=SourceCageGeometryOutcome::ResourceLimit;
            markBookkeepingResource(out);throw;
        }catch(const std::length_error&) {
            if(!returned)facts.geometryValueEvaluationCountComplete=false;
            a.outcome=SourceCageGeometryOutcome::ResourceLimit;out.outcome=SourceCageSlqpOutcome::ResourceLimit;
            if(event.status!=SourceCageSlqpEvaluationStatus::ResourceLimit) {
                a.checkCompleted=false;markBookkeepingResource(out);
            }throw;
        }
    }
};
}

SourceCageSlqpResources planSourceCageSlqpResources(const SourceCagePaperFunctional& f,const SourceCageSlqpSettings& s) {
    validateSettings(s);SourceCageSlqpResources r;
    const auto b=f.betaCount(),h=f.handleCount(),n=add(add(b,f.hasAngle()?1:0),mul(4,h)),k=f.constraintCount();
    const auto p=std::size_t(f.source().boundary().points().rows()),outer=f.source().boundary().outerPointCount();
    if(n!=f.stateCount()||n>s.budget.maximumStateValues||n>std::size_t(INT_MAX/4)||
       add(n,mul(2,k))>std::size_t(INT_MAX)||(k!=2&&k!=5&&k!=6))resource("state or native LP dimensions");
    if(mul(k,n)>std::size_t(INT_MAX))resource("native trace matrix coefficient count");
    const auto i=std::size_t(s.solver.maximumIterations),l=std::size_t(s.solver.maximumLineSearchSteps),m=std::size_t(s.solver.lbfgsMemory);
    r.stateValues=n;r.constraints=k;r.traceRows=add(i,1);r.callbackRecords=add(2,mul(mul(3,i),l));
    if(r.traceRows>s.budget.maximumTraceRows||r.callbackRecords>s.budget.maximumCallbackRecords)
        resource("trace or callback record cap");
    const auto tapeValues=mul(n,r.callbackRecords);
    if(tapeValues>std::size_t(std::numeric_limits<Eigen::Index>::max()))resource("Eigen tape dimensions");
    r.tapeBytes=add(mul(8,tapeValues),mul(r.callbackRecords,sizeof(SourceCageSlqpEvaluationRecord)));
    const auto sourceValues=sum({mul(3,p),mul(10,add(b,mul(4,h))),mul(add(11,k),n),k});
    // Two complete outputs cover transient/final coexistence and final return
    // duplication. Actual callback mapping moves gradient, J and constraints.
    r.sourceOutputBytes=mul(2,add(mul(8,sourceValues),add(sizeof(SourceCagePaperEvaluation),1024)));
    // Split beta/handles + fill + BBW weight rows + shape residual/vertex/handle
    // gradients. All upstream caps remain enforced independently as well.
    r.sourceWorkspaceBytes=add(mul(8,sum({n,b,mul(2,h),mul(6,outer),mul(4,h)})),16384);
    const auto rowValues=add(mul(add(k,5),n),mul(2,k));
    const auto trials=mul(3,l);
    const auto rowBytes=sum({sizeof(SlqpIteration),mul(8,rowValues),n,1024,
        mul(mul(2,trials),add(sizeof(RestorationTrial),128))});
    // Generic vector growth <=2T + finish() full-result copy <=T + local rows.
    r.traceBytes=mul(add(mul(3,r.traceRows),2),rowBytes);
    r.historyBytes=sum({mul(8,add(mul(mul(2,m),n),m)),mul(128,m),4096});
    // Current/result/candidate evaluations, inverseConstraint/freeA, Schur/COD,
    // GEMM packing, state/direction/KKT temporaries; no N x N buffer is used.
    r.solverWorkspaceBytes=add(mul(8,sum({mul(32,n),mul(mul(12,k),n),mul(16,mul(k,k)),mul(8,m)})),16384);
    const auto nnz=std::max(mul(mul(2,n),add(k,1)),mul(k,add(n,2)));
    const auto rows=std::max(mul(2,n),k),cols=std::max(add(k,1),add(n,mul(2,k)));
    // Explicit LP value/index arrays (8-byte index allowance), vector capacity,
    // row scratch and known passModel copy. HiGHS opaque internals are separate.
    r.explicitLpBytes=add(mul(4,sum({mul(16,nnz),mul(24,rows),mul(24,cols)})),16384);
    const auto serializedRow=sum({mul(32,add(rowValues,40)),1024,mul(trials,192)});
    const auto serialized=add(mul(add(r.traceRows,1),serializedRow),add(mul(64,n),4096));
    // Existing generic ostringstream and str() full serialization copies, with
    // growth reserve. Adapter SHA itself streams and needs fixed workspace only.
    r.hashWorkspaceBytes=add(mul(8,serialized),16384);
    r.requiredBytes=sum({r.tapeBytes,r.sourceOutputBytes,r.sourceWorkspaceBytes,r.traceBytes,r.historyBytes,
        r.solverWorkspaceBytes,r.explicitLpBytes,r.hashWorkspaceBytes,mul(8,n),sizeof(SourceCageSlqpResult),16384});
    if(r.requiredBytes>s.budget.maximumAdapterBytes)resource("simultaneous known buffer cap");
    return r;
}

static SourceCageSlqpResult runSourceCageSlqp(const SourceCagePaperFunctional& f,
    const Eigen::VectorXd& initial,const SourceCageSlqpSettings& settings,GuardContext* guard) {
    const auto resources=planSourceCageSlqpResources(f,settings);
    if(std::size_t(initial.size())!=f.stateCount()||!initial.allFinite()||
       (initial.head(f.betaCount()).array()<0).any()||(initial.head(f.betaCount()).array()>1).any())
        throw std::invalid_argument("invalid source cage SLQP full initial state");
    for(std::size_t h=0;h<f.handleCount();++h)if(initial[f.handleOffset()+4*h]<=0)
        throw std::invalid_argument("source cage SLQP initial handle scale must be positive");
    try {
        SourceCageSlqpResult out;out.settings=settings;out.resources=resources;out.initialState=initial;
        out.configSha256=configHash(f,initial,settings);
        out.records.resize(resources.callbackRecords);
        out.callbackStates=Eigen::MatrixXd::Zero(resources.stateValues,resources.callbackRecords);
        const auto evaluate=[&](const Eigen::VectorXd& x,SourceCageSlqpEvaluationPurpose purpose) {
            if(out.recordCount>=resources.callbackRecords)resource("internal callback bound exceeded");
            const auto index=out.recordCount++;auto& event=out.records[index];event.sequence=index;event.purpose=purpose;
            if(guard)guard->begin(event);
            if(purpose==SourceCageSlqpEvaluationPurpose::SolverCallback)++out.callbackCount;
            out.callbackStates.col(index)=x;Digest state;state.matrix(x);text(event.stateSha256,state.finish().c_str());
            if(guard)guard->check(x,out,event);
            if(out.outcome==SourceCageSlqpOutcome::ResourceLimit||out.sourceEvaluationCount>=settings.budget.maximumSourceEvaluations) {
                out.outcome=SourceCageSlqpOutcome::ResourceLimit;
                failure(event,SourceCageSlqpEvaluationStatus::ResourceLimit,"RESOURCE_LIMIT: source evaluation cap or prior resource latch");
                resource("source evaluation cap or prior resource latch");
            }
            ++out.sourceEvaluationCount;event.sourceInvoked=true;
            try {
                auto e=f.evaluate(x);text(event.outputSha256,e.outputSha256.c_str());
                event.status=SourceCageSlqpEvaluationStatus::Evaluated;return e;
            }catch(const std::length_error& e) {
                out.outcome=SourceCageSlqpOutcome::ResourceLimit;failure(event,SourceCageSlqpEvaluationStatus::ResourceLimit,e.what());throw;
            }catch(const std::bad_alloc&) {
                out.outcome=SourceCageSlqpOutcome::ResourceLimit;failure(event,SourceCageSlqpEvaluationStatus::ResourceLimit,"RESOURCE_LIMIT: source allocation failure");throw;
            }catch(const std::exception& e) {
                failure(event,SourceCageSlqpEvaluationStatus::InvalidNumericalTrial,e.what());throw;
            }catch(...) {
                failure(event,SourceCageSlqpEvaluationStatus::InvalidNumericalTrial,"NONSTANDARD_SOURCE_EXCEPTION");throw;
            }
        };
        EqualityBoxProblem problem;problem.variableCount=f.stateCount();problem.betaCount=f.betaCount();
        problem.evaluate=[&](const Eigen::VectorXd& x) {
          try {
            auto source=evaluate(x,SourceCageSlqpEvaluationPurpose::SolverCallback);
            AdaptiveFillEvaluation e;e.objective=source.totalObjective;e.paperObjective=source.paperObjective;
            e.regularizationEnergy=source.fillPenalty+source.shapeEnergy+source.displacementEnergy;
            e.gradient=std::move(source.jointGradient);e.constraints=std::move(source.constraints);
            e.constraintJacobian=std::move(source.jointConstraintJacobian);e.finalMoments=source.spinMoments;return e;
          }catch(const std::bad_alloc&) {
            // Also latch allocation failures in state/output hashing and mapping,
            // outside the source invocation itself: generic trial catches all.
            markBookkeepingResource(out);if(guard)guard->bookkeepingResource();throw;
          }catch(const std::length_error&) {
            out.outcome=SourceCageSlqpOutcome::ResourceLimit;
            if(out.recordCount&&out.records[out.recordCount-1].status==SourceCageSlqpEvaluationStatus::Unused)
                markBookkeepingResource(out);
            if(guard)guard->bookkeepingResource();
            throw;
          }
        };
        try {
            out.numerical=solveSlqpLbfgs(problem,initial,settings.solver);
            out.numericalResultAvailable=true;
        }catch(...) {
            // Generic initial evaluation is outside its trial catch/recovery.
            // Retain the actual attempted source record, without synthesizing
            // any generic trace, final state, KKT witness, or termination.
            if(out.callbackCount==1&&out.recordCount==1&&(out.records[0].sourceInvoked||guard)&&
               (out.records[0].status!=SourceCageSlqpEvaluationStatus::Evaluated||
                out.outcome==SourceCageSlqpOutcome::ResourceLimit)) {
                if(out.outcome!=SourceCageSlqpOutcome::ResourceLimit)
                    out.outcome=SourceCageSlqpOutcome::InitialEvaluationFailure;
                finishIdentity(out);return out;
            }
            throw;
        }
        out.numerical.formulation=Spin2014Formulation::PaperGivens;
        out.numerical.initialPhi=f.hasAngle()?initial[f.angleOffset()]:0.;
        out.numerical.finalPhi=f.hasAngle()?out.numerical.finalState[f.angleOffset()]:0.;
        for(auto& t:out.numerical.trace)t.phi=f.hasAngle()?t.state[f.angleOffset()]:0.;
        out.numerical.traceHash=hashAuditedSlqp(out.numerical,settings.solver);
        try {
            out.finalEvaluation=evaluate(out.numerical.finalState,SourceCageSlqpEvaluationPurpose::FinalAcceptedState);
            out.finalEvaluationAvailable=true;
        }catch(const std::length_error&) {out.outcome=SourceCageSlqpOutcome::ResourceLimit;}
        catch(const std::bad_alloc&) {markBookkeepingResource(out);if(guard)guard->bookkeepingResource();}
        catch(...) {if(!guard)throw;}
        finishIdentity(out);return out;
    }catch(const std::bad_alloc&) {resource("allocation failure outside recoverable callback");}
}
SourceCageSlqpResult solveSourceCageSlqp(const SourceCagePaperFunctional& f,
    const Eigen::VectorXd& initial,const SourceCageSlqpSettings& settings) {
    return runSourceCageSlqp(f,initial,settings,nullptr);
}

SourceCageGuardedSlqpResources planSourceCageGuardedSlqpResources(const SourceCagePaperFunctional& f,
    const SourceCageGuardedSlqpSettings& s) {
    const SourceCageGeometryChecker validate(f,s.geometry);
    SourceCageGuardedSlqpResources r;r.adapter=planSourceCageSlqpResources(f,s.source);
    for(auto cap:{s.budget.maximumGuardCalls,s.budget.maximumRetainedStageRecords})
        if(!cap||cap>std::size_t(INT_MAX))resource("zero or unrepresentable guard count cap");
    if(!s.budget.maximumGuardedBytes||s.budget.maximumGuardedBytes>std::size_t(PTRDIFF_MAX))resource("zero or unrepresentable guarded byte cap");
    const auto p=f.source().boundary().primitiveBoundaries().size();
    r.stagesPerGuard=add(3,add(mul(3,p),p?mul(p/2,p-1+p%2):0));
    r.stageCapacity=mul(r.stagesPerGuard,s.budget.maximumGuardCalls);
    if(r.stageCapacity>s.budget.maximumRetainedStageRecords)resource("retained stage count cap");
    r.attemptCapacity=r.adapter.callbackRecords;
    // Flat vectors reserve their complete declared capacities before callbacks.
    // Only clearance stages carry dynamic strings; allow 2x bounded payload
    // for std::string growth/copy while the current checker audit coexists.
    // No mesh or Nef object is retained per attempt. Return is move-only.
    r.retainedEvidenceBytes=sum({mul(r.stageCapacity,sizeof(SourceCageGeometryRecord)),
        mul(r.attemptCapacity,sizeof(SourceCageSlqpGuardAttempt)),
        mul(mul(s.budget.maximumGuardCalls,p),mul(2,add(s.geometry.maximumExactStringBytes,1))),16384});
    // One live checker plus its returned audit/diagnostic snapshot and metadata.
    // 2x covers checker vector capacity, retained audit return copies and Nef
    // diagnostic snapshots; opaque native allocation internals remain external.
    r.checkerScratchBytes=add(mul(2,s.geometry.maximumCheckerBytes),65536);
    r.requiredBytes=sum({r.adapter.requiredBytes,r.retainedEvidenceBytes,r.checkerScratchBytes,
        sizeof(SourceCageGuardedSlqpFacts),65536});
    if(r.requiredBytes>s.budget.maximumGuardedBytes)resource("simultaneous guarded buffer cap");
    return r;
}

SourceCageGuardedSlqpAudit solveSourceCageGuardedSlqp(const SourceCagePaperFunctional& f,
    const Eigen::VectorXd& initial,const SourceCageGuardedSlqpSettings& settings) {
    const auto resources=planSourceCageGuardedSlqpResources(f,settings);
    SourceCageGuardedSlqpAudit audit;auto& out=audit.facts_;
    out.settings=settings;out.resources=resources;
    out.attempts.reserve(resources.attemptCapacity);out.stages.reserve(resources.stageCapacity);
    GuardContext guard{SourceCageGeometryChecker(f,settings.geometry),out};
    out.source=runSourceCageSlqp(f,initial,settings.source,&guard);
    if(out.source.outcome==SourceCageSlqpOutcome::ResourceLimit)out.finalStatus=SourceCageGuardedFinalStatus::ResourceLimit;
    else if(out.source.finalEvaluationAvailable)out.finalStatus=SourceCageGuardedFinalStatus::ChecksPassed;
    else if(out.finalRecheckCount&&!out.attempts.empty())out.finalStatus=
        out.attempts.back().outcome==SourceCageGeometryOutcome::RejectedGeometry?
        SourceCageGuardedFinalStatus::RejectedGeometry:SourceCageGuardedFinalStatus::Unknown;
    Digest config;config.string("source-cage-guarded-slqp-v1 accounting-v1 move-only; feasibility UNKNOWN; no seal physics stable spin or original-paper grant");
    config.string(out.source.configSha256);config.real(settings.geometry.minimumOuterDomainClearance);
    const auto& g=settings.geometry;
    for(auto v:{g.maximumPoints,g.maximumTriangles,g.maximumPrimitives,g.maximumPairs,g.maximumRecords,
        g.maximumCheckerBytes,g.maximumNefComplexity,g.maximumConnectedShells,g.maximumExactStringBytes,
        g.maximumExactPairCacheEntries,g.maximumPrimitiveRegionCacheEntries,
        g.maximumPrimitiveRegionCacheBytes,g.maximumPrimitiveRegionCacheComplexity,
        settings.budget.maximumGuardCalls,settings.budget.maximumRetainedStageRecords,settings.budget.maximumGuardedBytes})config.integer(v);
    out.configSha256=config.finish();
    Digest h;h.string(out.configSha256);h.string(out.source.scientificSha256);
    h.integer(unsigned(out.finalStatus));h.integer(unsigned(out.feasibility));
    h.integer(out.checkerInvocationCount);h.integer(out.geometryValueEvaluationCount);
    h.integer(out.geometryValueEvaluationCountComplete);h.integer(out.finalRecheckCount);
    h.integer(out.attempts.size());
    for(const auto& a:out.attempts) {
        h.integer(a.callbackRecord);h.integer(a.firstStage);h.integer(a.stageCount);
        h.integer(a.checkerInvoked);h.integer(a.checkCompleted);h.integer(a.sourceValuesInvoked);h.integer(unsigned(a.outcome));
        h.string(a.inputSha256.data());h.string(a.sourceOutputSha256.data());h.string(a.scientificSha256.data());
        h.integer(a.plannedCheckerBytes);h.integer(a.pairsVisited);h.integer(a.exactPairRelations);
    }
    h.integer(out.stages.size());
    for(const auto& stage:out.stages) {
        h.integer(unsigned(stage.stage));h.integer(unsigned(stage.outcome));
        h.integer(stage.first);h.integer(stage.second);h.integer(stage.betaOwner);h.integer(stage.secondBetaOwner);
        h.string(stage.evidenceSha256.data());h.string(stage.reasonSha256.data());h.string(stage.reason.data());
        h.string(stage.squaredDistanceExact);h.integer(stage.witnessFaces.has_value());
        if(stage.witnessFaces)for(auto face:*stage.witnessFaces)h.integer(face);
    }
    out.scientificSha256=h.finish();return audit;
}
}
