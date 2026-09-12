#include "spin2014/SourceCageSlqp.h"
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
    const auto& s=settings.solver;
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
    if(b.maximumStateValues>std::size_t(INT_MAX/4)||b.maximumSourceEvaluations>std::size_t(INT_MAX)||
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
}

SourceCageSlqpResources planSourceCageSlqpResources(const SourceCagePaperFunctional& f,const SourceCageSlqpSettings& s) {
    validateSettings(s);SourceCageSlqpResources r;
    const auto b=f.betaCount(),h=f.handleCount(),n=add(add(b,f.hasAngle()?1:0),mul(4,h)),k=f.constraintCount();
    const auto p=std::size_t(f.source().boundary().points().rows()),outer=f.source().boundary().outerPointCount();
    if(n!=f.stateCount()||n>s.budget.maximumStateValues||n>std::size_t(INT_MAX/4)||
       add(n,mul(2,k))>std::size_t(INT_MAX)||(k!=2&&k!=5&&k!=6))resource("state or native LP dimensions");
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

SourceCageSlqpResult solveSourceCageSlqp(const SourceCagePaperFunctional& f,
    const Eigen::VectorXd& initial,const SourceCageSlqpSettings& settings) {
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
            if(purpose==SourceCageSlqpEvaluationPurpose::SolverCallback)++out.callbackCount;
            out.callbackStates.col(index)=x;Digest state;state.matrix(x);text(event.stateSha256,state.finish().c_str());
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
            markBookkeepingResource(out);throw;
          }catch(const std::length_error&) {
            out.outcome=SourceCageSlqpOutcome::ResourceLimit;
            if(out.recordCount&&out.records[out.recordCount-1].status==SourceCageSlqpEvaluationStatus::Unused)
                markBookkeepingResource(out);
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
            if(out.callbackCount==1&&out.recordCount==1&&out.records[0].sourceInvoked&&
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
        catch(const std::bad_alloc&) {markBookkeepingResource(out);}
        finishIdentity(out);return out;
    }catch(const std::bad_alloc&) {resource("allocation failure outside recoverable callback");}
}
}
