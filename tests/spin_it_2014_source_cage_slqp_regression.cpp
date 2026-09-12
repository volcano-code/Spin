#include <iostream>
#if !__has_include("spin2014/SourceCageSlqp.h")
int main() {std::cerr<<"FAIL: typed actual-source audited joint SLQP adapter is not implemented\n";return 1;}
#else
#include "spin2014/SourceCageSlqp.h"
#include "spin2014/SlqpAuditTrace.h"
#include "source_surface_test_utils.h"
#include <Eigen/Geometry>
#include <limits>
#include <type_traits>
#include <cstdlib>
#include <cstring>
#include <new>
#include <climits>
// Test-only, one-shot real allocator failure. No production callback or hook.
// Sweep actual 65-byte SHA result allocations rather than assuming an ordinal.
namespace allocator_probe {
thread_local bool armed=false,triggered=false;
thread_local std::size_t remaining=0;
}
void* operator new(std::size_t n) {
    if(allocator_probe::armed&&n==65&&allocator_probe::remaining--==0) {
        allocator_probe::armed=false;allocator_probe::triggered=true;throw std::bad_alloc();
    }
    if(void* p=std::malloc(n?n:1))return p;throw std::bad_alloc();
}
void operator delete(void* p) noexcept {std::free(p);}
void operator delete(void* p,std::size_t) noexcept {std::free(p);}
void* operator new[](std::size_t n) {return ::operator new(n);}
void operator delete[](void* p) noexcept {std::free(p);}
void operator delete[](void* p,std::size_t) noexcept {std::free(p);}
using namespace spin2014;
using namespace surface_test;
namespace {
// Independently specified synthetic boxes, never original paper experiments.
SourceConformingBoundary fixture(bool changed=false) {
    auto outer=cube(.3,.7);
    auto design=changed?box({.34,.36,.38},{.66,.64,.62}):cube(.34,.66);
    SourceSurfaceOctreeConfig sc;sc.maximumDepth=5;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
    SourceAdaptiveDomain domain(source,2,5);
    SourceSurfaceCageConfig cc;cc.cageLevel=2;
    SourceCageGridConfig gc;gc.fineLevel=3;
    SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),{}),{});
    require(bbw.success(),"actual FEM/BBW succeeds");
    return SourceConformingBoundary(SourceCageBinding(domain,bbw,{}),{});
}
template<class E,class F> void rejects(F f,const char* name) {
    bool caught=false;try {f();}catch(const E& e) {
        caught=true;
        if constexpr(std::is_same_v<E,std::length_error>)require(std::string(e.what()).find("RESOURCE_LIMIT")!=std::string::npos,"resource token");
    }require(caught,name);
}
void near(double a,double b,const char* name,double tol=2e-10) {
    require(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=tol*(1+std::abs(b)),name);
}
// Wrong multipliers, omitted handle columns, changed trace states, or a forged
// convergence flag must fail this independent scalar first-order algebra.
void checkKkt(const SourceCagePaperEvaluation& e,const Eigen::VectorXd& x,std::size_t b,const KktResidual& k) {
    require(k.valid,"actual source witness valid");
    require(k.lambda.size()==e.constraints.size()&&k.alpha.size()==x.size()&&k.eta.size()==x.size(),"full multiplier dimensions");
    double eq=0,bound=0,station=0,dual=0,comp=0;
    for(int j=0;j<e.constraints.size();++j)eq=std::max(eq,std::abs(e.constraints[j]));
    for(int i=0;i<x.size();++i) {
        double r=e.jointGradient[i]-k.alpha[i]+k.eta[i];
        for(int j=0;j<e.constraints.size();++j)r+=e.jointConstraintJacobian(j,i)*k.lambda[j];
        station=std::max(station,std::abs(r));dual=std::max({dual,-k.alpha[i],-k.eta[i]});
        if(std::size_t(i)<b) {
            bound=std::max({bound,-x[i],x[i]-1});
            comp=std::max({comp,std::abs(x[i]*k.alpha[i]),std::abs((1-x[i])*k.eta[i])});
        } else require(k.alpha[i]==0&&k.eta[i]==0,"unbounded full handle duals zero");
    }
    near(k.equalityPrimal,eq,"equality primal");near(k.boundPrimal,bound,"bound primal");
    near(k.stationarity,station,"stationarity");near(k.dualSign,dual,"dual sign");near(k.complementarity,comp,"complementarity");
    near(k.full,std::max({eq,bound,station,dual,comp}),"full KKT");
}
void check(const SourceCagePaperFunctional& f,const Eigen::VectorXd& initial,const SourceCageSlqpSettings& s,const SourceCageSlqpResult& r) {
    require(r.outcome==SourceCageSlqpOutcome::Completed&&r.numericalResultAvailable&&r.finalEvaluationAvailable,"completed source evaluation");
    require(r.geometryStatus==SourceCageSlqpGeometryStatus::NotEvaluated&&r.feasibility==SourceCageSlqpFeasibility::Unknown,"numeric status grants no geometry");
    require(r.numerical.algorithm==SlqpAlgorithm::AuditedKktRestorationV1&&r.numerical.boundedVariableCount==8,"audited beta-only prefix");
    require(r.numerical.finalState.size()==initial.size()&&r.numerical.beta.size()==initial.size(),"raw complete state preserved");
    require((r.initialState-initial).norm()==0,"initial state preserved");
    const auto direct=f.evaluate(r.numerical.finalState);
    require(r.finalEvaluation.outputSha256==direct.outputSha256,"final accepted state direct source digest");
    require((r.finalEvaluation.world.deformedPoints-direct.world.deformedPoints).norm()==0,"every original/material point preserved");
    near(r.numerical.evaluation.objective,direct.totalObjective,"total objective mapped");
    near(r.numerical.evaluation.paperObjective,direct.paperObjective,"paper objective mapped");
    near(r.numerical.evaluation.regularizationEnergy,direct.fillPenalty+direct.shapeEnergy+direct.displacementEnergy,"all three penalties retained");
    require((r.numerical.evaluation.finalMoments.values()-direct.spinMoments.values()).norm()==0,"generic moments use fixed spin frame");
    require((r.numerical.evaluation.gradient-direct.jointGradient).norm()==0,"all joint gradient columns");
    require((r.numerical.evaluation.constraintJacobian-direct.jointConstraintJacobian).norm()==0,"all joint constraint columns");
    require(r.numerical.trace.size()<=std::size_t(s.solver.maximumIterations+1)&&!r.numerical.trace.empty(),"bounded complete trace");
    require(r.numerical.trace.front().historyBefore==0,"fresh history");
    bool accepted=false;
    for(const auto& t:r.numerical.trace) {
        const auto e=f.evaluate(t.state);checkKkt(e,t.state,8,t.kkt);
        require((t.gradient-e.jointGradient).norm()==0&&(t.constraints-e.constraints).norm()==0&&(t.constraintJacobian-e.jointConstraintJacobian).norm()==0,"trace direct source algebra");
        near(t.objective,e.totalObjective,"trace objective");near(t.meritValue,e.totalObjective+s.solver.meritPenalty*e.constraints.lpNorm<1>(),"trace merit");
        require(t.phi==(f.hasAngle()?t.state[8]:0.),"static scale is not diagnostic phi");
        require(t.activeMembership.size()==f.stateCount(),"complete membership");
        for(std::size_t i=8;i<f.stateCount();++i)require(t.activeMembership[i]==0,"unbounded coordinates never active");
        require(t.membershipHash==slqpMembershipHash(t.activeMembership),"membership digest");
        for(const auto& trial:t.objectiveTrials)accepted|=trial.accepted;
        for(const auto& trial:t.restoration.trials)accepted|=trial.accepted;
    }
    require(accepted,"nonstationary actual source has accepted objective/restoration step");
    require((r.numerical.finalState.head(8)-initial.head(8)).norm()>1e-10,"beta participates in accepted solve");
    require((r.numerical.finalState.tail(104)-initial.tail(104)).norm()>1e-10,"full handles participate in accepted solve");
    require(r.numerical.finalPhi==(f.hasAngle()?r.numerical.finalState[8]:0.)&&r.numerical.initialPhi==(f.hasAngle()?initial[8]:0.),"correct initial/final phi");
    require(r.numerical.traceHash==hashAuditedSlqp(r.numerical,s.solver),"corrected trace rehashed");
    checkKkt(direct,r.numerical.finalState,8,r.numerical.finalKkt);
    if(r.numerical.status==SlqpStatus::Converged) {
        require(r.numerical.finalKkt.equalityPrimal<=s.solver.equalityTolerance,"convergence equality threshold");
        require(std::max({r.numerical.finalKkt.stationarity,r.numerical.finalKkt.boundPrimal,r.numerical.finalKkt.dualSign,r.numerical.finalKkt.complementarity})<=s.solver.projectedGradientTolerance,"convergence full first-order threshold");
    }
    require(r.callbackCount+1==r.recordCount&&r.records[r.recordCount-1].purpose==SourceCageSlqpEvaluationPurpose::FinalAcceptedState,"explicit final evaluation tape");
    std::size_t calls=0;
    for(std::size_t i=0;i<r.recordCount;++i) {
        const auto& event=r.records[i];require(event.sequence==i,"deterministic invocation sequence");
        if(event.sourceInvoked)++calls;
        if(event.status==SourceCageSlqpEvaluationStatus::Evaluated) {
            const auto e=f.evaluate(r.callbackStates.col(i));
            require(std::string(event.outputSha256.data())==e.outputSha256,"every actual callback output bound");
        } else require(event.reason[0]!='\0'&&event.reasonSha256[0]!='\0',"invalid callbacks have reason and digest");
        require(event.stateSha256[0]!='\0',"all attempted states bound");
    }
    require(calls==r.sourceEvaluationCount,"actual-source call count");
}
SourceCagePaperFunctionalConfig config(Spin2014Objective objective) {
    SourceCagePaperFunctionalConfig c;c.objective=objective;c.inertiaWeight=.8;c.centerWeight=.6;
    c.fillWeight=.03;c.shapeWeight=.04;c.displacementWeight=.05;return c;
}
Eigen::VectorXd initialState(const SourceCagePaperFunctional& f) {
    Eigen::VectorXd beta(8);beta<<.17,.24,.31,.38,.45,.52,.59,.66;
    auto x=f.initialState(beta,f.hasAngle()?.23:0.);
    for(std::size_t h=0;h<26;++h) {
        x[f.handleOffset()+4*h]=1.04+.0003*h;
        x[f.handleOffset()+4*h+1]=.007+.0001*h;
        x[f.handleOffset()+4*h+2]=-.011+.00007*h;
        x[f.handleOffset()+4*h+3]=.004-.00006*h;
    }return x;
}
}
int main() {try {
    const auto boundary=fixture();const SourceCageMomentFunctional moments(boundary,1.3,{});
    // Fixed settings declared before observations, not tuned to force a paper case.
    SourceCageSlqpSettings s;s.solver.maximumIterations=6;s.solver.maximumLineSearchSteps=12;s.solver.lbfgsMemory=4;
    SourceCagePaperFrame frame;frame.origin=Eigen::Vector3d(.48,.51,.3);
    frame.rotation=Eigen::AngleAxisd(.17,Eigen::Vector3d::UnitZ()).toRotationMatrix();
    for(const auto objective:{Spin2014Objective::StaticLowering,Spin2014Objective::TopMassAndLowering,Spin2014Objective::TopInertiaOnly,Spin2014Objective::TopLoweringOnly,Spin2014Objective::YoYo}) {
        SourceCagePaperFunctional f(moments,config(objective),frame);const auto x=initialState(f);
        require(x.size()==(f.hasAngle()?113:112)&&f.handleCount()==26,"all original handles and beta dimensions");
        const auto r=solveSourceCageSlqp(f,x,s);check(f,x,s,r);
        const auto repeat=solveSourceCageSlqp(f,x,s);check(f,x,s,repeat);
        require(r.scientificSha256==repeat.scientificSha256&&r.callbackTapeSha256==repeat.callbackTapeSha256&&(r.numerical.finalState-repeat.numerical.finalState).norm()==0,"deterministic two-run state and event identity");
        auto corrupt=r;corrupt.numerical.trace.front().gradient[8]+=.1;
        rejects<std::runtime_error>([&]{check(f,x,s,corrupt);},"independent checker rejects raw trace corruption");
        corrupt=r;corrupt.numerical.status=SlqpStatus::Converged;corrupt.numerical.finalKkt.stationarity=0;
        rejects<std::runtime_error>([&]{check(f,x,s,corrupt);},"independent checker rejects false convergence");
        std::size_t invalid=0;for(std::size_t i=0;i<r.recordCount;++i)invalid+=r.records[i].status==SourceCageSlqpEvaluationStatus::InvalidNumericalTrial;
        std::cout<<"objective "<<int(objective)<<" "<<slqpStatusToken(r.numerical.status)<<" records="<<r.recordCount<<" invalid="<<invalid<<" PASS\n";
    }
    SourceCagePaperFunctional f(moments,config(Spin2014Objective::StaticLowering),frame);const auto x=initialState(f);
    auto shortSettings=s;shortSettings.solver.maximumIterations=1;
    const auto shortRun=solveSourceCageSlqp(f,x,shortSettings);
    require(shortRun.numerical.status==SlqpStatus::MaximumIterations,"honest short iteration limit");
    const auto reference=solveSourceCageSlqp(f,x,s);
    for(int change=0;change<4;++change) {
        auto c=config(Spin2014Objective::StaticLowering);auto fr=frame;auto settings=s;
        if(change==0)c.fillWeight+=.01;if(change==1)fr.origin.x()+=.01;
        if(change==2)++settings.budget.maximumSourceEvaluations;
        const SourceCageMomentFunctional changedMoments(boundary,change==3?1.4:1.3,{});
        SourceCagePaperFunctional changed(changedMoments,c,fr);
        require(solveSourceCageSlqp(changed,x,settings).scientificSha256!=reference.scientificSha256,"frame/weight/density/budget identity");
    }
    const SourceCageMomentFunctional changedMoments(fixture(true),1.3,{});
    SourceCagePaperFunctional changed(changedMoments,config(Spin2014Objective::StaticLowering),frame);
    const auto changedResult=solveSourceCageSlqp(changed,initialState(changed),s);
    require(changedResult.scientificSha256!=reference.scientificSha256&&changedResult.numerical.trace.front().historyBefore==0,"source rebinding resets history and identity");
    const auto plan=planSourceCageSlqpResources(f,s);
    // The generic trace serializer visits K*N matrix coefficients using int.
    // Reject oversized public caps without allocating an oversized source.
    auto nativeCap=s;nativeCap.budget.maximumStateValues=std::size_t(INT_MAX)/6;
    require(planSourceCageSlqpResources(f,nativeCap).stateValues==112,
        "exact conservative native serialization state cap permits small source");
    ++nativeCap.budget.maximumStateValues;
    rejects<std::length_error>([&]{planSourceCageSlqpResources(f,nativeCap);},
        "one-above native serialization state cap rejected before allocation");
    auto exact=s;exact.budget.maximumStateValues=f.stateCount();exact.budget.maximumTraceRows=plan.traceRows;
    exact.budget.maximumCallbackRecords=plan.callbackRecords;exact.budget.maximumAdapterBytes=plan.requiredBytes;
    require(solveSourceCageSlqp(f,x,exact).finalEvaluationAvailable,"exact accounting caps permit solve");
    for(int cap=0;cap<5;++cap) {
        auto bad=exact;
        if(cap==0)--bad.budget.maximumStateValues;if(cap==1)--bad.budget.maximumTraceRows;
        if(cap==2)--bad.budget.maximumCallbackRecords;if(cap==3)--bad.budget.maximumAdapterBytes;
        if(cap==4)bad.budget.maximumSourceEvaluations=0;
        rejects<std::length_error>([&]{solveSourceCageSlqp(f,x,bad);},"one-below/zero resource caps");
    }
    auto exhausted=s;exhausted.budget.maximumSourceEvaluations=1;
    const auto limited=solveSourceCageSlqp(f,x,exhausted);
    require(limited.outcome==SourceCageSlqpOutcome::ResourceLimit&&!limited.finalEvaluationAvailable&&limited.sourceEvaluationCount==1,"resource latch overrides engine nonprogress");
    require(!limited.numerical.trace.empty()&&(limited.numerical.finalState-x).norm()==0,"resource failure keeps last numerical state/trace");
    bool latched=false;
    for(std::size_t i=0;i<limited.recordCount;++i) {
        if(limited.records[i].status==SourceCageSlqpEvaluationStatus::ResourceLimit)latched=true;
        if(latched)require(!limited.records[i].sourceInvoked,"no source work after callback cap");
    }
    require(latched&&limited.callbackCount>1,"actual swallowed trial exhaustion exercised");
    auto finalCap=s;finalCap.budget.maximumSourceEvaluations=reference.sourceEvaluationCount;
    require(solveSourceCageSlqp(f,x,finalCap).finalEvaluationAvailable,"exact actual-source call cap includes final");
    --finalCap.budget.maximumSourceEvaluations;
    const auto finalLimited=solveSourceCageSlqp(f,x,finalCap);
    require(finalLimited.outcome==SourceCageSlqpOutcome::ResourceLimit&&!finalLimited.finalEvaluationAvailable&&
        (finalLimited.numerical.finalState-reference.numerical.finalState).norm()==0,"one-below final cap preserves accepted numerical state");
    // Independently symmetric stationary input checks the narrow convergence
    // interpretation; the preceding five nonstationary runs exercise motion.
    SourceCagePaperFrame centered;centered.origin=Eigen::Vector3d::Constant(.5);
    SourceCagePaperFunctional stationary(moments,config(Spin2014Objective::StaticLowering),centered);
    const auto rest=stationary.initialState(Eigen::VectorXd::Constant(8,.3),0.);
    const auto converged=solveSourceCageSlqp(stationary,rest,s);
    require(converged.numerical.status==SlqpStatus::Converged&&converged.geometryStatus==SourceCageSlqpGeometryStatus::NotEvaluated&&
        converged.feasibility==SourceCageSlqpFeasibility::Unknown,"first-order convergence still gives no geometry/feasibility grant");
    const auto convergedLimited=solveSourceCageSlqp(stationary,rest,exhausted);
    require(convergedLimited.numerical.status==SlqpStatus::Converged&&convergedLimited.outcome==SourceCageSlqpOutcome::ResourceLimit&&
        !convergedLimited.finalEvaluationAvailable,"resource final failure cannot promote underlying convergence");
    auto lowDerivative=config(Spin2014Objective::StaticLowering);lowDerivative.maximumJointDerivativeValues=1;
    SourceCagePaperFunctional cannotStart(moments,lowDerivative,frame);
    const auto initialFailure=solveSourceCageSlqp(cannotStart,x,s);
    require(initialFailure.outcome==SourceCageSlqpOutcome::ResourceLimit&&initialFailure.recordCount==1&&
        initialFailure.sourceEvaluationCount==1&&initialFailure.records[0].sourceInvoked&&
        initialFailure.records[0].status==SourceCageSlqpEvaluationStatus::ResourceLimit&&
        initialFailure.records[0].reasonSha256[0]!='\0'&&!initialFailure.finalEvaluationAvailable,
        "initial actual-source resource failure retains attempted-state evidence");
    require(!initialFailure.numericalResultAvailable&&initialFailure.numerical.trace.empty()&&
        initialFailure.numerical.finalState.size()==0&&(initialFailure.initialState-x).norm()==0,
        "first callback failure invents no generic state/trace/KKT result");
    require(solveSourceCageSlqp(cannotStart,x,s).scientificSha256==initialFailure.scientificSha256,
        "initial resource failure deterministic provenance");
    auto huge=x;for(std::size_t h=0;h<26;++h)huge[f.handleOffset()+4*h]=1e200;
    const auto overflow=solveSourceCageSlqp(f,huge,s);
    require(overflow.outcome==SourceCageSlqpOutcome::InitialEvaluationFailure&&!overflow.numericalResultAvailable&&
        overflow.recordCount==1&&overflow.sourceEvaluationCount==1&&overflow.records[0].sourceInvoked&&
        overflow.records[0].status==SourceCageSlqpEvaluationStatus::InvalidNumericalTrial&&
        overflow.records[0].reasonSha256[0]!='\0'&&!overflow.finalEvaluationAvailable,
        "finite positive initial scale overflow retains actual numerical failure record");
    std::size_t injected=0,returnedResource=0;
    auto probeSettings=s;probeSettings.solver.maximumIterations=1;
    for(std::size_t ordinal=0;ordinal<40;++ordinal) {
        allocator_probe::remaining=ordinal;allocator_probe::triggered=false;allocator_probe::armed=true;
        SourceCageSlqpResult probe;bool explicitResourceException=false;
        try {probe=solveSourceCageSlqp(f,x,probeSettings);}
        catch(const std::length_error& e) {
            allocator_probe::armed=false;explicitResourceException=true;
            require(std::string(e.what()).find("RESOURCE_LIMIT")!=std::string::npos,"allocator exception category");
        }catch(...) {allocator_probe::armed=false;throw;}
        allocator_probe::armed=false;
        if(!allocator_probe::triggered)break;
        ++injected;
        if(explicitResourceException)continue;
        require(probe.outcome==SourceCageSlqpOutcome::ResourceLimit,"real allocation failure cannot return success");
        ++returnedResource;
        for(std::size_t i=0;i<probe.recordCount;++i) {
            const auto& event=probe.records[i];
            if(event.status==SourceCageSlqpEvaluationStatus::Unused)
                std::cerr<<"allocator ordinal="<<ordinal<<" event="<<i<<" remained Unused\n";
            require(event.status!=SourceCageSlqpEvaluationStatus::Unused,"allocator recovery records every attempted event");
            require(std::strlen(event.stateSha256.data())==64,"allocator recovery completes state digest");
            if(event.status!=SourceCageSlqpEvaluationStatus::Evaluated)
                require(event.reason[0]!='\0'&&std::strlen(event.reasonSha256.data())==64,"allocator recovery completes reason digest");
        }
    }
    require(injected>=3&&returnedResource>0,"real allocator calibration reached recoverable callback sites");
    std::cout<<"one-shot allocator faults="<<injected<<" returned resource evidence="<<returnedResource<<" PASS\n";
    for(int cap=0;cap<5;++cap) {
        auto bad=s;
        if(cap==0)bad.budget.maximumStateValues=0;if(cap==1)bad.budget.maximumTraceRows=0;
        if(cap==2)bad.budget.maximumCallbackRecords=0;if(cap==3)bad.budget.maximumAdapterBytes=0;if(cap==4)bad.budget.maximumSourceEvaluations=0;
        rejects<std::length_error>([&]{solveSourceCageSlqp(f,x,bad);},"all zero caps rejected");
    }
    for(int cap=0;cap<5;++cap) {
        auto bad=s;
        auto value=std::numeric_limits<std::size_t>::max();
        if(cap==0)bad.budget.maximumStateValues=value;if(cap==1)bad.budget.maximumTraceRows=value;
        if(cap==2)bad.budget.maximumCallbackRecords=value;if(cap==3)bad.budget.maximumAdapterBytes=value;if(cap==4)bad.budget.maximumSourceEvaluations=value;
        rejects<std::length_error>([&]{solveSourceCageSlqp(f,x,bad);},"overflow cap rejected");
    }
    for(int kind=0;kind<5;++kind) {
        auto bad=x;
        if(kind==0)bad.conservativeResize(x.size()-1);if(kind==1)bad[0]=-.1;if(kind==2)bad[0]=1.1;
        if(kind==3)bad[0]=std::numeric_limits<double>::quiet_NaN();if(kind==4)bad[8]=0.;
        rejects<std::invalid_argument>([&]{solveSourceCageSlqp(f,bad,s);},"bad full state rejected");
    }
    for(int kind=0;kind<6;++kind) {
        auto bad=s;if(kind==0)bad.solver.algorithm=SlqpAlgorithm::LegacyV3;
        if(kind==1)bad.solver.algorithm=static_cast<SlqpAlgorithm>(27);if(kind==2)bad.solver.maximumIterations=0;
        if(kind==3)bad.solver.armijo=std::numeric_limits<double>::quiet_NaN();if(kind==4)bad.solver.activeSetTolerance=.5;
        if(kind==5)bad.solver.audit.lpTimeLimit=0.;
        rejects<std::invalid_argument>([&]{solveSourceCageSlqp(f,x,bad);},"invalid audited settings rejected");
    }
    std::cout<<"actual source all-mode SLQP, KKT, provenance, caps and resource latch PASS\n";return 0;
}catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}}
#endif
