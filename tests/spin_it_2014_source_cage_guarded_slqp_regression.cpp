#include <iostream>
#if !__has_include("spin2014/SourceCageGuardedSlqp.h")
int main() {std::cerr<<"FAIL: actual-source geometry gate cannot reject an invalid initial state before SLQP evaluation\n";return 1;}
#else
#include "spin2014/SourceCageGuardedSlqp.h"
#include "source_surface_test_utils.h"
#include <type_traits>
#include <iomanip>
#include <cstdlib>
#include <new>
#include <limits>
namespace allocation_probe {
thread_local bool armed=false,triggered=false;
thread_local std::size_t remaining=0;
}
void* operator new(std::size_t n) {
    if(allocation_probe::armed&&n==65&&allocation_probe::remaining--==0) {
        allocation_probe::armed=false;allocation_probe::triggered=true;throw std::bad_alloc();
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
SourceConformingBoundary fixture() {
    const auto outer=cube(.3125,.6875),design=cube(.34,.66);
    SourceSurfaceOctreeConfig sc;sc.maximumDepth=3;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
    SourceAdaptiveDomain domain(source,2,3);
    SourceSurfaceCageConfig cc;cc.cageLevel=2;SourceCageGridConfig gc;gc.fineLevel=3;
    SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),{}),{});
    require(bbw.success(),"actual full26-handle BBW succeeds");
    return SourceConformingBoundary(SourceCageBinding(domain,bbw,{}),{});
}
template<class F> void resourceRejected(F f) {
    bool rejected=false;try {f();}catch(const std::length_error& e){rejected=std::string(e.what()).find("RESOURCE_LIMIT")!=std::string::npos;}
    require(rejected,"invalid resource cap rejected before solve");
}
SourceCagePaperFunctionalConfig config(Spin2014Objective objective) {
    SourceCagePaperFunctionalConfig c;c.objective=objective;c.inertiaWeight=.8;c.centerWeight=.6;
    c.fillWeight=.03;c.shapeWeight=.04;c.displacementWeight=.05;return c;
}
Eigen::VectorXd fullInitial(const SourceCagePaperFunctional& f) {
    Eigen::VectorXd b(8);b<<.17,.24,.31,.38,.45,.52,.59,.66;
    return f.initialState(b,f.hasAngle()?.23:0.);
}
void checkLinks(const SourceCagePaperFunctional& f,const SourceCageGuardedSlqpFacts& a) {
    const auto& r=a.source;
    require(a.attempts.size()==r.recordCount,"every attempted callback has a guard disposition, including bookkeeping resource failure");
    std::size_t checker=0,values=0,derivatives=0,stages=0,callback=1;
    std::size_t exactRelations=0,cacheHits=0,cacheMisses=0,booleanEvaluations=0;
    bool resource=false;
    for(std::size_t i=0;i<r.recordCount;++i) {
        const auto& attempt=a.attempts[i];const auto& event=r.records[i];
        require(attempt.callbackRecord==i&&event.sequence==i&&attempt.firstStage==stages,"ordered raw state/guard linkage");
        stages+=attempt.stageCount;require(stages<=a.stages.size(),"bounded complete stage span");
        if(resource)require(!attempt.checkerInvoked&&!event.sourceInvoked,"resource latch prevents geometry and derivative source work");
        checker+=attempt.checkerInvoked;values+=attempt.sourceValuesInvoked;derivatives+=event.sourceInvoked;
        if(attempt.sourceValuesInvoked) {
            require(attempt.exactPairCacheHits+attempt.exactPairCacheMisses==attempt.exactPairRelations,
                    "each completed guarded geometry call accounts for every exact pair result");
            require(attempt.exactPairBooleanEvaluations==attempt.exactPairCacheMisses,
                    "guard attempt distinguishes exact cache hits from Boolean work");
            exactRelations+=attempt.exactPairRelations;cacheHits+=attempt.exactPairCacheHits;
            cacheMisses+=attempt.exactPairCacheMisses;booleanEvaluations+=attempt.exactPairBooleanEvaluations;
        }
        if(event.status==SourceCageSlqpEvaluationStatus::ResourceLimit)resource=true;
        if(event.sourceInvoked)require(attempt.checkCompleted&&attempt.outcome==SourceCageGeometryOutcome::ChecksPassed,"actual completed geometry gate before every derivative call");
        if(attempt.sourceOutputSha256[0])require(f.evaluateValues(r.callbackStates.col(i)).outputSha256==attempt.sourceOutputSha256.data(),"guard source hash belongs to exact raw state");
        if(event.status==SourceCageSlqpEvaluationStatus::Evaluated)require(f.evaluate(r.callbackStates.col(i)).outputSha256==event.outputSha256.data(),"full derivative hash belongs to exact raw state");
    }
    require(stages==a.stages.size()&&checker==a.checkerInvocationCount&&values==a.geometryValueEvaluationCount&&derivatives==r.sourceEvaluationCount,"separate actual geometry and full evaluation counts");
    require(exactRelations==a.exactPairRelations&&cacheHits==a.exactPairCacheHits&&
            cacheMisses==a.exactPairCacheMisses&&booleanEvaluations==a.exactPairBooleanEvaluations,
            "guarded audit aggregates cache/work diagnostics without losing per-attempt counts");
    for(const auto& row:r.numerical.trace) {
        bool admitted=false;
        for(std::size_t i=0;i<r.recordCount;++i)if(r.records[i].status==SourceCageSlqpEvaluationStatus::Evaluated&&
            (r.callbackStates.col(i)-row.state).norm()==0)admitted=true;
        require(admitted,"each accepted solver state links to passing geometry and source output");
        auto trials=[&](const std::vector<RestorationTrial>& list,bool restoration) {
            for(const auto& trial:list) {
                require(callback<r.callbackCount,"solver trial has one corresponding actual callback");
                const auto& event=r.records[callback];const auto& attempt=a.attempts[callback];
                require(trial.valid==(event.status==SourceCageSlqpEvaluationStatus::Evaluated),"solver invalid trial retains actual rejection cause");
                if(trial.accepted)require(attempt.outcome==SourceCageGeometryOutcome::ChecksPassed,"accepted trial completed gate");
                if(restoration)require((r.callbackStates.col(callback)-row.state-trial.step*row.restoration.direction).norm()<1e-12,"actual restoration trial state from recorded direction and step");
                ++callback;
            }
        };
        trials(row.objectiveTrials,false);trials(row.restoration.trials,true);
    }
    if(r.numericalResultAvailable)require(callback==r.callbackCount,"no dropped or extra solver trials");
    if(a.finalStatus==SourceCageGuardedFinalStatus::ChecksPassed) {
        require(r.finalEvaluationAvailable&&a.finalRecheckCount==1&&a.attempts.back().checkCompleted,"fresh final checker and material output");
        require(r.records[r.recordCount-1].purpose==SourceCageSlqpEvaluationPurpose::FinalAcceptedState,"explicit final callback purpose");
        require((r.callbackStates.col(r.recordCount-1)-r.numerical.finalState).norm()==0,"final recheck uses final raw accepted state");
        const auto direct=f.evaluate(r.numerical.finalState);
        require(direct.outputSha256==r.finalEvaluation.outputSha256&&
            (direct.spinMoments.values()-r.numerical.evaluation.finalMoments.values()).norm()==0,"fresh material moments and source output");
    }
    require(a.feasibility==SourceCageSlqpFeasibility::Unknown&&r.geometryStatus==SourceCageSlqpGeometryStatus::NotEvaluated,"private guarded facts grant no feasibility or legacy promotion");
}
}
int main(int argc,char** argv) {try {
    static_assert(!std::is_default_constructible_v<SourceCageGuardedSlqpAudit>);
    static_assert(!std::is_copy_constructible_v<SourceCageGuardedSlqpAudit>);
    static_assert(std::is_const_v<std::remove_reference_t<decltype(std::declval<SourceCageGuardedSlqpAudit>().facts())>>);
    const SourceCageMomentFunctional moments(fixture(),1.3,{});
    SourceCagePaperFunctionalConfig c;c.objective=Spin2014Objective::StaticLowering;
    const SourceCagePaperFunctional f(moments,c);
    const auto x=f.initialState(Eigen::VectorXd::Constant(8,.25),0.);
    if(argc>1&&std::string(argv[1])=="--move") {
        SourceCageGuardedSlqpSettings settings;settings.source.solver.maximumIterations=1;settings.source.solver.maximumLineSearchSteps=2;
        auto original=solveSourceCageGuardedSlqp(f,x,settings);
        const auto& retained=original.facts();
        const auto digest=original.facts().scientificSha256;
        const auto moved=std::move(original);
        require(moved.facts().scientificSha256==digest,"move preserves trustworthy evidence");
        require(retained.finalStatus==SourceCageGuardedFinalStatus::Unavailable&&!retained.source.finalEvaluationAvailable&&
            !retained.source.numericalResultAvailable,"pre-move const reference cannot retain stale passing scalar or material availability");
        bool rejected=false;try {(void)original.facts();}catch(const std::logic_error&){rejected=true;}
        require(rejected,"moved-from audit cannot expose stale passing scalar with emptied evidence");
        std::cout<<"move-only audit validity PASS\n";return 0;
    }
    if(argc>1&&std::string(argv[1])=="--probe-exhaustion") {
        SourceCagePaperFrame frame;frame.origin=Eigen::Vector3d(.48,.51,.3);
        const SourceCagePaperFunctional actual(moments,config(Spin2014Objective::TopMassAndLowering),frame);
        for(int lines:{1,2,3,4}) {
            SourceCageGuardedSlqpSettings settings;settings.source.solver.maximumIterations=1;
            settings.source.solver.maximumLineSearchSteps=lines;settings.geometry.minimumOuterDomainClearance=.04;
            const auto run=solveSourceCageGuardedSlqp(actual,fullInitial(actual),settings);const auto& a=run.facts();checkLinks(actual,a);
            std::cout<<"EXHAUSTION lines="<<lines<<" status="<<slqpStatusToken(a.source.numerical.status)
                <<" displacement="<<(a.source.numerical.finalState-a.source.initialState).norm()<<" digest="<<a.scientificSha256<<"\n";
            for(const auto& attempt:a.attempts)std::cout<<" attempt="<<attempt.callbackRecord<<" outcome="<<sourceCageGeometryOutcomeToken(attempt.outcome)<<"\n";
            for(const auto& row:a.source.numerical.trace) {
                for(const auto& trial:row.objectiveTrials)std::cout<<" objective step="<<trial.step<<" valid="<<trial.valid<<" accepted="<<trial.accepted<<"\n";
                for(const auto& trial:row.restoration.trials)std::cout<<" restoration step="<<trial.step<<" valid="<<trial.valid<<" accepted="<<trial.accepted<<"\n";
            }std::cout.flush();
        }return 0;
    }
    if(argc>1&&std::string(argv[1])=="--allocation") {
        SourceCageGuardedSlqpSettings settings;settings.source.solver.maximumIterations=1;settings.source.solver.maximumLineSearchSteps=2;
        int injected=0,retained=0;
        for(std::size_t ordinal=0;ordinal<48;++ordinal) {
            allocation_probe::remaining=ordinal;allocation_probe::triggered=false;allocation_probe::armed=true;
            try {
                const auto run=solveSourceCageGuardedSlqp(f,x,settings);allocation_probe::armed=false;
                require(allocation_probe::triggered,"bounded allocator probe reached actual SHA allocation");++injected;++retained;
                const auto& a=run.facts();require(a.source.outcome==SourceCageSlqpOutcome::ResourceLimit&&a.finalStatus==SourceCageGuardedFinalStatus::ResourceLimit,"allocation failure cannot yield guarded pass");
                checkLinks(f,a);
                std::cout<<"allocation ordinal="<<ordinal<<" retained resource records="<<a.source.recordCount<<"\n";
            }catch(const std::length_error&) {allocation_probe::armed=false;require(allocation_probe::triggered,"propagated resource from real allocation");++injected;}
            catch(const std::bad_alloc&) {allocation_probe::armed=false;require(allocation_probe::triggered,"propagated actual allocation");++injected;}
        }
        require(injected==48&&retained>0,"allocation probe covers recoverable and propagated resource failures");
        std::cout<<"allocator resource gate PASS injected="<<injected<<" retained="<<retained<<"\n";return 0;
    }
    if(argc>1&&std::string(argv[1])=="--probe") {
        for(auto objective:{Spin2014Objective::StaticLowering,Spin2014Objective::TopMassAndLowering,
            Spin2014Objective::TopInertiaOnly,Spin2014Objective::TopLoweringOnly,Spin2014Objective::YoYo})
        for(double clearance:{0.,.005,.02,.031,.04}) {
            SourceCagePaperFunctionalConfig cfg;cfg.objective=objective;cfg.inertiaWeight=.8;cfg.centerWeight=.6;
            cfg.fillWeight=.03;cfg.shapeWeight=.04;cfg.displacementWeight=.05;
            SourceCagePaperFrame frame;frame.origin=Eigen::Vector3d(.48,.51,.3);
            const SourceCagePaperFunctional actual(moments,cfg,frame);
            Eigen::VectorXd beta(8);beta<<.17,.24,.31,.38,.45,.52,.59,.66;
            const auto initial=actual.initialState(beta,actual.hasAngle()?.23:0.);
            SourceCageGuardedSlqpSettings settings;settings.source.solver.maximumIterations=1;
            settings.source.solver.maximumLineSearchSteps=8;settings.geometry.minimumOuterDomainClearance=clearance;
            const auto run=solveSourceCageGuardedSlqp(actual,initial,settings);const auto& a=run.facts();
            std::cout<<std::setprecision(17)<<"PROBE objective="<<int(objective)<<" clearance="<<clearance
                <<" available="<<a.source.numericalResultAvailable<<" status="<<slqpStatusToken(a.source.numerical.status)
                <<" digest="<<a.scientificSha256<<"\n";
            for(const auto& attempt:a.attempts) {
                std::cout<<" attempt="<<attempt.callbackRecord<<" outcome="<<sourceCageGeometryOutcomeToken(attempt.outcome)
                    <<" derivatives="<<a.source.records[attempt.callbackRecord].sourceInvoked;
                if(attempt.stageCount)std::cout<<" stage="<<int(a.stages[attempt.firstStage+attempt.stageCount-1].stage)
                    <<" reason="<<a.stages[attempt.firstStage+attempt.stageCount-1].reason.data();
                std::cout<<"\n";
            }
            for(const auto& t:a.source.numerical.trace) {
                for(const auto& trial:t.objectiveTrials)std::cout<<" objective step="<<trial.step<<" valid="<<trial.valid<<" accepted="<<trial.accepted<<"\n";
                for(const auto& trial:t.restoration.trials)std::cout<<" restoration step="<<trial.step<<" valid="<<trial.valid<<" accepted="<<trial.accepted<<"\n";
            }
            std::cout.flush();
        }return 0;
    }
    SourceCageGuardedSlqpSettings s;s.source.solver.maximumIterations=1;s.source.solver.maximumLineSearchSteps=4;
    s.geometry.minimumOuterDomainClearance=.2;
    const auto audit=solveSourceCageGuardedSlqp(f,x,s);const auto& a=audit.facts();
    require(a.source.outcome==SourceCageSlqpOutcome::InitialEvaluationFailure,"initial rejected geometry is distinguished");
    require(!a.source.numericalResultAvailable&&a.source.numerical.trace.empty()&&a.source.numerical.finalState.size()==0,"no fabricated initial numerical result or KKT");
    require(a.source.sourceEvaluationCount==0&&a.geometryValueEvaluationCount==1,"rejected geometry never reaches full derivatives");
    require(a.attempts.size()==1&&a.attempts[0].outcome==SourceCageGeometryOutcome::RejectedGeometry,"actual initial geometry attempt retained");
    require((a.source.callbackStates.col(0)-x).norm()==0,"attempted full state retained");
    require(a.finalStatus==SourceCageGuardedFinalStatus::Unavailable,"no manufactured final recheck");
    checkLinks(f,a);
    auto base=s;base.geometry.minimumOuterDomainClearance=0.;base.source.solver.maximumLineSearchSteps=8;
    SourceCagePaperFrame frame;frame.origin=Eigen::Vector3d(.48,.51,.3);
    for(auto objective:{Spin2014Objective::StaticLowering,Spin2014Objective::TopMassAndLowering,
        Spin2014Objective::TopInertiaOnly,Spin2014Objective::TopLoweringOnly,Spin2014Objective::YoYo}) {
        const SourceCagePaperFunctional actual(moments,config(objective),frame);const auto initial=fullInitial(actual);
        require(initial.size()==(actual.hasAngle()?113:112)&&actual.handleCount()==26,"every actual beta/phi/full handle variable retained");
        const auto run=solveSourceCageGuardedSlqp(actual,initial,base);const auto& r=run.facts();checkLinks(actual,r);
        const auto repeat=solveSourceCageGuardedSlqp(actual,initial,base);checkLinks(actual,repeat.facts());
        require(r.finalStatus==SourceCageGuardedFinalStatus::ChecksPassed&&r.scientificSha256==repeat.facts().scientificSha256,"two identical all-mode runs with fresh final check");
        require((r.source.numerical.finalState-initial).norm()>1e-10,"actual solver moves full state");
        if(objective==Spin2014Objective::TopMassAndLowering) {
            const auto& row=r.source.numerical.trace.front();require(row.objectiveTrials.size()>=5,"actual backtracking retained");
            require(r.attempts[3].outcome==SourceCageGeometryOutcome::RejectedGeometry&&r.attempts[4].outcome==SourceCageGeometryOutcome::RejectedGeometry,"subset rejection distinct from earlier UNKNOWN source trials");
            require(row.objectiveTrials[4].accepted&&row.objectiveTrials[4].step<row.objectiveTrials[3].step,"geometry rejection shrinks to checked accepted state");
            Eigen::VectorXd direction=(r.source.callbackStates.col(1)-initial)/row.objectiveTrials.front().step;
            for(std::size_t j=0;j<row.objectiveTrials.size();++j) {
                Eigen::VectorXd expected=initial+row.objectiveTrials[j].step*direction;
                expected.head(8)=expected.head(8).cwiseMax(0.).cwiseMin(1.);
                require((expected-r.source.callbackStates.col(j+1)).norm()<1e-12,"full beta/phi/handle trial state equals actual solver direction times reduced step");
            }
        }
        std::cout<<"objective="<<int(objective)<<" final="<<slqpStatusToken(r.source.numerical.status)<<" guard_calls="<<r.checkerInvocationCount<<" values="<<r.geometryValueEvaluationCount<<" derivatives="<<r.source.sourceEvaluationCount<<" digest="<<r.scientificSha256<<" PASS\n";
    }
    const SourceCagePaperFunctional top(moments,config(Spin2014Objective::TopMassAndLowering),frame);
    const auto topInitial=fullInitial(top);
    auto rollbackSettings=base;rollbackSettings.geometry.minimumOuterDomainClearance=.04;rollbackSettings.source.solver.maximumLineSearchSteps=4;
    const auto rollback=solveSourceCageGuardedSlqp(top,topInitial,rollbackSettings);checkLinks(top,rollback.facts());
    const auto& rr=rollback.facts().source;
    require(rr.numerical.status==SlqpStatus::RestorationStalled&&(rr.numerical.finalState-topInitial).norm()==0,"actual exhausted geometry trials preserve last accepted state");
    const auto direct=top.evaluate(topInitial);
    require((rr.numerical.evaluation.finalMoments.values()-direct.spinMoments.values()).norm()==0&&rr.numerical.evaluation.objective==direct.totalObjective,"exhaustion preserves accepted moments and objective");
    require(rr.numerical.trace.front().restoration.attempted&&!rr.numerical.trace.front().restoration.accepted,"failed restoration traverses same actual gate");
    for(const auto& row:rr.numerical.trace)require(row.historyBefore==0&&row.historyAfter==0&&!row.curvaturePairAccepted,"failed trials never add curvature history");
    const SourceCagePaperFunctional yoyo(moments,config(Spin2014Objective::YoYo),frame);
    auto restorationSettings=base;restorationSettings.geometry.minimumOuterDomainClearance=.04;
    const auto restored=solveSourceCageGuardedSlqp(yoyo,fullInitial(yoyo),restorationSettings);checkLinks(yoyo,restored.facts());
    const auto& restoration=restored.facts().source.numerical.trace.front().restoration;
    require(restoration.accepted&&restoration.trials.size()>1&&restoration.trials.back().step<restoration.trials.front().step,"actual restoration shrinks rejected candidates to guarded acceptance");
    auto guardCap=base;guardCap.budget.maximumGuardCalls=1;
    const auto limited=solveSourceCageGuardedSlqp(top,topInitial,guardCap);checkLinks(top,limited.facts());
    require(limited.facts().source.outcome==SourceCageSlqpOutcome::ResourceLimit&&limited.facts().checkerInvocationCount==1&&limited.facts().source.sourceEvaluationCount==1,"permanent guard budget resource latch");
    require(limited.facts().source.callbackCount>2&&limited.facts().finalRecheckCount==0,"generic later callbacks recorded without checker or source work");
    auto sourceCap=base;sourceCap.source.budget.maximumSourceEvaluations=1;
    const auto sourceLimited=solveSourceCageGuardedSlqp(top,topInitial,sourceCap);checkLinks(top,sourceLimited.facts());
    require(sourceLimited.facts().source.outcome==SourceCageSlqpOutcome::ResourceLimit&&sourceLimited.facts().source.sourceEvaluationCount==1,"derivative budget resource latch");
    auto preflight=base;preflight.geometry.maximumPoints=1;
    const auto early=solveSourceCageGuardedSlqp(top,topInitial,preflight);checkLinks(top,early.facts());
    require(early.facts().checkerInvocationCount==1&&early.facts().geometryValueEvaluationCount==0&&!early.facts().source.numericalResultAvailable,"checker preflight failure does not count source values or invent numerical result");
    auto native=base;native.geometry.maximumNefComplexity=1;
    const auto nef=solveSourceCageGuardedSlqp(top,topInitial,native);checkLinks(top,nef.facts());
    require(nef.facts().source.outcome==SourceCageSlqpOutcome::ResourceLimit&&nef.facts().geometryValueEvaluationCount==1&&nef.facts().source.sourceEvaluationCount==0,"actual native geometry resource failure remains distinct");
    auto cfg=config(Spin2014Objective::TopMassAndLowering);cfg.maximumJointDerivativeValues=1;
    const SourceCagePaperFunctional lowDerivative(moments,cfg,frame);
    const auto derivative=solveSourceCageGuardedSlqp(lowDerivative,topInitial,base);checkLinks(lowDerivative,derivative.facts());
    require(derivative.facts().attempts.front().outcome==SourceCageGeometryOutcome::ChecksPassed&&derivative.facts().source.outcome==SourceCageSlqpOutcome::ResourceLimit&&!derivative.facts().source.numericalResultAvailable,"initial full-source derivative failure follows passing geometry without numerical manufacture");
    auto huge=topInitial;for(std::size_t h=0;h<26;++h)huge[top.handleOffset()+4*h]=1e200;
    const auto unknown=solveSourceCageGuardedSlqp(top,huge,base);checkLinks(top,unknown.facts());
    require(unknown.facts().source.outcome==SourceCageSlqpOutcome::InitialEvaluationFailure&&unknown.facts().attempts.front().outcome==SourceCageGeometryOutcome::Unknown&&unknown.facts().source.sourceEvaluationCount==0,"initial actual source value failure UNKNOWN distinguished from rejected geometry");
    const auto reference=solveSourceCageGuardedSlqp(f,x,base);checkLinks(f,reference.facts());
    std::cout<<"guarded exact-pair work relations="<<reference.facts().exactPairRelations
        <<" cache_hits="<<reference.facts().exactPairCacheHits
        <<" cache_misses="<<reference.facts().exactPairCacheMisses
        <<" boolean_evaluations="<<reference.facts().exactPairBooleanEvaluations<<'\n';
    auto exact=base;exact.budget.maximumGuardCalls=reference.facts().checkerInvocationCount;
    auto plan=planSourceCageGuardedSlqpResources(f,exact);
    exact.budget.maximumRetainedStageRecords=plan.stageCapacity;exact.budget.maximumGuardedBytes=plan.requiredBytes;
    require(solveSourceCageGuardedSlqp(f,x,exact).facts().finalStatus==SourceCageGuardedFinalStatus::ChecksPassed,"exact simultaneous bytes/stages/guard calls includes fresh final check");
    auto below=exact;--below.budget.maximumGuardedBytes;resourceRejected([&]{planSourceCageGuardedSlqpResources(f,below);});
    below=exact;--below.budget.maximumRetainedStageRecords;resourceRejected([&]{planSourceCageGuardedSlqpResources(f,below);});
    below=exact;--below.budget.maximumGuardCalls;
    const auto finalLimited=solveSourceCageGuardedSlqp(f,x,below);checkLinks(f,finalLimited.facts());
    require(finalLimited.facts().finalStatus==SourceCageGuardedFinalStatus::ResourceLimit&&!finalLimited.facts().source.finalEvaluationAvailable&&
        (finalLimited.facts().source.numerical.finalState-reference.facts().source.numerical.finalState).norm()==0,"one-below final guard cap preserves accepted state but never passes final check");
    auto finalSource=base;finalSource.source.budget.maximumSourceEvaluations=reference.facts().source.sourceEvaluationCount-1;
    const auto finalSourceLimited=solveSourceCageGuardedSlqp(f,x,finalSource);
    require(finalSourceLimited.facts().finalRecheckCount==1&&finalSourceLimited.facts().finalStatus==SourceCageGuardedFinalStatus::ResourceLimit&&!finalSourceLimited.facts().source.finalEvaluationAvailable,"final material recomputation resource failure cannot promote a passed final geometry check");
    for(auto member:{&SourceCageGuardedSlqpBudget::maximumGuardCalls,&SourceCageGuardedSlqpBudget::maximumRetainedStageRecords,&SourceCageGuardedSlqpBudget::maximumGuardedBytes})
    for(auto value:{std::size_t(0),std::numeric_limits<std::size_t>::max()}) {
        auto bad=base;bad.budget.*member=value;resourceRejected([&]{planSourceCageGuardedSlqpResources(f,bad);});
    }
    std::cout<<"actual-source guarded SLQP PASS\n";return 0;
}catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;}}
#endif
