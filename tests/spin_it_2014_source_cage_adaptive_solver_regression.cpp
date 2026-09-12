#include <iostream>
#if !__has_include("spin2014/SourceCageAdaptiveSolver.h")
int main() {std::cerr<<"FAIL: actual guarded adaptive solver cannot retain solved and pending endpoints across fresh passes\n";return 1;}
#else
#include "spin2014/SourceCageAdaptiveSolver.h"
#include "source_surface_test_utils.h"
#include <Eigen/Eigenvalues>
#include <cstring>
#include <limits>
#include <type_traits>
using namespace spin2014;
using namespace surface_test;
namespace {
bool identical(const Eigen::VectorXd& a,const Eigen::VectorXd& b) {
    return a.size()==b.size()&&std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0;
}
SourceCagePaperFunctional fixture(std::size_t leafCap=1000) {
    const auto outer=cube(.3125,.6875),design=cube(.36,.515);
    SourceSurfaceOctreeConfig sc;sc.maximumDepth=4;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
    SourceAdaptiveDomain domain(source,3,4,leafCap);
    SourceSurfaceCageConfig cc;cc.cageLevel=2;SourceCageGridConfig gc;gc.fineLevel=3;
    SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),{}),{});
    require(bbw.success(),"actual full BBW succeeds");
    SourceConformingBoundaryConfig bc;bc.maximumPatches=9000;
    SourceCageMomentFunctionalConfig mc;mc.maximumWeightValues=9000;
    SourceCagePaperFunctionalConfig pc;pc.objective=Spin2014Objective::TopMassAndLowering;
    pc.inertiaWeight=0.;pc.centerWeight=0.;
    return SourceCagePaperFunctional(SourceCageMomentFunctional(
        SourceConformingBoundary(SourceCageBinding(domain,bbw,{}),bc),1.3,mc),pc);
}
// KNOWN_SOLUTION_VALIDATION only: zero objective weights and principal frame
// make a stationary feasible relaxed state; this does not test optimization.
SourceCagePaperFunctional stationary(const SourceCagePaperFunctional& f,const Eigen::VectorXd& x) {
    const auto m=f.evaluate(x).world.finalMoments.values();
    SourceCagePaperFrame frame;frame.origin=m.segment<3>(1)/m[0];
    Eigen::Matrix3d second;second<<m[4],m[9],m[8],m[9],m[5],m[7],m[8],m[7],m[6];
    const Eigen::Matrix3d covariance=second-m[0]*frame.origin*frame.origin.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(covariance);
    // This cube/cuboid fixture has two equal transverse eigenvalues. Put that
    // plane in xy so a nonzero phi remains an actual stationary constraint state.
    frame.rotation.row(0)=eig.eigenvectors().col(1).transpose();
    frame.rotation.row(1)=eig.eigenvectors().col(2).transpose();
    frame.rotation.row(2)=eig.eigenvectors().col(0).transpose();
    if(frame.rotation.determinant()<0.)frame.rotation.row(0)*=-1.;
    return SourceCagePaperFunctional(f.source(),f.config(),frame);
}
SourceCageAdaptiveSolverSettings settings() {
    SourceCageAdaptiveSolverSettings s;s.maximumPasses=3;
    s.solve.source.solver.maximumIterations=1;s.solve.source.solver.maximumLineSearchSteps=2;
    s.solve.budget.maximumGuardCalls=8;
    s.budget.maximumKnownBufferBytes=std::size_t(8)*1024*1024*1024;
    return s;
}
void passLinks(const SourceCageAdaptiveSolverAudit& run) {
    const auto& facts=run.facts();std::size_t cumulative=0;
    for(std::size_t i=0;i<facts.passes.size();++i) {
        const auto& p=facts.passes[i];require(p.passNumber==i+1,"ordered pass numbering");
        if(p.initialState.size())require(p.source&&p.startingStateValues==std::size_t(p.initialState.size())&&p.startingBetaValues==p.source->betaCount(),"full starting dimensions");
        require(p.cumulativeKnownBufferBytes>=cumulative,"cumulative reservation never refunded");cumulative=p.cumulativeKnownBufferBytes;
        if(!p.solve)continue;
        const auto& g=p.solve->facts();const auto& raw=g.source;
        require(identical(raw.initialState,p.initialState),"fresh solve receives actual full state");
        require(p.guardedScientificSha256==g.scientificSha256,"raw guarded scientific digest retained");
        if(raw.numericalResultAvailable) {
            require(!raw.numerical.trace.empty()&&raw.numerical.trace.front().historyBefore==0&&p.freshHistory,"actual first trace has fresh history");
            require(raw.numerical.trace.front().state.size()==p.initialState.size(),"new trace has full changed dimension");
        }
        if(p.rebind) {
            require(raw.numericalResultAvailable&&raw.numerical.status==SlqpStatus::Converged&&raw.finalEvaluationAvailable&&
                g.finalStatus==SourceCageGuardedFinalStatus::ChecksPassed&&g.finalRecheckCount==1,"only converged freshly checked solve adapts");
            require(p.rebindScientificSha256==p.rebind->facts().scientificSha256,"actual transaction digest retained");
        }
        if(i&&facts.passes[i-1].rebind) {
            const auto checkpoint=facts.passes[i-1].rebind->checkpoint();
            require(checkpoint&&identical(checkpoint->state(),p.initialState),"private rebound checkpoint feeds next fresh solve");
            require(p.oldPartitionSha256==facts.passes[i-1].newPartitionSha256&&p.oldFunctionalSha256==facts.passes[i-1].newFunctionalSha256,"changed identity chain");
        }
    }
}
}
int main(int argc,char** argv) {try {
    static_assert(!std::is_default_constructible_v<SourceCageAdaptiveSolverAudit>);
    static_assert(!std::is_default_constructible_v<SourceCageAdaptiveSolvedEndpoint>);
    static_assert(std::is_const_v<std::remove_reference_t<decltype(std::declval<SourceCageAdaptiveSolvedEndpoint>().state())>>);
    auto original=fixture();auto x=original.initialState(Eigen::VectorXd::Constant(1,.25),.23);
    for(std::size_t h=0;h<original.handleCount();++h){x[original.handleOffset()+4*h]=1.01;x[original.handleOffset()+4*h+1]=.001;}
    const auto f=stationary(original,x);const auto s=settings();
    for(double threshold:{s.rebind.fractionalEpsilon,1.-s.rebind.fractionalEpsilon}) {
        auto edge=x;edge[0]=threshold;
        const auto edgeSource=stationary(original,edge);
        const auto boundary=solveSourceCageAdaptive(edgeSource,edge,s);
        require(boundary.facts().outcome==SourceCageAdaptiveSolverOutcome::RelaxedAdaptiveComplete&&boundary.facts().passes.size()==2&&
            boundary.facts().remainingFractionalLeaves==8&&boundary.facts().refinableFractionalLeaves==0,
            "inclusive exact-epsilon leaves refine and remain counted at maximum depth");
    }
    std::cout<<"KNOWN_SOLUTION_VALIDATION inclusive epsilon endpoints PASS\n";
    if(argc>1&&std::string(argv[1])=="--epsilon")return 0;
    const auto multi=solveSourceCageAdaptive(f,x,s);passLinks(multi);
    require(multi.facts().outcome==SourceCageAdaptiveSolverOutcome::RelaxedAdaptiveComplete&&multi.facts().passes.size()==2,"actual changed partition then unchanged completion");
    require(multi.lastSolved()&&multi.lastSolved()->functional().betaCount()==8&&!multi.pendingCheckpoint(),"completed endpoint is solved depth-four partition");
    require(multi.facts().passes[0].newPartitionSha256!=multi.facts().passes[0].oldPartitionSha256,"first pass really refines");
    require(identical(multi.lastSolved()->state().tail(4*f.handleCount()),x.tail(4*f.handleCount()))&&
        multi.lastSolved()->state()[multi.lastSolved()->functional().angleOffset()]==x[f.angleOffset()],"every handle and phi survive real remapping");
    require(multi.facts().remainingFractionalLeaves==8&&multi.facts().refinableFractionalLeaves==0,"depth-limited fractional completion is explicitly relaxed");
    std::cout<<"KNOWN_SOLUTION_VALIDATION multiple passes PASS\n";
    auto limit=s;limit.maximumPasses=1;
    const auto pending=solveSourceCageAdaptive(f,x,limit);passLinks(pending);
    require(pending.facts().outcome==SourceCageAdaptiveSolverOutcome::MaximumPasses&&pending.lastSolved()&&pending.pendingCheckpoint(),"pass limit keeps distinct solved and pending endpoints");
    require(pending.lastSolved()->functional().betaCount()==1&&pending.pendingCheckpoint()->functional().betaCount()==8,"unsolved rebound is never relabeled solved");
    const auto repeat=solveSourceCageAdaptive(f,x,limit);
    require(repeat.facts().scientificSha256==pending.facts().scientificSha256,"scientific digest reproducible without time");
    auto tiny=s;tiny.budget.maximumKnownBufferBytes=multi.facts().passes.front().cumulativeKnownBufferBytes;
    const auto aggregate=solveSourceCageAdaptive(f,x,tiny);passLinks(aggregate);
    require(aggregate.facts().outcome==SourceCageAdaptiveSolverOutcome::ResourceLimit&&aggregate.lastSolved()&&aggregate.pendingCheckpoint(),"aggregate next-pass exhaustion preserves endpoints");
    require(aggregate.facts().passes.size()==2&&!aggregate.facts().passes.back().solve,"next solve rejected before work");
    require(!aggregate.facts().passes.back().source&&aggregate.facts().passes.back().initialState.size()==0,"pre-admission pass has no fabricated source snapshot or solve");
    auto dimensionLimit=s;dimensionLimit.solve.source.budget.maximumStateValues=f.stateCount();
    const auto laterFailure=solveSourceCageAdaptive(f,x,dimensionLimit);passLinks(laterFailure);
    require(laterFailure.facts().outcome==SourceCageAdaptiveSolverOutcome::ResourceLimit&&laterFailure.lastSolved()&&laterFailure.pendingCheckpoint()&&
        laterFailure.facts().passes.size()==2&&!laterFailure.facts().passes.back().solve,"later actual adapter dimension failure retains solved and unsolved rebound");
    auto overflow=s;overflow.maximumPasses=std::numeric_limits<std::size_t>::max();
    const auto rejected=solveSourceCageAdaptive(f,x,overflow);
    require(rejected.facts().outcome==SourceCageAdaptiveSolverOutcome::ResourceLimit&&rejected.facts().passes.empty(),"checked aggregate metadata arithmetic overflow before work");
    overflow.maximumPasses=std::numeric_limits<std::size_t>::max()/2;
    require(solveSourceCageAdaptive(f,x,overflow).facts().outcome==SourceCageAdaptiveSolverOutcome::ResourceLimit,"checked metadata multiplication overflow before allocation");
    auto low=s;low.budget.maximumKnownBufferBytes=1;
    require(solveSourceCageAdaptive(f,x,low).facts().outcome==SourceCageAdaptiveSolverOutcome::ResourceLimit,"known-buffer cap before first work");
    auto graphs=s;graphs.budget.maximumRetainedSourcePairs=1;
    const auto graphLimit=solveSourceCageAdaptive(f,x,graphs);
    require(graphLimit.facts().outcome==SourceCageAdaptiveSolverOutcome::ResourceLimit&&graphLimit.lastSolved()&&!graphLimit.pendingCheckpoint(),"rebuilt graph count admitted before transaction");
    const auto unchanged=solveSourceCageAdaptive(multi.lastSolved()->functional(),multi.lastSolved()->state(),s);
    require(unchanged.facts().outcome==SourceCageAdaptiveSolverOutcome::RelaxedAdaptiveComplete&&unchanged.facts().passes.size()==1,"unchanged converged partition completes once");
    const auto& up=unchanged.facts().passes.front();
    require(up.oldPartitionSha256==up.newPartitionSha256&&up.solvedStateSha256==up.newStateSha256&&up.oldFunctionalSha256==up.newFunctionalSha256,"unchanged completion requires identities");
    const auto nontrivial=SourceCagePaperFunctional(f.source(),SourceCagePaperFunctionalConfig{},{});
    const auto nx=nontrivial.initialState(Eigen::VectorXd::Constant(1,.25),0.);
    const auto stalled=solveSourceCageAdaptive(nontrivial,nx,s);passLinks(stalled);
    require(stalled.facts().outcome==SourceCageAdaptiveSolverOutcome::SolveNotConverged&&!stalled.lastSolved()&&!stalled.pendingCheckpoint()&&
        stalled.facts().passes.front().solve&&!stalled.facts().passes.front().rebind,"actual nontrivial nonconvergence prevents adaptation");
    auto badGeometry=s;badGeometry.solve.geometry.minimumOuterDomainClearance=.2;
    const auto bad=solveSourceCageAdaptive(f,x,badGeometry);passLinks(bad);
    std::cout<<"initial geometry outcome="<<int(bad.facts().outcome)<<" reason="<<bad.facts().terminalReason.data()<<std::endl;
    if(bad.facts().passes.front().solve)std::cout<<"guard final="<<int(bad.facts().passes.front().solve->facts().finalStatus)<<" first attempt="<<int(bad.facts().passes.front().solve->facts().attempts.front().outcome)<<std::endl;
    require(bad.facts().outcome==SourceCageAdaptiveSolverOutcome::SolveGeometryRejected&&!bad.lastSolved()&&!bad.facts().passes.front().rebind,"actual initial geometry rejection prevents adaptation");
    auto rebindBad=s;rebindBad.rebind.geometry.minimumOuterDomainClearance=.2;
    const auto rb=solveSourceCageAdaptive(f,x,rebindBad);passLinks(rb);
    require(rb.facts().outcome==SourceCageAdaptiveSolverOutcome::RebindGeometryRejected&&rb.lastSolved()&&!rb.pendingCheckpoint()&&
        rb.facts().passes.front().rebind&&!rb.facts().passes.front().rebind->checkpoint(),"failed rebind cannot publish candidate");
    auto rebindResource=s;rebindResource.rebind.budget.maximumOwnedBytes=1;
    const auto rr=solveSourceCageAdaptive(f,x,rebindResource);
    require(rr.facts().outcome==SourceCageAdaptiveSolverOutcome::ResourceLimit&&rr.lastSolved()&&!rr.pendingCheckpoint()&&
        rr.facts().passes.front().rebind&&rr.facts().passes.front().rebind->facts().outcome==SourceCageAdaptiveRebindOutcome::ResourceLimit,
        "actual rebind resource failure retains solve and failed transaction evidence");
    auto move=solveSourceCageAdaptive(f,x,limit);const auto& beforeMove=move.facts();const auto keep=move.lastSolved();const auto checkpoint=move.pendingCheckpoint();
    const auto moved=std::move(move);require(beforeMove.scientificSha256==moved.facts().scientificSha256,"shared immutable evidence survives move");
    bool threw=false;try{(void)move.facts();}catch(const std::logic_error&){threw=true;}require(threw,"moved-from accessor rejected");
    auto copied=moved.facts();copied.outcome=SourceCageAdaptiveSolverOutcome::RelaxedAdaptiveComplete;
    require(moved.facts().outcome==SourceCageAdaptiveSolverOutcome::MaximumPasses&&moved.pendingCheckpoint(),"copied scalar facts cannot change owned conclusion or pending endpoint");
    auto lifetime=[&]{return solveSourceCageAdaptive(f,x,limit).lastSolved();}();
    require(lifetime->functional().evaluate(lifetime->state()).outputSha256==keep->functional().evaluate(keep->state()).outputSha256&&checkpoint->state().size()>keep->state().size(),"owned endpoint outlives audit");
    std::cout<<"adaptive solver regression PASS; no feasibility, binary, seal, physics or paper grant\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
#endif
