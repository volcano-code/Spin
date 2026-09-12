#include "spin2014/SourceCageAdaptiveSolver.h"
#include "source_surface_test_utils.h"
#include "task-5x-controller-guard-export.h"

const char* outcomeToken(spin2014::SourceCageAdaptiveSolverOutcome value) {
    using O=spin2014::SourceCageAdaptiveSolverOutcome;
    switch(value) {
    case O::RelaxedAdaptiveComplete:return "RELAXED_ADAPTIVE_COMPLETE";
    case O::MaximumPasses:return "MAXIMUM_PASSES";
    case O::SolveNotConverged:return "SOLVE_NOT_CONVERGED";
    case O::SolveGeometryRejected:return "SOLVE_GEOMETRY_REJECTED";
    case O::SolveUnknown:return "SOLVE_UNKNOWN";
    case O::RebindGeometryRejected:return "REBIND_GEOMETRY_REJECTED";
    case O::RebindUnknown:return "REBIND_UNKNOWN";
    case O::IdentityMismatch:return "IDENTITY_MISMATCH";
    case O::ResourceLimit:return "RESOURCE_LIMIT";
    }
    throw std::runtime_error("unrecognized adaptive termination");
}

int main() {
    using namespace spin2014;
    std::cout<<std::setprecision(17)<<std::boolalpha;
    // Fixed centered geometry. The zero objective has known stationary states;
    // this is an independent control-flow protocol, not a paper experiment.
    const auto outer=surface_test::cube(.3125,.6875),design=surface_test::cube(.425,.575);
    SourceSurfaceOctreeConfig sc;sc.maximumDepth=4;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
    const SourceAdaptiveDomain domain(source,2,4);
    SourceSurfaceCageConfig cc;cc.cageLevel=2;
    SourceCageGridConfig gc;gc.fineLevel=3;
    const SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),{}),{});
    surface_test::require(bbw.success(),"actual full BBW required");
    const SourceCageMomentFunctional material(SourceConformingBoundary(SourceCageBinding(domain,bbw,{}),{}),1.17,{});
    const char* names[]={"stationary_multi","stationary_pass_limit","nontrivial_stop",
                         "initial_geometry_reject","rebind_geometry_reject","aggregate_limit"};
    for(int test=0;test<6;++test) {
        SourceCagePaperFunctionalConfig pc;pc.objective=Spin2014Objective::StaticLowering;
        pc.centerWeight=test==2?1.:0.;pc.inertiaWeight=0.;
        SourceCagePaperFrame frame;frame.origin=Eigen::Vector3d(.5,.5,.3);
        if(test==2)frame.origin=Eigen::Vector3d(.49,.52,.3);
        const SourceCagePaperFunctional f(material,pc,frame);
        const Eigen::VectorXd x=f.initialState(Eigen::VectorXd::Constant(f.betaCount(),.31),0.);
        SourceCageAdaptiveSolverSettings s;s.maximumPasses=test==1?1:4;
        s.solve.source.solver.maximumIterations=1;s.solve.source.solver.maximumLineSearchSteps=2;
        s.solve.budget.maximumGuardCalls=8;
        if(test==3)s.solve.geometry.minimumOuterDomainClearance=.2;
        if(test==4)s.rebind.geometry.minimumOuterDomainClearance=.2;
        if(test==5)s.budget.maximumKnownBufferBytes=320U*1024U*1024U;
        const auto audit=solveSourceCageAdaptive(f,x,s);const auto& a=audit.facts();
        std::cout<<"{\"id\":"<<std::quoted(names[test])<<",\"scope\":"
            <<std::quoted(test==2?"SYNTHETIC_REGRESSION":"KNOWN_SOLUTION_VALIDATION")
            <<",\"sealed\":false,\"paper_reproduction\":false,\"outcome\":"<<std::quoted(outcomeToken(a.outcome))
            <<",\"hash\":"<<std::quoted(a.scientificSha256)<<",\"reason\":"<<std::quoted(a.terminalReason.data())
            <<",\"maximum_passes\":"<<s.maximumPasses<<",\"maximum_level\":4,\"epsilon\":"<<s.rebind.fractionalEpsilon
            <<",\"known_buffer_cap\":"<<s.budget.maximumKnownBufferBytes
            <<",\"known_buffer_bytes\":"<<a.resources.cumulativeKnownBufferBytes
            <<",\"guarded_reserved\":"<<a.resources.guardedReservedBytes
            <<",\"rebind_reserved\":"<<a.resources.rebindReservedBytes
            <<",\"admitted_pairs\":"<<a.resources.admittedSourcePairs
            <<",\"retained_pairs\":"<<a.resources.retainedSourcePairs
            <<",\"counts_available\":"<<a.fractionalCountsAvailable
            <<",\"fractional\":"<<a.remainingFractionalLeaves<<",\"refinable\":"<<a.refinableFractionalLeaves
            <<",\"passes\":[";
        for(std::size_t index=0;index<a.passes.size();++index) {
            if(index)std::cout<<',';const auto& p=a.passes[index];
            std::cout<<"{\"number\":"<<p.passNumber<<",\"old_partition\":"<<std::quoted(p.oldPartitionSha256)
                <<",\"new_partition\":"<<std::quoted(p.newPartitionSha256)
                <<",\"old_functional\":"<<std::quoted(p.oldFunctionalSha256)
                <<",\"new_functional\":"<<std::quoted(p.newFunctionalSha256)
                <<",\"old_state_hash\":"<<std::quoted(p.oldStateSha256)
                <<",\"solved_state_hash\":"<<std::quoted(p.solvedStateSha256)
                <<",\"new_state_hash\":"<<std::quoted(p.newStateSha256)
                <<",\"fresh_history\":"<<p.freshHistory
                <<",\"known_buffer_bytes\":"<<p.cumulativeKnownBufferBytes<<",\"source\":";
            dumpSource(*p.source,p.initialState);
            if(p.solve) {
                surface_test::require(p.guardedScientificSha256==p.solve->facts().scientificSha256,"pass must retain actual guarded digest");
                std::cout<<",\"guarded\":";dumpGuarded(*p.source,*p.solve);
            }
            if(p.rebind) {
                const auto& r=p.rebind->facts();
                surface_test::require(p.rebindScientificSha256==r.scientificSha256,"pass must retain actual rebind digest");
                std::cout<<",\"rebind\":{\"hash\":"<<std::quoted(r.scientificSha256)
                    <<",\"outcome\":"<<int(r.outcome)<<",\"last_stage\":"<<int(r.lastStage);
                const auto checkpoint=p.rebind->checkpoint();
                if(checkpoint) {std::cout<<",\"source\":";dumpSource(checkpoint->functional(),checkpoint->state());}
                std::cout<<'}';
            }
            std::cout<<'}';
        }
        std::cout<<']';
        if(const auto last=audit.lastSolved()) {std::cout<<",\"last_solved\":";dumpSource(last->functional(),last->state());}
        if(const auto pending=audit.pendingCheckpoint()) {std::cout<<",\"pending\":";dumpSource(pending->functional(),pending->state());}
        std::cout<<"}\n";
    }
}
