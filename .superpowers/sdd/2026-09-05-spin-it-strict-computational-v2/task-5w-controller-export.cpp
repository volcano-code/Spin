#include "spin2014/SourceCageAdaptiveRebind.h"
#include "spin2014/SourceCageGuardedSlqp.h"
#include "spin2014/SlqpAuditTrace.h"
#include "source_surface_test_utils.h"
#include "task-5w-controller-export-helpers.h"

int main() {
    using namespace spin2014;
    using surface_test::require;
    std::cout<<std::setprecision(17)<<std::boolalpha;
    // Predetermined synthetic geometry, not a reconstructed author case.
    const auto outer=surface_test::cube(.3125,.6875);
    const auto design=surface_test::cube(.355,.52);
    SourceSurfaceOctreeConfig sc;sc.maximumDepth=4;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
    SourceSurfaceCageConfig cc;cc.cageLevel=2;
    SourceCageGridConfig gc;gc.fineLevel=3;
    const SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),{}),{});
    require(bbw.success(),"actual coupled BBW run required");
    for(int test=0;test<6;++test) {
        const bool split=test<4;
        const SourceAdaptiveDomain domain(source,split?3:4,4);
        const SourceCageMomentFunctional material(SourceConformingBoundary(SourceCageBinding(domain,bbw,{}),{}),1.17,{});
        SourceCagePaperFunctionalConfig pc;
        pc.objective=test==3?Spin2014Objective::StaticLowering:Spin2014Objective::TopMassAndLowering;
        pc.inertiaWeight=.9;pc.centerWeight=.8;pc.fillWeight=.025;
        pc.shapeWeight=.03;pc.displacementWeight=.07;
        SourceCagePaperFrame frame;frame.origin=Eigen::Vector3d(.46,.48,.3);
        const SourceCagePaperFunctional f(material,pc,frame);
        Eigen::VectorXd beta=Eigen::VectorXd::Constant(f.betaCount(),split?.31:.01);
        if(test==5)for(Eigen::Index i=0;i<beta.size();++i)beta[i]=.011+.001*i;
        Eigen::VectorXd x=f.initialState(beta,f.hasAngle()?.27:0.);
        std::string id;
        if(test==0)id="split_rest";
        if(test==1) {
            id="split_affine";
            const auto& cage=bbw.fem().grid().cage();
            const Eigen::Vector3d offset(.03,-.02,.01);
            for(std::size_t h=0;h<f.handleCount();++h) {
                x[f.handleOffset()+4*h]=1.125;
                x.segment<3>(f.handleOffset()+4*h+1)=offset+.125*cage.boundaryVertices().row(cage.handleVertexIndices()[h]).transpose();
            }
        }
        if(test==2) {
            id="split_nonaffine";
            for(std::size_t h=0;h<f.handleCount();++h) {
                x[f.handleOffset()+4*h]=1.+.002*std::sin(double(h)+.2);
                x[f.handleOffset()+4*h+1]=.001*std::cos(2.*double(h)+.3);
                x[f.handleOffset()+4*h+2]=.0015*std::sin(3.*double(h)+.4);
                x[f.handleOffset()+4*h+3]=.0007*std::cos(4.*double(h)+.5);
            }
        }
        if(test==3)id="split_static";
        if(test==4)id="merge_equal";
        if(test==5)id="unequal_no_merge";
        const auto oldHash=f.evaluate(x).outputSha256;
        const auto oldPartition=domain.partitionHash();
        const auto audit=rebindSourceCageAdaptive(f,x,{});
        const auto checkpoint=audit.checkpoint();
        require(bool(checkpoint),"controller fixture must actually produce a checkpoint");
        require(f.evaluate(x).outputSha256==oldHash && domain.partitionHash()==oldPartition,"original input remains unchanged");
        const auto& next=checkpoint->functional();const auto& state=checkpoint->state();
        const auto& facts=audit.facts();
        require(facts.beforeMaterial.available && facts.afterMaterial.available && facts.deformedMomentDifferenceAvailable,"complete actual material records required");
        const auto newMaterial=next.evaluate(state);
        require(facts.beforeMaterial.outputSha256==oldHash && facts.afterMaterial.outputSha256==newMaterial.outputSha256,"retained derivative outputs match direct input/output source evaluations");
        require(facts.afterMaterial.worldMoments.values()==newMaterial.world.finalMoments.values(),"retained moments bind actual output");
        require(facts.inputGeometry && facts.outputGeometry,"both actual audits retained");
        require(facts.inputGeometry->facts().outcome==SourceCageGeometryOutcome::ChecksPassed && facts.outputGeometry->facts().outcome==SourceCageGeometryOutcome::ChecksPassed,"both actual geometry checks complete");
        require(facts.inputGeometry->facts().sourceOutputSha256==f.evaluateValues(x).outputSha256 && facts.outputGeometry->facts().sourceOutputSha256==next.evaluateValues(state).outputSha256,"geometry audits bind actual source states");
        std::cout<<"{\"id\":"<<std::quoted(id)<<",\"scope\":\"SYNTHETIC_REGRESSION\",\"sealed\":false,\"converged\":false,\"checkpoint\":true"
            <<",\"scientific_hash\":"<<std::quoted(audit.facts().scientificSha256)<<",\"before\":";
        dumpSource(f,x);std::cout<<",\"after\":";dumpSource(next,state);
        std::cout<<",\"delta\":";arrayJson(audit.facts().deformedMomentDifference);
        if(test==3) {
            SourceCageGuardedSlqpSettings settings;
            settings.source.solver.maximumIterations=1;
            settings.source.solver.maximumLineSearchSteps=2;
            const auto solve=solveSourceCageGuardedSlqp(next,state,settings);
            const auto& facts=solve.facts();const auto& run=facts.source;
            require(run.numericalResultAvailable && !run.numerical.trace.empty(),"actual rebound consumer must evaluate");
            require(run.initialState.isApprox(state,0.),"consumer uses remapped state");
            std::cout<<",\"consumer\":{\"status\":"<<std::quoted(slqpStatusToken(run.numerical.status))
                <<",\"guard_calls\":"<<facts.checkerInvocationCount<<",\"final_rechecks\":"<<facts.finalRecheckCount
                <<",\"final_status\":"<<int(facts.finalStatus)<<",\"initial\":";arrayJson(run.initialState);
            std::cout<<",\"trace0\":{\"iteration\":"<<run.numerical.trace[0].iteration;
            writeSlqpAuditJson(std::cout,run.numerical.trace[0]);std::cout<<"}}";
        }
        std::cout<<"}\n";
    }
}
