#include "spin2014/SourceCageGuardedSlqp.h"
#include "spin2014/SlqpAuditTrace.h"
#include "source_surface_test_utils.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>

template<class T> void vectorJson(const T& a) {
    std::cout << '[';
    for (Eigen::Index i=0;i<a.size();++i) { if(i)std::cout<<',';std::cout<<a(i); }
    std::cout << ']';
}
template<class T> void matrixJson(const T& a) {
    std::cout << '[';
    for(Eigen::Index i=0;i<a.rows();++i) {if(i)std::cout<<',';vectorJson(a.row(i));}
    std::cout << ']';
}
int main() {
    using namespace spin2014;
    std::cout<<std::setprecision(17)<<std::boolalpha;
    const auto outer=surface_test::cube(5./16,11./16),design=surface_test::cube(.34,.66);
    SourceSurfaceOctreeConfig sc;sc.maximumDepth=3;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
    SourceAdaptiveDomain domain(source,2,3);
    SourceSurfaceCageConfig cc;cc.cageLevel=2;
    SourceCageGridConfig gc;gc.fineLevel=3;
    SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),SourceCageFemConfig{}),SourceCageBbwConfig{});
    if(!bbw.success())throw std::runtime_error(bbw.terminalReason());
    SourceCageBinding binding(domain,bbw,SourceCageBindingConfig{});
    SourceConformingBoundary boundary(binding,SourceConformingBoundaryConfig{});
    SourceCageMomentFunctional material(boundary,1.25,SourceCageMomentFunctionalConfig{});
    const Spin2014Objective modes[]={Spin2014Objective::YoYo,Spin2014Objective::TopInertiaOnly,
        Spin2014Objective::TopLoweringOnly,Spin2014Objective::TopMassAndLowering,Spin2014Objective::StaticLowering};
    for(int test=0;test<10;++test) {
        const int mode=test<5?test:(test==8?0:(test==9?3:2));
        SourceCagePaperFunctionalConfig pc;pc.objective=modes[mode];
        pc.inertiaWeight=1.1;pc.centerWeight=.7;pc.fillWeight=.02;pc.shapeWeight=.06;pc.displacementWeight=.09;
        SourceCagePaperFrame frame;frame.origin=Eigen::Vector3d(.49,.52,.3);
        if(test>=8) {
            pc.inertiaWeight=.8;pc.centerWeight=.6;pc.fillWeight=.03;pc.shapeWeight=.04;pc.displacementWeight=.05;
            frame.origin=Eigen::Vector3d(.48,.51,.3);
        }
        SourceCagePaperFunctional paper(material,pc,frame);
        if(paper.betaCount()!=8||paper.handleCount()!=26)throw std::runtime_error("full source dimensions changed");
        Eigen::VectorXd beta(8);for(int i=0;i<8;++i)beta[i]=.13+.017*i;
        if(test>=8)beta<<.17,.24,.31,.38,.45,.52,.59,.66;
        auto x=paper.initialState(beta,mode==4?0.:(test>=8?.23:.29));
        SourceCageGuardedSlqpSettings settings;
        settings.source.solver.maximumIterations=2;settings.source.solver.maximumLineSearchSteps=8;
        std::string id="mode"+std::to_string(mode);
        if(test==5) {id="initial_reject";x.head(8).setZero();x[paper.handleOffset()+4*3+1]+=1.;}
        if(test==6) {id="guard_cap";settings.budget.maximumGuardCalls=1;}
        if(test==7) {id="source_cap";settings.source.budget.maximumSourceEvaluations=1;}
        if(test>=8) {
            id=test==8?"restoration":"exhaustion";
            settings.geometry.minimumOuterDomainClearance=.04;
            settings.source.solver.maximumIterations=1;
            if(test==9)settings.source.solver.maximumLineSearchSteps=4;
        }
        auto audit=solveSourceCageGuardedSlqp(paper,x,settings);
        const auto& a=audit.facts();const auto& s=a.source;
        std::cout<<"{\"id\":"<<std::quoted(id)<<",\"mode\":"<<mode<<",\"scope\":\"SYNTHETIC_REGRESSION\",\"sealed\":false"
            <<",\"feasibility\":"<<int(a.feasibility)<<",\"final_status\":"<<int(a.finalStatus)
            <<",\"source_outcome\":"<<int(s.outcome)<<",\"numerical_available\":"<<s.numericalResultAvailable
            <<",\"final_evaluation_available\":"<<s.finalEvaluationAvailable<<",\"legacy_geometry\":"<<int(s.geometryStatus)
            <<",\"source_evaluations\":"<<s.sourceEvaluationCount<<",\"callbacks\":"<<s.callbackCount
            <<",\"checker_calls\":"<<a.checkerInvocationCount<<",\"geometry_values\":"<<a.geometryValueEvaluationCount
            <<",\"geometry_values_complete\":"<<a.geometryValueEvaluationCountComplete<<",\"final_rechecks\":"<<a.finalRecheckCount
            <<",\"input_hash\":"<<std::quoted(a.configSha256)<<",\"result_hash\":"<<std::quoted(a.scientificSha256)
            <<",\"maximum_guard_calls\":"<<settings.budget.maximumGuardCalls
            <<",\"minimum_clearance\":"<<settings.geometry.minimumOuterDomainClearance<<",\"initial\":";vectorJson(x);
        std::cout<<",\"controls\":";writeSlqpAuditControlsJson(std::cout,settings.source.solver);
        std::cout<<",\"calls\":[";
        for(std::size_t i=0;i<s.recordCount;++i) {
            if(i)std::cout<<',';const auto& r=s.records[i];
            const Eigen::VectorXd state=s.callbackStates.col(i);
            if(r.status==SourceCageSlqpEvaluationStatus::Evaluated && paper.evaluate(state).outputSha256!=r.outputSha256.data())
                throw std::runtime_error("gradient output not bound to actual state");
            std::cout<<"{\"sequence\":"<<r.sequence<<",\"purpose\":"<<int(r.purpose)<<",\"status\":"<<int(r.status)
                <<",\"source_invoked\":"<<r.sourceInvoked<<",\"state_hash\":"<<std::quoted(r.stateSha256.data())
                <<",\"output_hash\":"<<std::quoted(r.outputSha256.data())<<",\"reason_hash\":"<<std::quoted(r.reasonSha256.data())
                <<",\"reason\":"<<std::quoted(r.reason.data())<<",\"state\":";vectorJson(state);std::cout<<'}';
        }
        std::cout<<"],\"guards\":[";
        for(std::size_t i=0;i<a.attempts.size();++i) {
            if(i)std::cout<<',';const auto& r=a.attempts[i];
            if(r.sourceOutputSha256[0] && paper.evaluateValues(s.callbackStates.col(r.callbackRecord)).outputSha256!=r.sourceOutputSha256.data())
                throw std::runtime_error("guard value output not bound to actual state");
            std::cout<<"{\"callback\":"<<r.callbackRecord<<",\"first_stage\":"<<r.firstStage<<",\"stage_count\":"<<r.stageCount
                <<",\"invoked\":"<<r.checkerInvoked<<",\"completed\":"<<r.checkCompleted<<",\"values_invoked\":"<<r.sourceValuesInvoked
                <<",\"outcome\":"<<int(r.outcome)<<",\"input_hash\":"<<std::quoted(r.inputSha256.data())
                <<",\"value_hash\":"<<std::quoted(r.sourceOutputSha256.data())<<",\"result_hash\":"<<std::quoted(r.scientificSha256.data())<<'}';
        }
        std::cout<<"],\"stages\":[";
        for(std::size_t i=0;i<a.stages.size();++i) {
            if(i)std::cout<<',';const auto& r=a.stages[i];
            std::cout<<"{\"stage\":"<<int(r.stage)<<",\"outcome\":"<<int(r.outcome)<<",\"first\":"<<r.first<<",\"second\":"<<r.second
                <<",\"reason\":"<<std::quoted(r.reason.data())<<",\"evidence_hash\":"<<std::quoted(r.evidenceSha256.data())
                <<",\"distance\":"<<std::quoted(r.squaredDistanceExact)<<'}';
        }
        std::cout<<"],\"trace\":[";
        for(std::size_t i=0;i<s.numerical.trace.size();++i) {
            if(i)std::cout<<',';const auto& t=s.numerical.trace[i];
            const auto e=paper.evaluate(t.state);
            if(e.totalObjective!=t.objective||!e.jointGradient.isApprox(t.gradient,0.)||!e.constraints.isApprox(t.constraints,0.)||!e.jointConstraintJacobian.isApprox(t.constraintJacobian,0.))
                throw std::runtime_error("accepted solver state not actual source output");
            std::cout<<"{\"iteration\":"<<t.iteration<<",\"phi\":"<<t.phi<<",\"objective\":"<<t.objective<<",\"merit\":"<<t.meritValue
                <<",\"accepted_step\":"<<t.acceptedStep<<",\"line_search_evaluations\":"<<t.lineSearchEvaluations;
            writeSlqpAuditJson(std::cout,t);std::cout<<'}';
        }
        std::cout<<"],\"numerical_status\":"<<std::quoted(slqpStatusToken(s.numerical.status))<<",\"final\":";vectorJson(s.numerical.finalState);
        if(s.finalEvaluationAvailable) {
            if(paper.evaluate(s.numerical.finalState).outputSha256!=s.finalEvaluation.outputSha256)
                throw std::runtime_error("final evaluation not bound to actual final state");
            std::cout<<",\"material\":{\"density\":"<<material.density()<<",\"points\":";matrixJson(s.finalEvaluation.world.deformedPoints);
            std::cout<<",\"outer_faces\":";matrixJson(outer.f);std::cout<<",\"moments\":";vectorJson(s.finalEvaluation.world.finalMoments.values());
            std::cout<<",\"primitives\":[";
            for(std::size_t p=0;p<boundary.primitiveBoundaries().size();++p) {
                if(p)std::cout<<',';std::cout<<"{\"owner\":"<<binding.primitives()[p].betaIndex<<",\"faces\":[";
                const auto& refs=boundary.primitiveBoundaries()[p];
                for(std::size_t k=0;k<refs.size();++k) {
                    if(k)std::cout<<',';auto t=boundary.triangles()[refs[k].triangle];if(refs[k].orientation<0)std::swap(t[1],t[2]);
                    std::cout<<'['<<t[0]<<','<<t[1]<<','<<t[2]<<']';
                }std::cout<<"]}";
            }std::cout<<"]}";
        }
        std::cout<<"}\n";
    }
}
