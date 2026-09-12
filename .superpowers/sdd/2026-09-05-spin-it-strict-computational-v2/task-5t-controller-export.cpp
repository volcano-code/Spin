#include "spin2014/SourceCageSlqp.h"
#include "spin2014/SlqpAuditTrace.h"
#include <fstream>
#include "source_surface_test_utils.h"
#include <Eigen/Geometry>
#include <iomanip>
#include <iostream>

template<class T> void row(const T& v) {
    for(Eigen::Index i=0;i<v.size();++i) std::cout << ' ' << v(i);
}
template<class T> void seq(const T& v) {for(auto x:v)std::cout << ' ' << x;}
template<class T> void jsonVector(std::ostream& out,const T& a) {
    out<<'[';for(Eigen::Index i=0;i<a.size();++i){if(i)out<<',';out<<a(i);}out<<']';
}
void traceJson(std::ostream& out,int fixture,int mode,
    const spin2014::SourceCagePaperFunctional& paper,const Eigen::VectorXd& initial,
    const spin2014::SlqpLbfgsSettings& controls,const spin2014::SlqpLbfgsResult& result,
    const std::string& configHash,const std::string& resultHash,
    const char* geometry,const char* feasibility,const spin2014::SourceCageSlqpResult& sourceResult) {
    using namespace spin2014;
    out<<std::setprecision(17)<<"{\"fixture\":"<<fixture<<",\"mode\":"<<mode
       <<",\"beta_count\":"<<paper.betaCount()<<",\"handle_count\":"<<paper.handleCount()
       <<",\"geometry_status\":"<<std::quoted(geometry)<<",\"feasibility_status\":"<<std::quoted(feasibility)
       <<",\"algorithm\":"<<std::quoted(slqpAlgorithmToken(result.algorithm))
       <<",\"status\":"<<std::quoted(slqpStatusToken(result.status))
       <<",\"source_hash\":"<<std::quoted(paper.upstreamSha256())
       <<",\"config_hash\":"<<std::quoted(configHash)<<",\"result_hash\":"<<std::quoted(resultHash)
       <<",\"trace_hash\":"<<std::quoted(result.traceHash)<<",\"initial_phi\":"<<result.initialPhi
       <<",\"final_phi\":"<<result.finalPhi
       <<",\"numerical_only\":"<<((std::string(geometry)=="NOT_EVALUATED"&&std::string(feasibility)=="UNKNOWN")?"true":"false")
       <<",\"initial\":";jsonVector(out,initial);out<<",\"final\":";jsonVector(out,result.finalState);
    out<<",\"controls\":";writeSlqpAuditControlsJson(out,controls);out<<",\"trace\":[";
    for(std::size_t i=0;i<result.trace.size();++i) {
        if(i)out<<',';const auto& t=result.trace[i];
        const auto actual=paper.evaluate(t.state);
        if(actual.totalObjective!=t.objective || !actual.jointGradient.isApprox(t.gradient,0.) ||
           !actual.constraints.isApprox(t.constraints,0.) || !actual.jointConstraintJacobian.isApprox(t.constraintJacobian,0.))
            throw std::runtime_error("solver trace is not the actual source functional at its recorded state");
        out<<"{\"iteration\":"<<t.iteration<<",\"phi\":"<<t.phi<<",\"objective\":"<<t.objective
           <<",\"kkt_infinity\":"<<t.kktInfinityNorm<<",\"merit\":"<<t.meritValue
           <<",\"accepted_step\":"<<t.acceptedStep<<",\"line_search_evaluations\":"<<t.lineSearchEvaluations;
        writeSlqpAuditJson(out,t);out<<'}';
    }out<<"],\"callback_count\":"<<sourceResult.callbackCount
         <<",\"source_evaluation_count\":"<<sourceResult.sourceEvaluationCount<<",\"calls\":[";
    for(std::size_t i=0;i<sourceResult.recordCount;++i) {
        if(i)out<<',';const auto& event=sourceResult.records[i];
        const Eigen::VectorXd state=sourceResult.callbackStates.col(i);
        if(event.status==SourceCageSlqpEvaluationStatus::Evaluated &&
           paper.evaluate(state).outputSha256!=event.outputSha256.data())
            throw std::runtime_error("source call output digest not bound to actual recorded state");
        out<<"{\"sequence\":"<<event.sequence<<",\"purpose\":"<<int(event.purpose)
           <<",\"status\":"<<int(event.status)<<",\"source_invoked\":"<<(event.sourceInvoked?"true":"false")
           <<",\"state_hash\":"<<std::quoted(event.stateSha256.data())
           <<",\"output_hash\":"<<std::quoted(event.outputSha256.data())
           <<",\"reason_hash\":"<<std::quoted(event.reasonSha256.data())
           <<",\"reason\":"<<std::quoted(event.reason.data())<<",\"state\":";
        jsonVector(out,state);out<<'}';
    }out<<"]}\n";
}
int main(int argc,char** argv) {
    if(argc!=2)throw std::runtime_error("expected new trace output path");
    std::ofstream traceStream(argv[1]);if(!traceStream)throw std::runtime_error("trace output unavailable");
    using namespace spin2014;
    const Spin2014Objective modes[]={Spin2014Objective::YoYo,
        Spin2014Objective::TopInertiaOnly,Spin2014Objective::TopLoweringOnly,
        Spin2014Objective::TopMassAndLowering,Spin2014Objective::StaticLowering};
    std::cout << std::setprecision(17);
    for(int f=0;f<2;++f) {
        const double scale=f?1.7:1.;
        const Eigen::Vector3d offset=f?Eigen::Vector3d(-.2,.15,.6):Eigen::Vector3d::Zero();
        auto outer=f?surface_test::box({.31,.32,.3},{.69,.65,.7}):surface_test::cube(.3,.7);
        auto design=f?surface_test::box({.35,.36,.38},{.61,.6,.65}):surface_test::cube(.34,.66);
        outer.v=(scale*outer.v).rowwise()+offset.transpose();
        design.v=(scale*design.v).rowwise()+offset.transpose();
        SourceSurfaceOctreeConfig sc;sc.maximumDepth=5;sc.rootOrigin=offset;sc.rootSide=scale;
        auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
        SourceAdaptiveDomain domain(source,2,5);
        SourceSurfaceCageConfig cc;cc.cageLevel=2;cc.rootOrigin=offset;cc.rootSide=scale;
        SourceCageGridConfig gc;gc.fineLevel=3;
        const SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),SourceCageFemConfig{}),SourceCageBbwConfig{});
        if(!bbw.success())throw std::runtime_error(bbw.terminalReason());
        const SourceCageBinding binding(domain,bbw,SourceCageBindingConfig{});
        const SourceConformingBoundary boundary(binding,SourceConformingBoundaryConfig{});
        const double density=f?2.3:.8;
        const SourceCageMomentFunctional material(boundary,density,SourceCageMomentFunctionalConfig{});
        Eigen::VectorXd beta(material.betaCount());
        for(Eigen::Index i=0;i<beta.size();++i)beta(i)=.13+.017*i;
        auto handles=material.initialHandleParameters();
        for(std::size_t h=0;h<material.handleCount();++h) {
            handles(4*h)=1+.01*(int(h%5)-2);
            handles(4*h+1)=scale*.006*(int(h%3)-1);
            handles(4*h+2)=scale*.004*(int(h%5)-2);
            handles(4*h+3)=scale*.005*(int(h%7)-3);
        }
        SourceCagePaperFrame frame;frame.origin=offset+scale*Eigen::Vector3d(.41,.43,.19);
        frame.rotation=(Eigen::AngleAxisd(.31,Eigen::Vector3d::UnitX())*
            Eigen::AngleAxisd(-.23,Eigen::Vector3d::UnitY())*
            Eigen::AngleAxisd(.17,Eigen::Vector3d::UnitZ())).toRotationMatrix();
        const auto& grid=bbw.fem().grid();const auto& cage=grid.cage();
        const auto initialBeta=beta,initialHandles=handles;
        for(int mode=0;mode<5;++mode) {
            auto beta=initialBeta,handles=initialHandles;
            SourceCagePaperFunctionalConfig cfg;cfg.objective=modes[mode];
            cfg.inertiaWeight=1.7;cfg.centerWeight=.6;cfg.fillWeight=.23;
            cfg.shapeWeight=.37;cfg.displacementWeight=.19;
            const SourceCagePaperFunctional paper(material,cfg,frame);
            double phi=mode==4?0.:.37;
            auto initial=paper.initialState(beta,phi);initial.tail(handles.size())=handles;
            SourceCageSlqpSettings settings;settings.solver.maximumIterations=6;
            if(f==0 && mode==0) {
                auto low=cfg;low.maximumJointDerivativeValues=1;
                const SourceCagePaperFunctional capped(material,low,frame);
                const auto failed=solveSourceCageSlqp(capped,initial,settings);
                if(failed.outcome!=SourceCageSlqpOutcome::ResourceLimit || failed.numericalResultAvailable ||
                   failed.finalEvaluationAvailable || failed.recordCount!=1 || !failed.records[0].sourceInvoked ||
                   failed.records[0].status!=SourceCageSlqpEvaluationStatus::ResourceLimit ||
                   !failed.numerical.trace.empty())
                    throw std::runtime_error("first actual resource failure lost or fabricated numerical result");
                auto overflow=initial;
                for(std::size_t h=0;h<paper.handleCount();++h)overflow[paper.handleOffset()+4*h]=1e200;
                const auto invalid=solveSourceCageSlqp(paper,overflow,settings);
                if(invalid.outcome!=SourceCageSlqpOutcome::InitialEvaluationFailure || invalid.numericalResultAvailable ||
                   invalid.finalEvaluationAvailable || invalid.recordCount!=1 || !invalid.records[0].sourceInvoked ||
                   invalid.records[0].status!=SourceCageSlqpEvaluationStatus::InvalidNumericalTrial ||
                   !invalid.numerical.trace.empty())
                    throw std::runtime_error("first actual numerical failure lost or fabricated numerical result");
                std::cerr<<"PASS actual initial resource/numerical failures retain records; no invented KKT/termination\n";
            }
            const auto solved=solveSourceCageSlqp(paper,initial,settings);
            if(!solved.numericalResultAvailable || !solved.finalEvaluationAvailable || solved.outcome!=SourceCageSlqpOutcome::Completed)
                throw std::runtime_error("final source evaluation unavailable or resource-limited");
            const auto& final=solved.numerical.finalState;
            beta=final.head(paper.betaCount());handles=final.tail(4*paper.handleCount());
            phi=paper.hasAngle()?final[paper.angleOffset()]:0.;
            const auto& r=solved.finalEvaluation;const auto v=paper.evaluateValues(final);
            traceJson(traceStream,f,mode,paper,initial,settings.solver,solved.numerical,
                solved.configSha256,solved.scientificSha256,
                solved.geometryStatus==SourceCageSlqpGeometryStatus::NotEvaluated?"NOT_EVALUATED":"UNEXPECTED",
                solved.feasibility==SourceCageSlqpFeasibility::Unknown?"UNKNOWN":"UNEXPECTED",solved);
            if(!paper.shapeFunctional())throw std::runtime_error("missing actual shape functional");
            const auto shape=paper.shapeFunctional()->evaluate(r.world.deformedPoints.topRows(binding.outerPointCount()));
            const auto shapeValues=paper.shapeFunctional()->evaluateValues(r.world.deformedPoints.topRows(binding.outerPointCount()));
            auto disabledConfig=cfg;disabledConfig.shapeWeight=0;disabledConfig.displacementWeight=0;
            const SourceCagePaperFunctional disabled(material,disabledConfig,frame);
            const auto baseline=disabled.evaluate(beta,phi,handles);
            std::cout << "CASE " << f << ' ' << mode << ' ' << density << ' ' << beta.size() << ' ' << material.handleCount() << '\n';
            std::cout << "IDENT " << paper.configSha256() << ' ' << r.outputSha256 << ' ' << v.outputSha256 << '\n';
            std::cout << "FLAGS " << r.hasJacobians << ' ' << v.hasJacobians << ' ' << v.jointGradient.size() << ' ' << v.jointConstraintJacobian.size() << ' ' << v.spinMomentJacobian.size() << ' ' << v.world.betaJacobian.size() << ' ' << v.world.handleJacobian.size() << '\n';
            std::cout << "ORIGIN";row(frame.origin);std::cout << '\n';
            for(int i=0;i<3;++i){std::cout << "R " << i;row(frame.rotation.row(i));std::cout << '\n';}
            std::cout << "WEIGHTS " << cfg.inertiaWeight << ' ' << cfg.centerWeight << ' ' << cfg.fillWeight << ' ' << cfg.shapeWeight << ' ' << cfg.displacementWeight << '\n';
            std::cout << "OUTERCOUNT " << binding.outerPointCount() << '\n';
            std::cout << "SHAPE " << r.shapeEnergy << ' ' << r.displacementEnergy << '\n';
            std::cout << "VSHAPE " << v.shapeEnergy << ' ' << v.displacementEnergy << '\n';
            std::cout << "SHAPEFLAGS " << shape.hasGradient << ' ' << shapeValues.hasGradient << ' ' << shapeValues.vertexGradient.size() << ' ' << shapeValues.handleGradient.size() << '\n';
            for(Eigen::Index i=0;i<shape.vertexGradient.rows();++i){
                std::cout << "SHAPEVGRAD " << i;row(shape.vertexGradient.row(i));std::cout << '\n';
                std::cout << "SHAPERES " << i;row(shape.residual.row(i));std::cout << '\n';
            }
            std::cout << "SHAPEHGRAD";row(shape.handleGradient);std::cout << '\n';
            std::cout << "BASEGRAD";row(baseline.jointGradient);std::cout << '\n';
            if(!r.world.finalMoments.values().isApprox(baseline.world.finalMoments.values(),0) ||
               !r.constraints.isApprox(baseline.constraints,0) ||
               !r.spinMomentJacobian.isApprox(baseline.spinMomentJacobian,0) ||
               !r.jointConstraintJacobian.isApprox(baseline.jointConstraintJacobian,0))
                throw std::runtime_error("regularization mutated moments or constraints");
            std::cout << "PHI " << phi << '\n';
            std::cout << "BETA";row(beta);std::cout << '\n';
            std::cout << "PARAM";row(handles);std::cout << '\n';
            for(Eigen::Index i=0;i<grid.nodePositions().rows();++i){std::cout << "NODE " << i;row(grid.nodePositions().row(i));row(bbw.finalWeights().row(i));std::cout << '\n';}
            for(std::size_t i=0;i<grid.cells().size();++i){std::cout << "CELL " << i;seq(grid.cells()[i]);seq(grid.cellNodes()[i]);std::cout << '\n';}
            for(std::size_t h=0;h<material.handleCount();++h){std::cout << "HANDLE " << h;row(cage.boundaryVertices().row(cage.handleVertexIndices()[h]));std::cout << '\n';}
            for(Eigen::Index i=0;i<boundary.points().rows();++i){
                std::cout << "POINT " << i;row(boundary.points().row(i));std::cout << '\n';
                std::cout << "MOVED " << i;row(r.world.deformedPoints.row(i));std::cout << '\n';
                std::cout << "VMOVED " << i;row(v.world.deformedPoints.row(i));std::cout << '\n';
            }
            for(Eigen::Index i=0;i<outer.f.rows();++i){std::cout << "FACE " << i;row(outer.f.row(i));std::cout << '\n';}
            for(std::size_t i=0;i<boundary.triangles().size();++i){std::cout << "TRI " << i;seq(boundary.triangles()[i]);std::cout << '\n';}
            for(const auto& p:boundary.patches()){
                auto owner=[&](std::size_t id)->long long{return id==SourceConformingBoundary::noOwner?-1:binding.primitives()[id].betaIndex;};
                std::cout << "PATCH " << owner(p.negativeOwner) << ' ' << owner(p.positiveOwner) << ' ' << p.triangleBegin << ' ' << p.triangleCount << '\n';
            }
            for(auto e:domain.edges())std::cout << "EDGE " << e.first << ' ' << e.second << '\n';
            std::cout << "OBJ " << r.totalObjective << ' ' << r.paperObjective << ' ' << r.fillPenalty << '\n';
            std::cout << "VOBJ " << v.totalObjective << ' ' << v.paperObjective << ' ' << v.fillPenalty << '\n';
            std::cout << "WORLD";row(r.world.finalMoments.values());std::cout << '\n';
            std::cout << "SPIN";row(r.spinMoments.values());std::cout << '\n';
            std::cout << "VSPIN";row(v.spinMoments.values());std::cout << '\n';
            std::cout << "CONSTRAINT";row(r.constraints);std::cout << '\n';
            std::cout << "VCONSTRAINT";row(v.constraints);std::cout << '\n';
            std::cout << "GRAD";row(r.jointGradient);std::cout << '\n';
            for(Eigen::Index i=0;i<r.spinMomentJacobian.cols();++i){std::cout << "JM " << i;row(r.spinMomentJacobian.col(i));std::cout << '\n';}
            for(Eigen::Index i=0;i<r.jointConstraintJacobian.cols();++i){std::cout << "JC " << i;row(r.jointConstraintJacobian.col(i));std::cout << '\n';}
            std::cout << "END\n";
        }
    }
}
