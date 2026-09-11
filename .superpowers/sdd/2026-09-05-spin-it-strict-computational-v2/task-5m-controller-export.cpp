#include "spin2014/CageSpinDeformation.h"
#include "spin2014/SourceCageBbw.h"
#include "source_surface_test_utils.h"
#include <Eigen/Geometry>
#include <cmath>
#include <iomanip>
#include <iostream>

template<class Matrix> void matrix(const Matrix& m) {
    std::cout<<'[';
    for(int i=0;i<m.rows();++i) {
        if(i)std::cout<<',';std::cout<<'[';
        for(int j=0;j<m.cols();++j) {if(j)std::cout<<',';std::cout<<m(i,j);}
        std::cout<<']';
    }
    std::cout<<']';
}
template<class Vector> void vector(const Vector& v) {
    std::cout<<'[';for(int i=0;i<v.size();++i) {if(i)std::cout<<',';std::cout<<v(i);}std::cout<<']';
}
void mesh(const spin::TriangleMesh& m,const Eigen::MatrixXd& weights) {
    std::cout<<"{\"vertices\":";matrix(m.vertices);
    std::cout<<",\"faces\":";matrix(m.faces);
    std::cout<<",\"weights\":";matrix(weights);std::cout<<'}';
}
int main() {
    try {
        using namespace spin2014;
        auto original=surface_test::cube(.3,.7);
        original.v.col(1)=.91*original.v.col(1).array()+.037;
        original.v.col(2)=.87*original.v.col(2).array()+.051;
        SourceSurfaceCageConfig cc;cc.cageLevel=2;
        SourceCageGridConfig gc;gc.fineLevel=3;
        const SourceCageGrid grid(SourceSurfaceCage(original.v,original.f,cc),gc);
        const SourceCageFem fem(grid,{});
        const SourceCageBbw bbw(fem);
        if(!bbw.success())throw std::runtime_error(bbw.terminalReason());
        CageSpinDeformationRequest request;
        request.outerSurface.vertices=original.v;request.outerSurface.faces=original.f;
        request.outerWeights=grid.interpolate(bbw.finalWeights(),original.v);
        for(auto bounds:{std::pair<double,double>{.38,.44},{.54,.61}}) {
            const auto box=surface_test::cube(bounds.first,bounds.second);
            spin::TriangleMesh cell;cell.vertices=box.v;cell.faces=box.f;
            request.voidCells.push_back(cell);
            request.voidCellWeights.push_back(grid.interpolate(bbw.finalWeights(),cell.vertices));
        }
        for(int n:grid.handleNodeIndices()) {
            deformation::Handle h;h.restPivot=grid.nodePositions().row(n).transpose();
            if(request.handles.size()%7==0)
                h.fixedRotation=Eigen::AngleAxisd(.05,Eigen::Vector3d(1,2,-1).normalized()).toRotationMatrix();
            h.constrainedNodes={n};request.handles.push_back(h);
        }
        if(request.handles.size()<=16)throw std::runtime_error("fixture lost complete cage");
        request.fillLaplacian.resize(2,2);request.fillLaplacian<<1,-1,-1,1;
        request.fillLaplacianWeight=.03;request.biharmonicWeight=.005;request.displacementWeight=.009;
        Eigen::VectorXd beta(2);beta<<.27,.61;
        Eigen::VectorXd handles=deformation::HandleSet(request.handles).optimizationParameters();
        for(int h=0;h<int(request.handles.size());++h) {
            for(int a=0;a<3;++a)handles(4*h+a)=.012*std::sin(.4+h*.31+a*.7);
            handles(4*h+3)=1+.018*std::cos(.1+h*.29);
        }
        std::cout<<std::setprecision(17)<<"{\"evidence_class\":\"SYNTHETIC_REGRESSION\",\"bbw_success\":true,\"bbw_output_hash\":\""
                 <<bbw.outputSha256()<<"\",\"grid_output_hash\":\""<<grid.outputSha256()
                 <<"\",\"gamma_i\":0.6,\"gamma_c\":0.7,\"gamma_l\":0.03,\"gamma_dl\":0.005,\"gamma_dv\":0.009,\"outer\":";
        mesh(request.outerSurface,request.outerWeights);
        std::cout<<",\"cells\":[";
        for(int i=0;i<2;++i) {if(i)std::cout<<',';mesh(request.voidCells[i],request.voidCellWeights[i]);}
        std::cout<<"],\"pivots\":[";
        for(int i=0;i<int(request.handles.size());++i) {if(i)std::cout<<',';vector(request.handles[i].restPivot);}
        std::cout<<"],\"rotations\":[";
        for(int i=0;i<int(request.handles.size());++i) {if(i)std::cout<<',';matrix(request.handles[i].fixedRotation);}
        std::cout<<"],\"fill_laplacian\":";matrix(request.fillLaplacian);
        std::cout<<",\"cases\":[";bool first=true;
        for(auto mode:{Spin2014Objective::YoYo,Spin2014Objective::TopInertiaOnly,Spin2014Objective::TopLoweringOnly,
                      Spin2014Objective::TopMassAndLowering,Spin2014Objective::StaticLowering}) {
            if(!first)std::cout<<',';first=false;
            request.problem=SpinIt2014MomentProblem::fromObjective(mode,.6,.7);
            const CageSpinDeformationFunctional functional(request);
            const bool angle=functional.paperHasAngle();
            const auto value=functional.evaluatePaper(beta,angle?.31:0,handles);
            Eigen::VectorXd state(functional.jointParameterCount());state.head(2)=beta;
            if(angle)state(2)=.31;state.tail(handles.size())=handles;
            Eigen::VectorXd values(1+value.deformation.constraints.size());
            values(0)=value.deformation.objective;values.tail(value.deformation.constraints.size())=value.deformation.constraints;
            Eigen::MatrixXd derivatives(values.size(),state.size());derivatives.row(0)=value.jointGradient.transpose();
            derivatives.bottomRows(value.jointConstraintJacobian.rows())=value.jointConstraintJacobian;
            std::cout<<"{\"mode\":\""<<spin2014ObjectiveToken(mode)<<"\",\"state\":";vector(state);
            std::cout<<",\"values\":";vector(values);
            std::cout<<",\"moments\":";vector(value.deformation.finalMoments.values());
            std::cout<<",\"deformed_outer\":";matrix(value.deformation.deformedOuterVertices);
            std::cout<<",\"moment_derivatives\":";matrix(value.rawMomentJacobian);
            std::cout<<",\"derivatives\":";matrix(derivatives);std::cout<<'}';
        }
        std::cout<<"]}\n";
    }catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
