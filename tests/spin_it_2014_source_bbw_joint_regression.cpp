#include "spin2014/CageSpinDeformation.h"
#include "spin2014/SourceCageBbw.h"
#include "source_surface_test_utils.h"
#include <Eigen/Geometry>
#include <iostream>
#include <set>

using namespace spin2014;
using surface_test::require;

namespace {
void near(double a,double b) {
    require(std::isfinite(a) && std::abs(a-b)<=2e-9+1e-6*std::abs(b),"joint derivative or value mismatch");
}
// Oriented tetrahedra (origin, face vertices), with exact polynomial simplex
// averages; this does not call the production raw-moment differential.
surface_test::Ten moments(const Eigen::MatrixXd& v,const Eigen::MatrixXi& f) {
    surface_test::Ten out=surface_test::Ten::Zero();
    for(int i=0;i<f.rows();++i) {
        Eigen::Matrix3d t;
        for(int j=0;j<3;++j) t.row(j)=v.row(f(i,j));
        const double volume=t.determinant()/6.;
        const Eigen::Vector3d sum=t.colwise().sum();
        const Eigen::Matrix3d q=volume*(t.transpose()*t+sum*sum.transpose())/20.;
        surface_test::Ten r;
        r<<volume,volume*sum/4.,q(0,0),q(1,1),q(2,2),q(1,2),q(0,2),q(0,1);
        out+=r;
    }
    return out;
}
spin::TriangleMesh mesh(const surface_test::Mesh& m) { return {m.v,m.f}; }
}
int main() {
    try {
        const auto original=surface_test::box({.30,.31,.32},{.69,.66,.70});
        SourceSurfaceCageConfig cc; cc.cageLevel=2;
        const SourceSurfaceCage cage(original.v,original.f,cc);
        SourceCageGridConfig gc; gc.fineLevel=3;
        const SourceCageGrid grid(cage,gc);
        const SourceCageFem fem(grid,{});
        const SourceCageBbw bbw(fem);
        require(bbw.success() && bbw.finalAssessment().primalPassed && bbw.finalAssessment().optimalityPassed,"actual coupled source BBW solve succeeds");
        require(!bbw.attempts().empty() && bbw.attempts()[0].runInvoked && bbw.attempts()[0].qpIterations>0,"real backend solve");
        CageSpinDeformationRequest r;
        r.outerSurface=mesh(original);
        r.voidCells={mesh(surface_test::box({.36,.37,.38},{.40,.42,.43})),
                     mesh(surface_test::box({.53,.51,.54},{.59,.57,.60}))};
        r.outerWeights=grid.interpolate(bbw.finalWeights(),r.outerSurface.vertices);
        for(const auto& cell:r.voidCells) r.voidCellWeights.push_back(grid.interpolate(bbw.finalWeights(),cell.vertices));
        for(int index:cage.handleVertexIndices()) {
            deformation::Handle h; h.restPivot=cage.boundaryVertices().row(index);
            h.constrainedNodes={grid.handleNodeIndices().at(r.handles.size())};
            r.handles.push_back(h);
        }
        const int count=r.handles.size();
        require(count>16 && count==int(grid.handleNodeIndices().size()) && r.outerWeights.cols()==count,"all cage handles, no sampling");
        r.fillLaplacian.resize(2,2); r.fillLaplacian<<1,-1,-1,1;
        r.fillLaplacianWeight=.11; r.biharmonicWeight=.007; r.displacementWeight=.013;
        r.problem=SpinIt2014MomentProblem::fromObjective(Spin2014Objective::TopMassAndLowering,.83,.29);
        CageSpinDeformationFunctional f(r);
        Eigen::VectorXd x(f.jointParameterCount()); x.head<3>()<<.23,.64,.31;
        x.tail(4*count)=f.initialHandleParameters();
        for(int h=0;h<count;++h) {
            x(3+4*h)+=.002*std::sin(h+1.);
            x(4+4*h)+=.001*std::cos(2*h+1.);
            x(5+4*h)+=.0015*std::sin(3*h+1.);
            x(6+4*h)+=.01*std::cos(h+1.);
        }
        const auto p=f.evaluatePaperJoint(x); const auto& d=p.deformation;
        require(x.size()==3+4*count && p.jointGradient.size()==x.size() && p.jointConstraintJacobian.cols()==x.size() && p.rawMomentJacobian.cols()==x.size(),"every free handle column present");
        auto lbs=[&](const Eigen::MatrixXd& v,const Eigen::MatrixXd& w) {
            Eigen::MatrixXd out=Eigen::MatrixXd::Zero(v.rows(),3);
            for(int i=0;i<v.rows();++i) for(int h=0;h<count;++h)
                out.row(i)+=w(i,h)*(x(6+4*h)*(v.row(i).transpose()-r.handles[h].restPivot)+r.handles[h].restPivot+x.segment<3>(3+4*h)).transpose();
            return out;
        };
        const auto outer=lbs(r.outerSurface.vertices,r.outerWeights);
        require((outer-d.deformedOuterVertices).norm()<1e-12,"independent outer LBS geometry");
        auto raw=moments(outer,r.outerSurface.faces);
        for(int k=0;k<2;++k) {
            const auto v=lbs(r.voidCells[k].vertices,r.voidCellWeights[k]);
            require((v-d.deformedVoidCellVertices[k]).norm()<1e-12,"independent void LBS geometry");
            raw-=x(k)*moments(v,r.voidCells[k].faces);
        }
        require((raw-d.finalMoments.values()).norm()<1e-12,"independent signed tetrahedra moments");
        const double c=std::cos(x(2)),s=std::sin(x(2)),delta=raw(3)*raw(3)/raw(0);
        const double a=c*c*(raw(5)+raw(6))+s*s*(raw(4)+raw(6))+2*c*s*raw(9)-delta;
        const double b=s*s*(raw(5)+raw(6))+c*c*(raw(4)+raw(6))-2*c*s*raw(9)-delta;
        near(d.spinObjective,.83*(a*a+b*b)/std::pow(raw(4)+raw(5),2)+.29*raw(3)*raw(3));
        Eigen::VectorXd expected(5); expected<<raw(1),raw(2),raw(8),raw(7),c*s*(raw(4)-raw(5))+(c*c-s*s)*raw(9);
        require((expected-d.constraints).norm()<1e-12,"independent all explicit equalities");
        Eigen::MatrixXd L=Eigen::MatrixXd::Zero(8,8);
        std::set<std::pair<int,int>> edges;
        for(int i=0;i<original.f.rows();++i) for(int j=0;j<3;++j) {
            int a=original.f(i,j),b=original.f(i,(j+1)%3); if(a>b)std::swap(a,b); edges.emplace(a,b);
        }
        for(auto e:edges) { L(e.first,e.first)++; L(e.second,e.second)++; L(e.first,e.second)--; L(e.second,e.first)--; }
        near(d.shapeEnergy,.007*(L*outer).squaredNorm());
        near(d.displacementEnergy,.013*(outer-original.v).squaredNorm());
        near(d.fillRegularizationEnergy,.5*.11*std::pow(x(0)-x(1),2));
        near(d.objective,d.spinObjective+d.shapeEnergy+d.displacementEnergy+d.fillRegularizationEnergy);
        std::vector<double> errors;
        for(double step:{1e-3,1e-4,1e-5}) {
            double error=0.;
            for(int j=0;j<x.size();++j) {
                auto xp=x,xm=x; xp(j)+=step; xm(j)-=step;
                const auto ep=f.evaluatePaperJoint(xp),em=f.evaluatePaperJoint(xm);
                const double fd=(ep.deformation.objective-em.deformation.objective)/(2*step);
                error=std::max(error,std::abs(fd-p.jointGradient(j))/(1+std::abs(p.jointGradient(j))));
                if(step==1e-5) near(fd,p.jointGradient(j));
                for(int row=0;row<d.constraints.size();++row) {
                    const double dc=(ep.deformation.constraints(row)-em.deformation.constraints(row))/(2*step);
                    error=std::max(error,std::abs(dc-p.jointConstraintJacobian(row,j))/(1+std::abs(p.jointConstraintJacobian(row,j))));
                    if(step==1e-5) near(dc,p.jointConstraintJacobian(row,j));
                }
                const auto dm=(ep.deformation.finalMoments.values()-em.deformation.finalMoments.values())/(2*step);
                if(step==1e-5) for(int row=0;row<10;++row) near(dm(row),p.rawMomentJacobian(row,j));
            }
            errors.push_back(error);
            std::cout<<"source BBW joint step="<<step<<" relative error="<<error<<'\n';
        }
        require(errors[1]<.1*errors[0] && errors[2]<errors[1],"all-column central difference convergence");
        std::cout<<"PASS mesh-based integration handles="<<count<<" jointColumns="<<x.size()<<" BBW="<<bbw.outputSha256()<<'\n';
        return 0;
    } catch(const std::exception& e) { std::cerr<<"FAIL source BBW joint: "<<e.what()<<'\n'; return 1; }
}
