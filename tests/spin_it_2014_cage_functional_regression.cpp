#include "spin2014/CageSpinDeformation.h"

#include "deformation/HandleSet.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>


namespace
{


spin::TriangleMesh tetra(double scale, const Eigen::Vector3d& offset)
{
    spin::TriangleMesh mesh;
    mesh.vertices.resize(4, 3);
    mesh.vertices <<
        0.0, 0.0, 0.0,
        scale, 0.0, 0.0,
        0.0, scale, 0.0,
        0.0, 0.0, scale;
    mesh.vertices.rowwise() += offset.transpose();
    mesh.faces.resize(4, 3);
    mesh.faces << 1, 2, 3, 0, 2, 1, 0, 1, 3, 0, 3, 2;
    return mesh;
}


deformation::HandleVector handles()
{
    deformation::Handle first;
    first.restPivot = Eigen::Vector3d::Zero();
    first.uniformScale = 1.0;
    first.constrainedNodes = {0};
    deformation::Handle second;
    second.restPivot = Eigen::Vector3d(2.0, 2.0, 2.0);
    second.uniformScale = 1.0;
    second.constrainedNodes = {1};
    return {first, second};
}


Eigen::MatrixXd weights(const Eigen::MatrixXd& vertices)
{
    Eigen::MatrixXd result(vertices.rows(), 2);
    for (Eigen::Index row = 0; row < vertices.rows(); ++row)
    {
        const double blend = std::clamp(
            0.15 + 0.1 * vertices.row(row).sum(), 0.05, 0.95
        );
        result(row, 0) = 1.0 - blend;
        result(row, 1) = blend;
    }
    return result;
}


} // namespace

namespace {
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
void near(double actual, double expected, double absolute=2e-10, double relative=1e-10) {
    if (!(std::isfinite(actual) && std::abs(actual-expected)<=absolute+relative*std::abs(expected)))
        std::cerr<<"actual="<<actual<<" expected="<<expected<<" difference="<<actual-expected<<'\n';
    require(std::isfinite(actual) && std::abs(actual-expected)<=absolute+relative*std::abs(expected),
            "independent value or derivative mismatch");
}
// Uniform tetrahedron expectations, independent of the surface differential.
Eigen::Matrix<double,10,1> tetraMoments(const Eigen::MatrixXd& v) {
    const Eigen::Vector3d e1=v.row(1)-v.row(0),e2=v.row(2)-v.row(0),e3=v.row(3)-v.row(0);
    const double volume=e1.dot(e2.cross(e3))/6.;
    const Eigen::Vector3d sum=v.colwise().sum();
    const Eigen::Matrix3d second=volume*(v.transpose()*v+sum*sum.transpose())/20.;
    Eigen::Matrix<double,10,1> r;
    r << volume, (volume*sum/4.).eval(), second(0,0),second(1,1),second(2,2),
        second(1,2),second(0,2),second(0,1);
    return r;
}
spin2014::CageSpinDeformationRequest paperFixture() {
    spin2014::CageSpinDeformationRequest r;
    r.outerSurface=tetra(3.,{-1.,-.8,.2});
    r.outerSurface.vertices(1,1)+=.35;
    r.outerSurface.vertices(2,2)+=.2;
    r.voidCells={tetra(.25,{-.5,-.3,.7}),tetra(.32,{.05,-.2,.8})};
    r.handles=handles();
    r.handles[0].fixedRotation=Eigen::AngleAxisd(.13,Eigen::Vector3d(1,2,3).normalized()).toRotationMatrix();
    deformation::Handle locked; locked.locked=true; locked.restPivot={.2,.3,.4};
    locked.translation={.01,-.02,.03}; locked.constrainedNodes={2};
    r.handles.push_back(locked);
    auto bind=[](const Eigen::MatrixXd& v) {
        Eigen::MatrixXd w(v.rows(),3); w.leftCols(2)=.9*weights(v); w.col(2).setConstant(.1); return w;
    };
    r.outerWeights=bind(r.outerSurface.vertices);
    for (const auto& cell:r.voidCells) r.voidCellWeights.push_back(bind(cell.vertices));
    r.fillLaplacian.resize(2,2); r.fillLaplacian<<1,-1,-1,1;
    r.fillLaplacianWeight=.17; r.biharmonicWeight=.023; r.displacementWeight=.047;
    return r;
}
void paperRegression() {
    using O=spin2014::Spin2014Objective;
    std::vector<double> topValues;
    for (auto mode:{O::YoYo,O::TopInertiaOnly,O::TopLoweringOnly,O::TopMassAndLowering,O::StaticLowering}) {
        auto r=paperFixture(); r.problem=spin2014::SpinIt2014MomentProblem::fromObjective(mode,.73,.41);
        spin2014::CageSpinDeformationFunctional f(r);
        Eigen::Vector2d beta(.27,.68); auto q=f.initialHandleParameters();
        q(0)=.024; q(3)=.96; q(5)=-.031; q(7)=1.035;
        const double phi=f.paperHasAngle()?.37:0.;
        const auto p=f.evaluatePaper(beta,phi,q); const auto& d=p.deformation;
        require(q.size()==8 && f.betaCount()==2 && f.jointParameterCount()==10+int(f.paperHasAngle()),"all free columns / locked handle contract");
        Eigen::VectorXd x(f.jointParameterCount()); x.head(2)=beta;
        if(f.paperHasAngle()) x(2)=phi;
        x.tail(q.size())=q;
        near(f.evaluatePaperJoint(x).deformation.objective,d.objective);
        // Reconstruct LBS without the production free-parameter model.
        auto lbs=[&](const Eigen::MatrixXd& v,const Eigen::MatrixXd& w) {
            Eigen::MatrixXd out=Eigen::MatrixXd::Zero(v.rows(),3);
            for(int i=0;i<v.rows();++i) for(int h=0;h<3;++h) {
                const auto& handle=r.handles[h];
                const Eigen::Vector3d t=h==2?handle.translation:Eigen::Vector3d(q.segment<3>(4*h));
                const double scale=h==2?handle.uniformScale:q(4*h+3);
                out.row(i)+=w(i,h)*(scale*handle.fixedRotation*(v.row(i).transpose()-handle.restPivot)+handle.restPivot+t).transpose();
            }
            return out;
        };
        const auto outer=lbs(r.outerSurface.vertices,r.outerWeights);
        require((outer-d.deformedOuterVertices).norm()<1e-12,"actual returned outer LBS geometry");
        auto raw=tetraMoments(outer);
        for(int k=0;k<2;++k) {
            const auto cell=lbs(r.voidCells[k].vertices,r.voidCellWeights[k]);
            require((cell-d.deformedVoidCellVertices[k]).norm()<1e-12,"actual returned cell LBS geometry");
            raw-=beta(k)*tetraMoments(cell);
            require((d.momentFillJacobian.col(k)+tetraMoments(cell)).norm()<1e-12,"negative cell raw beta column");
        }
        require((raw-d.finalMoments.values()).norm()<1e-12 && std::abs(raw(9))>1e-3,"independent nonzero cross raw moments");
        const double c=std::cos(phi),s=std::sin(phi),delta=mode==O::YoYo?0:raw(3)*raw(3)/raw(0);
        const double a=c*c*(raw(5)+raw(6))+s*s*(raw(4)+raw(6))+2*c*s*raw(9)-delta;
        const double b=s*s*(raw(5)+raw(6))+c*c*(raw(4)+raw(6))-2*c*s*raw(9)-delta;
        double expected=f.paperHasAngle()?.73*(a*a+b*b)/std::pow(raw(4)+raw(5),2):0.;
        if(mode==O::TopMassAndLowering) expected+=.41*raw(3)*raw(3);
        if(mode==O::TopLoweringOnly || mode==O::StaticLowering) expected+=.41*std::pow(raw(3)/raw(0),2);
        near(d.spinObjective,expected);
        if(mode==O::TopInertiaOnly || mode==O::TopLoweringOnly || mode==O::TopMassAndLowering) topValues.push_back(d.spinObjective);
        std::vector<double> constraints={raw(1),raw(2)};
        if(f.paperHasAngle()) {
            if(mode==O::YoYo) constraints.push_back(raw(3));
            constraints.push_back(raw(8)); constraints.push_back(raw(7));
            constraints.push_back(c*s*(raw(4)-raw(5))+(c*c-s*s)*raw(9));
        }
        require(d.constraints.size()==int(constraints.size()),"all mode-specific equalities");
        for(int i=0;i<d.constraints.size();++i) near(d.constraints(i),constraints[i]);
        Eigen::Matrix4d L=4*Eigen::Matrix4d::Identity()-Eigen::Matrix4d::Ones();
        near(d.shapeEnergy,.023*(L*outer).squaredNorm());
        near(d.displacementEnergy,.047*(outer-r.outerSurface.vertices).squaredNorm());
        near(d.fillRegularizationEnergy,.5*.17*std::pow(beta(0)-beta(1),2));
        near(d.objective,expected+d.shapeEnergy+d.displacementEnergy+d.fillRegularizationEnergy);
        if(f.paperHasAngle()) {
            require(p.rawMomentJacobian.col(2).norm()==0,"phi independent of raw moments");
            // Legacy centers all inertia entries; paper uses the centroid
            // equalities. Compare the objectives on that common manifold.
            auto centered=r;
            centered.outerSurface.vertices=d.deformedOuterVertices;
            for(int k=0;k<2;++k) centered.voidCells[k].vertices=d.deformedVoidCellVertices[k];
            Eigen::Vector3d shift=-raw.segment<3>(1)/raw(0);
            if(mode!=O::YoYo) shift.z()=0.;
            for(auto& h:centered.handles) {
                h.fixedRotation.setIdentity(); h.uniformScale=1.; h.translation=shift;
            }
            spin2014::CageSpinDeformationFunctional cf(centered);
            const auto cq=cf.initialHandleParameters();
            const auto cr=cf.evaluatePaper(beta,0.,cq).deformation.finalMoments.values();
            near(cr(1),0.); near(cr(2),0.); if(mode==O::YoYo) near(cr(3),0.);
            const double diagonalPhi=.5*std::atan2(-2*cr(9),cr(4)-cr(5));
            near(cf.evaluatePaper(beta,diagonalPhi,cq).deformation.objective,cf.evaluate(beta,cq).objective);
        } else near(d.objective,f.evaluate(beta,q).objective);
        std::vector<double> errors;
        for(double step:{1e-3,1e-4,1e-5}) {
            double error=0.;
            for(int j=0;j<x.size();++j) {
                auto plus=x,minus=x; plus(j)+=step; minus(j)-=step;
                auto ep=f.evaluatePaperJoint(plus),em=f.evaluatePaperJoint(minus);
                const double fd=(ep.deformation.objective-em.deformation.objective)/(2*step);
                error=std::max(error,std::abs(fd-p.jointGradient(j))/(1+std::abs(p.jointGradient(j))));
                if(step==1e-5) near(fd,p.jointGradient(j),2e-8,1e-6);
                for(int row=0;row<d.constraints.size();++row) {
                    const double df=(ep.deformation.constraints(row)-em.deformation.constraints(row))/(2*step);
                    error=std::max(error,std::abs(df-p.jointConstraintJacobian(row,j))/(1+std::abs(p.jointConstraintJacobian(row,j))));
                    if(step==1e-5) near(df,p.jointConstraintJacobian(row,j),2e-8,1e-6);
                }
                const auto dm=(ep.deformation.finalMoments.values()-em.deformation.finalMoments.values())/(2*step);
                if(step==1e-5) for(int row=0;row<10;++row) near(dm(row),p.rawMomentJacobian(row,j),2e-8,1e-6);
            }
            errors.push_back(error);
            std::cout<<"paper mode="<<int(mode)<<" step="<<step<<" relative error="<<error<<'\n';
        }
        require(errors[1]<.1*errors[0] && errors[2]<errors[1],"central differences converge before noise");
        auto rejects=[](const auto& call) { bool caught=false; try{call();}catch(const std::invalid_argument&){caught=true;} require(caught,"malformed state explicitly rejected"); };
        rejects([&]{f.evaluatePaper(Eigen::VectorXd::Zero(1),phi,q);});
        rejects([&]{f.evaluatePaper(beta,phi,Eigen::VectorXd::Zero(7));});
        rejects([&]{f.evaluatePaperJoint(Eigen::VectorXd::Zero(x.size()+1));});
        for(double bad:{-.01,1.01}) {
            auto badBeta=beta; badBeta(0)=bad;
            rejects([&]{f.evaluatePaper(badBeta,phi,q);});
        }
        for(double bad:{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}) {
            rejects([&]{f.evaluatePaper(beta,bad,q);});
            auto badBeta=beta; badBeta(0)=bad; rejects([&]{f.evaluatePaper(badBeta,phi,q);});
            auto badQ=q; badQ(1)=bad; rejects([&]{f.evaluatePaper(beta,phi,badQ);});
            auto badX=x; badX(0)=bad; rejects([&]{f.evaluatePaperJoint(badX);});
        }
        if(!f.paperHasAngle()) rejects([&]{f.evaluatePaper(beta,.01,q);});
    }
    require(std::abs(topValues[0]-topValues[1])>.01 && std::abs(topValues[1]-topValues[2])>.01,"tea-pot objectives have distinct values");
    auto locked=paperFixture();
    for(auto& h:locked.handles) h.locked=true;
    spin2014::CageSpinDeformationFunctional lf(locked);
    Eigen::VectorXd x(3); x<<.27,.68,.37;
    const auto lp=lf.evaluatePaperJoint(x);
    require(lf.initialHandleParameters().size()==0 && lp.deformation.handleGradient.size()==0 &&
            lp.jointGradient.size()==3 && lp.jointConstraintJacobian.cols()==3,"zero free handles retain beta/phi layout");
    auto huge=paperFixture(); huge.biharmonicWeight=std::numeric_limits<double>::max();
    spin2014::CageSpinDeformationFunctional hf(huge);
    bool overflowRejected=false;
    try { hf.evaluatePaper(Eigen::Vector2d(.27,.68),.37,hf.initialHandleParameters()); }
    catch(const std::invalid_argument&) { overflowRejected=true; }
    require(overflowRejected,"finite inputs with nonfinite computed output rejected");
}
}

int main()
{
    try { paperRegression(); }
    catch(const std::exception& error) { std::cerr<<"[FAIL] explicit cage: "<<error.what()<<'\n'; return EXIT_FAILURE; }
    spin2014::CageSpinDeformationRequest request;
    request.outerSurface = tetra(3.0, Eigen::Vector3d(-1.0, -1.0, 0.0));
    request.voidCells = {tetra(0.3, Eigen::Vector3d(-0.1, -0.1, 0.5))};
    request.outerWeights = weights(request.outerSurface.vertices);
    request.voidCellWeights = {weights(request.voidCells.front().vertices)};
    request.handles = handles();
    request.problem = spin2014::SpinIt2014MomentProblem::staticBalance(0.7);
    request.fillLaplacian = Eigen::MatrixXd::Zero(1, 1);
    request.fillLaplacianWeight = 0.1;
    request.biharmonicWeight = 0.03;
    request.displacementWeight = 0.04;
    const spin2014::CageSpinDeformationFunctional functional(request);

    Eigen::VectorXd beta(1);
    beta << 0.37;
    Eigen::VectorXd q = functional.initialHandleParameters();
    q(0) = 0.02;
    q(5) = -0.03;
    q(3) = 0.97;
    q(7) = 1.04;
    const auto value = functional.evaluate(beta, q);
    const double step = 1e-6;
    double gradientError = 0.0;
    double jacobianError = 0.0;

    for (Eigen::Index index = 0; index < beta.size(); ++index)
    {
        Eigen::VectorXd plus = beta;
        Eigen::VectorXd minus = beta;
        plus(index) += step;
        minus(index) -= step;
        const auto p = functional.evaluate(plus, q);
        const auto m = functional.evaluate(minus, q);
        gradientError = std::max(
            gradientError,
            std::abs((p.objective - m.objective) / (2.0 * step) -
                     value.fillGradient(index))
        );
        jacobianError = std::max(
            jacobianError,
            ((p.constraints - m.constraints) / (2.0 * step) -
             value.fillConstraintJacobian.col(index)).lpNorm<Eigen::Infinity>()
        );
    }
    for (Eigen::Index index = 0; index < q.size(); ++index)
    {
        Eigen::VectorXd plus = q;
        Eigen::VectorXd minus = q;
        plus(index) += step;
        minus(index) -= step;
        const auto p = functional.evaluate(beta, plus);
        const auto m = functional.evaluate(beta, minus);
        gradientError = std::max(
            gradientError,
            std::abs((p.objective - m.objective) / (2.0 * step) -
                     value.handleGradient(index))
        );
        jacobianError = std::max(
            jacobianError,
            ((p.constraints - m.constraints) / (2.0 * step) -
             value.handleConstraintJacobian.col(index)).lpNorm<Eigen::Infinity>()
        );
    }
    if (
        gradientError > 2e-7 || jacobianError > 2e-8 ||
        value.shapeEnergy <= 0.0 || value.displacementEnergy <= 0.0
    )
    {
        std::cerr << "[FAIL] scenario=spin-it-2014-cage-functional gradient="
                  << gradientError << " jacobian=" << jacobianError
                  << " shape=" << value.shapeEnergy << '\n';
        return EXIT_FAILURE;
    }
    std::cout << "[PASS] scenario=spin-it-2014-cage-functional gradient="
              << gradientError << " jacobian=" << jacobianError << '\n';
    return EXIT_SUCCESS;
}
