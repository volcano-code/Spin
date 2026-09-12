#include "spin2014/SourceCageMomentFunctional.h"
#include "source_surface_test_utils.h"
#include <Eigen/Geometry>
#include <iostream>
#include <limits>
using namespace spin2014;
using namespace surface_test;

SourceConformingBoundary fixture(bool transformed=false,bool empty=false) {
    auto outer=cube(.3,.7);
    auto design=empty?cube(.499,.501):(transformed?box({.34,.36,.38},{.66,.64,.62}):cube(.34,.66));
    const double scale=transformed?2.:1.;
    const Eigen::Vector3d shift=transformed?Eigen::Vector3d(1.25,-.75,.125):Eigen::Vector3d::Zero();
    outer.v=(outer.v*scale).rowwise()+shift.transpose();
    design.v=(design.v*scale).rowwise()+shift.transpose();
    SourceSurfaceOctreeConfig sc; sc.maximumDepth=5; sc.rootOrigin=shift; sc.rootSide=scale;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
    SourceAdaptiveDomain domain(source,2,5);
    SourceSurfaceCageConfig cc; cc.cageLevel=2; cc.rootOrigin=shift; cc.rootSide=scale;
    SourceCageGridConfig gc; gc.fineLevel=3;
    SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),{}),{});
    require(bbw.success(),"actual coupled BBW succeeds");
    return SourceConformingBoundary(SourceCageBinding(domain,bbw,{}),{});
}

Ten signedTriangle(const Eigen::MatrixXd& points,SourceConformingBoundary::Triangle t) {
    Eigen::Matrix<double,3,3> v;
    for(int i=0;i<3;++i) v.row(i)=points.row(t[i]);
    const double volume=v.determinant()/6;
    const Eigen::Vector3d sum=v.colwise().sum();
    const Eigen::Matrix3d second=volume*(sum*sum.transpose()+v.transpose()*v)/20;
    Ten result;
    result<<volume,volume*sum/4,second(0,0),second(1,1),second(2,2),second(1,2),second(0,2),second(0,1);
    return result;
}
Ten independent(const SourceConformingBoundary& b,const Eigen::MatrixXd& points,const Eigen::VectorXd& beta,double density) {
    Ten result=Ten::Zero();
    const auto& faces=b.binding().source().outerFaces();
    for(Eigen::Index i=0;i<faces.rows();++i)
        result+=signedTriangle(points,{{std::size_t(faces(i,0)),std::size_t(faces(i,1)),std::size_t(faces(i,2))}});
    for(std::size_t owner=0;owner<b.primitiveBoundaries().size();++owner)
        for(const auto& ref:b.primitiveBoundaries()[owner])
            result-=beta[b.binding().primitives()[owner].betaIndex]*ref.orientation*signedTriangle(points,b.triangles()[ref.triangle]);
    return density*result;
}
template<class Exception,class F> void throwsAs(F action,const char* label) {
    bool caught=false;
    try {action();} catch(const Exception& e) {caught=true;}
    require(caught,label);
}
void derivatives(const SourceCageMomentFunctional& functional,Eigen::VectorXd beta,Eigen::VectorXd parameters) {
    // Catches wrong beta ownership/sign or a missing full-BBW coordinate contraction.
    for(Eigen::Index h=0;h<parameters.size()/4;++h) {
        parameters[4*h]=1.+.015*std::sin(h+1);
        for(int a=0;a<3;++a) parameters[4*h+1+a]=.006*std::cos((h+1)*(a+1));
    }
    for(Eigen::Index i=0;i<beta.size();++i) beta[i]=.13+.017*i;
    const auto result=functional.evaluate(beta,parameters);
    require(result.betaJacobian.rows()==10 && result.betaJacobian.cols()==Eigen::Index(functional.betaCount()) &&
            result.handleJacobian.rows()==10 && result.handleJacobian.cols()==Eigen::Index(4*functional.handleCount()),"all derivative output columns retained");
    close(result.finalMoments.values(),independent(functional.boundary(),result.deformedPoints,beta,functional.density()),"independent real primitive tetra moments");
    const auto values=functional.evaluateValues(beta,parameters);
    require(!values.hasJacobians && values.betaJacobian.size()==0 && values.handleJacobian.size()==0,"explicit value-only absence");
    require(result.hasJacobians && result.finalMoments.values()==values.finalMoments.values() && result.deformedPoints==values.deformedPoints,"same value polynomial path");
    require(functional.evaluate(beta,parameters).outputSha256==result.outputSha256,"deterministic output identity");
    const double steps[]={3e-4,1e-4,3e-5};
    for(int kind=0;kind<2;++kind) {
        const Eigen::MatrixXd& jac=kind?result.handleJacobian:result.betaJacobian;
        for(double step:steps) {
            double worstRelative=0,worstAbsolute=0;
            for(Eigen::Index col=0;col<jac.cols();++col) {
                auto plusBeta=beta,minusBeta=beta,plusParameters=parameters,minusParameters=parameters;
                (kind?plusParameters:plusBeta)[col]+=step;
                (kind?minusParameters:minusBeta)[col]-=step;
                const Ten fd=(functional.evaluateValues(plusBeta,plusParameters).finalMoments.values()-functional.evaluateValues(minusBeta,minusParameters).finalMoments.values())/(2*step);
                for(int row=0;row<10;++row) {
                    const double error=std::abs(fd[row]-jac(row,col));
                    worstAbsolute=std::max(worstAbsolute,error);
                    // Fixed before observations: native fixture units, near-zero threshold 1e-7.
                    if(std::abs(jac(row,col))>1e-7) {
                        worstRelative=std::max(worstRelative,error/std::abs(jac(row,col)));
                        require(error/std::abs(jac(row,col))<=1e-6,"all-column FD relative tolerance");
                    } else require(error<=1e-11,"all-column FD near-zero absolute tolerance");
                }
            }
            std::cout<<"FD "<<(kind?"handles":"beta")<<" step="<<step<<" worst_relative="<<worstRelative<<" worst_absolute="<<worstAbsolute<<'\n';
        }
    }
    // Catches density omission and wrong handle origin in a joint affine motion.
    auto uniform=functional.initialHandleParameters();
    const auto& cage=functional.boundary().binding().bbw().fem().grid().cage();
    const double scale=1.07;
    const Eigen::Vector3d translation(.01,-.02,.03);
    for(std::size_t h=0;h<functional.handleCount();++h) {
        const Eigen::Vector3d center=cage.boundaryVertices().row(cage.handleVertexIndices()[h]);
        uniform[4*h]=scale;
        uniform.segment<3>(4*h+1)=(scale-1)*center+translation;
    }
    const auto rest=functional.evaluateValues(beta,functional.initialHandleParameters());
    close(functional.evaluateValues(beta,uniform).finalMoments.values(),rest.finalMoments.transformed(scale*Eigen::Matrix3d::Identity(),translation).values(),"uniform affine covariance");
    SourceCageMomentFunctional unit(functional.boundary(),1.,{});
    close(result.finalMoments.values(),functional.density()*unit.evaluateValues(beta,parameters).finalMoments.values(),"density scales all ten");
}

void invalidInputs(const SourceCageMomentFunctional& f) {
    auto beta=Eigen::VectorXd::Constant(f.betaCount(),.2).eval();
    auto p=f.initialHandleParameters();
    const double nan=std::numeric_limits<double>::quiet_NaN(),inf=std::numeric_limits<double>::infinity();
    for(double bad:{-1.,1.01,nan}) {
        auto b=beta; b[0]=bad;
        throwsAs<std::invalid_argument>([&]{f.evaluateValues(b,p);},"invalid beta rejected");
    }
    throwsAs<std::invalid_argument>([&]{f.evaluateValues(Eigen::VectorXd::Zero(1),p);},"beta dimension rejected");
    throwsAs<std::invalid_argument>([&]{f.evaluateValues(beta,Eigen::VectorXd::Zero(1));},"handle dimension rejected");
    for(int slot:{0,1}) for(double bad:{nan,inf}) {
        auto q=p; q[slot]=bad;
        throwsAs<std::invalid_argument>([&]{f.evaluateValues(beta,q);},"nonfinite handle rejected");
    }
    for(double bad:{0.,-1.}) {
        auto q=p; q[0]=bad;
        throwsAs<std::invalid_argument>([&]{f.evaluateValues(beta,q);},"nonpositive scale rejected");
    }
    for(double bad:{0.,-1.,nan,inf})
        throwsAs<std::invalid_argument>([&]{SourceCageMomentFunctional x(f.boundary(),bad,{});},"invalid density rejected");
    for(int budget=0;budget<2;++budget) for(auto bad:{std::size_t(0),std::numeric_limits<std::size_t>::max()}) {
        SourceCageMomentFunctionalConfig c;
        (budget?c.maximumJacobianValues:c.maximumWeightValues)=bad;
        throwsAs<std::length_error>([&]{SourceCageMomentFunctional x(f.boundary(),1.,c);},"malformed budget rejected");
    }
    SourceCageMomentFunctionalConfig c; c.maximumWeightValues=f.handleCount(); c.maximumJacobianValues=1;
    SourceCageMomentFunctional limited(f.boundary(),f.density(),c);
    close(limited.evaluateValues(beta,p).finalMoments.values(),f.evaluateValues(beta,p).finalMoments.values(),"one-row weight budget, no dense global cache");
    throwsAs<std::length_error>([&]{limited.evaluate(beta,p);},"Jacobian budget enforced only for derivatives");
    c.maximumWeightValues=f.handleCount()-1;
    throwsAs<std::length_error>([&]{SourceCageMomentFunctional x(f.boundary(),1.,c);},"all handles row budget cannot be reduced");
    auto huge=p; huge[1]=std::numeric_limits<double>::max();
    throwsAs<std::overflow_error>([&]{f.evaluateValues(beta,huge);},"finite input arithmetic overflow");
    auto collapsed=p;
    const auto& cage=f.boundary().binding().bbw().fem().grid().cage();
    // Source API arithmetic categories must not inherit the legacy surface
    // normal-overflow classification. Check normal and later polynomial cases.
    for(double scale:{1e200,1e100}) {
        auto scaled=p;
        for(std::size_t h=0;h<f.handleCount();++h) {
            scaled[4*h]=scale;
            scaled.segment<3>(4*h+1)=(scale-1)*cage.boundaryVertices().row(cage.handleVertexIndices()[h]).transpose();
        }
        throwsAs<std::overflow_error>([&]{f.evaluateValues(beta,scaled);},"source value normal/polynomial overflow stays arithmetic");
        throwsAs<std::overflow_error>([&]{f.evaluate(beta,scaled);},"source derivative normal/polynomial overflow stays arithmetic");
    }
    for(std::size_t h=0;h<f.handleCount();++h) {
        collapsed[4*h]=std::numeric_limits<double>::denorm_min();
        collapsed.segment<3>(4*h+1)=-cage.boundaryVertices().row(cage.handleVertexIndices()[h]).transpose();
    }
    throwsAs<std::invalid_argument>([&]{f.evaluateValues(Eigen::VectorXd::Zero(beta.size()),collapsed);},"finite collapsed triangles rejected even at zero beta");
    auto reversed=p;
    for(std::size_t h=0;h<f.handleCount();++h)
        reversed.segment<3>(4*h+1)=-4*cage.boundaryVertices().row(cage.handleVertexIndices()[h]).transpose();
    bool nonpositiveMass=false;
    try {f.evaluateValues(Eigen::VectorXd::Zero(beta.size()),reversed);}
    catch(const std::invalid_argument& e) {nonpositiveMass=std::string(e.what()).find("final mass")!=std::string::npos;}
    require(nonpositiveMass,"actual full-BBW orientation reversal rejects nonpositive final mass");
    const auto baseline=f.evaluateValues(beta,p);
    p[1]+=1e-4;
    require(f.evaluateValues(beta,p).outputSha256!=baseline.outputSha256,"input state bound in identity");
    require(limited.evaluateValues(beta,f.initialHandleParameters()).outputSha256!=baseline.outputSha256,"budgets bound in identity");
}

int main() {
    try {
        // Catches substituting an adaptive bounding box or dropping actual handles.
        const auto boundary=fixture();
        SourceCageMomentFunctional functional(boundary,2.3,{});
        require(functional.betaCount()==8 && functional.handleCount()==26 && boundary.binding().primitives().size()==496,"actual fixture counts");
        auto beta=Eigen::VectorXd::Zero(functional.betaCount()).eval();
        const auto parameters=functional.initialHandleParameters();
        auto result=functional.evaluate(beta,parameters);
        close(result.finalMoments.values(),2.3*boxMoments(Eigen::Vector3d::Constant(.3),Eigen::Vector3d::Constant(.7)),"rest outer all ten moments");
        beta.setOnes();
        result=functional.evaluate(beta,parameters);
        close(result.finalMoments.values(),2.3*(boxMoments(Eigen::Vector3d::Constant(.3),Eigen::Vector3d::Constant(.7))-boxMoments(Eigen::Vector3d::Constant(11./32),Eigen::Vector3d::Constant(21./32))),"actual accepted interior subtraction");
        derivatives(functional,beta,parameters);
        invalidInputs(functional);
        const SourceCageMomentFunctional transformed(fixture(true),.8,{});
        derivatives(transformed,Eigen::VectorXd::Zero(transformed.betaCount()),transformed.initialHandleParameters());
        const SourceCageMomentFunctional empty(fixture(false,true),1.,{});
        require(empty.boundary().binding().primitives().empty(),"actual empty carveable binding");
        const auto emptyResult=empty.evaluate(Eigen::VectorXd::Zero(empty.betaCount()),empty.initialHandleParameters());
        close(emptyResult.finalMoments.values(),boxMoments(Eigen::Vector3d::Constant(.3),Eigen::Vector3d::Constant(.7)),"empty source outer material");
        require(emptyResult.betaJacobian.isZero(),"empty source no invented cavity derivatives");
        std::cout<<"source cage moments PASS\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n'; return 1;
    }
}
