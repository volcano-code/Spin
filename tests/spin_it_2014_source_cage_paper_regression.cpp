#include "spin2014/SourceCagePaperFunctional.h"
#include "source_surface_test_utils.h"
#include <Eigen/Geometry>
#include <complex>
#include <iostream>
#include <limits>
#include <type_traits>
using namespace spin2014;
using namespace surface_test;

// Actual Task5q source/BBW fixtures; no replacement mass, weights or Jacobians.
SourceConformingBoundary fixture(bool transformed) {
    auto outer=cube(.3,.7);
    auto design=transformed?box({.34,.36,.38},{.66,.64,.62}):cube(.34,.66);
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
    require(bbw.success(),"actual BBW fixture");
    return SourceConformingBoundary(SourceCageBinding(domain,bbw,{}),{});
}

Ten tetra(const Eigen::MatrixXd& points,SourceConformingBoundary::Triangle t) {
    Eigen::Matrix3d v;
    for(int i=0;i<3;++i) v.row(i)=points.row(t[i]);
    const double volume=v.determinant()/6.;
    const Eigen::Vector3d sum=v.colwise().sum();
    const Eigen::Matrix3d second=volume*(sum*sum.transpose()+v.transpose()*v)/20.;
    Ten r; r<<volume,volume*sum/4.,second(0,0),second(1,1),second(2,2),second(1,2),second(0,2),second(0,1);
    return r;
}
Ten integrated(const SourceConformingBoundary& b,const Eigen::MatrixXd& points,const Eigen::VectorXd& beta,double density) {
    Ten r=Ten::Zero();
    const auto& faces=b.binding().source().outerFaces();
    for(Eigen::Index i=0;i<faces.rows();++i)
        r+=tetra(points,{{std::size_t(faces(i,0)),std::size_t(faces(i,1)),std::size_t(faces(i,2))}});
    for(std::size_t owner=0;owner<b.primitiveBoundaries().size();++owner)
        for(const auto& ref:b.primitiveBoundaries()[owner])
            r-=beta[b.binding().primitives()[owner].betaIndex]*ref.orientation*tetra(points,b.triangles()[ref.triangle]);
    return density*r;
}
template<class S> struct Expected {S objective; Eigen::Matrix<S,Eigen::Dynamic,1> constraints;};
template<class S> Expected<S> independentPaper(const Eigen::Matrix<S,10,1>& m,S phi,const SourceCagePaperFunctionalConfig& c) {
    // Build and rotate the inertia tensor, independently of production Givens.
    Eigen::Matrix<S,3,3> q;
    q<<m[4],m[9],m[8],m[9],m[5],m[7],m[8],m[7],m[6];
    Eigen::Matrix<S,3,3> inertia=q.trace()*Eigen::Matrix<S,3,3>::Identity()-q;
    const Eigen::Matrix<S,3,1> u(std::cos(phi),-std::sin(phi),0),v(std::sin(phi),std::cos(phi),0);
    Expected<S> r; r.objective=0;
    if(c.objective==Spin2014Objective::StaticLowering) {
        r.constraints.resize(2); r.constraints<<m[1],m[2];
    } else {
        const S shift=c.objective==Spin2014Objective::YoYo?S(0):m[3]*m[3]/m[0];
        const S a=(u.transpose()*inertia*u)[0]-shift,b=(v.transpose()*inertia*v)[0]-shift,z=inertia(2,2);
        r.objective=c.inertiaWeight*(a*a+b*b)/(z*z);
        if(c.objective==Spin2014Objective::YoYo) {
            r.constraints.resize(6);r.constraints<<m[1],m[2],m[3],m[8],m[7],-(u.transpose()*inertia*v)[0];
        } else {r.constraints.resize(5);r.constraints<<m[1],m[2],m[8],m[7],-(u.transpose()*inertia*v)[0];}
    }
    if(c.objective==Spin2014Objective::TopMassAndLowering) r.objective+=c.centerWeight*m[3]*m[3];
    if(c.objective==Spin2014Objective::StaticLowering || c.objective==Spin2014Objective::TopLoweringOnly)
        r.objective+=c.centerWeight*std::pow(m[3]/m[0],2);
    return r;
}
template<class Exception,class F> void throwsAs(F f,const char* label) {
    bool caught=false;try {f();}catch(const Exception& e) {
        caught=true;
        if constexpr(std::is_same_v<Exception,std::length_error>)require(std::string(e.what()).find("RESOURCE_LIMIT")!=std::string::npos,"resource failure category token");
    }require(caught,label);
}
void near(double a,double b,const char* label,double tolerance=2e-10) {
    require(std::abs(a-b)<=tolerance*(1+std::abs(b)),label);
}

using Complex=std::complex<double>;
using ComplexTen=Eigen::Matrix<Complex,10,1>;
using ComplexPoints=Eigen::Matrix<Complex,Eigen::Dynamic,3>;
template<class S> Eigen::Matrix<S,10,1> referenceTetra(const Eigen::Matrix<S,Eigen::Dynamic,3>& points,SourceConformingBoundary::Triangle t) {
    Eigen::Matrix<S,3,3> v;for(int i=0;i<3;++i)v.row(i)=points.row(t[i]);
    const S volume=v.determinant()/6.;
    const Eigen::Matrix<S,3,1> sum=v.colwise().sum();
    const Eigen::Matrix<S,3,3> second=volume*(sum*sum.transpose()+v.transpose()*v)/20.;
    Eigen::Matrix<S,10,1> r;r<<volume,volume*sum/4.,second(0,0),second(1,1),second(2,2),second(1,2),second(0,2),second(0,1);return r;
}
Eigen::MatrixXd complexDerivatives(const SourceCagePaperFunctional& f,const Eigen::VectorXd& x,const SourceCagePaperEvaluation& result) {
    // Independent complex step on transformed POINTS and signed tetra integrals.
    // Only actual rest embedding and final BBW nodal data are reused; no source
    // moment Jacobian, production frame map or production Givens calls enter.
    const auto& boundary=f.source().boundary();const auto& bbw=boundary.binding().bbw();
    const auto& grid=bbw.fem().grid();const auto& cage=grid.cage();
    constexpr double step=1e-30;
    Eigen::MatrixXd jac(11+f.constraintCount(),x.size());
    for(Eigen::Index col=0;col<x.size();++col) {
        auto beta=x.head(f.betaCount()).cast<Complex>().eval();
        Complex phi=f.hasAngle()?x[f.angleOffset()]:0.;
        ComplexPoints points=result.world.deformedPoints.cast<Complex>();
        if(col<Eigen::Index(f.betaCount()))beta[col]+=Complex(0,step);
        else if(f.hasAngle() && col==Eigen::Index(f.angleOffset()))phi+=Complex(0,step);
        else {
            const std::size_t h=(col-f.handleOffset())/4,slot=(col-f.handleOffset())%4;
            const Eigen::Vector3d center=cage.boundaryVertices().row(cage.handleVertexIndices()[h]);
            for(Eigen::Index p=0;p<points.rows();++p) {
                const Eigen::Vector3d rest=boundary.points().row(p);
                const auto embedding=grid.embed(rest);double w=0;
                for(int node=0;node<8;++node)w+=embedding.coefficients[node]*bbw.finalWeights()(embedding.nodeIndices[node],h);
                const Eigen::Vector3d velocity=slot?Eigen::Vector3d::Unit(slot-1).eval():(rest-center).eval();
                for(int a=0;a<3;++a)points(p,a)+=Complex(0,step*w*velocity[a]);
            }
        }
        points=((points.rowwise()-f.frame().origin.cast<Complex>().transpose())*f.frame().rotation.transpose().cast<Complex>()).eval();
        ComplexTen moments=ComplexTen::Zero();
        const auto& faces=boundary.binding().source().outerFaces();
        for(Eigen::Index i=0;i<faces.rows();++i)
            moments+=referenceTetra(points,{{std::size_t(faces(i,0)),std::size_t(faces(i,1)),std::size_t(faces(i,2))}});
        for(std::size_t owner=0;owner<boundary.primitiveBoundaries().size();++owner)
            for(const auto& ref:boundary.primitiveBoundaries()[owner])
                moments-=beta[boundary.binding().primitives()[owner].betaIndex]*double(ref.orientation)*referenceTetra(points,boundary.triangles()[ref.triangle]);
        moments*=f.source().density();
        auto expected=independentPaper(moments,phi,f.config());
        for(auto e:boundary.binding().domain().edges()) {
            const Complex delta=beta[e.first]-beta[e.second];expected.objective+=.5*f.config().fillWeight*delta*delta;
        }
        jac(0,col)=expected.objective.imag()/step;
        jac.block(1,col,f.constraintCount(),1)=expected.constraints.imag()/step;
        jac.block(1+f.constraintCount(),col,10,1)=moments.imag()/step;
    }
    return jac;
}
using LongPoints=Eigen::Matrix<long double,Eigen::Dynamic,3>;
using LongVector=Eigen::Matrix<long double,Eigen::Dynamic,1>;
LongVector longReference(const SourceCagePaperFunctional& f,const LongVector& beta,long double phi,LongPoints points) {
    const auto& boundary=f.source().boundary();
    points=((points.rowwise()-f.frame().origin.cast<long double>().transpose())*f.frame().rotation.transpose().cast<long double>()).eval();
    Eigen::Matrix<long double,10,1> moments=Eigen::Matrix<long double,10,1>::Zero();
    const auto& faces=boundary.binding().source().outerFaces();
    for(Eigen::Index i=0;i<faces.rows();++i)
        moments+=referenceTetra(points,{{std::size_t(faces(i,0)),std::size_t(faces(i,1)),std::size_t(faces(i,2))}});
    for(std::size_t owner=0;owner<boundary.primitiveBoundaries().size();++owner)
        for(const auto& ref:boundary.primitiveBoundaries()[owner])
            moments-=beta[boundary.binding().primitives()[owner].betaIndex]*ref.orientation*referenceTetra(points,boundary.triangles()[ref.triangle]);
    moments*=static_cast<long double>(f.source().density());
    auto expected=independentPaper(moments,phi,f.config());
    for(auto e:boundary.binding().domain().edges()) {
        const long double delta=beta[e.first]-beta[e.second];expected.objective+=.5L*f.config().fillWeight*delta*delta;
    }
    LongVector result(11+f.constraintCount());result<<expected.objective,expected.constraints,moments;return result;
}
void longDoubleDifferences(const SourceCagePaperFunctional& f,const Eigen::VectorXd& x,const SourceCagePaperEvaluation& result) {
    // Precision repair for the independent numerical checker, not for the
    // production algorithm: long double from BBW accumulation onward and a
    // five-point centered stencil at exactly the original three fixed steps.
    require(std::numeric_limits<long double>::digits>=64,"reference requires at least 64 mantissa bits");
    const auto& boundary=f.source().boundary();const auto& bbw=boundary.binding().bbw();
    const auto& grid=bbw.fem().grid();const auto& cage=grid.cage();
    const auto weights=[&](Eigen::Index p) {
        const auto embedding=grid.embed(boundary.points().row(p));
        LongVector row=LongVector::Zero(f.handleCount());
        for(std::size_t h=0;h<f.handleCount();++h)for(int node=0;node<8;++node)
            row[h]+=static_cast<long double>(embedding.coefficients[node])*static_cast<long double>(bbw.finalWeights()(embedding.nodeIndices[node],h));
        return row;
    };
    LongPoints base=LongPoints::Zero(boundary.points().rows(),3);
    for(Eigen::Index p=0;p<base.rows();++p) {
        const Eigen::Matrix<long double,3,1> rest=boundary.points().row(p).transpose().cast<long double>();
        const auto row=weights(p);
        for(std::size_t h=0;h<f.handleCount();++h) {
            const Eigen::Matrix<long double,3,1> center=cage.boundaryVertices().row(cage.handleVertexIndices()[h]).transpose().cast<long double>();
            base.row(p)+=(row[h]*(center+static_cast<long double>(x[f.handleOffset()+4*h])*(rest-center)+x.segment<3>(f.handleOffset()+4*h+1).cast<long double>())).transpose();
        }
    }
    double worst[3]={0,0,0},absolute[3]={0,0,0};
    const long double steps[]={3e-4L,1e-4L,3e-5L};
    for(Eigen::Index col=0;col<x.size();++col) {
        LongPoints velocity=LongPoints::Zero(base.rows(),3);
        if(col>=Eigen::Index(f.handleOffset())) {
            const std::size_t h=(col-f.handleOffset())/4,slot=(col-f.handleOffset())%4;
            const Eigen::Matrix<long double,1,3> center=cage.boundaryVertices().row(cage.handleVertexIndices()[h]).cast<long double>();
            for(Eigen::Index p=0;p<base.rows();++p) {
                const long double w=weights(p)[h];
                if(slot)velocity(p,slot-1)=w;
                else velocity.row(p)=w*(boundary.points().row(p).cast<long double>()-center);
            }
        }
        for(int step=0;step<3;++step) {
            const long double h=steps[step];
            const auto value=[&](long double delta) {
                LongVector beta=x.head(f.betaCount()).cast<long double>();
                long double phi=f.hasAngle()?static_cast<long double>(x[f.angleOffset()]):0;
                if(col<Eigen::Index(f.betaCount()))beta[col]+=delta;
                else if(f.hasAngle() && col==Eigen::Index(f.angleOffset()))phi+=delta;
                return longReference(f,beta,phi,(base+delta*velocity).eval());
            };
            const LongVector fd=(-value(2*h)+8*value(h)-8*value(-h)+value(-2*h))/(12*h);
            for(Eigen::Index row=0;row<fd.size();++row) {
                const double actual=row==0?result.jointGradient[col]:(row<=result.constraints.size()?result.jointConstraintJacobian(row-1,col):result.spinMomentJacobian(row-1-result.constraints.size(),col));
                const double error=std::abs(fd[row]-static_cast<long double>(actual));absolute[step]=std::max(absolute[step],error);
                if(std::abs(actual)>1e-7) {worst[step]=std::max(worst[step],error/std::abs(actual));require(error/std::abs(actual)<=1e-6,"independent long-double five-point FD relative target");}
                else require(error<=2e-9,"independent long-double five-point FD absolute rule");
            }
        }
    }
    for(int step=0;step<3;++step)std::cout<<"LD5 "<<spin2014ObjectiveToken(f.config().objective)<<" columns="<<x.size()<<" h="<<double(steps[step])<<" relative="<<worst[step]<<" absolute="<<absolute[step]<<std::endl;
}
std::size_t doubleFinestDiscrepancies=0;
void allColumns(const SourceCagePaperFunctional& f,const Eigen::VectorXd& x,const SourceCagePaperEvaluation& r) {
    // Fixed input-unit scales, declared before observations. Near-zero magnitude
    // threshold 1e-7; absolute rule 2e-9. All other entries target 1e-6 relative.
    // Three ordinary double FD steps remain diagnostic, including discrepancies.
    // Accuracy gates above use every long-double five-point and complex-step
    // entry with these same thresholds, per the task's verification addendum.
    const double steps[]={3e-4,1e-4,3e-5};
    const auto cs=complexDerivatives(f,x,r);
    double csWorst=0,csAbsolute=0;
    for(Eigen::Index col=0;col<x.size();++col)for(Eigen::Index row=0;row<cs.rows();++row) {
        const double actual=row==0?r.jointGradient[col]:(row<=r.constraints.size()?r.jointConstraintJacobian(row-1,col):r.spinMomentJacobian(row-1-r.constraints.size(),col));
        const double error=std::abs(cs(row,col)-actual);csAbsolute=std::max(csAbsolute,error);
        if(std::abs(actual)>1e-7) {csWorst=std::max(csWorst,error/std::abs(actual));require(error/std::abs(actual)<=1e-6,"independent complex-step all-column relative target");}
        else require(error<=2e-9,"independent complex-step near-zero absolute rule");
    }
    std::cout<<"CS "<<spin2014ObjectiveToken(f.config().objective)<<" columns="<<x.size()<<" relative="<<csWorst<<" absolute="<<csAbsolute<<std::endl;
    longDoubleDifferences(f,x,r);
    for(double h:steps) {
        double worst=0,absolute=0;
        for(Eigen::Index col=0;col<x.size();++col) {
            auto plus=x,minus=x;plus[col]+=h;minus[col]-=h;
            const auto a=f.evaluateValues(plus),b=f.evaluateValues(minus);
            const auto check=[&](double fd,double actual) {
                const double error=std::abs(fd-actual);absolute=std::max(absolute,error);
                if(std::abs(actual)>1e-7) {
                    worst=std::max(worst,error/std::abs(actual));
                    if(error/std::abs(actual)>1e-6) {
                        std::cerr<<"FD discrepancy mode="<<spin2014ObjectiveToken(f.config().objective)<<" column="<<col<<" h="<<h<<" actual="<<actual<<" fd="<<fd<<" relative="<<error/std::abs(actual)<<'\n';
                        if(h==steps[2])++doubleFinestDiscrepancies;
                    }
                } else require(error<=2e-9,"all joint columns near-zero absolute tolerance");
            };
            check((a.totalObjective-b.totalObjective)/(2*h),r.jointGradient[col]);
            for(Eigen::Index row=0;row<r.constraints.size();++row)
                check((a.constraints[row]-b.constraints[row])/(2*h),r.jointConstraintJacobian(row,col));
            for(int row=0;row<10;++row)
                check((a.spinMoments.values()[row]-b.spinMoments.values()[row])/(2*h),r.spinMomentJacobian(row,col));
        }
        std::cout<<"FD "<<spin2014ObjectiveToken(f.config().objective)<<" columns="<<x.size()<<" h="<<h<<" relative="<<worst<<" absolute="<<absolute<<std::endl;
    }
}

void cases(bool transformed) {
    SourceCageMomentFunctional source(fixture(transformed),transformed?.8:2.3,{});
    require(source.betaCount()==8 && source.handleCount()==26,"real source dimensions");
    Eigen::VectorXd beta(8);for(int i=0;i<8;++i)beta[i]=.13+.017*i;
    auto handles=source.initialHandleParameters();
    for(int h=0;h<26;++h) {
        handles[4*h]=1.+.015*std::sin(h+1);
        for(int a=0;a<3;++a)handles[4*h+1+a]=.006*std::cos((h+1)*(a+1));
    }
    SourceCagePaperFrame frame;frame.origin={.19,-.23,.11};
    frame.rotation=Eigen::AngleAxisd(.41,Eigen::Vector3d(2,-1,3).normalized()).toRotationMatrix();
    for(auto objective:{Spin2014Objective::YoYo,Spin2014Objective::TopInertiaOnly,Spin2014Objective::TopLoweringOnly,Spin2014Objective::TopMassAndLowering,Spin2014Objective::StaticLowering}) {
        SourceCagePaperFunctionalConfig c;c.objective=objective;c.inertiaWeight=.7;c.centerWeight=1.3;c.fillWeight=.19;
        SourceCagePaperFunctional f(source,c,frame);
        const double phi=f.hasAngle()?.23:0;
        auto x=f.initialState(beta,phi);x.tail(handles.size())=handles;
        require(x.size()==(f.hasAngle()?113:112) && f.handleOffset()==std::size_t(f.hasAngle()?9:8),"no redundant static angle");
        require(f.angleOffset()==(f.hasAngle()?8:std::numeric_limits<std::size_t>::max()),"angle offset sentinel");
        const auto r=f.evaluate(x),values=f.evaluateValues(beta,phi,handles);
        require(r.hasJacobians && !values.hasJacobians && values.jointGradient.size()==0 && values.jointConstraintJacobian.size()==0 && values.spinMomentJacobian.size()==0 && !values.world.hasJacobians,"value-only global derivatives absent");
        require(r.spinMoments.values()==values.spinMoments.values() && r.world.deformedPoints==values.world.deformedPoints && r.totalObjective==values.totalObjective && r.constraints==values.constraints,"value and derivative results identical");
        close(r.world.finalMoments.values(),integrated(source.boundary(),r.world.deformedPoints,beta,source.density()),"independent world tetra moments");
        const Eigen::MatrixXd points=(r.world.deformedPoints.rowwise()-frame.origin.transpose())*frame.rotation.transpose();
        const Ten spin=integrated(source.boundary(),points,beta,source.density());
        close(r.spinMoments.values(),spin,"fixed-frame point integration all ten mixed moments");
        close(f.worldToSpinMomentMap()*r.world.finalMoments.values(),spin,"exact linear map");
        const auto expected=independentPaper(spin,phi,c);
        near(r.paperObjective,expected.objective,"independent inertia tensor paper objective");
        require((r.constraints-expected.constraints).norm()<2e-10,"independent constraints and ordering");
        double edgeEnergy=0;Eigen::VectorXd edgeGradient=Eigen::VectorXd::Zero(8);
        for(auto e:source.boundary().binding().domain().edges()) {
            const double delta=beta[e.first]-beta[e.second];edgeEnergy+=.5*c.fillWeight*delta*delta;
            edgeGradient[e.first]+=c.fillWeight*delta;edgeGradient[e.second]-=c.fillWeight*delta;
        }
        near(r.fillPenalty,edgeEnergy,"sparse fill equals actual edge energy");
        near(r.totalObjective,r.paperObjective+edgeEnergy,"independent weights sum");
        auto nofill=c;nofill.fillWeight=0;const auto unregularized=SourceCagePaperFunctional(source,nofill,frame).evaluate(x);
        require((r.jointGradient.head(8)-unregularized.jointGradient.head(8)-edgeGradient).norm()<2e-12,"edge fill gradient");
        require(r.jointGradient.tail(x.size()-8)==unregularized.jointGradient.tail(x.size()-8),"fill handle and angle derivatives zero");
        near(f.evaluateValues(Eigen::VectorXd::Constant(8,.4),phi,handles).fillPenalty,0,"constant beta fill nullspace");
        if(f.hasAngle()) {
            require(r.spinMomentJacobian.col(f.angleOffset()).isZero(0),"exact zero raw phi column");
            auto changed=x;changed[f.angleOffset()]+=.1;
            require((f.evaluateValues(changed).constraints-r.constraints).norm()>1e-8,"phi affects transverse equality");
        }
        require(f.evaluate(x).outputSha256==r.outputSha256 && values.outputSha256!=r.outputSha256,"deterministic mode-bound digest");
        auto changed=x;changed[0]+=.001;require(f.evaluateValues(changed).outputSha256!=values.outputSha256,"state-bound digest");
        auto capped=c;capped.maximumJointDerivativeValues=1;
        SourceCagePaperFunctional limited(source,capped,frame);
        require(limited.evaluateValues(x).totalObjective==r.totalObjective,"value-only ignores derivative budget");
        throwsAs<std::length_error>([&]{limited.evaluate(x);},"joint derivative resource guard");
        require(limited.configSha256()!=f.configSha256() && limited.evaluateValues(x).outputSha256!=values.outputSha256,"budget-bound identity");
        auto exact=c;exact.maximumJointDerivativeValues=(11+f.constraintCount())*f.stateCount()+f.betaCount();
        require(SourceCagePaperFunctional(source,exact,frame).evaluate(x).jointGradient==r.jointGradient,"exact combined derivative capacity accepted");
        --exact.maximumJointDerivativeValues;
        throwsAs<std::length_error>([&]{SourceCagePaperFunctional(source,exact,frame).evaluate(x);},"one below combined capacity rejected");
        SourceCageMomentFunctionalConfig upstreamCap;upstreamCap.maximumJacobianValues=1;
        SourceCagePaperFunctional upstreamLimited(SourceCageMomentFunctional(source.boundary(),source.density(),upstreamCap),c,frame);
        require(upstreamLimited.evaluateValues(x).totalObjective==r.totalObjective,"value-only upstream derivative budget independence");
        throwsAs<std::length_error>([&]{upstreamLimited.evaluate(x);},"upstream derivative budget preserved");
        if(f.hasAngle())for(double bad:{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()})
            throwsAs<std::invalid_argument>([&]{f.initialState(beta,bad);},"nonstatic initial phi finite");
        auto doubled=SourceCagePaperFunctional(SourceCageMomentFunctional(source.boundary(),2*source.density(),{}),c,frame).evaluateValues(x);
        near(doubled.fillPenalty,r.fillPenalty,"density-independent fill");
        const double extra=objective==Spin2014Objective::TopMassAndLowering?3*c.centerWeight*spin[3]*spin[3]:0;
        near(doubled.paperObjective,r.paperObjective+extra,"density covariance five objectives");
        allColumns(f,x,r);
    }
    SourceCagePaperFunctionalConfig c;SourceCagePaperFunctional f(source,c,frame);
    auto x=f.initialState(beta,0);
    const double nan=std::numeric_limits<double>::quiet_NaN(),inf=std::numeric_limits<double>::infinity();
    for(double phi:{.1,nan,inf})throwsAs<std::invalid_argument>([&]{f.evaluateValues(beta,phi,handles);},"static/nonfinite phi rejected");
    throwsAs<std::invalid_argument>([&]{f.initialState(beta,.1);},"static initial builder no redundant phi");
    throwsAs<std::invalid_argument>([&]{f.evaluateValues(x.head(x.size()-1));},"joint size rejected");
    for(double bad:{-.1,1.1,nan,inf}) {auto q=x;q[0]=bad;throwsAs<std::invalid_argument>([&]{f.evaluateValues(q);},"invalid beta");}
    for(double bad:{0.,-1.,nan,inf}) {auto q=x;q[f.handleOffset()]=bad;throwsAs<std::invalid_argument>([&]{f.evaluateValues(q);},"invalid handle scale");}
    throwsAs<std::invalid_argument>([&]{f.evaluateValues(beta,0,handles.head(1));},"split handle dimensions");
    for(int field=0;field<3;++field)for(double bad:{-1.,nan,inf}) {
        auto b=c;(field==0?b.inertiaWeight:(field==1?b.centerWeight:b.fillWeight))=bad;
        throwsAs<std::invalid_argument>([&]{SourceCagePaperFunctional invalid(source,b,frame);},"invalid independent weight");
    }
    auto invalid=c;invalid.objective=static_cast<Spin2014Objective>(-1);
    throwsAs<std::invalid_argument>([&]{SourceCagePaperFunctional bad(source,invalid,frame);},"invalid objective");
    for(auto cap:{std::size_t(0),std::numeric_limits<std::size_t>::max()}) {auto b=c;b.maximumJointDerivativeValues=cap;throwsAs<std::length_error>([&]{SourceCagePaperFunctional bad(source,b,frame);},"invalid cap");}
    for(int kind=0;kind<5;++kind) {
        auto b=frame;
        if(kind==0)b.origin[0]=nan;if(kind==1)b.rotation(0,0)=inf;
        if(kind==2)b.rotation.row(0)*=-1;if(kind==3)b.rotation(0,0)+=.01;if(kind==4)b.rotation*=2;
        throwsAs<std::invalid_argument>([&]{SourceCagePaperFunctional bad(source,c,b);},"invalid fixed frame");
    }
    auto huge=frame;huge.origin[0]=std::numeric_limits<double>::max();
    throwsAs<std::overflow_error>([&]{SourceCagePaperFunctional bad(source,c,huge);},"frame map arithmetic overflow");
    auto hugeHandles=handles;hugeHandles[1]=std::numeric_limits<double>::max();
    throwsAs<std::overflow_error>([&]{f.evaluateValues(beta,0,hugeHandles);},"upstream arithmetic failure preserved");
    auto hugeFill=c;hugeFill.fillWeight=std::numeric_limits<double>::max();
    Eigen::VectorXd alternating(8);for(int i=0;i<8;++i)alternating[i]=i%2;
    throwsAs<std::overflow_error>([&]{SourceCagePaperFunctional(source,hugeFill,frame).evaluateValues(alternating,0,handles);},"fill arithmetic fails closed");
    auto moved=frame;moved.origin[0]+=.01;
    require(SourceCagePaperFunctional(source,c,moved).configSha256()!=f.configSha256(),"frame identity sensitivity");
    for(int field=0;field<3;++field) {auto b=c;(field==0?b.inertiaWeight:(field==1?b.centerWeight:b.fillWeight))+=.1;require(SourceCagePaperFunctional(source,b,frame).configSha256()!=f.configSha256(),"weight identity sensitivity");}
}
int main() {
    try {
        std::cout<<"independent reference long-double mantissa bits="<<std::numeric_limits<long double>::digits<<std::endl;
        cases(false);cases(true);
        std::cout<<"ordinary double finest-step discrepancies="<<doubleFinestDiscrepancies<<" (retained diagnostic, not claimed passing FD)\n";
        std::cout<<"source cage paper SYNTHETIC_REGRESSION PASS: independent LD5 and complex-step accuracy gates\n";return 0;
    }
    catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
