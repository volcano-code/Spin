#include <iostream>
#if !__has_include("spin2014/LinearizedRotationLaplacian.h")
int main() {std::cerr<<"FAIL: source linearized-rotation shape operator is not implemented\n";return 1;}
#else
#include "spin2014/SourceCagePaperFunctional.h"
#include "spin2014/SourceCageShapeFunctional.h"
#include "source_surface_test_utils.h"
#include <Eigen/QR>
#include <Eigen/Geometry>
#include <set>
#include <type_traits>
using namespace spin2014;
using namespace surface_test;
namespace {
template<class Exception,class F> void throwsAs(F f,const char* label) {
    bool caught=false;try {f();}catch(const Exception& e) {
        caught=true;
        if constexpr(std::is_same_v<Exception,std::length_error>)
            require(std::string(e.what()).find("RESOURCE_LIMIT")!=std::string::npos,"resource category");
    }require(caught,label);
}
void near(double a,double b,const char* label,double tol=1e-11) {require(std::abs(a-b)<=tol*(1+std::abs(b)),label);}
using Long=long double;
using Points=Eigen::Matrix<Long,Eigen::Dynamic,3>;
// Independent reference: raw rest coordinates, QR fit on DEFORMED positions,
// then transform delta. Never reads or multiplies the production operator.
Points referenceResidual(const Eigen::MatrixXd& rest,const Eigen::MatrixXi& faces,const Points& moved) {
    std::vector<std::set<int>> rings(rest.rows());
    for(Eigen::Index f=0;f<faces.rows();++f)for(int a=0;a<3;++a)for(int b=0;b<3;++b)
        if(a!=b)rings[faces(f,a)].insert(faces(f,b));
    Points residual(rest.rows(),3);
    for(Eigen::Index i=0;i<rest.rows();++i) {
        std::vector<int> local{int(i)};local.insert(local.end(),rings[i].begin(),rings[i].end());
        Eigen::Matrix<Long,Eigen::Dynamic,Eigen::Dynamic> a=Eigen::Matrix<Long,Eigen::Dynamic,Eigen::Dynamic>::Zero(3*local.size(),7);
        Eigen::Matrix<Long,Eigen::Dynamic,1> b(3*local.size());
        Eigen::Matrix<Long,3,1> delta=rest.row(i).transpose().cast<Long>(),lap=moved.row(i).transpose();
        for(int j:rings[i]) {delta-=rest.row(j).transpose().cast<Long>()/Long(rings[i].size());lap-=moved.row(j).transpose()/Long(rings[i].size());}
        for(std::size_t j=0;j<local.size();++j) {
            const Long x=rest(local[j],0),y=rest(local[j],1),z=rest(local[j],2);
            a.block(3*j,0,3,7)<<x,0,z,-y,1,0,0,y,-z,0,x,0,1,0,z,y,-x,0,0,0,1;
            b.segment(3*j,3)=moved.row(local[j]).transpose();
        }
        const Eigen::Matrix<Long,7,1> p=a.colPivHouseholderQr().solve(b);
        const Eigen::Matrix<Long,3,1> h=p.segment<3>(1);
        residual.row(i)=(p[0]*delta+h.cross(delta)-lap).transpose();
    }
    return residual;
}
Long referenceEnergy(const Eigen::MatrixXd& rest,const Eigen::MatrixXi& faces,const Points& moved,double sw,double dw) {
    return Long(sw)*referenceResidual(rest,faces,moved).squaredNorm()+Long(dw)*(moved-rest.cast<Long>()).squaredNorm();
}
void derivativeGate(Long actual,Long expected,Long inputScale,const char* label) {
    const Long error=std::abs(actual-expected);
    if(std::abs(expected)>1e-7L*inputScale)require(error/std::abs(expected)<=1e-6L,label);
    else require(error<=1e-10L*inputScale,label);
}
void mathematical() {
    auto mesh=cube(.1,.9);
    mesh.v(6,0)+=.17;mesh.v(6,1)-=.04;mesh.v(3,2)+=.07;
    const auto rest=mesh.v;
    LinearizedRotationLaplacian op(rest,mesh.f,{});
    auto moved=rest;
    moved(0,0)+=.12;moved(4,1)-=.06;moved(7,2)+=.08;
    const auto result=op.evaluate(moved,.7,.2);
    const auto reference=referenceResidual(rest,mesh.f,moved.cast<Long>());
    require((result.residual-reference.cast<double>()).norm()<1e-12,"signed Sorkine residual vs independent raw QR");
    near(result.shapeEnergy,.7*double(reference.squaredNorm()),"shape without half factor");
    near(result.displacementEnergy,.2*(moved-rest).squaredNorm(),"displacement without half factor");
    require(result.shapeEnergy>1e-5,"local deformation penalized");
    const Eigen::Vector3d shift(.4,-.7,.2),h(.12,-.05,.07);
    Eigen::Matrix3d cross;cross<<0,-h.z(),h.y(),h.z(),0,-h.x(),-h.y(),h.x(),0;
    for(const Eigen::MatrixXd& field:std::vector<Eigen::MatrixXd>{rest,(rest.rowwise()+shift.transpose()).eval(),(1.3*rest).eval(),(rest*(Eigen::Matrix3d::Identity()+cross).transpose()).eval()})
        require(op.evaluateValues(field,1.,0.).shapeEnergy<1e-25,"rest translation uniform scale infinitesimal rotation nullspace");
    Eigen::Matrix3d shear=Eigen::Matrix3d::Identity();shear(0,1)=.3;
    require(op.evaluateValues(rest*shear.transpose(),1.,0.).shapeEnergy>1e-5,"shear not nullspace");
    const Eigen::Matrix3d rotation=Eigen::AngleAxisd(.9,Eigen::Vector3d(1,2,3).normalized()).toRotationMatrix();
    require(op.evaluateValues(rest*rotation.transpose(),1.,0.).shapeEnergy>1e-5,"finite rotation is not infinitesimal nullspace");
    const double scale=2.3;
    const Eigen::MatrixXd transformedRest=(scale*rest*rotation.transpose()).rowwise()+shift.transpose();
    const Eigen::MatrixXd transformedMoved=(scale*moved*rotation.transpose()).rowwise()+shift.transpose();
    LinearizedRotationLaplacian transformed(transformedRest,mesh.f,{});
    const auto tr=transformed.evaluate(transformedMoved,.7,.2);
    require((tr.residual-scale*result.residual*rotation.transpose()).norm()<1e-12,"covariant residual");
    near(tr.shapeEnergy,scale*scale*result.shapeEnergy,"scale squared shape units");
    near(tr.displacementEnergy,scale*scale*result.displacementEnergy,"scale squared displacement units");
    require((tr.vertexGradient-scale*result.vertexGradient*rotation.transpose()).norm()<1e-11,"covariant gradient");
    const Long inputScale=std::max(1.,rest.squaredNorm());
    for(Long step:{3e-4L,1e-4L,3e-5L})for(Eigen::Index v=0;v<rest.rows();++v)for(int axis=0;axis<3;++axis) {
        Points plus=moved.cast<Long>(),minus=plus;plus(v,axis)+=step;minus(v,axis)-=step;
        const Long fd=(referenceEnergy(rest,mesh.f,plus,.7,.2)-referenceEnergy(rest,mesh.f,minus,.7,.2))/(2*step);
        derivativeGate(result.vertexGradient(v,axis),fd,inputScale,"vertex gradient three-step raw-QR FD");
        Eigen::MatrixXd p=moved,m=moved;p(v,axis)+=double(step);m(v,axis)-=double(step);
        const auto rp=op.evaluateValues(p,.7,.2),rm=op.evaluateValues(m,.7,.2);
        derivativeGate(result.vertexGradient(v,axis),((rp.shapeEnergy+rp.displacementEnergy)-(rm.shapeEnergy+rm.displacementEnergy))/(2*step),inputScale,"vertex gradient production FD");
    }
    auto values=op.evaluateValues(moved,.7,.2);
    require(!values.hasGradient && values.vertexGradient.size()==0,"values have no vertex derivatives");
    require(values.shapeEnergy==result.shapeEnergy && values.displacementEnergy==result.displacementEnergy,"value identical energy");
    require(op.evaluate(moved,.7,.2).outputSha256==result.outputSha256,"exact repeated output hash");
    require(LinearizedRotationLaplacian(rest,mesh.f,{}).operatorSha256()==op.operatorSha256(),"exact repeated operator");
    require(transformed.operatorSha256()!=op.operatorSha256(),"geometry identity");
    require((mesh.v-rest).norm()==0,"immutable constructor input");
    require(op.sparseValueCount()==9*op.coefficientCount(),"full block accounting");
    LinearizedRotationLaplacianConfig constructionCap;
    constructionCap.maximumConstructionWorkspaceDoubles=op.requiredConstructionWorkspaceDoubles();
    LinearizedRotationLaplacian constructionBoundary(rest,mesh.f,constructionCap);
    --constructionCap.maximumConstructionWorkspaceDoubles;
    throwsAs<std::length_error>([&]{LinearizedRotationLaplacian bad(rest,mesh.f,constructionCap);},"constructor exact-boundary one below");
    LinearizedRotationLaplacianConfig cap;cap.maximumSparseValues=op.sparseValueCount();
    LinearizedRotationLaplacian boundary(rest,mesh.f,cap);
    --cap.maximumSparseValues;
    throwsAs<std::length_error>([&]{LinearizedRotationLaplacian bad(rest,mesh.f,cap);},"sparse boundary cap");
    cap={};cap.maximumEvaluationWorkspaceDoubles=3*rest.rows();
    require(op.requiredEvaluationWorkspaceDoubles(false)==3*std::size_t(rest.rows()) && op.requiredEvaluationWorkspaceDoubles(true)==6*std::size_t(rest.rows()),"actual required evaluation arrays");
    LinearizedRotationLaplacian valuesOnly(rest,mesh.f,cap);
    require(valuesOnly.evaluateValues(moved,1.,1.).vertexGradient.size()==0,"values under derivative shortage");
    throwsAs<std::length_error>([&]{valuesOnly.evaluate(moved,1.,1.);},"gradient workspace cap");
    cap.maximumEvaluationWorkspaceDoubles=3*rest.rows()-1;
    LinearizedRotationLaplacian tooSmall(rest,mesh.f,cap);
    throwsAs<std::length_error>([&]{tooSmall.evaluateValues(moved,1.,1.);},"value workspace cap");
    cap={};cap.maximumConstructionWorkspaceDoubles=1;
    throwsAs<std::length_error>([&]{LinearizedRotationLaplacian bad(rest,mesh.f,cap);},"constructor cap");
    cap={};cap.maximumVertices=rest.rows()-1;
    throwsAs<std::length_error>([&]{LinearizedRotationLaplacian bad(rest,mesh.f,cap);},"vertex cap");
    cap={};cap.maximumOneRingVertices=2;
    throwsAs<std::length_error>([&]{LinearizedRotationLaplacian bad(rest,mesh.f,cap);},"high valence cap");
    for(double tolerance:{0.,-1.,1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
        cap={};cap.relativeRankTolerance=tolerance;
        throwsAs<std::invalid_argument>([&]{LinearizedRotationLaplacian bad(rest,mesh.f,cap);},"invalid rank tolerance");
    }
    cap={};cap.maximumSparseValues=std::numeric_limits<std::size_t>::max();
    throwsAs<std::length_error>([&]{LinearizedRotationLaplacian bad(rest,mesh.f,cap);},"byte overflow budget");
    for(double invalid:{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
        auto bad=rest;bad(0,0)=invalid;
        throwsAs<std::invalid_argument>([&]{LinearizedRotationLaplacian x(bad,mesh.f,{});},"nonfinite rest");
        throwsAs<std::invalid_argument>([&]{op.evaluateValues(bad,1.,0.);},"nonfinite deformation");
        throwsAs<std::invalid_argument>([&]{op.evaluateValues(rest,invalid,0.);},"nonfinite weight");
    }
    auto faces=mesh.f;faces(0,0)=faces(0,1);
    throwsAs<std::invalid_argument>([&]{LinearizedRotationLaplacian bad(rest,faces,{});},"repeated face index");
    faces=mesh.f;faces(0,0)=rest.rows();
    throwsAs<std::invalid_argument>([&]{LinearizedRotationLaplacian bad(rest,faces,{});},"out of range face");
    auto bad=rest;bad.row(0)=bad.row(1);
    throwsAs<std::invalid_argument>([&]{LinearizedRotationLaplacian x(bad,mesh.f,{});},"zero area");
    bad.conservativeResize(rest.rows()+1,3);bad.row(rest.rows()).setZero();
    throwsAs<std::invalid_argument>([&]{LinearizedRotationLaplacian x(bad,mesh.f,{});},"isolated vertex");
    throwsAs<std::invalid_argument>([&]{op.evaluateValues(Eigen::MatrixXd::Zero(3,2),1.,0.);},"malformed deformed shape");
    Eigen::MatrixXd planar(4,3);planar<<0,0,0,1,0,0,1,1,0,0,1,0;
    Eigen::MatrixXi triangles(2,3);triangles<<0,1,2,0,2,3;
    LinearizedRotationLaplacian patch(planar,triangles,{});
    require(patch.evaluateValues(planar,1.,0.).shapeEnergy<1e-25,"planar noncollinear valid");
    planar(2,1)=1e-15;planar(3,1)=2e-15;
    throwsAs<std::invalid_argument>([&]{LinearizedRotationLaplacian x(planar,triangles,{});},"numerically rank deficient thin patch");
    auto overflow=moved;overflow(0,0)=1e308;
    throwsAs<std::overflow_error>([&]{op.evaluate(overflow,1e308,1e308);},"energy overflow");
    for(int n:{16,32}) {
        Eigen::MatrixXd fan=Eigen::MatrixXd::Zero(n+1,3);Eigen::MatrixXi fanFaces(n,3);
        for(int i=0;i<n;++i) {
            const double angle=2*std::acos(-1.)*i/n;fan.row(i+1)<<std::cos(angle),std::sin(angle),0.;
            fanFaces.row(i)<<0,i+1,(i+1)%n+1;
        }
        LinearizedRotationLaplacian grown(fan,fanFaces,{});
        require(grown.coefficientCount()==std::size_t(5*n+1) && grown.maximumOneRingSize()==std::size_t(n),"linear sparse fan growth with valence");
        cap={};cap.maximumOneRingVertices=n;
        LinearizedRotationLaplacian ringBoundary(fan,fanFaces,cap);
        --cap.maximumOneRingVertices;
        throwsAs<std::length_error>([&]{LinearizedRotationLaplacian x(fan,fanFaces,cap);},"high-valence exact-boundary cap");
    }
    for(int setting=0;setting<5;++setting) {
        cap={};
        if(setting==0)cap.maximumVertices=0;
        if(setting==1)cap.maximumOneRingVertices=0;
        if(setting==2)cap.maximumSparseValues=0;
        if(setting==3)cap.maximumConstructionWorkspaceDoubles=0;
        if(setting==4)cap.maximumEvaluationWorkspaceDoubles=0;
        throwsAs<std::length_error>([&]{LinearizedRotationLaplacian x(rest,mesh.f,cap);},"all zero budgets rejected");
    }
    std::cout<<"mathematical nullspaces, independent QR/FD, covariance, caps, invalid inputs PASS\n";
}
SourceConformingBoundary fixture(bool transformed) {
    auto outer=cube(.3,.7);auto design=transformed?box({.34,.36,.38},{.66,.64,.62}):cube(.34,.66);
    const double scale=transformed?2.:1.;const Eigen::Vector3d shift=transformed?Eigen::Vector3d(1.25,-.75,.125):Eigen::Vector3d::Zero();
    outer.v=(outer.v*scale).rowwise()+shift.transpose();design.v=(design.v*scale).rowwise()+shift.transpose();
    SourceSurfaceOctreeConfig sc;sc.maximumDepth=5;sc.rootOrigin=shift;sc.rootSide=scale;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
    SourceAdaptiveDomain domain(source,2,5);
    SourceSurfaceCageConfig cc;cc.cageLevel=2;cc.rootOrigin=shift;cc.rootSide=scale;
    SourceCageGridConfig gc;gc.fineLevel=3;
    SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),{}),{});
    require(bbw.success(),"actual BBW fixture");
    return SourceConformingBoundary(SourceCageBinding(domain,bbw,{}),{});
}
void actualSource(bool transformed) {
    const auto boundary=fixture(transformed);const auto& binding=boundary.binding();
    const auto& rest=binding.source().outerVertices();const auto& faces=binding.source().outerFaces();
    const auto& bbw=binding.bbw();const auto& grid=bbw.fem().grid();const auto& cage=grid.cage();
    const SourceCageMomentFunctional moment(boundary,2.7,{});
    require(moment.betaCount()==8 && moment.handleCount()==26,"all actual fixture variables");
    const Long inputScale=std::max(1.,rest.squaredNorm());
    for(auto mode:{Spin2014Objective::StaticLowering,Spin2014Objective::TopInertiaOnly,Spin2014Objective::TopMassAndLowering,Spin2014Objective::TopLoweringOnly,Spin2014Objective::YoYo}) {
        SourceCagePaperFunctionalConfig c;c.objective=mode;c.fillWeight=.31;c.inertiaWeight=.63;c.centerWeight=.42;
        SourceCagePaperFrame frame;frame.origin=transformed?Eigen::Vector3d(2.25,.25,.725):Eigen::Vector3d(.5,.5,.3);
        const SourceCagePaperFunctional disabled(moment,c,frame);
        c.shapeWeight=.7;c.displacementWeight=.2;
        const SourceCagePaperFunctional enabled(moment,c,frame);
        Eigen::VectorXd beta=Eigen::VectorXd::Constant(8,.2);
        auto state=enabled.initialState(beta,enabled.hasAngle()?.13:0.);
        for(std::size_t h=0;h<26;++h) {
            state[enabled.handleOffset()+4*h]=1.+.03*std::sin(double(h+1));
            for(int a=0;a<3;++a)state[enabled.handleOffset()+4*h+1+a]=.02*std::cos(double(2*h+a));
        }
        const auto result=enabled.evaluate(state),base=disabled.evaluate(state),values=enabled.evaluateValues(state);
        require(enabled.shapeFunctional()!=nullptr && disabled.shapeFunctional()==nullptr,"honest enabled identity");
        const auto shape=enabled.shapeFunctional()->evaluate(result.world.deformedPoints.topRows(rest.rows()));
        const auto reference=referenceResidual(rest,faces,result.world.deformedPoints.topRows(rest.rows()).cast<Long>());
        near(result.shapeEnergy,.7*double(reference.squaredNorm()),"actual source independent shape energy");
        near(result.displacementEnergy,.2*(result.world.deformedPoints.topRows(rest.rows())-rest).squaredNorm(),"actual displacement energy");
        near(result.totalObjective,base.totalObjective+result.shapeEnergy+result.displacementEnergy,"additive objective");
        require((result.world.finalMoments.values()-base.world.finalMoments.values()).norm()==0 && (result.constraints-base.constraints).norm()==0,"unchanged moments and constraints");
        require((result.spinMomentJacobian-base.spinMomentJacobian).norm()==0 && (result.jointConstraintJacobian-base.jointConstraintJacobian).norm()==0,"unchanged moment and constraint Jacobians");
        require((result.jointGradient.head(enabled.handleOffset())-base.jointGradient.head(enabled.handleOffset())).norm()==0,"shape beta and phi derivatives exactly zero");
        require((result.jointGradient.tail(104)-base.jointGradient.tail(104)-shape.handleGradient).norm()<1e-13,"shape chain integrated all handles");
        require(values.shapeEnergy==result.shapeEnergy && values.displacementEnergy==result.displacementEnergy && values.totalObjective==result.totalObjective,"actual values identical");
        require(values.jointGradient.size()==0 && values.spinMomentJacobian.size()==0,"value only no derivatives");
        require(enabled.evaluate(state).outputSha256==result.outputSha256,"actual repeat twice hash");
        require(result.outputSha256!=base.outputSha256 && enabled.configSha256()!=disabled.configSha256(),"new mode hashes");
        for(Eigen::Index col=0;col<state.size();++col) {
            Points velocity=Points::Zero(rest.rows(),3);
            if(col>=Eigen::Index(enabled.handleOffset())) {
                const std::size_t h=(col-enabled.handleOffset())/4,slot=(col-enabled.handleOffset())%4;
                const Eigen::Matrix<Long,1,3> center=cage.boundaryVertices().row(cage.handleVertexIndices()[h]).cast<Long>();
                for(Eigen::Index p=0;p<rest.rows();++p) {
                    const auto embedding=grid.embed(rest.row(p));Long w=0;
                    for(int n=0;n<8;++n)w+=Long(embedding.coefficients[n])*Long(bbw.finalWeights()(embedding.nodeIndices[n],h));
                    if(slot)velocity(p,slot-1)=w;else velocity.row(p)=w*(rest.row(p).cast<Long>()-center);
                }
            }
            const Long actual=result.jointGradient[col]-base.jointGradient[col];
            for(Long step:{3e-4L,1e-4L,3e-5L}) {
                const Points points=result.world.deformedPoints.topRows(rest.rows()).cast<Long>();
                const Long fd=(referenceEnergy(rest,faces,(points+step*velocity).eval(),.7,.2)-referenceEnergy(rest,faces,(points-step*velocity).eval(),.7,.2))/(2*step);
                derivativeGate(actual,fd,inputScale,"all joint columns three-step independent FD");
            }
        }
        c.shapeWeight=0.;const SourceCagePaperFunctional displacement(moment,c,frame);
        const auto displacementResult=displacement.evaluate(state);
        require(displacement.shapeFunctional()!=nullptr && displacementResult.shapeEnergy==0.,"displacement-only honest enabled mode");
        require(displacementResult.displacementEnergy==result.displacementEnergy,"shape weight independent");
        c.shapeWeight=.7;c.displacementWeight=0.;const SourceCagePaperFunctional shapeOnly(moment,c,frame);
        const auto shapeResult=shapeOnly.evaluate(state);
        require(shapeResult.shapeEnergy==result.shapeEnergy && shapeResult.displacementEnergy==0.,"displacement weight independent");
        c.maximumShapeHandleGradientValues=1;
        const SourceCagePaperFunctional limited(moment,c,frame);
        require(limited.evaluateValues(state).jointGradient.size()==0,"shape all-handle budget value-only");
        throwsAs<std::length_error>([&]{limited.evaluate(state);},"shape all-handle gradient cap");
        std::cout<<"actual fixture="<<transformed<<" mode="<<spin2014ObjectiveToken(mode)<<" columns="<<state.size()<<" three-step QR FD PASS\n";
    }
    SourceCageShapeFunctionalConfig sc;sc.shapeWeight=.7;sc.displacementWeight=.2;
    SourceCageShapeFunctional shape(binding,sc);
    const auto original=shape.evaluate(rest);
    SourceCageShapeFunctional copy=shape;
    require(copy.configSha256()==shape.configSha256() && copy.evaluate(rest).outputSha256==original.outputSha256,"immutable wrapper copies");
    SourceCageMomentFunctional otherDensity(boundary,5.4,{});
    SourceCagePaperFunctionalConfig pc;pc.shapeWeight=.7;pc.displacementWeight=.2;
    SourceCagePaperFunctional f1(moment,pc),f2(otherDensity,pc);
    auto state=f1.initialState(Eigen::VectorXd::Constant(8,.2),0.);state[f1.handleOffset()+1]=.01;
    const auto r1=f1.evaluateValues(state),r2=f2.evaluateValues(state);
    require(r1.shapeEnergy==r2.shapeEnergy && r1.displacementEnergy==r2.displacementEnergy,"no density factor");
    sc.maximumWeightValues=25;SourceCageShapeFunctional insufficient(binding,sc);
    require(insufficient.evaluateValues(rest).handleGradient.size()==0,"values require no weight row");
    throwsAs<std::length_error>([&]{insufficient.evaluate(rest);},"full BBW weight-row cap");
    const auto retained=binding.outputSha256();
    require(retained==shape.binding().outputSha256(),"source binding retained identity");
    SourceCageBbwConfig changedBbwConfig=bbw.config();++changedBbwConfig.maximumDenseValues;
    SourceCageBbw changedBbw(bbw.fem(),changedBbwConfig);require(changedBbw.success(),"actual changed BBW identity fixture");
    const SourceCageBinding changedBinding(binding.domain(),changedBbw,binding.config());
    sc.maximumWeightValues=26;sc.maximumHandleGradientValues=104;
    const SourceCageShapeFunctional changed(changedBinding,sc),exactCaps(binding,sc);
    require(exactCaps.evaluate(rest).handleGradient.size()==104,"all handle exact caps");
    require(changedBinding.bbw().outputSha256()!=binding.bbw().outputSha256() && changed.configSha256()!=exactCaps.configSha256(),"changed actual BBW input identity");
    require(changed.evaluateValues(rest).outputSha256!=exactCaps.evaluateValues(rest).outputSha256,"changed BBW distinguished under zero motion");
    sc.maximumHandleGradientValues=103;const SourceCageShapeFunctional below(binding,sc);
    throwsAs<std::length_error>([&]{below.evaluate(rest);},"one below actual handle derivative cap");
    static std::string previousSourceIdentity;
    if(!previousSourceIdentity.empty())require(previousSourceIdentity!=shape.configSha256(),"different source fixture identity under zero motion");
    previousSourceIdentity=shape.configSha256();
    pc.shapeWeight=0;pc.displacementWeight=0;pc.shapeOperatorConfig.maximumVertices=1;
    const SourceCagePaperFunctional disabledTiny(moment,pc);
    require(disabledTiny.shapeFunctional()==nullptr && disabledTiny.evaluateValues(state).shapeEnergy==0.,"disabled weights never construct LD under insufficient vertex cap");
    pc.shapeWeight=-1;
    throwsAs<std::invalid_argument>([&]{SourceCagePaperFunctional invalid(moment,pc);},"invalid paper regularizer weight");
    std::cout<<"actual source/BBW identity and exact wrapper caps PASS\n";
}
}
int main() {try {mathematical();actualSource(false);actualSource(true);std::cout<<"SYNTHETIC_REGRESSION source cage shape PASS\n";return 0;}catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}}
#endif
