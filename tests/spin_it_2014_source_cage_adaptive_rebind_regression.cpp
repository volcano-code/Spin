#include <iostream>
#if !__has_include("spin2014/SourceCageAdaptiveRebind.h")
int main() {
    std::cerr << "FAIL: actual source adaptive rebind cannot publish a freshly checked immutable checkpoint\n";
    return 1;
}
#else
#include "spin2014/SourceCageAdaptiveRebind.h"
#include "spin2014/SourceCageGuardedSlqp.h"
#include "source_surface_test_utils.h"
#include <Eigen/Geometry>
#include <cstring>
#include <iomanip>
#include <limits>
#include <type_traits>
using namespace spin2014;
using namespace surface_test;
namespace {
// Breaks caught: source/BBW replacement, wrong partition remap, stale geometry,
// mutable authorization, absent derivative reconstruction, inherited solver history.
SourceCagePaperFunctional fixture(Spin2014Objective objective=Spin2014Objective::TopMassAndLowering,
    std::size_t leaves=1000,std::size_t primitives=1000) {
    const auto outer=cube(.3125,.6875),design=cube(.36,.515);
    SourceSurfaceOctreeConfig sc;sc.maximumDepth=4;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
    SourceAdaptiveDomain domain(source,3,4,leaves);
    SourceSurfaceCageConfig cc;cc.cageLevel=2;SourceCageGridConfig gc;gc.fineLevel=3;
    SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),{}),{});
    require(bbw.success(),"actual full BBW succeeds");
    SourceCageBindingConfig bc;bc.maximumPrimitives=primitives;
    SourceConformingBoundaryConfig boundary;boundary.maximumPatches=9000;
    SourceCageMomentFunctionalConfig mc;mc.maximumWeightValues=9000;
    SourceCagePaperFunctionalConfig pc;pc.objective=objective;pc.fillWeight=.03;
    pc.shapeWeight=.02;pc.displacementWeight=.01;pc.centerWeight=.6;pc.inertiaWeight=.8;
    SourceCagePaperFrame frame;frame.origin=Eigen::Vector3d(.48,.51,.3);
    frame.rotation=Eigen::AngleAxisd(.13,Eigen::Vector3d::UnitZ()).toRotationMatrix();
    return SourceCagePaperFunctional(SourceCageMomentFunctional(
        SourceConformingBoundary(SourceCageBinding(domain,bbw,bc),boundary),1.3,mc),pc,frame);
}
bool identical(const Eigen::VectorXd& a,const Eigen::VectorXd& b) {
    return a.size()==b.size()&&std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0;
}
Ten tetra(const Eigen::MatrixXd& points,SourceConformingBoundary::Triangle t) {
    Eigen::Matrix3d v;for(int i=0;i<3;++i)v.row(i)=points.row(t[i]);
    const double volume=v.determinant()/6.;const Eigen::Vector3d sum=v.colwise().sum();
    const Eigen::Matrix3d q=volume*(sum*sum.transpose()+v.transpose()*v)/20.;
    Ten r;r<<volume,volume*sum/4.,q(0,0),q(1,1),q(2,2),q(1,2),q(0,2),q(0,1);return r;
}
Ten integrated(const SourceCagePaperFunctional& f,const Eigen::VectorXd& x) {
    const auto& b=f.source().boundary();const auto e=f.evaluateValues(x);Ten r=Ten::Zero();
    const auto& faces=b.binding().source().outerFaces();
    for(Eigen::Index i=0;i<faces.rows();++i)
        r+=tetra(e.world.deformedPoints,{{std::size_t(faces(i,0)),std::size_t(faces(i,1)),std::size_t(faces(i,2))}});
    for(std::size_t p=0;p<b.primitiveBoundaries().size();++p)
        for(const auto& t:b.primitiveBoundaries()[p])
            r-=x[b.binding().primitives()[p].betaIndex]*t.orientation*tetra(e.world.deformedPoints,b.triangles()[t.triangle]);
    return f.source().density()*r;
}
void failed(const SourceCageAdaptiveRebindAudit& a,SourceCageAdaptiveRebindOutcome outcome,
            SourceCageAdaptiveRebindStage stage) {
    require(!a.checkpoint()&&a.facts().outcome==outcome,"failure cannot publish usable checkpoint");
    require(a.facts().lastStage==stage,"failure retains exact transaction stage");
}
void verify(const SourceCagePaperFunctional& old,const Eigen::VectorXd& x,const SourceCageAdaptiveRebindAudit& a) {
    require(a.checkpoint()!=nullptr,"actual rebind succeeds");
    const auto& facts=a.facts();const auto& f=a.checkpoint()->functional();const auto& y=a.checkpoint()->state();
    require(facts.inputGeometry&&facts.outputGeometry&&
        facts.inputGeometry->facts().outcome==SourceCageGeometryOutcome::ChecksPassed&&
        facts.outputGeometry->facts().outcome==SourceCageGeometryOutcome::ChecksPassed,"both actual geometry gates completed");
    require(old.handleCount()==f.handleCount()&&identical(x.tail(4*old.handleCount()),y.tail(4*f.handleCount())),"all handle bytes preserved");
    require(old.hasAngle()==f.hasAngle()&&(!old.hasAngle()||x[old.angleOffset()]==y[f.angleOffset()]),"full phi uses actual offsets");
    const auto& oc=old.config();const auto& nc=f.config();
    require(oc.objective==nc.objective&&oc.inertiaWeight==nc.inertiaWeight&&oc.centerWeight==nc.centerWeight&&
        oc.fillWeight==nc.fillWeight&&oc.shapeWeight==nc.shapeWeight&&oc.displacementWeight==nc.displacementWeight&&
        oc.maximumJointDerivativeValues==nc.maximumJointDerivativeValues&&
        oc.maximumShapeHandleGradientValues==nc.maximumShapeHandleGradientValues&&oc.maximumShapeWeightValues==nc.maximumShapeWeightValues&&
        oc.shapeOperatorConfig.maximumVertices==nc.shapeOperatorConfig.maximumVertices&&
        oc.shapeOperatorConfig.maximumOneRingVertices==nc.shapeOperatorConfig.maximumOneRingVertices&&
        oc.shapeOperatorConfig.maximumSparseValues==nc.shapeOperatorConfig.maximumSparseValues&&
        oc.shapeOperatorConfig.maximumConstructionWorkspaceDoubles==nc.shapeOperatorConfig.maximumConstructionWorkspaceDoubles&&
        oc.shapeOperatorConfig.maximumEvaluationWorkspaceDoubles==nc.shapeOperatorConfig.maximumEvaluationWorkspaceDoubles&&
        oc.shapeOperatorConfig.relativeRankTolerance==nc.shapeOperatorConfig.relativeRankTolerance&&
        old.source().config().maximumWeightValues==f.source().config().maximumWeightValues&&
        old.source().config().maximumJacobianValues==f.source().config().maximumJacobianValues&&
        old.source().boundary().config().maximumFaceTreeNodes==f.source().boundary().config().maximumFaceTreeNodes&&
        old.source().boundary().config().maximumPatches==f.source().boundary().config().maximumPatches&&
        old.source().boundary().config().maximumPoints==f.source().boundary().config().maximumPoints&&
        old.source().boundary().config().maximumTriangles==f.source().boundary().config().maximumTriangles&&
        old.source().boundary().binding().configSha256()==f.source().boundary().binding().configSha256(),"immutable configs retained");
    require(old.source().density()==f.source().density()&&old.frame().origin==f.frame().origin&&old.frame().rotation==f.frame().rotation,"density and frame retained");
    require(&old.source().boundary().binding().bbw().finalWeights()==&f.source().boundary().binding().bbw().finalWeights(),"exact old BBW shared");
    require(&old.source().boundary().binding().source()==&f.source().boundary().binding().source(),"exact source retained");
    const auto e=f.evaluate(y);require(e.hasJacobians&&e.jointGradient.size()==y.size()&&e.spinMomentJacobian.cols()==y.size(),"rebuilt derivatives have new dimension");
    require(e.outputSha256==facts.afterMaterial.outputSha256,"retained material output belongs to rebuilt full evaluation");
    close(facts.beforeMaterial.worldMoments.values(),integrated(old,x),"independent before deformed integration");
    close(facts.afterMaterial.worldMoments.values(),integrated(f,y),"independent after deformed integration");
    require(facts.deformedMomentDifference==facts.afterMaterial.worldMoments.values()-facts.beforeMaterial.worldMoments.values(),"honest deformed difference");
    require(facts.inputGeometry->facts().sourceOutputSha256==old.evaluateValues(x).outputSha256&&
        facts.outputGeometry->facts().sourceOutputSha256==f.evaluateValues(y).outputSha256,"fresh geometry hashes bind actual values");
}
}
int main(int argc,char** argv) {try {
    static_assert(!std::is_default_constructible_v<SourceCageAdaptiveRebindAudit>);
    static_assert(!std::is_default_constructible_v<SourceCageAdaptiveCheckpoint>);
    static_assert(std::is_const_v<std::remove_reference_t<decltype(std::declval<SourceCageAdaptiveCheckpoint>().state())>>);
    const auto old=fixture();const auto x=old.initialState(Eigen::VectorXd::Constant(1,.25),.23);
    const auto saved=x;
    const auto before=old.evaluate(x).outputSha256;
    const auto a=rebindSourceCageAdaptive(old,x);verify(old,x,a);
    const auto& f=a.checkpoint()->functional();const auto& y=a.checkpoint()->state();
    if(argc>1&&std::string(argv[1])=="--probe-clearance") {
        auto clearance=[](const SourceCageGeometryAudit& audit) {
            long double minimum=1;
            for(const auto& r:audit.facts().records)if(!r.squaredDistanceExact.empty()) {
                const auto slash=r.squaredDistanceExact.find('/');
                const auto value=std::stold(r.squaredDistanceExact.substr(0,slash))/
                    (slash==std::string::npos?1:std::stold(r.squaredDistanceExact.substr(slash+1)));
                minimum=std::min(minimum,value);
            }return std::sqrt(minimum);
        };
        for(int mode=1;mode<=44;++mode) {
            auto px=x;const double amplitude=mode<=12?.02*(1+(mode-1)/4):.1;
            for(std::size_t h=0;h<old.handleCount();++h)for(int k=0;k<3;++k)
                px[old.handleOffset()+4*h+k+1]=amplitude*std::cos((h+1)*(k+1)*mode);
            const auto probe=rebindSourceCageAdaptive(old,px);
            if(probe.checkpoint())std::cout<<std::setprecision(20)<<"CLEARANCE mode="<<mode<<" amp="<<amplitude
                <<" before="<<clearance(*probe.facts().inputGeometry)<<" after="<<clearance(*probe.facts().outputGeometry)<<std::endl;
            else std::cout<<"CLEARANCE mode="<<mode<<" failure="<<int(probe.facts().outcome)<<" stage="<<int(probe.facts().lastStage)<<std::endl;
        }return 0;
    }
    require(f.betaCount()==8&&old.betaCount()==1&&f.source().boundary().binding().primitives().size()==8,"one full cube splits into eight actual primitives");
    require(a.facts().split.changes.size()==1&&a.facts().merge.changes.empty(),"fractional split without merge");
    require(a.facts().oldLeaves.size()==1&&a.facts().oldLeaves[0].level==3&&a.facts().oldLeaves[0].coordinates==std::array<std::uint32_t,3>{3,3,3},"old source leaf identity");
    for(std::size_t i=0;i<8;++i) {
        const auto& key=a.facts().newLeaves[i];require(key.level==4&&y[i]==.25,"inherited beta and refined level");
        for(auto c:key.coordinates)require(c==6||c==7,"child occupies actual parent octant");
        const auto& p=f.source().boundary().binding().primitives()[i];
        require(p.betaIndex==i&&p.sourceKey==key,"rebuilt primitive ownership matches ordered adaptive leaf");
    }
    require(f.source().boundary().points().rows()>old.source().boundary().points().rows()&&
        f.source().boundary().binding().outputSha256()!=old.source().boundary().binding().outputSha256(),"rebuild samples changed represented geometry");
    const Ten expected=1.3*(boxMoments(Eigen::Vector3d::Constant(.3125),Eigen::Vector3d::Constant(.6875))-
        .25*boxMoments(Eigen::Vector3d::Constant(.375),Eigen::Vector3d::Constant(.5)));
    close(a.facts().beforeMaterial.worldMoments.values(),expected,"independent analytic rest moments before");
    close(a.facts().afterMaterial.worldMoments.values(),expected,"independent analytic rest moments after");
    const Ten removed=.25*boxMoments(Eigen::Vector3d::Constant(.375),Eigen::Vector3d::Constant(.5));
    close(a.facts().split.removedBefore.values(),removed,"transition before contains rest geometric moments without density");
    close(a.facts().split.removedAfter.values(),removed,"transition after conserves independent rest moments");
    require(identical(x,saved)&&old.evaluate(x).outputSha256==before,"old state and derived objects intact");
    const auto repeat=rebindSourceCageAdaptive(old,x);require(repeat.facts().scientificSha256==a.facts().scientificSha256,"deterministic transaction hash");
    auto altered=x;altered[0]=.3;require(rebindSourceCageAdaptive(old,altered).facts().scientificSha256!=a.facts().scientificSha256,"input hash changes");
    auto moved=a;const auto& retained=moved.facts();const auto movedTo=std::move(moved);
    require(retained.scientificSha256==a.facts().scientificSha256&&movedTo.checkpoint(),"shared immutable evidence survives move without stale emptied pass");
    rejects([&]{moved.facts();},"moved-from accessor cannot read stale empty evidence");
    auto binary=y;binary.head(8).setConstant(.01);
    const auto merged=rebindSourceCageAdaptive(f,binary);verify(f,binary,merged);
    require(merged.checkpoint()->functional().betaCount()==1&&merged.facts().merge.changes.size()==4&&
        merged.facts().newLeaves[0].level==0&&merged.facts().newBeta[0]==.01,"equal near-bound siblings merge recursively with exact beta");
    auto unequal=binary;unequal[0]=.02;
    const auto unchanged=rebindSourceCageAdaptive(f,unequal);verify(f,unequal,unchanged);
    require(unchanged.facts().merge.changes.empty()&&unchanged.facts().split.changes.empty()&&identical(unchanged.checkpoint()->state(),unequal),"unequal near-bound siblings never merge; no-change still freshly checked");
    const auto& rebuiltDomain=unchanged.checkpoint()->functional().source().boundary().binding().domain();
    require(rebuiltDomain.edges().size()==12&&rebuiltDomain.uniformLaplacian().rows()==8,"new cubical adjacency has twelve independent face neighbors");
    require(std::abs(unchanged.checkpoint()->functional().evaluate(unequal).fillPenalty-4.5e-6)<1e-18,"rebuilt fill penalty from three unequal incident edges");
    const auto stat=fixture(Spin2014Objective::StaticLowering);
    const auto sx=stat.initialState(Eigen::VectorXd::Constant(1,.25),0.);
    const auto sa=rebindSourceCageAdaptive(stat,sx);verify(stat,sx,sa);
    require(!sa.checkpoint()->functional().hasAngle()&&sa.checkpoint()->state().size()==8+4*stat.handleCount(),"static layout never copies imaginary angle");
    auto nonaffine=x;
    for(std::size_t h=0;h<old.handleCount();++h) {
        nonaffine[old.handleOffset()+4*h]=1.+.001*std::sin(h+1);
        for(int k=0;k<3;++k)nonaffine[old.handleOffset()+4*h+k+1]=.001*std::cos((h+1)*(k+1));
    }
    const auto na=rebindSourceCageAdaptive(old,nonaffine);verify(old,nonaffine,na);
    auto affine=x;const auto& cage=old.source().boundary().binding().bbw().fem().grid().cage();
    for(std::size_t h=0;h<old.handleCount();++h) {
        affine[old.handleOffset()+4*h]=1.07;
        affine.segment<3>(old.handleOffset()+4*h+1)=.07*cage.boundaryVertices().row(cage.handleVertexIndices()[h]).transpose()+Eigen::Vector3d(.01,-.02,.03);
    }
    const auto aa=rebindSourceCageAdaptive(old,affine);verify(old,affine,aa);
    close(aa.facts().beforeMaterial.worldMoments.values(),a.facts().beforeMaterial.worldMoments.transformed(
        1.07*Eigen::Matrix3d::Identity(),Eigen::Vector3d(.01,-.02,.03)).values(),"global affine analytic covariance before");
    close(aa.facts().afterMaterial.worldMoments.values(),aa.facts().beforeMaterial.worldMoments.values(),"globally affine ten-moment conservation");
    std::cout<<std::setprecision(17)<<"NONAFFINE before "<<na.facts().beforeMaterial.worldMoments.values().transpose()
        <<"\nNONAFFINE after "<<na.facts().afterMaterial.worldMoments.values().transpose()
        <<"\nNONAFFINE delta "<<na.facts().deformedMomentDifference.transpose()<<'\n';
    SourceCageAdaptiveRebindSettings settings;settings.geometry.minimumOuterDomainClearance=.2;
    failed(rebindSourceCageAdaptive(old,x,settings),SourceCageAdaptiveRebindOutcome::RejectedGeometry,SourceCageAdaptiveRebindStage::InputGeometry);
    settings={};settings.geometry.maximumPrimitives=1;
    failed(rebindSourceCageAdaptive(old,x,settings),SourceCageAdaptiveRebindOutcome::ResourceLimit,SourceCageAdaptiveRebindStage::OutputGeometry);
    settings={};settings.budget.maximumOwnedBytes=1;
    failed(rebindSourceCageAdaptive(old,x,settings),SourceCageAdaptiveRebindOutcome::ResourceLimit,SourceCageAdaptiveRebindStage::Preflight);
    const auto leafLimited=fixture(Spin2014Objective::StaticLowering,1);
    failed(rebindSourceCageAdaptive(leafLimited,leafLimited.initialState(Eigen::VectorXd::Constant(1,.25),0.)),
        SourceCageAdaptiveRebindOutcome::ResourceLimit,SourceCageAdaptiveRebindStage::Split);
    const auto bindingLimited=fixture(Spin2014Objective::StaticLowering,1000,1);
    failed(rebindSourceCageAdaptive(bindingLimited,bindingLimited.initialState(Eigen::VectorXd::Constant(1,.25),0.)),
        SourceCageAdaptiveRebindOutcome::ResourceLimit,SourceCageAdaptiveRebindStage::Rebuild);
    auto mc=old.source().config();mc.maximumJacobianValues=1055;
    const SourceCagePaperFunctional outputMaterialLimited(SourceCageMomentFunctional(old.source().boundary(),old.source().density(),mc),old.config(),old.frame());
    failed(rebindSourceCageAdaptive(outputMaterialLimited,x),SourceCageAdaptiveRebindOutcome::ResourceLimit,SourceCageAdaptiveRebindStage::OutputMaterial);
    mc.maximumJacobianValues=1;
    const SourceCagePaperFunctional inputMaterialLimited(SourceCageMomentFunctional(old.source().boundary(),old.source().density(),mc),old.config(),old.frame());
    failed(rebindSourceCageAdaptive(inputMaterialLimited,x),SourceCageAdaptiveRebindOutcome::ResourceLimit,SourceCageAdaptiveRebindStage::InputMaterial);
    auto huge=x;for(std::size_t h=0;h<old.handleCount();++h)huge[old.handleOffset()+4*h]=1e200;
    failed(rebindSourceCageAdaptive(old,huge),SourceCageAdaptiveRebindOutcome::Unknown,SourceCageAdaptiveRebindStage::InputGeometry);
    rejects([&]{rebindSourceCageAdaptive(old,Eigen::VectorXd::Zero(1));},"invalid dimension");
    for(auto offset:{std::size_t(0),old.angleOffset(),old.handleOffset(),old.handleOffset()+1}) {
        auto nan=x;nan[offset]=std::numeric_limits<double>::quiet_NaN();rejects([&]{rebindSourceCageAdaptive(old,nan);},"NaN full state field");
    }
    for(double epsilon:{0.,.5,-.1,std::numeric_limits<double>::quiet_NaN()}) {
        settings={};settings.fractionalEpsilon=epsilon;rejects([&]{rebindSourceCageAdaptive(old,x,settings);},"invalid epsilon");
    }
    settings={};settings.fractionalEpsilon=.1;
    require(rebindSourceCageAdaptive(old,x,settings).facts().scientificSha256!=a.facts().scientificSha256,"changed settings bind hash even when same partition");
    SourceCageGuardedSlqpSettings solver;solver.source.solver.maximumIterations=1;solver.source.solver.maximumLineSearchSteps=2;
    const auto solved=solveSourceCageGuardedSlqp(f,y,solver);const auto& sf=solved.facts();
    require(identical(sf.source.initialState,y)&&sf.source.numericalResultAvailable&&!sf.source.numerical.trace.empty()&&
        sf.source.numerical.trace.front().historyBefore==0,"new solve consumes remapped dimension with empty initial history");
    require(sf.attempts.front().checkCompleted&&sf.attempts.front().outcome==SourceCageGeometryOutcome::ChecksPassed&&
        sf.finalRecheckCount==1&&sf.finalStatus==SourceCageGuardedFinalStatus::ChecksPassed,"solve independently rechecks initial and final actual geometry");
    require(sf.feasibility==SourceCageSlqpFeasibility::Unknown&&sf.source.finalEvaluation.outputSha256==f.evaluate(sf.source.numerical.finalState).outputSha256,"fresh material hash grants no numerical feasibility");
    std::cout<<"REBINDED_SOLVE status="<<slqpStatusToken(sf.source.numerical.status)<<" hash="<<sf.scientificSha256<<'\n';
    require(identical(x,saved)&&old.evaluate(x).outputSha256==before,"all success/failure calls leave original intact");
    std::cout<<"atomic actual-source adaptive rebind PASS\n";return 0;
} catch(const std::exception& e) {std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}}
#endif
