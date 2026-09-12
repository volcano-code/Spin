#include <iostream>
#if !__has_include("spin2014/SourceCageGeometry.h")
int main() {std::cerr<<"FAIL: actual-source represented geometry admission is not implemented\n";return 1;}
#else
#include "spin2014/SourceCageGeometry.h"
#include "verification/ExactBoundaryAudit.h"
#include "source_surface_test_utils.h"
#include <type_traits>
#include <limits>
#include <cstring>
#include <map>
#include <Eigen/QR>
using namespace spin2014;
using namespace surface_test;
namespace {
SourceConformingBoundary fixture(bool hollow=false,bool disconnected=false) {
    auto outer=cube(.3125,.6875),design=cube(.34,.66);
    if(hollow) {outer=join(outer,cube(.48,.52),true);design=box({.34,.34,.34},{.46,.66,.66});}
    if(disconnected)design=join(cube(.34,.46),box({.54,.34,.34},{.66,.46,.46}));
    SourceSurfaceOctreeConfig sc;sc.maximumDepth=(hollow||disconnected)?4:3;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
    SourceAdaptiveDomain domain(source,2,sc.maximumDepth);
    SourceSurfaceCageConfig cc;cc.cageLevel=2;SourceCageGridConfig gc;gc.fineLevel=3;
    SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),{}),{});
    require(bbw.success(),"actual full BBW succeeds");
    return SourceConformingBoundary(SourceCageBinding(domain,bbw,{}),{});
}
void near(double a,double b,const char* label) {require(std::abs(a-b)<1e-11,label);}
enum class MotionFailure {SelfIntersection,NegativeOrientation,Protrusion};
struct InvalidMotion {
    int handle;
    double translation;
    SourceCageGeometryStage stage;
    MotionFailure failure;
};
// The displayed reason is truncated. Derive complete exact boundary facts from
// the same actual output, preserving all primitive references and their winding.
// First-encounter local indexing is independent of the checker's sorted mapping.
void checkMotionFailure(const SourceCagePaperFunctional& source,const Eigen::VectorXd& state,
                        const SourceCageGeometryAudit& rejected,MotionFailure expected) {
    const auto& record=rejected.facts().records.back();
    const std::string reason=record.reason.data();
    if(expected==MotionFailure::Protrusion) {
        require(reason=="primitive not subset of actual outer material","protrusion retains exact subset-failure reason");
        return;
    }
    const auto evaluation=source.evaluateValues(state);
    require(evaluation.outputSha256==rejected.facts().sourceOutputSha256,"failure-category audit uses the same actual source output");
    const auto& boundary=source.source().boundary();
    Mesh mesh;
    if(record.stage==SourceCageGeometryStage::OuterRegion) {
        mesh.v=evaluation.world.deformedPoints.topRows(boundary.outerPointCount());
        mesh.f=boundary.binding().source().outerFaces();
    } else {
        const auto& refs=boundary.primitiveBoundaries().at(record.first);
        std::map<std::size_t,int> local;
        mesh.f.resize(refs.size(),3);
        for(std::size_t face=0;face<refs.size();++face) {
            auto triangle=boundary.triangles().at(refs[face].triangle);
            if(refs[face].orientation<0)std::swap(triangle[1],triangle[2]);
            for(int corner=0;corner<3;++corner) {
                const auto entry=local.emplace(triangle[corner],int(local.size()));
                mesh.f(face,corner)=entry.first->second;
            }
        }
        mesh.v.resize(local.size(),3);
        for(const auto& entry:local)mesh.v.row(entry.second)=evaluation.world.deformedPoints.row(entry.first);
    }
    const auto exact=verification::auditExactTriangleBoundary(mesh.v,mesh.f,verification::ExactBoundaryOrientation::AnyNonzero);
    const auto& facts=exact.facts();
    require(facts.intersectionCompleted,"full exact failure-category intersection scan completed");
    if(expected==MotionFailure::SelfIntersection) {
        require(!facts.intersectingPairs.empty(),"self-intersection fixture must have actual exact intersecting triangle pairs");
        require(reason.rfind("invalid exact boundary:",0)==0,"self-intersection retains boundary-failure reason");
    } else {
        require(facts.boundaryValid&&facts.intersectionFree&&facts.componentVolumeSigns==std::vector<int>{-1},
                "orientation fixture must be a valid intersection-free negative-volume boundary");
        require(reason=="component 0 orientation violates positive roots and alternating nesting",
                "negative orientation retains signed-nesting failure reason");
    }
}
}
int main(int argc,char** argv) {try {
    static_assert(!std::is_default_constructible_v<SourceCageGeometryAudit>);
    static_assert(!std::is_constructible_v<SourceCageGeometryAudit,SourceCageGeometryFacts>);
    static_assert(std::is_copy_constructible_v<SourceCageGeometryAudit>);
    const bool overlapProbe=argc==2&&!std::strcmp(argv[1],"--probe-overlap");
    const auto boundary=fixture(false,overlapProbe);const SourceCageMomentFunctional moments(boundary,1.,{});
    SourceCagePaperFunctionalConfig pc;pc.objective=Spin2014Objective::TopInertiaOnly;
    SourceCagePaperFunctional f(moments,pc);SourceCageGeometrySettings settings;
    const SourceCageGeometryChecker checker(f,settings);
    auto x=f.initialState(Eigen::VectorXd::Constant(f.betaCount(),.25),.1);
    if(argc>1) {
        require(argc==2&&(!std::strcmp(argv[1],"--probe-motion")||!std::strcmp(argv[1],"--probe-collapse")||overlapProbe),"unknown probe argument");
        if(overlapProbe) {
            const auto w=boundary.binding().pointWeights(0,boundary.outerPointCount());
            for(int selected:{0,4,10,12,13,15,21,25})for(double amplitude:{1.,2.,-1.,-2.}) {
                Eigen::VectorXd direction=Eigen::VectorXd::Zero(26);direction[selected]=1.;direction-=w.completeOrthogonalDecomposition().solve(w.col(selected));
                auto z=x;z.head(f.betaCount()).setZero();for(int h=0;h<26;++h)z[f.handleOffset()+4*h+1]+=amplitude*direction[h];
                const auto r=checker.check(z);const auto& last=r.facts().records.back();
                std::cout<<"overlap selected="<<selected<<" amplitude="<<amplitude<<" primitives="<<r.facts().primitiveCount<<" outcome="<<sourceCageGeometryOutcomeToken(r.facts().outcome)<<" stage="<<int(last.stage)<<" id="<<last.first<<" reason="<<last.reason.data()<<" digest="<<r.facts().scientificSha256<<std::endl;
            }return 0;
        }
        if(std::strcmp(argv[1],"--probe-collapse")==0) {
            auto cfg=pc;cfg.objective=Spin2014Objective::StaticLowering;SourceCagePaperFunctional sf(moments,cfg);
            for(int power=45;power<=55;++power) {
                auto z=sf.initialState(Eigen::VectorXd::Zero(sf.betaCount()),0.);const auto& cage=boundary.binding().bbw().fem().grid().cage();
                const double scale=std::ldexp(1.,-power);
                for(std::size_t h=0;h<26;++h) {
                    const Eigen::Vector3d center=cage.boundaryVertices().row(cage.handleVertexIndices()[h]);
                    z[sf.handleOffset()+4*h]=scale;z.segment<3>(sf.handleOffset()+4*h+1)=(scale-1)*center+Eigen::Vector3d::Constant(.5);
                }
                const auto r=SourceCageGeometryChecker(sf,settings).check(z);const auto& last=r.facts().records.back();
                std::cout<<"collapse power="<<power<<" outcome="<<sourceCageGeometryOutcomeToken(r.facts().outcome)<<" stage="<<int(last.stage)<<" id="<<last.first<<" reason="<<last.reason.data()<<" digest="<<r.facts().scientificSha256<<std::endl;
            }return 0;
        }
        for(std::size_t h=0;h<26;++h)for(double d:{.5,1.,-1.}) {
            auto z=x;z.head(f.betaCount()).setZero();z[f.handleOffset()+4*h+1]+=d;
            const auto r=checker.check(z);const auto& last=r.facts().records.back();
            std::cout<<"probe h="<<h<<" d="<<d<<" outcome="<<sourceCageGeometryOutcomeToken(r.facts().outcome)<<" stage="<<int(last.stage)<<" id="<<last.first<<" reason="<<last.reason.data()<<std::endl;
        }return 0;
    }
    const auto a=checker.check(x);const auto& facts=a.facts();
    require(facts.outcome==SourceCageGeometryOutcome::ChecksPassed,"identity geometry passes all represented checks");
    require(facts.sourceOutputSha256==f.evaluateValues(x).outputSha256,"actual source output identity");
    require(facts.primitiveCount==8&&facts.pairCount==28&&facts.pairsVisited==28,"complete ordered pair coverage");
    require(facts.scientificSha256==checker.check(x).facts().scientificSha256,"deterministic repeated scientific evidence");
    const auto copied=a;require(copied.facts().scientificSha256==facts.scientificSha256,"typed audit copy retains immutable evidence");
    std::size_t minimumDiagnosticBytes=0;
    for(const auto& refs:boundary.primitiveBoundaries())minimumDiagnosticBytes+=4*(refs.size()*(refs.size()-1)/2)*sizeof(std::array<int,2>);
    require(facts.plannedCheckerBytes>=minimumDiagnosticBytes,"account simultaneous exact diagnostics and vector capacity, including all possible intersection pairs");
    near(f.evaluateValues(x).world.finalMoments.values()[0],.052734375-.25*.015625,"analytic material volume");
    auto beta=x;beta[0]=.4;beta[f.angleOffset()]=.3;
    const auto b=checker.check(beta);
    require(b.facts().outcome==facts.outcome&&b.facts().scientificSha256!=facts.scientificSha256,"beta/phi preserve geometry and bind full state");
    auto moved=x;for(std::size_t h=0;h<f.handleCount();++h) {
        const auto& cage=boundary.binding().bbw().fem().grid().cage();
        const Eigen::Vector3d center=cage.boundaryVertices().row(cage.handleVertexIndices()[h]);
        moved[f.handleOffset()+4*h]=2.;moved.segment<3>(f.handleOffset()+4*h+1)=center+Eigen::Vector3d(.125,-.25,.375);
    }
    require(checker.check(moved).facts().outcome==SourceCageGeometryOutcome::ChecksPassed,"common scale and translation preserve validity");
    near(f.evaluateValues(moved).world.finalMoments.values()[0],8*(.052734375-.25*.015625),"scaled analytic material volume");
    auto nonaffine=x;nonaffine[f.handleOffset()+1]+=.001;
    require(checker.check(nonaffine).facts().outcome==SourceCageGeometryOutcome::ChecksPassed,"real nonaffine full control motion valid");
    auto low=settings;low.maximumPoints=boundary.points().rows()-1;
    require(SourceCageGeometryChecker(f,low).check(x).facts().outcome==SourceCageGeometryOutcome::ResourceLimit,"point cap cannot admit");
    low=settings;low.maximumPairs=27;
    require(SourceCageGeometryChecker(f,low).check(x).facts().outcome==SourceCageGeometryOutcome::ResourceLimit,"pair cap cannot admit");
    low=settings;low.maximumNefComplexity=1;
    const auto exhausted=SourceCageGeometryChecker(f,low).check(x);
    require(exhausted.facts().outcome==SourceCageGeometryOutcome::ResourceLimit,"native Nef resource exception remains resource");
    low=settings;low.minimumOuterDomainClearance=.2;
    require(SourceCageGeometryChecker(f,low).check(x).facts().outcome==SourceCageGeometryOutcome::RejectedGeometry,"requested conservative clearance rejects");
    low=settings;low.maximumNefComplexity=std::size_t(std::numeric_limits<int>::max())+1;rejects([&]{SourceCageGeometryChecker c(f,low);},"unrepresentable exact backend cap rejected early");
    // Controller's independent Fraction triangle-distance oracle derived this
    // literal from actual output, and proved its square equals the minimum.
    low=settings;low.minimumOuterDomainClearance=0x1.fffffffffffc0p-5;
    require(SourceCageGeometryChecker(f,low).check(x).facts().outcome==SourceCageGeometryOutcome::ChecksPassed,"exact requested clearance equality admitted");
    low.minimumOuterDomainClearance=0x1.fffffffffffc1p-5;
    require(SourceCageGeometryChecker(f,low).check(x).facts().outcome==SourceCageGeometryOutcome::RejectedGeometry,"one binary64 above exact clearance rejected");
    for(const auto& record:facts.records)if(record.stage==SourceCageGeometryStage::OuterDomainClearance) {
        require(!record.squaredDistanceExact.empty()&&record.witnessFaces.has_value(),"exact domain clearance exported");
        if(record.first>=4)require(record.squaredDistanceExact=="19807040628565802923409276929/5070602400912917605986812821504","independent exact minimum boundary distance");
    }
    // Missing preflight resource checks, truncated pair loops, or off-by-one
    // limits must fail these boundary cases. Caps bind actual full dimensions.
    low=settings;low.maximumPoints=facts.pointCount;low.maximumTriangles=facts.triangleCount;
    low.maximumPrimitives=facts.primitiveCount;low.maximumPairs=facts.pairCount;
    low.maximumRecords=facts.plannedRecords;low.maximumCheckerBytes=facts.plannedCheckerBytes;
    low.maximumNefComplexity=facts.nefBudget.peakComplexity;
    require(SourceCageGeometryChecker(f,low).check(x).facts().outcome==SourceCageGeometryOutcome::ChecksPassed,"all exact positive caps pass");
    for(auto member:{&SourceCageGeometrySettings::maximumTriangles,&SourceCageGeometrySettings::maximumPrimitives,
        &SourceCageGeometrySettings::maximumRecords,&SourceCageGeometrySettings::maximumCheckerBytes,&SourceCageGeometrySettings::maximumNefComplexity}) {
        auto tooSmall=low;--(tooSmall.*member);const auto limited=SourceCageGeometryChecker(f,tooSmall).check(x);
        require(limited.facts().outcome==SourceCageGeometryOutcome::ResourceLimit,"one-below resource cap never admits");
        require(limited.facts().records.back().outcome==SourceCageGeometryStageOutcome::ResourceLimit,"resource terminal stage retained");
    }
    low=settings;low.maximumExactStringBytes=1;
    require(SourceCageGeometryChecker(f,low).check(x).facts().outcome==SourceCageGeometryOutcome::ResourceLimit,"retained exact-string cap cannot admit");
    for(const auto& test:std::array<InvalidMotion,4>{{
        {2,1,SourceCageGeometryStage::OuterRegion,MotionFailure::SelfIntersection},
        {3,1,SourceCageGeometryStage::PrimitiveRegion,MotionFailure::SelfIntersection},
        {4,1,SourceCageGeometryStage::PrimitiveRegion,MotionFailure::NegativeOrientation},
        {1,-1,SourceCageGeometryStage::PrimitiveSubset,MotionFailure::Protrusion}}}) {
        auto z=x;z.head(f.betaCount()).setZero();z[f.handleOffset()+4*test.handle+1]+=test.translation;
        const auto rejected=checker.check(z);
        require(rejected.facts().outcome==SourceCageGeometryOutcome::RejectedGeometry,"actual outer intersection / primitive intersection / inversion / protrusion rejected");
        require(rejected.facts().records.back().stage==test.stage,"actual rejection stage retained");
        checkMotionFailure(f,z,rejected,test.failure);
    }
    auto invalid=x;invalid[f.handleOffset()]=0;
    const auto unknown=checker.check(invalid);
    require(unknown.facts().outcome==SourceCageGeometryOutcome::Unknown&&unknown.facts().records.back().stage==SourceCageGeometryStage::SourceEvaluation,"source failure never fabricates geometry diagnostics");
    auto staticConfig=pc;staticConfig.objective=Spin2014Objective::StaticLowering;SourceCagePaperFunctional sf(moments,staticConfig);
    auto collapsed=sf.initialState(Eigen::VectorXd::Zero(sf.betaCount()),0.);
    const auto& cage=boundary.binding().bbw().fem().grid().cage();const double tiny=std::ldexp(1.,-49);
    for(std::size_t h=0;h<26;++h) {
        const Eigen::Vector3d center=cage.boundaryVertices().row(cage.handleVertexIndices()[h]);
        collapsed[sf.handleOffset()+4*h]=tiny;collapsed.segment<3>(sf.handleOffset()+4*h+1)=(tiny-1)*center+Eigen::Vector3d::Constant(.5);
    }
    const auto degeneracy=SourceCageGeometryChecker(sf,settings).check(collapsed);
    require(degeneracy.facts().outcome==SourceCageGeometryOutcome::Unknown&&degeneracy.facts().records.back().stage==SourceCageGeometryStage::SourceEvaluation,"actual rounded degeneracy caught by source retains UNKNOWN before geometry");
    for(const auto objective:{Spin2014Objective::StaticLowering,Spin2014Objective::TopMassAndLowering,Spin2014Objective::TopLoweringOnly,Spin2014Objective::YoYo}) {
        auto cfg=pc;cfg.objective=objective;SourceCagePaperFunctional other(moments,cfg);
        const auto z=other.initialState(Eigen::VectorXd::Constant(other.betaCount(),.25),other.hasAngle()?.1:0.);
        require(SourceCageGeometryChecker(other,settings).check(z).facts().outcome==SourceCageGeometryOutcome::ChecksPassed,"all five objectives use world material geometry");
    }
    SourceCagePaperFrame frame;frame.origin=Eigen::Vector3d(.5,.5,.3125);
    auto changedConfig=pc;changedConfig.shapeWeight=.03;changedConfig.displacementWeight=.02;
    SourceCagePaperFunctional framed(moments,changedConfig,frame);
    const auto framedAudit=SourceCageGeometryChecker(framed,settings).check(x);
    require(framedAudit.facts().outcome==SourceCageGeometryOutcome::ChecksPassed&&framedAudit.facts().inputSha256!=facts.inputSha256,"frame and optional full-source regularizers bind input without changing world geometry");
    const auto hole=fixture(true);require(hole.outerPointCount()==16&&hole.binding().handleCount()==26,"original void rows and full handles preserved");
    SourceCagePaperFunctional holeFunctional(SourceCageMomentFunctional(hole,1.,{}),pc);
    const auto holeState=holeFunctional.initialState(Eigen::VectorXd::Constant(holeFunctional.betaCount(),.25),.1);
    const auto holeAudit=SourceCageGeometryChecker(holeFunctional,settings).check(holeState);
    require(holeAudit.facts().outcome==SourceCageGeometryOutcome::ChecksPassed,"actual depth4 source preserves signed nested outer void");
    for(auto member:{&SourceCageGeometrySettings::maximumPoints,&SourceCageGeometrySettings::maximumTriangles,&SourceCageGeometrySettings::maximumPrimitives,
        &SourceCageGeometrySettings::maximumPairs,&SourceCageGeometrySettings::maximumRecords,&SourceCageGeometrySettings::maximumCheckerBytes,
        &SourceCageGeometrySettings::maximumNefComplexity,&SourceCageGeometrySettings::maximumConnectedShells,&SourceCageGeometrySettings::maximumExactStringBytes}) {
        low=settings;low.*member=0;rejects([&]{SourceCageGeometryChecker c(f,low);},"every zero resource cap rejected");
    }
    low=settings;low.maximumTriangles=std::numeric_limits<std::size_t>::max();rejects([&]{SourceCageGeometryChecker c(f,low);},"unrepresentable cap rejected");
    low=settings;low.minimumOuterDomainClearance=-1;rejects([&]{SourceCageGeometryChecker c(f,low);},"negative clearance rejected");
    rejects([&]{checker.check(x.head(3));},"bad state dimensions rejected");
    std::cout<<"actual-source geometry admission PASS digest="<<facts.scientificSha256<<"\n";return 0;
}catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;}}
#endif
