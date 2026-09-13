#include <iostream>
// Admission and geometry assertions below exercise the now-present layer.
#if !__has_include("spin2014/SourceCageMaterialRegions.h")
int main(){std::cerr<<"FAIL: actual binary exact region/partition/design membership/outer clearance layer absent\n";return 1;}
#else
#include "spin2014/SourceCageMaterialRegions.h"
#include "source_surface_test_utils.h"
#include <Eigen/Eigenvalues>
#include <cstring>
#include <limits>
#include <type_traits>
using namespace spin2014;
using namespace surface_test;
namespace {
bool identical(const Eigen::VectorXd& a,const Eigen::VectorXd& b) {
    return a.size()==b.size()&&std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0;
}
SourceCagePaperFunctional fixture(int level=4,bool wider=false) {
    const auto outer=cube(.3125,.6875),design=wider?box({.3625,.425,.425},{.6375,.575,.575}):cube(.425,.575);
    SourceSurfaceOctreeConfig sc;sc.maximumDepth=4;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
    SourceAdaptiveDomain domain(source,level,4,1000);
    SourceSurfaceCageConfig cc;cc.cageLevel=2;SourceCageGridConfig gc;gc.fineLevel=3;
    SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),{}),{});
    require(bbw.success(),"actual full BBW succeeds");
    SourceConformingBoundaryConfig bc;bc.maximumPatches=9000;
    SourceCageMomentFunctionalConfig mc;mc.maximumWeightValues=9000;
    SourceCagePaperFunctionalConfig pc;pc.objective=Spin2014Objective::TopMassAndLowering;
    pc.inertiaWeight=0.;pc.centerWeight=0.;
    return SourceCagePaperFunctional(SourceCageMomentFunctional(
        SourceConformingBoundary(SourceCageBinding(domain,bbw,{}),bc),1.17,mc),pc);
}
Eigen::VectorXd state(const SourceCagePaperFunctional& f,double beta) {
    auto x=f.initialState(Eigen::VectorXd::Constant(f.betaCount(),beta),.23);
    for(std::size_t h=0;h<f.handleCount();++h) {
        x[f.handleOffset()+4*h]=1.02;
        x.segment<3>(f.handleOffset()+4*h+1)<<.011,-.007,.019;
    }
    return x;
}
// KNOWN_SOLUTION_VALIDATION: fixed principal frame and zero objective weights.
// A nonzero Givens angle is retained, with the transverse basis counter-rotated.
SourceCagePaperFunctional stationary(const SourceCagePaperFunctional& f,const Eigen::VectorXd& x) {
    const auto m=f.evaluateValues(x).world.finalMoments.values();
    SourceCagePaperFrame frame;frame.origin=m.segment<3>(1)/m[0];
    Eigen::Matrix3d q;q<<m[4],m[9],m[8],m[9],m[5],m[7],m[8],m[7],m[6];
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(q-m[0]*frame.origin*frame.origin.transpose());
    Eigen::Matrix3d base=eig.eigenvectors().transpose();if(base.determinant()<0)base.row(0)*=-1;
    const double c=std::cos(x[f.angleOffset()]),s=std::sin(x[f.angleOffset()]);
    Eigen::Matrix3d g;g<<c,s,0,-s,c,0,0,0,1;frame.rotation=g*base;
    return SourceCagePaperFunctional(f.source(),f.config(),frame);
}
SourceCageAdaptiveSolverSettings solverSettings() {
    SourceCageAdaptiveSolverSettings s;s.maximumPasses=3;
    s.solve.source.solver.maximumIterations=1;s.solve.source.solver.maximumLineSearchSteps=2;
    s.solve.budget.maximumGuardCalls=8;return s;
}
SourceCageAdaptiveSolverAudit endpoint(const SourceCagePaperFunctional& f,const Eigen::VectorXd& x) {
    auto audit=solveSourceCageAdaptive(stationary(f,x),x,solverSettings());
    require(audit.facts().outcome==SourceCageAdaptiveSolverOutcome::RelaxedAdaptiveComplete,"fixture must actually converge and complete");
    return audit;
}

}
int main(){try {
    static_assert(!std::is_default_constructible_v<SourceCageMaterialRegionsAudit>);
    static_assert(!std::is_default_constructible_v<SourceCageMaterialRegions>);
    static_assert(!std::is_constructible_v<SourceCageMaterialRegionsAudit,SourceCageMaterialRegionsFacts>);
    const auto f=fixture();
    const auto ep=endpoint(f,state(f,1.));
    const auto binary=validateSourceCageBinaryMaterial(ep);
    require(bool(binary.validatedMaterial()),"real upstream validated material");
    auto settings=SourceCageMaterialRegionsSettings{};
    settings.minimumOuterWall=.1;
    auto audit=validateSourceCageMaterialRegions(binary,settings);
    const auto& a=audit.facts();
    require(a.outcome==SourceCageMaterialRegionsOutcome::Validated&&audit.validatedRegions(),"exact region chain passes");
    require(a.selectedCount==8&&a.selected.size()==8,"all eight immutable actual selected witnesses");
    require(a.materialPartition.completePartition==true&&a.selectedPartition.completePartition==true,"both exact partitions completed");
    require(a.outerClearance.applicable==true&&a.outerClearance.contact==false&&a.outerClearance.atLeastMinimum==true,"positive outer clearance");
    require(a.material.region.materialComponents==1&&a.material.region.boundedVoidComponents==1&&a.carved.region.materialComponents==1,"native material and void components");
    require(a.representedDesignMembershipEstablished&&a.evaluationSha256==binary.facts().binaryFunctionalOutputSha256,"constructive represented design membership and actual evaluation binding");
    require(audit.validatedRegions()->binaryMaterial().get()==binary.validatedMaterial().get(),"owned link preserves actual arrays and moments");
    require(&audit.validatedRegions()->binaryAudit().facts()==&binary.facts(),"product preserves immutable upstream integrated moments and all original evidence");
    std::cout<<"actual_nonidentity_squared_clearance="<<*a.outerClearance.squaredDistanceExact<<std::endl;
    for(const auto& p:a.selected)require(p.subsetOuter.aSubsetB==true&&p.meshSha256.size()==64&&p.region.region.valid,"selected actual primitives each subset outer");
    const auto repeat=validateSourceCageMaterialRegions(binary,settings);
    require(a.scientificSha256==repeat.facts().scientificSha256,"deterministic complete hash");
    settings.minimumOuterWall=.2;
    const auto rejected=validateSourceCageMaterialRegions(binary,settings);
    require(rejected.facts().outcome==SourceCageMaterialRegionsOutcome::Rejected&&!rejected.validatedRegions()&&rejected.facts().outerClearance.atLeastMinimum==false,"too-large outer wall retains measured rejection");
    settings.minimumOuterWall=0.;
    require(validateSourceCageMaterialRegions(binary,settings).facts().outerClearance.contact==false,"zero threshold still requires no contact");
    const auto emptyBinary=validateSourceCageBinaryMaterial(endpoint(f,state(f,0.)));
    const auto empty=validateSourceCageMaterialRegions(emptyBinary,settings);
    require(empty.validatedRegions()&&empty.facts().selectedCount==0&&empty.facts().selected.empty(),"real empty selected cavity admitted");
    require(empty.facts().selectedPartition.partExplicitEmpty==std::vector<bool>{true},"one explicit empty algebra witness");
    const auto& clearance=empty.facts().outerClearance;
    require(clearance.completed&&clearance.state=="NOT_APPLICABLE"&&clearance.applicable==false&&!clearance.contact&&!clearance.squaredDistanceExact&&!clearance.atLeastMinimum,"empty clearance retains actual absent measurements");
    require(empty.facts().carved.region.materialComponents==0&&empty.facts().material.region.boundedVoidComponents==0,"empty native components");
    const auto lifetime=[&]{
        const auto localBinary=validateSourceCageBinaryMaterial(endpoint(f,state(f,0.)));
        const auto localRegions=validateSourceCageMaterialRegions(localBinary);
        return localRegions.validatedRegions();
    }(); // Both originating audits are destroyed here.
    require(lifetime&&identical(lifetime->binaryAudit().facts().worldIntegratedMoments,emptyBinary.facts().worldIntegratedMoments)&&
        identical(lifetime->binaryAudit().facts().spinIntegratedMoments,emptyBinary.facts().spinIntegratedMoments),
        "original world and spin moments survive both originating audits");
    const auto wide=fixture(4,true);auto splitState=state(wide,0.);
    const auto& wideLeaves=wide.source().boundary().binding().domain().leaves();
    for(std::size_t i=0;i<wideLeaves.size();++i)if(wideLeaves[i].coordinates[0]==6||wideLeaves[i].coordinates[0]==9)splitState[i]=1.;
    const auto splitBinary=validateSourceCageBinaryMaterial(endpoint(wide,splitState));
    require(bool(splitBinary.validatedMaterial()),"two separated actual binary cavities validate upstream");
    const auto split=validateSourceCageMaterialRegions(splitBinary);
    require(split.validatedRegions()&&split.facts().carved.region.materialComponents==2&&split.facts().material.region.boundedVoidComponents==2&&split.facts().material.region.materialComponents==1,"two actual disconnected cavities have native component counts");
    // Bounded probe: retain exact rational distance from actual identity-handle
    // BBW evaluation; mathematical affine intent alone does not guarantee that
    // its binary64 output has exactly planar faces or representable distance.
    auto identityState=f.initialState(Eigen::VectorXd::Ones(f.betaCount()),.23);
    const auto identityBinary=validateSourceCageBinaryMaterial(endpoint(f,identityState));
    const auto identityRegions=validateSourceCageMaterialRegions(identityBinary);
    require(bool(identityRegions.validatedRegions()),"actual identity-handle source also passes exact chain");
    std::cout<<"actual_identity_squared_clearance="<<*identityRegions.facts().outerClearance.squaredDistanceExact<<std::endl;
    // Exact actual output distance: (2^48-1)/2^51, independently factored
    // numerator (2^48-1)^2 and denominator 2^102; no tolerance or sqrt rounding.
    require(identityRegions.facts().outerClearance.squaredDistanceExact==std::string("79228162514263774643590529025/5070602400912917605986812821504"),"literal exact actual identity-handle squared distance");
    settings=SourceCageMaterialRegionsSettings{};settings.minimumOuterWall=std::ldexp(281474976710655.,-51);
    const auto equality=validateSourceCageMaterialRegions(identityBinary,settings);
    require(equality.validatedRegions()&&equality.facts().outerClearance.atLeastMinimum==true,"exactly measured representable wall threshold equality passes");
    settings.minimumOuterWall=std::nextafter(settings.minimumOuterWall,std::numeric_limits<double>::infinity());
    const auto above=validateSourceCageMaterialRegions(identityBinary,settings);
    require(!above.validatedRegions()&&above.facts().outcome==SourceCageMaterialRegionsOutcome::Rejected&&above.facts().outerClearance.atLeastMinimum==false,"nextafter above actual exact distance rejects");
    auto mixed=state(f,.4);
    const auto& leaves=f.source().boundary().binding().domain().leaves();
    for(std::size_t i=0;i<leaves.size();++i)if(leaves[i].coordinates[0]==8)mixed[i]=.6;
    mixed[0]=.5;
    const auto failed=validateSourceCageBinaryMaterial(endpoint(f,mixed));
    require(failed.candidate()&&!failed.validatedMaterial(),"actual binary rounding failure has candidate");
    const auto denied=validateSourceCageMaterialRegions(failed);
    require(denied.facts().outcome==SourceCageMaterialRegionsOutcome::NotAdmitted&&!denied.validatedRegions(),"candidate cannot admit");
    auto copied=denied.facts();copied.outcome=SourceCageMaterialRegionsOutcome::Validated;
    require(!denied.validatedRegions(),"copied facts cannot grant");
    for(double t:{-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
        settings.minimumOuterWall=t;
        rejects([&]{(void)validateSourceCageMaterialRegions(binary,settings);},"invalid thresholds reject");
    }
    for(int field=0;field<6;++field)for(std::size_t cap:{std::size_t(0),std::size_t(1),SIZE_MAX}) {
        settings=SourceCageMaterialRegionsSettings{};
        std::size_t* caps[]={&settings.maximumPrimitives,&settings.maximumPoints,&settings.maximumTriangles,&settings.maximumOwnedBytes,&settings.maximumExactDiagnosticBytes,&settings.maximumNefComplexity};
        *caps[field]=cap;
        const auto limited=validateSourceCageMaterialRegions(binary,settings);
        require(limited.facts().outcome==SourceCageMaterialRegionsOutcome::ResourceLimit&&!limited.validatedRegions(),"caps and overflow never publish");
    }
    const auto retained=audit.validatedRegions();const auto hash=a.scientificSha256;
    const auto moved=std::move(audit);
    require(moved.facts().scientificSha256==hash&&retained->binaryMaterial()->binaryState()[0]==1.,"immutable product lifetime across move");
    rejects([&]{(void)audit.facts();},"moved-from evidence rejects");
    // Generic backend negatives: these are NOT corruptible immutable upstream
    // admission fixtures and confer no source-binary-product authority.
    const auto unit=cube(0.,1.),left=cube(0.,.75),right=cube(.25,1.);
    const auto whole=verification::makeExactMaterialRegion(unit.v,unit.f);
    const auto l=verification::makeExactMaterialRegion(left.v,left.f);
    const auto r=verification::makeExactMaterialRegion(right.v,right.f);
    const auto overlap=verification::checkExactMaterialPartition(whole,{l,r});
    require(overlap.completed&&overlap.noPositiveVolumeOverlap==false&&overlap.completePartition==false,"raw overlapping generic parts reject");
    const auto missing=verification::checkExactMaterialPartition(whole,{l});
    require(missing.completed&&missing.missingWholeContent==true&&missing.completePartition==false,"raw missing generic content rejects");
    const auto badNesting=join(unit,cube(.25,.75));
    const auto nested=verification::makeExactMaterialRegion(badNesting.v,badNesting.f);
    require(nested.facts().completed&&!nested.facts().valid,"same-positive nested generic shells reject without reorientation");
    std::cout<<"KNOWN_SOLUTION_VALIDATION exact regions PASS; represented design and outer wall only; no seal, full thickness, feasibility, physics or paper grant\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
#endif
