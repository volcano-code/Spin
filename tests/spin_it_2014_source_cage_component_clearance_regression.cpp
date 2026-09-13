#include <iostream>
// Full actual geometry tests become active once the missing-layer RED is resolved.
#if !__has_include("spin2014/SourceCageComponentClearance.h")
int main(){std::cerr<<"FAIL: actual final boundary-component separation layer absent\n";return 1;}
#else
#include "spin2014/SourceCageComponentClearance.h"
#include "source_surface_test_utils.h"
#include <Eigen/Eigenvalues>
#include <cstring>
#include <limits>
#include <type_traits>
using namespace spin2014;
using namespace surface_test;
namespace {
SourceCagePaperFunctional fixture(bool wide=false) {
    const auto outer=cube(.3125,.6875);
    const auto design=wide?box({.4,.45625,.45625},{.56,.54375,.54375}):cube(.425,.575);
    const int level=wide?5:4;
    SourceSurfaceOctreeConfig sc;sc.maximumDepth=level;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,sc);
    SourceAdaptiveDomain domain(source,level,level,1000);
    SourceSurfaceCageConfig cc;cc.cageLevel=2;SourceCageGridConfig gc;gc.fineLevel=3;
    SourceCageBbw bbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),{}),{});
    require(bbw.success(),("actual full BBW succeeds: "+bbw.terminalReason()).c_str());
    SourceConformingBoundaryConfig bc;bc.maximumPatches=9000;
    SourceCageMomentFunctionalConfig mc;mc.maximumWeightValues=9000;
    SourceCagePaperFunctionalConfig pc;pc.objective=Spin2014Objective::TopMassAndLowering;pc.inertiaWeight=0.;pc.centerWeight=0.;
    return SourceCagePaperFunctional(SourceCageMomentFunctional(SourceConformingBoundary(SourceCageBinding(domain,bbw,{}),bc),1.17,mc),pc);
}
Eigen::VectorXd state(const SourceCagePaperFunctional& f,double beta) {
    auto x=f.initialState(Eigen::VectorXd::Constant(f.betaCount(),beta),.23);
    for(std::size_t h=0;h<f.handleCount();++h){x[f.handleOffset()+4*h]=1.02;x.segment<3>(f.handleOffset()+4*h+1)<<.011,-.007,.019;}
    return x;
}
// KNOWN_SOLUTION_VALIDATION: actual full handles, nonzero phi, stationary zero objective.
SourceCageBinaryMaterialAudit binary(const SourceCagePaperFunctional& f,const Eigen::VectorXd& x) {
    const auto m=f.evaluateValues(x).world.finalMoments.values();SourceCagePaperFrame frame;frame.origin=m.segment<3>(1)/m[0];
    Eigen::Matrix3d q;q<<m[4],m[9],m[8],m[9],m[5],m[7],m[8],m[7],m[6];
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(q-m[0]*frame.origin*frame.origin.transpose());
    Eigen::Matrix3d base=eig.eigenvectors().transpose();if(base.determinant()<0)base.row(0)*=-1;
    const double c=std::cos(x[f.angleOffset()]),s=std::sin(x[f.angleOffset()]);Eigen::Matrix3d g;g<<c,s,0,-s,c,0,0,0,1;frame.rotation=g*base;
    SourceCageAdaptiveSolverSettings settings;settings.maximumPasses=3;settings.solve.source.solver.maximumIterations=1;
    settings.solve.source.solver.maximumLineSearchSteps=2;settings.solve.budget.maximumGuardCalls=8;
    const auto ep=solveSourceCageAdaptive(SourceCagePaperFunctional(f.source(),f.config(),frame),x,settings);
    require(ep.facts().outcome==SourceCageAdaptiveSolverOutcome::RelaxedAdaptiveComplete,"real adaptive fixture completed");
    return validateSourceCageBinaryMaterial(ep);
}
void provenance(const SourceCageComponentClearanceAudit& a,const SourceCageMaterialRegionsAudit& upstream,std::size_t n) {
    const auto& f=a.facts();const auto& original=upstream.validatedRegions()->binaryMaterial()->materialBoundary();
    require(f.components.size()==n&&f.componentCount==n&&f.identitiesChecked&&f.coverageChecked,"complete actual extraction");
    require(f.materialInputSha256==upstream.facts().material.region.inputDigest,"actual mesh hash matches upstream");
    std::vector<int> vertices(original.vertices.rows(),0),faces(original.triangles.rows(),0);
    for(std::size_t i=0;i<n;++i){const auto& c=f.components[i];
        require(c.componentId==i&&c.originalSign==(i? -1:1)&&c.checkerReversed==(i!=0),"canonical component sign and reversal");
        require(c.parent==upstream.facts().material.region.parents[i]&&c.depth==upstream.facts().material.region.depths[i],"original nesting retained");
        require(std::is_sorted(c.sourceVertexIds.begin(),c.sourceVertexIds.end())&&std::is_sorted(c.sourceFaceIds.begin(),c.sourceFaceIds.end()),"canonical maps");
        for(auto v:c.sourceVertexIds)++vertices.at(v);
        for(auto face:c.sourceFaceIds){++faces.at(face);require(upstream.facts().material.boundary->faceComponents[face]==int(i),"map follows exact upstream component");}
        require(c.region.region.valid&&c.region.region.materialComponents==1&&c.region.region.boundedVoidComponents==0&&c.region.boundary->closed,"standalone positive whole component");
        require(c.region.region.inputDigest==c.checkerInputSha256&&c.region.boundary->inputDigest==c.checkerInputSha256,"checker exact identity");
    }
    for(auto k:vertices)require(k==1,"every original vertex exactly once");for(auto k:faces)require(k==1,"every original face exactly once");
    require(f.plannedPairCount==n*(n-1)/2&&f.pairs.size()==f.plannedPairCount,"all unordered pairs retained");
    std::size_t k=0;for(std::size_t i=0;i<n;++i)for(std::size_t j=i+1;j<n;++j){const auto& p=f.pairs[k++];
        require(p.componentIds==std::array<std::size_t,2>{i,j}&&p.clearance.completed&&p.clearance.squaredDistanceExact&&p.clearance.witnessFaces&&p.parentWitnessFaces,"complete pair measurement");
        require((*p.parentWitnessFaces)[0]==f.components[i].sourceFaceIds.at((*p.clearance.witnessFaces)[0])&&(*p.parentWitnessFaces)[1]==f.components[j].sourceFaceIds.at((*p.clearance.witnessFaces)[1]),"local to parent witness");
    }
    require(f.fullManufacturingThickness=="NOT_ESTABLISHED"&&f.scope=="BOUNDARY_COMPONENT_SEPARATION","bounded scope only");
    require(f.nefBudget.currentComplexity==0&&f.nefBudget.peakComplexity>0,"only released local Nef snapshot escapes");
}
}
int main(){try{
    static_assert(!std::is_default_constructible_v<SourceCageComponentClearanceAudit>);
    static_assert(!std::is_default_constructible_v<SourceCageComponentClearance>);
    static_assert(!std::is_constructible_v<SourceCageComponentClearanceAudit,SourceCageComponentClearanceFacts>);
    const auto f=fixture();const auto b=binary(f,state(f,1.));const auto r=validateSourceCageMaterialRegions(b);
    require(bool(r.validatedRegions()),"one-cavity upstream validated");
    SourceCageComponentClearanceSettings settings;settings.minimumComponentSeparation=.1;
    const auto a=validateSourceCageComponentClearance(r,settings);provenance(a,r,2);
    std::cout<<"actual one-cavity extraction measured"<<std::endl;
    require(a.validatedClearance()&&a.facts().outcome==SourceCageComponentClearanceOutcome::Validated&&a.facts().state=="MEASURED"&&a.facts().atLeastMinimum==true,"actual one cavity passes");
    require(a.facts().minimumSquaredDistanceExact==r.facts().outerClearance.squaredDistanceExact,"actual component gap agrees with outer-carved for one cavity");
    require(a.validatedClearance()->regions().get()==r.validatedRegions().get()&&&a.validatedClearance()->regionsAudit().facts()==&r.facts(),"entire original audit retained");
    require(a.facts().scientificSha256==validateSourceCageComponentClearance(r,settings).facts().scientificSha256,"deterministic repeat");
    settings.minimumComponentSeparation=.2;const auto rejected=validateSourceCageComponentClearance(r,settings);
    require(!rejected.validatedClearance()&&rejected.facts().outcome==SourceCageComponentClearanceOutcome::Rejected&&rejected.facts().minimumSquaredDistanceExact&&rejected.facts().atLeastMinimum==false,"measured rejection retained");
    const auto empty=validateSourceCageMaterialRegions(binary(f,state(f,0.)));const auto emptyA=validateSourceCageComponentClearance(empty);provenance(emptyA,empty,1);
    require(emptyA.validatedClearance()&&emptyA.facts().state=="NOT_APPLICABLE"&&emptyA.facts().applicable==false&&!emptyA.facts().minimumSquaredDistanceExact&&!emptyA.facts().minimumPairIndex&&!emptyA.facts().minimumParentWitnessFaces&&!emptyA.facts().atLeastMinimum&&!emptyA.facts().contact,"one component is explicit N/A with no invented measurement");
    std::cout<<"actual empty-cavity N/A checked"<<std::endl;
    const auto wide=fixture(true);auto x=state(wide,0.);const auto& leaves=wide.source().boundary().binding().domain().leaves();
    for(std::size_t i=0;i<leaves.size();++i)if(leaves[i].coordinates[0]==14||leaves[i].coordinates[0]==16)x[i]=1.;
    const auto wr=validateSourceCageMaterialRegions(binary(wide,x));require(bool(wr.validatedRegions()),"two-cavity wider actual fixture validated");
    settings.minimumComponentSeparation=.04;const auto split=validateSourceCageComponentClearance(wr,settings);provenance(split,wr,3);
    require(split.facts().pairs[0].clearance.atLeastMinimum==true&&split.facts().pairs[1].clearance.atLeastMinimum==true&&split.facts().pairs[2].clearance.atLeastMinimum==false,"outer gaps pass, thin actual inter-cavity gap rejects");
    require(split.facts().outcome==SourceCageComponentClearanceOutcome::Rejected&&!split.validatedClearance()&&split.facts().minimumPairIndex==2,"all-pair global minimum catches thin web");
    std::cout<<"actual thin inter-cavity rejection checked; squared_gap="<<*split.facts().minimumSquaredDistanceExact<<std::endl;
    settings.minimumComponentSeparation=.3;const auto early=validateSourceCageComponentClearance(wr,settings);
    require(early.facts().pairs.size()==3&&early.facts().pairs[0].clearance.atLeastMinimum==false&&early.facts().pairs[2].clearance.completed&&early.facts().minimumPairIndex==2,"first-pair rejection still measures every later pair");
    settings.minimumComponentSeparation=0.;const auto zero=validateSourceCageComponentClearance(wr,settings);require(zero.validatedClearance()&&zero.facts().contact==false,"zero threshold still measures all distinct gaps");
    const auto ir=validateSourceCageMaterialRegions(binary(f,f.initialState(Eigen::VectorXd::Ones(f.betaCount()),.23)));
    settings.minimumComponentSeparation=std::ldexp(281474976710655.,-51);const auto eq=validateSourceCageComponentClearance(ir,settings);
    require(eq.validatedClearance()&&eq.facts().minimumSquaredDistanceExact==std::string("79228162514263774643590529025/5070602400912917605986812821504"),"literal exact identity threshold equality passes");
    settings.minimumComponentSeparation=std::nextafter(settings.minimumComponentSeparation,std::numeric_limits<double>::infinity());
    require(validateSourceCageComponentClearance(ir,settings).facts().outcome==SourceCageComponentClearanceOutcome::Rejected,"one ULP above literal exact gap rejects");
    std::cout<<"literal exact equality/nextafter checked"<<std::endl;
    SourceCageMaterialRegionsSettings rs;rs.minimumOuterWall=.2;const auto denied=validateSourceCageComponentClearance(validateSourceCageMaterialRegions(b,rs));
    require(denied.facts().outcome==SourceCageComponentClearanceOutcome::NotAdmitted&&!denied.validatedClearance(),"rejected upstream cannot admit");auto copied=denied.facts();copied.outcome=SourceCageComponentClearanceOutcome::Validated;require(!denied.validatedClearance(),"copied facts cannot grant");
    for(double t:{-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}){settings.minimumComponentSeparation=t;rejects([&]{(void)validateSourceCageComponentClearance(r,settings);},"invalid threshold rejects");}
    for(int field=0;field<7;++field)for(std::size_t cap:{std::size_t(0),std::size_t(1),SIZE_MAX}){settings={};std::size_t* caps[]={&settings.maximumComponents,&settings.maximumPoints,&settings.maximumTriangles,&settings.maximumPairs,&settings.maximumOwnedBytes,&settings.maximumExactDiagnosticBytes,&settings.maximumNefComplexity};*caps[field]=cap;
        const auto limit=validateSourceCageComponentClearance(wr,settings);require(limit.facts().outcome==SourceCageComponentClearanceOutcome::ResourceLimit&&!limit.validatedClearance()&&limit.facts().state!="NOT_APPLICABLE","budgets/overflow cannot masquerade as N/A or success");}
    const auto lifetime=[&]{const auto local=validateSourceCageMaterialRegions(binary(f,state(f,1.)));return validateSourceCageComponentClearance(local).validatedClearance();}();
    const auto& actual=lifetime->regions()->binaryAudit().facts().worldIntegratedMoments;const auto& expected=b.facts().worldIntegratedMoments;
    require(actual.size()==expected.size()&&std::memcmp(actual.data(),expected.data(),actual.size()*sizeof(double))==0,"original moment bytes survive originating audits");
    const auto& spinActual=lifetime->regions()->binaryAudit().facts().spinIntegratedMoments;const auto& spinExpected=b.facts().spinIntegratedMoments;
    require(spinActual.size()==spinExpected.size()&&std::memcmp(spinActual.data(),spinExpected.data(),spinActual.size()*sizeof(double))==0,"original spin moment bytes survive originating audits");
    const auto& retainedMesh=lifetime->regions()->binaryMaterial()->materialBoundary();const auto& originalMesh=b.validatedMaterial()->materialBoundary();
    require(retainedMesh.vertices.rows()==originalMesh.vertices.rows()&&retainedMesh.triangles.rows()==originalMesh.triangles.rows()&&
        std::memcmp(retainedMesh.vertices.data(),originalMesh.vertices.data(),originalMesh.vertices.size()*sizeof(double))==0&&
        std::memcmp(retainedMesh.triangles.data(),originalMesh.triangles.data(),originalMesh.triangles.size()*sizeof(int))==0,"original actual material mesh bytes retained unchanged");
    auto movable=validateSourceCageComponentClearance(r);const auto kept=movable.validatedClearance();const auto moved=std::move(movable);require(kept&&moved.validatedClearance(),"immutable product move lifetime");rejects([&]{(void)movable.facts();},"moved-from missing authority rejects");
    std::cout<<"KNOWN_SOLUTION_VALIDATION actual component clearance PASS; no full thickness or seal\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
#endif
