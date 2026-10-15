// Missing-operation branch is preserved as the initial TDD RED evidence.
#if __has_include("spin2014/SourceMaterialBoundary.h")
#include "spin2014/SourceMaterialBoundary.h"
#include "spin2014/AdaptiveFillFunctional.h"
#include "spin2014/SlqpAuditTrace.h"
#include "source_surface_test_utils.h"
#include <cstring>
#include <iomanip>
#include <iostream>
#include <type_traits>
#if __has_include("verification/ExactNefBudget.h") // Recompile after the preserved missing-API RED.
#include "verification/ExactNefBudget.h"
#define SOURCE_MATERIAL_HAS_SCOPED_BUDGET 1
#endif
using namespace spin2014;
using namespace verification;
using namespace surface_test;

std::shared_ptr<const SourceSurfaceOctree> source(const Mesh& outer,const Mesh& design,
    int depth=4,Eigen::Vector3d origin={0,0,0},double side=1) {
    SourceSurfaceOctreeConfig cfg; cfg.maximumDepth=depth; cfg.rootOrigin=origin; cfg.rootSide=side;
    return std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,cfg);
}
ExactCuboid cuboid(Eigen::Vector3d lo,Eigen::Vector3d hi) {return {lo,hi};}
Ten integrated(const ExactCuboidUnion& u) {
    require(u.facts().status=="EXACT_UNION_VALIDATED","union must produce accepted mesh");
    require(u.facts().roundtripEqual==true,"exact region roundtrip equality");
    return physics::RawMassMoments10::fromMassMoments(computeSurfaceMassMoments(u.vertices(),u.faces(),1)).values();
}
void unionTests() {
    // Catch concatenation/shared faces, double-counted overlap, omitted shells,
    // filled facet holes, lossy coordinates and acceptance after a resource cap.
    const auto a=cuboid({0,0,0},{1,1,1});
    auto adjacent=makeExactCuboidUnion({a,cuboid({1,0,0},{2,1,1})});
    close(integrated(adjacent),boxMoments({0,0,0},{2,1,1}),"adjacent exact union box");
    require(adjacent.region()->facts().materialComponents==1,"adjacent shared faces cancelled");
    auto overlap=makeExactCuboidUnion({a,cuboid({.5,0,0},{1.5,1,1})});
    close(integrated(overlap),boxMoments({0,0,0},{1.5,1,1}),"overlap counted exactly once");
    auto disconnected=makeExactCuboidUnion({a,cuboid({2,0,0},{3,1,1})});
    close(integrated(disconnected),boxMoments({0,0,0},{1,1,1})+boxMoments({2,0,0},{3,1,1}),"disconnected union");
    require(disconnected.region()->facts().materialComponents==2,"retain both components");
    std::vector<ExactCuboid> t={a};
    for(int y=0;y<2;++y) for(int z=0;z<2;++z)
        t.push_back(cuboid({1,y*.5,z*.5},{2,(y+1)*.5,(z+1)*.5}));
    close(integrated(makeExactCuboidUnion(t)),boxMoments({0,0,0},{2,1,1}),"coarse fine T junction conforming");
    for(auto b:{cuboid({1,1,0},{2,2,1}),cuboid({1,1,1},{2,2,2})}) {
        auto bad=makeExactCuboidUnion({a,b});
        require(bad.facts().status=="NONMANIFOLD_GEOMETRY" && bad.vertices().rows()==0 && !bad.region(),"edge/point contact rejected without partial mesh");
    }
    std::vector<ExactCuboid> hollow;
    for(int x=0;x<3;++x) for(int y=0;y<3;++y) for(int z=0;z<3;++z)
        if(x!=1||y!=1||z!=1) hollow.push_back(cuboid({double(x),double(y),double(z)},{double(x+1),double(y+1),double(z+1)}));
    auto shell=makeExactCuboidUnion(hollow); Ten expected;
    expected<<26,39,39,39,236./3,236./3,236./3,117./2,117./2,117./2;
    close(integrated(shell),expected,"hollow26 literal all ten");
    require(shell.region()->facts().boundedVoidComponents==1,"hollow26 inner void retained");
    const auto inertia=computeSurfaceMassMoments(shell.vertices(),shell.faces(),1).inertiaAboutCenter();
    require((inertia.diagonal().array()-121./3).abs().maxCoeff()<1e-12,"hollow26 central inertia");
    // Through-hole creates polygonal facets with holes, not just nested shells.
    std::vector<ExactCuboid> ring;
    for(int x=0;x<3;++x) for(int y=0;y<3;++y) if(x!=1||y!=1)
        ring.push_back(cuboid({double(x),double(y),0},{double(x+1),double(y+1),1}));
    close(integrated(makeExactCuboidUnion(ring)),boxMoments({0,0,0},{3,3,1})-boxMoments({1,1,0},{2,2,1}),"facet holes remain open");
    auto repeated=makeExactCuboidUnion(hollow);
    require(shell.facts().outputDigest==repeated.facts().outputDigest && shell.vertices()==repeated.vertices() && shell.faces()==repeated.faces(),"deterministic canonical mesh");
    std::reverse(hollow.begin(),hollow.end());
    require(makeExactCuboidUnion(hollow).facts().inputDigest!=shell.facts().inputDigest,"ordered input identity");
    for(auto b:{cuboid({0,0,0},{0,1,1}),cuboid({0,0,0},{-1,1,1}),cuboid({0,0,0},{INFINITY,1,1})})
        rejects([&]{makeExactCuboidUnion({b});},"malformed cuboid rejected");
    ExactCuboidUnionSettings caps;
    for(int field=0;field<4;++field) {
        auto invalid=caps;
        std::size_t* p[]={&invalid.maximumInputCuboids,&invalid.maximumOutputVertices,&invalid.maximumOutputTriangles,&invalid.maximumIntermediateComplexity};
        *p[field]=0; rejects([&]{makeExactCuboidUnion({},invalid);},"invalid budget rejected even for empty");
        auto limited=caps;
        std::size_t* q[]={&limited.maximumInputCuboids,&limited.maximumOutputVertices,&limited.maximumOutputTriangles,&limited.maximumIntermediateComplexity};
        *q[field]=1;auto failed=makeExactCuboidUnion({a,cuboid({1,0,0},{2,1,1})},limited);
        require(failed.facts().status=="RESOURCE_LIMIT_UNKNOWN" && failed.vertices().rows()==0 && !failed.region(),"cap fails atomically");
    }
    require(shell.facts().booleanOperations==25 && shell.facts().peakComplexity>0 && shell.facts().finalComplexity.vertices>0,"actual resource counters");
    auto empty=makeExactCuboidUnion({});
    require(empty.facts().status=="EXACT_UNION_VALIDATED" && empty.region()->facts().explicitEmpty && empty.vertices().rows()==0,"explicit empty algebra");
}
SourceAdaptiveSpinRequest goodRequest() {
    SourceAdaptiveSpinRequest r;
    r.source=source(box({-1,-1,0},{1,1,2}),box({-.5,-.5,1},{.5,.5,1.5}),4,{-1,-1,0},2);
    r.initialLevel=2;r.maximumOptimizationLevel=4;r.initialFill=.5;r.initialPhi=0;
    r.laplacianWeight=1;r.centerWeight=1;r.density=1;
    r.slqp.maximumIterations=300;r.slqp.equalityTolerance=1e-10;r.slqp.projectedGradientTolerance=1e-8;
    return r;
}
void sourceDualDensityFunctionalTests() {
    auto r=goodRequest();
    SourceAdaptiveDomain domain(r.source,r.initialLevel,r.maximumOptimizationLevel);
    Eigen::VectorXd beta=Eigen::VectorXd::Constant(domain.leaves().size(),.375);
    const double shellDensity=1.25,interiorDensity=6.5;
    AdaptiveFillFunctional dual(domain,
        SpinIt2014MomentProblem::fromObjective(r.objective,r.inertiaWeight,r.centerWeight),
        r.laplacianWeight,shellDensity,interiorDensity);
    const auto expected=(r.source->outerSurfaceMoments()-r.source->designSurfaceMoments()).scaledDensity(shellDensity)
        +r.source->designSurfaceMoments().scaledDensity(interiorDensity)
        -domain.removedMoments(beta).scaledDensity(interiorDensity);
    const auto actual=dual.evaluatePaper(beta,0).finalMoments;
    close(actual.values(),expected.values(),"source dual density applies shell/interior separately");

    AdaptiveFillFunctional legacy(domain,
        SpinIt2014MomentProblem::fromObjective(r.objective,r.inertiaWeight,r.centerWeight),
        r.laplacianWeight,r.density);
    AdaptiveFillFunctional equal(domain,
        SpinIt2014MomentProblem::fromObjective(r.objective,r.inertiaWeight,r.centerWeight),
        r.laplacianWeight,r.density,r.density);
    close(legacy.evaluatePaper(beta,0).finalMoments.values(),equal.evaluatePaper(beta,0).finalMoments.values(),
        "source dual equal densities preserve legacy moment semantics");

    auto request=r;request.boundaryDensity=shellDensity;request.interiorDensity=interiorDensity;
    validateSourceAdaptiveSpinRequest(request);
    require(request.resolvedBoundaryDensity()==shellDensity && request.resolvedInteriorDensity()==interiorDensity,
        "source request resolves explicit dual densities");
    const auto dualHash=computeSourceAdaptiveSpinRequestHash(request);
    request.boundaryDensity.reset();request.interiorDensity.reset();
    require(computeSourceAdaptiveSpinRequestHash(request)!=dualHash,"dual densities bind source request identity");
    request.boundaryDensity=shellDensity;
    rejects([&]{validateSourceAdaptiveSpinRequest(request);},"partial dual-density request rejected");

    auto product=makeSourceMaterialBoundary(domain,Eigen::VectorXd::Ones(domain.leaves().size()),shellDensity,interiorDensity);
    require(bool(product.product()),product.facts().diagnostic.c_str());
    require(product.facts().materialComponentsExported && product.facts().dualMaterialPartition.completePartition==true,
        "shell/interior/carve geometry is an exact partition of original solid");
    require(product.product()->boundaryMaterialFaces().rows()>0 && product.product()->interiorMaterialFaces().rows()>0,
        "both source surface material boundaries are exported");
    const auto expectedProduct=(r.source->outerSurfaceMoments()-r.source->designSurfaceMoments()).scaledDensity(shellDensity)
        +r.source->designSurfaceMoments().scaledDensity(interiorDensity)
        -domain.removedMoments(Eigen::VectorXd::Ones(domain.leaves().size())).scaledDensity(interiorDensity);
    close(product.facts().integratedMoments.values(),expectedProduct.values(),
        "two exported source material meshes independently integrate to weighted total");
    const auto singleDensity=makeSourceMaterialBoundary(domain,Eigen::VectorXd::Ones(domain.leaves().size()),1.);
    const auto equalDensity=makeSourceMaterialBoundary(domain,Eigen::VectorXd::Ones(domain.leaves().size()),1.,1.);
    require(bool(singleDensity.product()) && bool(equalDensity.product()) &&
        equalDensity.facts().materialComponentsExported,
        "equal-density dual material still exports both component surfaces");
    close(singleDensity.facts().integratedMoments.values(),equalDensity.facts().integratedMoments.values(),
        "equal-density dual material preserves single-material total moments");
    const auto noCarve=makeSourceMaterialBoundary(domain,Eigen::VectorXd::Zero(domain.leaves().size()),shellDensity,interiorDensity);
    require(bool(noCarve.product()) && noCarve.facts().dualMaterialPartition.completePartition==true &&
        noCarve.product()->interiorMaterialFaces().rows()==r.source->designFaces().rows(),
        "empty carve retains the full independently exported interior material");
    rejects([&]{makeSourceMaterialBoundary(domain,Eigen::VectorXd::Zero(domain.leaves().size()),0.,interiorDensity);},
        "invalid boundary density is rejected before geometry work");

    auto solveRequest=goodRequest();solveRequest.boundaryDensity=shellDensity;solveRequest.interiorDensity=interiorDensity;
    SourceMaterialBoundarySettings settings;settings.minimumOuterWall=.5;
    const auto solved=solveSourceAdaptiveSpinMaterial(solveRequest,settings);
    require(solved.solve().status==SourceAdaptiveSpinStatus::BinaryMomentValidated &&
        solved.boundary().product() && solved.boundary().facts().materialComponentsExported,
        "source dual-density request reaches optimizer, binary rounding and independent material export");
    require(solved.boundary().facts().normalizedMomentError<=1e-10 &&
        solved.solve().binaryCheck && solved.solve().binaryCheck->passed,
        "source dual-density final geometry and constraints independently pass");
}
void exactOutputCapTests() {
    // Catch using a loose growth estimate as if it were actual output size.
    ExactCuboidUnionSettings caps;caps.maximumOutputVertices=8;caps.maximumOutputTriangles=12;
    auto cube=makeExactCuboidUnion({cuboid({0,0,0},{1,1,1})},caps);
    require(bool(cube.region()),"exact 8-vertex/12-triangle cube budget must be accepted");
    require(cube.vertices().rows()==8 && cube.faces().rows()==12,"actual cube output counts");
    caps.maximumOutputVertices=7;
    require(!makeExactCuboidUnion({cuboid({0,0,0},{1,1,1})},caps).region(),"one fewer output vertex fails");
    caps.maximumOutputVertices=8;
    caps.maximumOutputTriangles=11;
    require(!makeExactCuboidUnion({cuboid({0,0,0},{1,1,1})},caps).region(),"one fewer output triangle fails");
}
void checkOuter(const SourceMaterialBoundary& p,const Mesh& outer) {
    require(p.outerVertices().rows()==outer.v.rows() && p.outerFaces().rows()==outer.f.rows(),"original array dimensions");
    require(std::memcmp(p.outerVertices().data(),outer.v.data(),outer.v.size()*sizeof(double))==0 &&
        std::memcmp(p.outerFaces().data(),outer.f.data(),outer.f.size()*sizeof(int))==0,"original outer byte identical");
    require(p.materialVertices().topRows(outer.v.rows())==outer.v && p.materialFaces().topRows(outer.f.rows())==outer.f,"untouched outer prefix");
}
void materialTests() {
    auto r=goodRequest(); SourceAdaptiveDomain domain(r.source,2,4);
    Eigen::VectorXd beta=Eigen::VectorXd::Ones(domain.leaves().size());
    SourceMaterialBoundarySettings cfg;cfg.minimumOuterWall=.5;
    auto result=makeSourceMaterialBoundary(domain,beta,1,cfg);
    require(bool(result.product()),result.facts().diagnostic.c_str());
    const auto& p=*result.product();const auto& f=p.facts();
    std::cout<<std::setprecision(17)<<"material product="<<f.productDigest<<" mesh="<<f.materialMeshHash
        <<" union="<<f.unionFacts.outputDigest<<" cuboids="<<f.unionFacts.inputCuboids
        <<" operations="<<f.unionFacts.booleanOperations<<" peak_complexity="<<f.unionFacts.peakComplexity
        <<" material_budget_peak="<<f.resourceBudget.peakComplexity<<" material_budget_live_snapshot="<<f.resourceBudget.currentComplexity
        <<" final_complexity="<<f.unionFacts.finalComplexity.total()<<" residual="<<f.normalizedMomentError
        <<" min_inertia="<<f.minimumCentralInertiaEigenvalue<<" moments="<<f.integratedMoments.values().transpose()<<'\n';
    checkOuter(p,box({-1,-1,0},{1,1,2}));
    close(f.integratedMoments.values(),boxMoments({-1,-1,0},{1,1,2})-boxMoments({-3./8,-3./8,9./8},{3./8,3./8,11./8}),"independent final mesh all ten");
    require(std::abs(f.integratedMoments.mass()-503./64)<1e-12 && std::abs(f.integratedMoments.values()(3)-2003./256)<1e-12,"literal final mass and firstz");
    require(f.partition.completePartition==true && f.carvedSubsetDesign.aSubsetB==true,"independent exact partition and design subset");
    require(f.clearance.squaredDistanceExact=="25/64" && f.clearance.contact==false && f.clearance.witnessFaces.has_value(),"actual 5/8 clearance witness");
    require(f.materialRegion.materialComponents==1 && f.materialRegion.boundedVoidComponents==1 && f.carvedRegion.materialComponents==1,"one bounded new cavity");
    cfg.minimumOuterWall=5./8;auto equal=makeSourceMaterialBoundary(domain,beta,1,cfg);require(bool(equal.product()),"clearance threshold equality accepted");
    require(equal.facts().productDigest!=result.facts().productDigest,"wall request binds product identity");
    cfg.minimumOuterWall=std::nextafter(5./8,INFINITY);auto above=makeSourceMaterialBoundary(domain,beta,1,cfg);
    require(!above.product() && above.facts().clearance.atLeastMinimum==false && above.facts().clearance.squaredDistanceExact=="25/64","nextafter threshold fails without changing geometry");
    cfg.minimumOuterWall=.5;
    auto twice=makeSourceMaterialBoundary(domain,beta,2,cfg);require(bool(twice.product()),"density2 accepted");
    close(twice.product()->facts().integratedMoments.values(),2*f.integratedMoments.values(),"density doubles all ten");
    require(twice.product()->facts().clearance.squaredDistanceExact==f.clearance.squaredDistanceExact && twice.facts().productDigest!=result.facts().productDigest,"density binds product but not clearance");
    auto same=makeSourceMaterialBoundary(domain,beta,1,cfg);require(same.facts().productDigest==result.facts().productDigest,"repeat product digest");
    auto capChange=cfg;++capChange.unionSettings.maximumOutputVertices;
    auto changed=makeSourceMaterialBoundary(domain,beta,1,capChange);require(changed.facts().productDigest!=result.facts().productDigest,"caps bind product");
    auto exactCaps=cfg;exactCaps.unionSettings.maximumOutputVertices=16;exactCaps.unionSettings.maximumOutputTriangles=24;
    require(bool(makeSourceMaterialBoundary(domain,beta,1,exactCaps).product()),"combined original-plus-cavity cap equality");
    exactCaps.unionSettings.maximumOutputVertices=15;
    auto vertexCap=makeSourceMaterialBoundary(domain,beta,1,exactCaps);
    require(!vertexCap.product() && vertexCap.facts().status=="RESOURCE_LIMIT_UNKNOWN","combined material vertex cap excludes partial product");
    exactCaps.unionSettings.maximumOutputVertices=16;exactCaps.unionSettings.maximumOutputTriangles=23;
    auto triangleCap=makeSourceMaterialBoundary(domain,beta,1,exactCaps);
    require(!triangleCap.product() && triangleCap.facts().status=="RESOURCE_LIMIT_UNKNOWN","combined material triangle cap excludes partial product");
    beta.setZero();auto zero=makeSourceMaterialBoundary(domain,beta,1,cfg);require(bool(zero.product()),"zero carve geometry accepted");
    require(zero.product()->materialVertices()==r.source->outerVertices() && zero.product()->materialFaces()==r.source->outerFaces() && zero.product()->facts().carvedRegion.explicitEmpty,"zero carve unchanged outer");
    require(zero.product()->facts().clearance.applicable==false && !zero.product()->facts().clearance.squaredDistanceExact,"zero carve no fabricated distance");
    require(zero.facts().productDigest!=result.facts().productDigest && result.product()->facts().integratedMoments.mass()>7,"input mutation cannot alter product");
    beta(0)=.5;rejects([&]{makeSourceMaterialBoundary(domain,beta,1,cfg);},"fractional input rejected");
    beta(0)=NAN;rejects([&]{makeSourceMaterialBoundary(domain,beta,1,cfg);},"nonfinite input rejected");
    rejects([&]{makeSourceMaterialBoundary(domain,Eigen::VectorXd::Zero(beta.size()),0,cfg);},"density validation");
    cfg.minimumOuterWall=-1;rejects([&]{makeSourceMaterialBoundary(domain,Eigen::VectorXd::Zero(beta.size()),1,cfg);},"wall validation");
    // Translation/asymmetry makes xy,xz,yz nonzero; independent box formulas.
    auto shifted=source(box({1,2,3},{3,4,5}),box({1.5,2.5,4},{2.5,3.5,4.5}),4,{1,2,3},2);
    SourceAdaptiveDomain sd(shifted,2,4);auto asym=makeSourceMaterialBoundary(sd,Eigen::VectorXd::Ones(sd.leaves().size()),1);
    require(bool(asym.product()),"asymmetric source accepted");
    close(asym.product()->facts().integratedMoments.values(),boxMoments({1,2,3},{3,4,5})-boxMoments({13./8,21./8,33./8},{19./8,27./8,35./8}),"shifted asymmetric all ten mixed moments");
    // Canonical non-dyadic endpoints survive exact union conversion unchanged.
    auto nonDyadic=source(cube(0,1),cube(.125,.875),4,{-.3,-.2,-.1},1.4);
    SourceAdaptiveDomain nd(nonDyadic,3,4);auto ndOut=makeSourceMaterialBoundary(nd,Eigen::VectorXd::Ones(nd.leaves().size()),1);
    require(bool(ndOut.product()),"non-dyadic source exact union and tetrahedra agree");
    for(bool existingHole:{false,true}) {
        const auto outer=existingHole?join(cube(0,1),cube(.375,.625),true):cube(0,1);
        auto nested=source(outer,join(cube(.125,.875),cube(.375,.625),true));
        SourceAdaptiveDomain d(nested,2,4);auto island=makeSourceMaterialBoundary(d,Eigen::VectorXd::Ones(d.leaves().size()),1);
        require(bool(island.product()),island.facts().diagnostic.c_str());
        const auto& facts=island.product()->facts();
        require(facts.materialRegion.materialComponents==2 && facts.materialRegion.boundedVoidComponents==(existingHole?2:1) && facts.carvedRegion.boundedVoidComponents==1,"retained island nesting and existing void count");
        checkOuter(*island.product(),outer);
    }
    auto ds=source(cube(0,1),join(box({.125,.125,.125},{.375,.375,.375}),box({.625,.625,.625},{.875,.875,.875})));
    SourceAdaptiveDomain dd(ds,2,4);auto dis=makeSourceMaterialBoundary(dd,Eigen::VectorXd::Ones(dd.leaves().size()),1);
    require(bool(dis.product()) && dis.product()->facts().carvedRegion.materialComponents==2 && dis.product()->facts().materialRegion.boundedVoidComponents==2,"disconnected newly carved components");
}
void printSolve(const SourceAdaptiveSpinResult& out) {
    std::cout<<std::setprecision(17)<<"status="<<sourceAdaptiveSpinStatusToken(out.status)<<" input="<<out.inputHash<<" trace="<<out.traceHash<<" output="<<out.outputHash<<" beta="<<out.relaxedFill.transpose()<<'\n';
    for(const auto& p:out.passes) {
        std::cout<<"pass="<<p.index<<" before="<<p.beforeHash<<" after="<<p.afterHash<<" slqp="<<p.solve.traceHash<<'\n';
        for(const auto& t:p.solve.trace) {std::cout<<"{\"iteration\":"<<t.iteration<<",\"objective\":"<<t.objective;writeSlqpAuditJson(std::cout,t);std::cout<<"}\n";}
    }
}
void wrapperTests() {
    auto r=goodRequest(); SourceMaterialBoundarySettings cfg;cfg.minimumOuterWall=.5;
    auto result=solveSourceAdaptiveSpinMaterial(r,cfg);printSolve(result.solve());
    require(result.solve().status==SourceAdaptiveSpinStatus::BinaryMomentValidated && bool(result.boundary().product()),"STOP: fixed actual solve must reach full-carve material; no retuning allowed");
    require((result.solve().binaryFill.array()==1).all(),"actual fixed optimum full removal");
    std::cout<<"wrapper="<<result.outputDigest()<<" material_product="<<result.boundary().facts().productDigest<<'\n';
    auto repeated=solveSourceAdaptiveSpinMaterial(r,cfg);
    require(repeated.outputDigest()==result.outputDigest() && repeated.solve().traceHash==result.solve().traceHash,"wrapper repeat complete identity");
    require(std::abs(result.boundary().product()->facts().integratedMoments.mass()-503./64)<1e-12,"real solver mesh mass");
    auto limited=cfg;limited.unionSettings.maximumInputCuboids=1;
    auto cap=solveSourceAdaptiveSpinMaterial(r,limited);
    require(cap.solve().status==SourceAdaptiveSpinStatus::BinaryMomentValidated && !cap.boundary().product() && cap.boundary().facts().status=="RESOURCE_LIMIT_UNKNOWN","numerical status separate from geometry resource failure");
    r.source=source(box({-.5+1./8192,-.5,0},{.5+1./8192,.5,1}),box({.125,-.125,.625},{.375,.125,.875}),5,{-1,-1,0},2);
    r.maximumOptimizationLevel=5;
    auto quarter=solveSourceAdaptiveSpinMaterial(r,cfg);printSolve(quarter.solve());
    require(quarter.solve().status==SourceAdaptiveSpinStatus::BinaryValidationFailed && !quarter.boundary().product(),"quarter-fill counterexample never material success");
    r.maximumAdaptivePasses=1;auto unknown=solveSourceAdaptiveSpinMaterial(r,cfg);
    require(unknown.solve().status==SourceAdaptiveSpinStatus::AdaptivePassLimit && !unknown.boundary().product(),"pass UNKNOWN no material product");
    limited=cfg;limited.unionSettings.maximumOutputVertices=0;
    rejects([&]{solveSourceAdaptiveSpinMaterial(r,limited);},"boundary caps validated before unsuccessful solve");
    r=goodRequest();r.source=source(cube(-.5,.5),cube(.39,.40),3,{-1,-1,-1},2);r.maximumOptimizationLevel=3;
    for(auto obj:{Spin2014Objective::StaticLowering,Spin2014Objective::YoYo,Spin2014Objective::TopInertiaOnly,Spin2014Objective::TopLoweringOnly,Spin2014Objective::TopMassAndLowering}) {
        r.objective=obj;r.initialPhi=obj==Spin2014Objective::StaticLowering?0:.37;
        auto empty=solveSourceAdaptiveSpinMaterial(r,cfg);
        require(empty.solve().status==SourceAdaptiveSpinStatus::RelaxationUnknown && !empty.boundary().product(),"all five empty-source UNKNOWN states preserved");
    }
    static_assert(!std::is_default_constructible<SourceMaterialBoundary>::value,"no caller fabricated accepted product");
    static_assert(std::is_const<std::remove_reference_t<decltype(result.solve())>>::value,"immutable wrapper solve access");
}
void zeroCarveBudgetTest() {
    // Catch an empty union bypassing all subsequent original/design/final Nefs.
    auto r=goodRequest();SourceAdaptiveDomain domain(r.source,2,4);
    SourceMaterialBoundarySettings limits;limits.unionSettings.maximumIntermediateComplexity=4;
    auto result=makeSourceMaterialBoundary(domain,Eigen::VectorXd::Zero(domain.leaves().size()),1,limits);
    std::cout<<"zero-carve cap4 status="<<result.facts().status<<" product="<<bool(result.product())<<'\n';
    require(!result.product() && result.facts().status=="RESOURCE_LIMIT_UNKNOWN","zero carve cannot bypass cap4 through ordinary-cube material audits");
    require(result.facts().resourceBudget.exceeded && result.facts().resourceBudget.currentComplexity==0,
        "failed zero-carve chain retains measured rejection and releases all Nefs");
    auto baseline=makeSourceMaterialBoundary(domain,Eigen::VectorXd::Zero(domain.leaves().size()),1);
    const auto peak=baseline.facts().resourceBudget.peakComplexity;
    limits.unionSettings.maximumIntermediateComplexity=peak;
    auto equal=makeSourceMaterialBoundary(domain,Eigen::VectorXd::Zero(domain.leaves().size()),1,limits);
    require(bool(equal.product()),"zero-carve complete original/design/final/partition exact peak equality");
    limits.unionSettings.maximumIntermediateComplexity=peak-1;
    auto below=makeSourceMaterialBoundary(domain,Eigen::VectorXd::Zero(domain.leaves().size()),1,limits);
    require(!below.product() && below.facts().status=="RESOURCE_LIMIT_UNKNOWN","zero-carve full-chain one-below peak fails closed");
    std::cout<<"zero-carve full-chain peak="<<peak<<'\n';
}
void auditorBudgetTests() {
#ifdef SOURCE_MATERIAL_HAS_SCOPED_BUDGET
    static_assert(!std::is_copy_constructible<ExactNefBudget>::value && !std::is_copy_assignable<ExactNefBudget>::value,
        "live operation contexts cannot be copied or reset under existing leases");
    const auto shape=join(cube(0,3),cube(1,2),true);
    auto unlimited=std::make_shared<ExactNefBudget>(100000);
    auto region=makeExactMaterialRegionWithBudget(shape.v,shape.f,unlimited);
    require(region.facts().valid,"budgeted signed nesting build");
    const auto builderPeak=unlimited->facts().peakComplexity;
    require(builderPeak>138,"construction/nesting/assembly count simultaneous explicit solids");
    auto at=std::make_shared<ExactNefBudget>(builderPeak);
    auto same=makeExactMaterialRegionWithBudget(shape.v,shape.f,at);
    require(same.facts().valid && same.facts().regionDigest==makeExactMaterialRegion(shape.v,shape.f).facts().regionDigest,"budget equality preserves exact legacy region identity");
    auto below=std::make_shared<ExactNefBudget>(builderPeak-1);
    auto denied=makeExactMaterialRegionWithBudget(shape.v,shape.f,below);
    require(!denied.facts().valid && below->facts().exceeded,"signed auditor must reject one below actual peak");
    auto meshA=cube(0,2),meshB=cube(1,3);
    const auto a=makeExactMaterialRegion(meshA.v,meshA.f),b=makeExactMaterialRegion(meshB.v,meshB.f);
    auto relationBudget=std::make_shared<ExactNefBudget>(100000);
    const auto relation=relateExactMaterialRegionsWithBudget(a,b,relationBudget);
    require(relation.completed && relation.positiveVolumeOverlap==true,"budgeted overlap relation");
    const auto unitMesh=box({0,0,0},{1,1,1}),shiftedMesh=box({0.5,0,0},{1.5,1,1});
    const auto unit=makeExactMaterialRegion(unitMesh.v,unitMesh.f);
    const auto shifted=makeExactMaterialRegion(shiftedMesh.v,shiftedMesh.f);
    const auto exactDifference=relateExactMaterialRegions(unit,shifted);
    require(exactDifference.completed && exactDifference.equal==false &&
        exactDifference.symmetricDifferenceVolumeExact==std::optional<std::string>("1/1"),
        "exact relation reports rational symmetric-difference volume for shifted unit cubes");
    const auto exactEqual=relateExactMaterialRegions(unit,unit);
    require(exactEqual.completed && exactEqual.equal==true &&
        exactEqual.symmetricDifferenceVolumeExact==std::optional<std::string>("0/1"),
        "exact equality has zero regularized symmetric-difference volume");
    const auto relationPeak=relationBudget->facts().peakComplexity;
    require(relationPeak>92 && relationBudget->facts().currentComplexity==0,"relation counts both existing operands and live Boolean results, then releases leases");
    at=std::make_shared<ExactNefBudget>(relationPeak);
    require(relateExactMaterialRegionsWithBudget(a,b,at).completed,"relation exact peak equality");
    below=std::make_shared<ExactNefBudget>(relationPeak-1);
    auto failedRelation=relateExactMaterialRegionsWithBudget(a,b,below);
    require(!failedRelation.completed && !failedRelation.equal && below->facts().exceeded,"relation cap cannot publish partial facts");
    const auto left=box({0,0,0},{1,2,2}),right=box({1,0,0},{2,2,2});
    const std::vector<ExactMaterialRegion> parts={makeExactMaterialRegion(left.v,left.f),makeExactMaterialRegion(right.v,right.f)};
    auto partitionBudget=std::make_shared<ExactNefBudget>(100000);
    auto partition=checkExactMaterialPartitionWithBudget(a,parts,partitionBudget);
    require(partition.completePartition==true,"budgeted actual inclusive partition");
    const auto partitionPeak=partitionBudget->facts().peakComplexity;
    require(partitionPeak>138 && partitionBudget->facts().currentComplexity==0,"partition counts whole/parts/union/results then releases");
    at=std::make_shared<ExactNefBudget>(partitionPeak);
    require(checkExactMaterialPartitionWithBudget(a,parts,at).completePartition==true,"partition exact peak equality");
    below=std::make_shared<ExactNefBudget>(partitionPeak-1);
    auto failedPartition=checkExactMaterialPartitionWithBudget(a,parts,below);
    require(!failedPartition.completed && !failedPartition.completePartition && below->facts().exceeded,"partition cap cannot publish partial facts");
    // Persisted region solids remain charged until their last owner releases.
    auto lifetime=std::make_shared<ExactNefBudget>(100000);
    {auto owned=makeExactMaterialRegionWithBudget(meshA.v,meshA.f,lifetime);
        require(owned.facts().valid && lifetime->facts().currentComplexity==46,"actual cube Nef retained by region");}
    require(lifetime->facts().currentComplexity==0,"region destruction releases actual solid lease");
    auto aliases=std::make_shared<ExactNefBudget>(140);
    require(relateExactMaterialRegionsWithBudget(a,a,aliases).equal==true &&
        aliases->facts().peakComplexity==140 && aliases->facts().observationsByOperation.at("existing operand")==1,
        "same existing operand is admitted once, not twice");
    auto ownedBudget=std::make_shared<ExactNefBudget>(100000);
    auto owned=makeExactMaterialRegionWithBudget(meshA.v,meshA.f,ownedBudget);
    const auto admissions=ownedBudget->facts().observationsByOperation.count("existing operand");
    require(relateExactMaterialRegionsWithBudget(owned,owned,ownedBudget).equal==true &&
        ownedBudget->facts().observationsByOperation.count("existing operand")==admissions,
        "already charged context operands are not double-counted");
    const auto ownedSnapshot=hashExactNefBudgetFacts(ownedBudget->facts());
    require(relateExactMaterialRegions(owned,owned).equal==true &&
        hashExactNefBudgetFacts(ownedBudget->facts())==ownedSnapshot,
        "absent-budget legacy operation never mutates an operand context");
    ExactCuboidUnionSettings settings;
    auto shared=std::make_shared<ExactNefBudget>(settings.maximumIntermediateComplexity);
    auto output=makeExactCuboidUnionWithBudget({cuboid({0,0,0},{1,1,1})},settings,shared);
    require(bool(output.region()),"shared union success");
    const auto frozen=hashExactNefBudgetFacts(output.facts().budget),identity=output.facts().outputDigest;
    {auto later=makeExactMaterialRegionWithBudget(meshA.v,meshA.f,shared);
        require(later.facts().valid,"later context operation");}
    require(output.facts().outputDigest==identity && hashExactNefBudgetFacts(output.facts().budget)==frozen &&
        hashExactNefBudgetFacts(shared->facts())!=frozen,"published budget snapshot and identity cannot mutate across later operations");
    rejects([]{ExactNefBudget invalid(0);},"zero shared budget rejected");
    std::cout<<"scoped budgets builder_peak="<<builderPeak<<" relation_peak="<<relationPeak<<" partition_peak="<<partitionPeak<<'\n';
#else
    require(false,"RED: scoped shared-auditor Nef budget operations are unavailable");
#endif
}
int main(int argc,char** argv) {try {
    if(argc==1 || std::string(argv[1])=="--budget-zero") zeroCarveBudgetTest();
    if(argc==1 || std::string(argv[1])=="--auditor-budget") auditorBudgetTests();
    if(argc==1 || std::string(argv[1])=="--union") unionTests();
    if(argc==1 || std::string(argv[1])=="--cap") exactOutputCapTests();
    if(argc==1 || std::string(argv[1])=="--material") materialTests();
    if(argc==1 || std::string(argv[1])=="--dual-density") sourceDualDensityFunctionalTests();
    if(argc==1 || std::string(argv[1])=="--solve") wrapperTests();
    std::cout<<"SYNTHETIC_REGRESSION source material tests passed; faster="<<MAKE_IT_STAND_HAS_SPIN_IT_FASTER<<'\n';return 0;
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}}
#else
#include <iostream>
int main(){std::cerr<<"RED: exact source cuboid union/material boundary operation is unavailable\n";return 1;}
#endif
