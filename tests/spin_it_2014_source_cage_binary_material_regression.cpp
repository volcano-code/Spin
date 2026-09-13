#include <iostream>
#if !__has_include("spin2014/SourceCageBinaryMaterial.h")
int main(){std::cerr<<"FAIL: actual deformed binary material admission and validation layer is absent\n";return 1;}
#else
#include "spin2014/SourceCageBinaryMaterial.h"
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
SourceCagePaperFunctional fixture(int level=4) {
    const auto outer=cube(.3125,.6875),design=cube(.425,.575);
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
void candidateChecks(const SourceCageBinaryMaterialAudit& audit) {
    const auto m=audit.candidate();require(bool(m),"candidate geometry retained");
    const auto& f=m->functional();const auto values=f.evaluateValues(m->binaryState());
    require(identical(m->relaxedState().tail(m->relaxedState().size()-f.betaCount()),
        m->binaryState().tail(m->binaryState().size()-f.betaCount())),"rounding preserves phi and all handles bit-for-bit");
    require((m->outerBoundary().triangles.array()==f.source().boundary().binding().source().outerFaces().array()).all(),"actual original outer connectivity retained");
    require((m->outerBoundary().vertices-values.world.deformedPoints.topRows(m->outerBoundary().vertices.rows())).norm()==0.,"emitted outer coordinates are actual binary deformation");
    require((m->outerBoundary().vertices-f.source().boundary().binding().source().outerVertices()).norm()>1e-3,"nonidentity deformation exercised");
    require((audit.facts().worldIntegratedMoments-values.world.finalMoments.values()).cwiseAbs().maxCoeff()<1e-12,"all ten independently emitted moments agree with actual functional");
    require((audit.facts().spinIntegratedMoments-values.spinMoments.values()).cwiseAbs().maxCoeff()<1e-12,"fixed frame independent transform agrees");
    require((audit.facts().rawConstraints-values.constraints).norm()<1e-12,"Givens constraints recomputed on actual binary output");
    require((audit.facts().centralInertia-values.spinMoments.inertiaAboutCenter()).norm()<1e-12,"central inertia uses independent binary moments in fixed nonidentity paper frame");
    require(audit.facts().centralInertiaPsd&&audit.facts().centralInertiaEigenvalues.minCoeff()>0.,"positive material has checked central inertia spectrum");
    require(identical(audit.facts().stateDelta,m->binaryState()-m->relaxedState()),"full raw rounding delta retained");
    Eigen::VectorXd normalized=values.constraints;
    for(Eigen::Index i=0;i<normalized.size();++i)normalized[i]/=f.source().density()*std::pow(f.source().boundary().binding().source().config().rootSide,i<2?4:5);
    require((audit.facts().normalizedConstraints-normalized).norm()<1e-12,"binary Givens residuals use fixed INPUT scales");
    require(audit.facts().geometry&&audit.facts().exactBoundary,"both actual checker and final boundary audit run");
}
}
int main(){try {
    static_assert(!std::is_default_constructible_v<SourceCageBinaryMaterialAudit>);
    static_assert(!std::is_default_constructible_v<SourceCageBinaryMaterial>);
    static_assert(std::is_const_v<std::remove_reference_t<decltype(std::declval<SourceCageBinaryMaterial>().binaryState())>>);
    const auto f=fixture();require(f.betaCount()==8,"eight actual source owners");
    for(double beta:{0.,.25,.5,.75,1.}) {
        const auto x=state(f,beta);const auto ep=endpoint(f,x);
        const auto audit=validateSourceCageBinaryMaterial(ep);candidateChecks(audit);
        require(audit.facts().outcome==SourceCageBinaryMaterialOutcome::Validated&&audit.validatedMaterial(),"stationary symmetric binary validates");
        require((audit.candidate()->binaryState().head(audit.candidate()->functional().betaCount()).array()==(beta>=.5?1.:0.)).all(),"zero/one and exact .5 tie round only beta");
        require(audit.facts().emptyCavity==(beta<.5),"empty cavity explicit");
        require(audit.candidate()->carvedBoundary().triangles.rows()==(beta<.5?0:96),"internal faces cancelled in carved union");
        require(audit.facts().exactBoundary->facts().boundaryValid,"actual emitted material exact audit valid");
        require(audit.facts().normalizedMomentResiduals.cwiseAbs().maxCoeff()<1e-10,"fixed input normalization gate");
        require(audit.facts().inputSha256.size()==64&&audit.facts().functionalSha256.size()==64&&audit.facts().partitionSha256.size()==64&&
            audit.facts().configSha256.size()==64&&audit.facts().outputSha256.size()==64,"complete immutable identity digests");
        if(beta==.5) {
            const auto repeated=validateSourceCageBinaryMaterial(ep);
            require(repeated.facts().scientificSha256==audit.facts().scientificSha256,"deterministic complete evidence digest");
        }
        std::cout<<"KNOWN_SOLUTION_VALIDATION beta="<<beta<<" PASS\n";
    }
    auto mixed=state(f,.4);
    const auto& leaves=f.source().boundary().binding().domain().leaves();
    for(std::size_t i=0;i<leaves.size();++i)if(leaves[i].coordinates[0]==8)mixed[i]=.6;
    mixed[0]=.5;
    const auto mixedEndpoint=endpoint(f,mixed);auto mixedAudit=validateSourceCageBinaryMaterial(mixedEndpoint);candidateChecks(mixedAudit);
    require(mixedAudit.facts().outcome==SourceCageBinaryMaterialOutcome::BinaryValidationFailed&&!mixedAudit.validatedMaterial(),"genuine rounding-induced fixed-frame failure stays BINARY_VALIDATION_FAILED");
    require(mixedAudit.facts().rawConstraints.cwiseAbs().maxCoeff()>1e-6,"raw rounding constraint failure visible");
    require(mixedAudit.candidate()->binaryState()[0]==1.,"mixed exact .5 tie hollowed");
    // The same geometry at low density has small RAW first moments, while
    // dimensionless fixed-input residuals retain the genuine axis violation.
    const SourceCagePaperFunctional dilute(SourceCageMomentFunctional(f.source().boundary(),.001,f.source().config()),f.config());
    const auto diluteEndpoint=endpoint(dilute,mixed);
    const auto diluteAudit=validateSourceCageBinaryMaterial(diluteEndpoint);candidateChecks(diluteAudit);
    require(diluteAudit.facts().rawConstraints.cwiseAbs().maxCoeff()<1e-6,"low-density raw residual alone would pass");
    require(diluteAudit.facts().outcome==SourceCageBinaryMaterialOutcome::BinaryValidationFailed&&!diluteAudit.validatedMaterial(),
        "fixed input density normalization must reject despite small raw residual");
    require(diluteAudit.facts().normalizedConstraints.cwiseAbs().maxCoeff()>1e-6,"dimensionless rounding violation retained");
    SourceCageBinaryMaterialSettings smallGeometry;smallGeometry.geometry.maximumPoints=1;
    const auto bothFail=validateSourceCageBinaryMaterial(diluteEndpoint,smallGeometry);
    require(bothFail.facts().geometry&&bothFail.facts().geometry->facts().outcome==SourceCageGeometryOutcome::ResourceLimit&&
        bothFail.facts().outcome==SourceCageBinaryMaterialOutcome::BinaryValidationFailed&&!bothFail.validatedMaterial(),
        "established binary constraint failure stays failed when checker cannot complete");
    auto relaxedCopy=mixedAudit.facts();relaxedCopy.outcome=SourceCageBinaryMaterialOutcome::Validated;
    require(!mixedAudit.validatedMaterial(),"copied facts cannot authorize material");
    const auto& retained=mixedAudit.facts();const auto candidate=mixedAudit.candidate();const auto moved=std::move(mixedAudit);
    require(retained.scientificSha256==moved.facts().scientificSha256,"immutable facts survive move");
    rejects([&]{(void)mixedAudit.facts();},"moved-from access rejected");
    require(candidate->functional().evaluateValues(candidate->binaryState()).outputSha256==moved.facts().binaryFunctionalOutputSha256,"candidate retains owned source lifetime");
    auto diagonal=state(f,.2);
    for(std::size_t i=0;i<leaves.size();++i)if(leaves[i].coordinates[0]==leaves[i].coordinates[1])diagonal[i]=.8;
    const auto nonmanifold=validateSourceCageBinaryMaterial(endpoint(f,diagonal));candidateChecks(nonmanifold);
    require(nonmanifold.facts().outcome==SourceCageBinaryMaterialOutcome::BinaryValidationFailed&&!nonmanifold.validatedMaterial()&&
        nonmanifold.facts().geometry->facts().outcome==SourceCageGeometryOutcome::ChecksPassed&&
        !nonmanifold.facts().exactBoundary->facts().boundaryValid,"edge contact cavity union honestly rejected by final boundary audit without repair");
    const auto complete=endpoint(f,state(f,1.));
    SourceCageBinaryMaterialSettings badGeometry;badGeometry.geometry.minimumOuterDomainClearance=.2;
    const auto geometryRejected=validateSourceCageBinaryMaterial(complete,badGeometry);
    require(geometryRejected.facts().outcome==SourceCageBinaryMaterialOutcome::BinaryValidationFailed&&!geometryRejected.validatedMaterial()&&
        geometryRejected.facts().geometry->facts().outcome==SourceCageGeometryOutcome::RejectedGeometry,"fresh explicit binary geometry settings can reject");
    for(double t:{0.,-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
        SourceCageBinaryMaterialSettings bad;bad.constraintTolerance=t;
        rejects([&]{(void)validateSourceCageBinaryMaterial(complete,bad);},"finite positive tolerance required");
    }
    for(std::size_t cap:{std::size_t(0),std::size_t(1),std::numeric_limits<std::size_t>::max()}) {
        SourceCageBinaryMaterialSettings bad;bad.maximumTriangles=cap;
        const auto limited=validateSourceCageBinaryMaterial(complete,bad);
        require(limited.facts().outcome==SourceCageBinaryMaterialOutcome::ResourceLimit&&!limited.candidate()&&!limited.validatedMaterial(),"invalid or overflowing mesh budgets rejected before emission");
    }
    SourceCageBinaryMaterialSettings tiny;tiny.maximumOwnedBytes=1;
    require(validateSourceCageBinaryMaterial(complete,tiny).facts().outcome==SourceCageBinaryMaterialOutcome::ResourceLimit,"owned buffer cap enforced");
    tiny=SourceCageBinaryMaterialSettings{};tiny.maximumExactDiagnosticBytes=1;
    require(validateSourceCageBinaryMaterial(complete,tiny).facts().outcome==SourceCageBinaryMaterialOutcome::ResourceLimit,"final exact diagnostic cap enforced before allocation");
    const auto coarse=fixture(3);const auto cx=state(coarse,.25);const auto cf=stationary(coarse,cx);
    auto limit=solverSettings();limit.maximumPasses=1;
    const auto pending=solveSourceCageAdaptive(cf,cx,limit);
    require(bool(pending.pendingCheckpoint()),"actual changed checkpoint exists");
    const auto notAdmitted=validateSourceCageBinaryMaterial(pending);
    require(notAdmitted.facts().outcome==SourceCageBinaryMaterialOutcome::NotAdmitted&&!notAdmitted.candidate(),"pending checkpoint never admitted as terminal");
    auto limitedSettings=solverSettings();limitedSettings.budget.maximumKnownBufferBytes=1;
    const auto resource=solveSourceCageAdaptive(cf,cx,limitedSettings);
    require(validateSourceCageBinaryMaterial(resource).facts().outcome==SourceCageBinaryMaterialOutcome::NotAdmitted,"resource endpoint never admitted");
    SourceCagePaperFunctional nonstationary(f.source(),SourceCagePaperFunctionalConfig{},{});
    auto nx=nonstationary.initialState(Eigen::VectorXd::Constant(8,.25),0.);
    const auto stalled=solveSourceCageAdaptive(nonstationary,nx,solverSettings());
    require(stalled.facts().outcome==SourceCageAdaptiveSolverOutcome::SolveNotConverged,"actual nonconverged fixture");
    require(validateSourceCageBinaryMaterial(stalled).facts().outcome==SourceCageBinaryMaterialOutcome::NotAdmitted,"nonconverged endpoint never admitted");
    std::cout<<"binary material regression PASS; bounded pre-seal object only; no sealing, feasibility, physics or paper grant\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
#endif
