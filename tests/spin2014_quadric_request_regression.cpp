#include "spin2014/QuadricHollowing.h"
#include "spin/CanonicalFullRankHollowing.h"
#include "manufacturing/ManufacturingGeometry.h"
#include <Eigen/Geometry>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
#include "verification/ExactRationalCoordinates.h"
#endif

namespace {
void require(bool ok,const char* why) {if(!ok)throw std::runtime_error(why);}
spin::TriangleMesh box(double a,double b) {
    spin::TriangleMesh m;m.vertices.resize(8,3);m.faces.resize(12,3);
    m.vertices<<a,a,a,b,a,a,b,b,a,a,b,a,a,a,b,b,a,b,b,b,b,a,b,b;
    m.faces<<0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;
    return m;
}
spin::SpinDesignRequest request() {
    spin::SpinDesignRequest r;r.outerSurface=box(-2,2);r.designDomain=box(-1,1);
    r.mode=spin::SpinMode::SpinningTop;r.materialDensity=5;
    r.settings.equationProfile=spin::SpinEquationProfile::DualMonomialRawV1;
    r.settings.maximumRestarts=1;r.settings.maximumFunctionEvaluationsPerRestart=1;
    return r;
}
template<class F> void rejected(F f) {
    bool ok=false;try {f();} catch(const std::invalid_argument&) {ok=true;}
    require(ok,"invalid geometry/objective input was admitted");
}
void admission() {
    auto r=request();spin2014::QuadricObjectiveSettings objective;
    spin::QuadricMomentEvaluatorSettings legacyAudit,exactAudit;
    exactAudit.exactBinary64SelfIntersection=true;
    require(spin::computeQuadricEvaluatorRequestHash(r,legacyAudit)!=
            spin::computeQuadricEvaluatorRequestHash(r,exactAudit),"explicit exact audit policy is missing from effective identity");
    require(spin::auditQuadricMaterialBoundary(r.designDomain,true).strictBoundaryValid(),
            "exact intersection policy rejects a legal box");
    auto overlapping=box(-1,1);auto other=overlapping;
    other.vertices.col(0).array()+=.5;
    overlapping.vertices.conservativeResize(16,3);overlapping.vertices.bottomRows(8)=other.vertices;
    overlapping.faces.conservativeResize(24,3);overlapping.faces.bottomRows(12)=other.faces.array()+8;
    require(!spin::auditQuadricMaterialBoundary(overlapping,true).strictBoundaryValid(),
            "exact policy bypassed a genuine crossing between closed components");
    auto small=box(0,1e-8);small.vertices.col(0).array()+=3;
    auto mixed=box(-1,1);mixed.vertices.conservativeResize(16,3);mixed.vertices.bottomRows(8)=small.vertices;
    mixed.faces.conservativeResize(24,3);mixed.faces.bottomRows(12)=small.faces.array()+8;
    require(!spin::auditQuadricMaterialBoundary(mixed,true).nondegenerateFaces,
            "exact intersection policy relaxed the existing scaled-area gate");
    auto open=r.designDomain;open.faces.conservativeResize(11,3);
    require(!spin::auditQuadricMaterialBoundary(open,true).strictBoundaryValid(),
            "exact intersection policy admitted an open boundary");
    const auto owned=spin2014::prepareSpin2014QuadricRequest(r,objective,3);
    require(std::abs(owned.fixedShellMoments.mass()-168)<1e-12,"shell density not applied exactly once to volume 64-8");
    require(owned.materialDensity==5,"shell density overwrote design density");
    require(r.fixedShellMoments.mass()==0 && owned.outerSurface.vertices==r.outerSurface.vertices &&
            owned.outerSurface.faces==r.outerSurface.faces,"request mutation or outer-shape change");
    const auto equal=spin2014::prepareSpin2014QuadricRequest(r,objective,5);
    require(std::abs(equal.fixedShellMoments.mass()+40-320)<1e-12,"equal density is not equivalent to the full single-material box");
    auto translated=r;translated.targetFrame.origin=Eigen::Vector3d(.2,-.1,.3);
    translated.targetFrame.axes=Eigen::AngleAxisd(.3,Eigen::Vector3d::UnitZ()).toRotationMatrix();
    const auto moved=spin2014::prepareSpin2014QuadricRequest(translated,objective,3);
    const auto target=moved.fixedShellMoments.transformed(moved.targetFrame.axes.transpose(),
        -moved.targetFrame.axes.transpose()*moved.targetFrame.origin);
    require((target.centerOfMass()+moved.targetFrame.axes.transpose()*moved.targetFrame.origin).norm()<1e-14,
            "shell target frame transformed twice or not at all");
    rejected([&] {spin2014::prepareSpin2014QuadricRequest(r,objective,-1);});
    auto bad=r;bad.mode=spin::SpinMode::YoYo;
    rejected([&] {spin2014::prepareSpin2014QuadricRequest(bad,objective,3);});
    bad=r;bad.fixedShellMoments=owned.fixedShellMoments.scaledDensity(2);
    rejected([&] {spin2014::prepareSpin2014QuadricRequest(bad,objective,3);});
    bad=r;bad.designDomain.vertices.col(0).array()+=3;
    rejected([&] {spin2014::prepareSpin2014QuadricRequest(bad,objective,3);});
    bad=r;bad.designDomain.faces.row(0).swap(bad.designDomain.faces.row(1));
    bad.designDomain.faces(0,1)=bad.designDomain.faces(0,0);
    rejected([&] {spin2014::prepareSpin2014QuadricRequest(bad,objective,3);});
    std::cout<<"[PASS] owned real shell, unequal densities, target frame, domain and invalid-input admission\n";
}
void realCall() {
    auto r=request();spin::QuadricMomentEvaluatorSettings es;
    es.derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
    es.triangulation.icosphereSubdivisions=2;es.triangulation.truncationRadius=8;
    spin::QuadricCoefficients10::Vector p;
    p<<-.31,.07,-.09,.11,.83,1.17,-.72,.06,-.04,.03;
    const auto q=spin::QuadricCoefficients10::fromUnnormalized(p);
    for(auto objective : {spin2014::Spin2014Objective::YoYo,spin2014::Spin2014Objective::TopInertiaOnly,
                         spin2014::Spin2014Objective::TopLoweringOnly,spin2014::Spin2014Objective::TopMassAndLowering}) {
        spin2014::QuadricObjectiveSettings settings;settings.objective=objective;
        r.mode=spin2014::spin2014ObjectiveMode(objective);
        const auto run=spin2014::solveSpin2014QuadricHollowing(r,settings,3,es,{q});
        require(run.objectiveProfile=="SPIN2014_OBJECTIVES_QUADRIC_V1" &&
                run.objectiveIdentity==spin2014::spin2014QuadricObjectiveHash(settings),"wrong objective identity in material result");
        require(run.solve.numerical.functionEvaluations==1,"one-evaluation diagnostic exceeded budget");
        const auto value=spin::evaluateCanonicalFullRankMaterialMoments(run.measuredRequest,q,es);
        const auto expected=spin2014::makeSpin2014QuadricMomentProblem(settings).evaluate(value.totalMomentsTarget.values());
        bool observed=false;
        for(const auto& event:run.solve.numerical.events) if(event.objective) {
            require(std::abs(*event.objective-expected.objective)<1e-11,"Newton factory ignored the actual 2014 objective on Boolean geometry");
            observed=true;
        }
        require(observed,"no actual objective evaluation recorded");
        require(!run.solve.finalGeometryVerified,"budget-only diagnostic falsely claimed a qualified final material");
        std::cout<<"actual-boolean objective="<<spin2014::spin2014ObjectiveToken(objective)
                 <<" status="<<spin::spinSolverStatusToken(run.solve.numerical.status)<<'\n';
    }
    std::cout<<"[PASS] real CGAL Boolean maps consume four original rotational objectives (diagnostic, not solved case)\n";
}
void rationalComplement() {
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    // Geometry-only production call. No converged state or solver pass is
    // fabricated to test a saved material's same-retained complement.
    auto r=request();spin::QuadricMomentEvaluatorSettings es;
    es.derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
    es.triangulation.rationalSourceCoordinates=true;
    es.triangulation.icosphereSubdivisions=1;
    spin::QuadricCoefficients10::Vector p;p<<-.25,0,0,0,1,1,1,0,0,0;
    const auto retained=spin::evaluateCanonicalFullRankMaterialMoments(r,
        spin::QuadricCoefficients10::fromUnnormalized(p),es);
    const auto cavity=spin2014::exportSpin2014RationalCavity(r,retained);
    const auto whole=spin::integrateCanonicalFullRankBoundaryMoments(
        retained.designExactVerticesTarget,r.designDomain.faces,r.materialDensity);
    for(int i=0;i<10;++i)require(
        verification::parseExactRationalCoordinate(whole[i])==
        verification::parseExactRationalCoordinate(retained.variableMomentsExactTarget[i])+
        verification::parseExactRationalCoordinate(cavity.variableMomentsExactTarget[i]),
        "same-retained rational complement moment partition differs");
    require(!cavity.coefficientJacobian && cavity.sourceScientificHash==retained.scientificHash,
        "complement falsely claims a quadric Jacobian or loses its retained-material identity");
    bool caught=false;
    auto bad=retained;bad.variableMomentsExactTarget[0]="1/1";
    try{spin2014::exportSpin2014RationalCavity(r,bad);}catch(const std::runtime_error&){caught=true;}
    require(caught,"complement admitted stale same-call material moments");
    caught=false;auto wrongDensity=r;wrongDensity.materialDensity=2.;
    try{spin2014::exportSpin2014RationalCavity(wrongDensity,retained);}
    catch(const std::runtime_error&){caught=true;}
    require(caught,"complement reused authority at the wrong density");
    caught=false;bad=retained;bad.sourceConstructionIdentity.clear();
    try{spin2014::exportSpin2014RationalCavity(r,bad);}catch(const std::runtime_error&){caught=true;}
    require(caught,"complement admitted missing source identity");
    std::cout<<"[PASS] same-retained rational cavity and exact independent partition\n";
#else
    throw std::runtime_error("rational complement test requires exact geometry");
#endif
}
}
int main(int argc,char**argv) {
    try {
        std::string arg=argc>1?argv[1]:"--all";
        if(arg=="--admission"||arg=="--all")admission();
        if(arg=="--real-call"||arg=="--all")realCall();
        if(arg=="--rational-complement")rationalComplement();
        return EXIT_SUCCESS;
    }catch(const std::exception&e){std::cerr<<"[FAIL] "<<e.what()<<'\n';return EXIT_FAILURE;}
}
