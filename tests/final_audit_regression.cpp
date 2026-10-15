#include "spin/CanonicalFullRankHollowing.h"
#include <Eigen/Geometry>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace spin {
// Only the separately compiled test subject refers to these hooks.
FullRankHollowingResult auditTestHollowingSolve(const SpinDesignRequest&,
    const QuadricMomentEvaluatorSettings&, const std::vector<QuadricCoefficients10>&);
}
namespace {
using V = spin::QuadricCoefficients10::Vector;
std::string fault;
bool numericalFinished = false;
int injections = 0;
void require(bool value, const std::string& diagnostic) {
    if (!value) throw std::runtime_error(diagnostic);
}
physics::RawMassMoments10 raw(const V& v) {
    MassMoments m; m.mass=v(0); m.firstMoment=v.segment<3>(1);
    m.secondMoment << v(4),v(9),v(8),v(9),v(5),v(7),v(8),v(7),v(6);
    return physics::RawMassMoments10::fromMassMoments(m);
}
spin::TriangleMesh box(double a, double b) {
    spin::TriangleMesh m; m.vertices.resize(8,3); m.faces.resize(12,3);
    m.vertices << a,a,a,b,a,a,b,b,a,a,b,a,a,a,b,b,a,b,b,b,b,a,b,b;
    m.faces << 0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,
        3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;
    return m;
}
// Independent volume formula: outer cube minus inner cube, uniform density 1.
// No production moment integration/transform helper computes this expectation.
physics::RawMassMoments10 analyticShell(const Eigen::Vector3d& c) {
    const double mass = 8.0 - 8.0*.8*.8*.8;
    const double central = (8.0 - 8.0*.8*.8*.8*.8*.8)/3.0;
    V v;
    v << mass, mass*c.x(), mass*c.y(), mass*c.z(),
        central+mass*c.x()*c.x(), central+mass*c.y()*c.y(),
        central+mass*c.z()*c.z(), mass*c.y()*c.z(),
        mass*c.x()*c.z(), mass*c.x()*c.y();
    return raw(v);
}
struct Fixture {
    spin::SpinDesignRequest request;
    spin::QuadricMomentEvaluatorSettings settings;
    spin::QuadricCoefficients10 seed;
};
Fixture fixture() {
    spin::SpinDesignRequest r;
    r.designDomain=box(-.8,.8); r.outerSurface=box(-1,1);
    r.designDomain.vertices.col(2).array()+=2;
    r.outerSurface.vertices.col(2).array()+=2;
    r.materialDensity=1; r.minimumWallThickness=.15;
    r.mode=spin::SpinMode::SpinningTop; r.centerOfMassWeight=1;
    r.settings.equationProfile=spin::SpinEquationProfile::DualMonomialRawV1;
    r.settings.maximumRestarts=1; r.settings.maximumFunctionEvaluationsPerRestart=20;
    r.fixedShellMoments=analyticShell(Eigen::Vector3d(0,0,2));
    spin::QuadricMomentEvaluatorSettings s;
    s.triangulation.icosphereSubdivisions=1; s.triangulation.truncationRadius=8;
    V p; p << -.31,.07,-.09,.11,.83,1.17,-.72,.06,-.04,.03;
    const auto base=spin::QuadricCoefficients10::fromUnnormalized(p);
    const Eigen::Vector3d c(0,0,2);
    p=base.values();
    p(0)=base.constant()-base.linearVector().dot(c)+c.dot(base.quadraticMatrix()*c);
    p.segment<3>(1)=base.linearVector()-2*base.quadraticMatrix()*c;
    return {r,s,spin::QuadricCoefficients10::fromUnnormalized(p)};
}
void requireRealConvergence(const spin::FullRankHollowingResult& result) {
    require(result.numerical.status==spin::SpinSolverStatus::Converged,
        "real numerical solver failed to converge");
    require(result.numerical.functionEvaluations>1 && result.numerical.trace.size()>1,
        "not a nontrivial real numerical solve");
}
void requireVerified(const spin::FullRankHollowingResult& result) {
    requireRealConvergence(result);
    require(result.finalGeometryVerified && result.finalEvaluation &&
        result.independentTotalMoments && result.independentlyRecomputedRawResidual &&
        result.theoremOneNumericalAudit,
        "final audit failed: "+result.verificationDiagnostic);
    require(*result.independentlyRecomputedRawResidual<=1e-12,
        "original raw residual tolerance was not met");
    require(
        result.theoremOneNumericalAudit->status ==
            spin::FasterTheoremConditionStatus::SatisfiedNumerically,
        "final material passed without the narrow Eq. 7 condition audit"
    );
    require(
        result.theoremOneNumericalAudit->proofScope ==
            "NUMERICAL_EQ7_AND_EQUALITY_CONDITIONS_ONLY_NOT_A_THEOREM_PROOF",
        "final material widened the theorem-audit claim"
    );
}
void frameTest() {
    auto f=fixture();
    fault="control"; numericalFinished=false; injections=0;
    const auto base=spin::solveCanonicalFullRankQuadricHollowing(f.request,f.settings,{f.seed});
    requireVerified(base);
    Eigen::Matrix3d axes; axes << 0,0,1,1,0,0,0,1,0;
    const Eigen::Vector3d origin(3,-2,1);
    for (auto* mesh : {&f.request.designDomain,&f.request.outerSurface})
        for (Eigen::Index i=0;i<mesh->vertices.rows();++i)
            mesh->vertices.row(i)=(axes*mesh->vertices.row(i).transpose()+origin).transpose();
    f.request.targetFrame.axes=axes; f.request.targetFrame.origin=origin;
    f.request.fixedShellMoments=analyticShell(Eigen::Vector3d(5,-2,1));
    const auto moved=spin::auditTestHollowingSolve(f.request,f.settings,{f.seed});
    requireVerified(moved);
    const auto publicMoved=spin::solveCanonicalFullRankQuadricHollowing(f.request,f.settings,{f.seed});
    requireVerified(publicMoved);
    require(publicMoved.independentTotalMoments->values()==moved.independentTotalMoments->values(),
        "instrumented frame audit differs from the public canonical path");
    const auto& a=base.independentTotalMoments->values();
    const auto& b=moved.independentTotalMoments->values();
    const double error=(a-b).norm()/a.norm();
    require(error<=1e-10,"nonidentity source frame changed final target moments");
    const auto shell=moved.finalEvaluation->totalMomentsTarget.values()-
        moved.finalEvaluation->variableMomentsTarget.values();
    const double shellError=(shell-analyticShell(Eigen::Vector3d(0,0,2)).values()).norm();
    require(shellError<=1e-12,"fixed shell was not transformed to the target frame");
    require(base.inputHash!=publicMoved.inputHash,
        "changed source representation omitted from canonical input identity");
    // Compare the same public policy namespace and change only the frame origin.
    // These budget-limited identity probes are not successful design evidence.
    auto hashRequest=f.request;
    hashRequest.settings.maximumFunctionEvaluationsPerRestart=1;
    const auto hashBase=spin::solveCanonicalFullRankQuadricHollowing(hashRequest,f.settings,{f.seed});
    hashRequest.targetFrame.origin.x()+=.125;
    const auto hashChanged=spin::solveCanonicalFullRankQuadricHollowing(hashRequest,f.settings,{f.seed});
    require(hashBase.inputHash!=hashChanged.inputHash,"frame origin omitted from input identity");
    std::cout << "target_moment_relative_error=" << error << " shell_absolute_error="
        << shellError << " raw=" << *moved.independentlyRecomputedRawResidual << '\n';
}
void auditTest(const std::string& mode) {
    const auto f=fixture(); fault=mode; numericalFinished=false; injections=0;
    const auto out=spin::auditTestHollowingSolve(f.request,f.settings,{f.seed});
    requireRealConvergence(out);
    if (mode=="control") {
        requireVerified(out);
        const auto actual=spin::solveCanonicalFullRankQuadricHollowing(f.request,f.settings,{f.seed});
        requireVerified(actual);
        require(actual.numerical.eventHash==out.numerical.eventHash &&
            actual.numerical.traceHash==out.numerical.traceHash &&
            actual.independentTotalMoments->values()==out.independentTotalMoments->values(),
            "test instrumentation changed unmodified numerical/final result");
        require(injections==0,"no-fault control injected a fault");
    } else {
        require(injections==1,"fault was not injected exactly once after real convergence");
        require(!out.finalGeometryVerified,"damaged converged result was incorrectly certified");
        require(!out.verificationDiagnostic.empty() && !out.outputHash.empty(),
            "final rejection lost its diagnostic/identity");
        require(out.finalVerificationSeconds>0,"post-convergence audit did not execute");
    }
    std::cout << "numerical=" << spin::spinSolverStatusToken(out.numerical.status)
        << " verified=" << out.finalGeometryVerified << " evaluations="
        << out.numerical.functionEvaluations << " diagnostic=" << out.verificationDiagnostic << '\n';
}
}
namespace spin {
SpinNewtonResult auditTestNumericalSolve(const SpinMomentProblem& problem,
    const SpinSolverSettings& settings, const SpinMomentEvaluator& evaluator,
    const std::vector<QuadricCoefficients10>& guesses) {
    auto result=solveSpinMomentSystem(problem,settings,evaluator,guesses);
    numericalFinished=true;
    if (result.status!=SpinSolverStatus::Converged) return result;
    if (fault=="stale-lambda") { require(result.multipliers.size()>0,"missing fixture multiplier");
        result.multipliers(0)+=.5; ++injections; }
    else if (fault=="invalid-mu") { result.multiplierMu=-1; ++injections; }
    else if (fault=="missing-quadric") { result.quadric.reset(); ++injections; }
    else if (fault=="residual") { result.rawEquationResidual=.1; ++injections; }
    return result;
}
FullRankMomentEvaluation auditTestMaterialEvaluation(const SpinDesignRequest& r,
    const QuadricCoefficients10& q, const QuadricMomentEvaluatorSettings& s) {
    auto value=evaluateCanonicalFullRankMaterialMoments(r,q,s);
    if (!numericalFinished) return value;
    if (fault=="missing-j") { value.coefficientJacobian.reset(); ++injections; }
    else if (fault=="nonfinite-j") {
        require(value.coefficientJacobian.has_value(),"fixture derivative missing");
        (*value.coefficientJacobian)(0,0)=std::numeric_limits<double>::quiet_NaN(); ++injections;
    } else if (fault=="moments") {
        V changed=value.totalMomentsTarget.values(); changed(0)+=.125;
        value.totalMomentsTarget=raw(changed); ++injections;
    } else if (fault=="inward") {
        value.materialBoundaryTarget.faces.col(0).swap(value.materialBoundaryTarget.faces.col(2));
        ++injections;
    }
    return value;
}
}
int main(int argc,char** argv) {
    try {
        require(argc==2,"one registered test mode required");
        const std::string mode=argv[1]; std::cout<<std::setprecision(17);
        if (mode=="frame") frameTest(); else auditTest(mode);
        std::cout<<"[PASS] "<<mode<<'\n'; return 0;
    } catch (const std::exception& e) { std::cerr<<"[FAIL] "<<e.what()<<'\n'; return 1; }
}
