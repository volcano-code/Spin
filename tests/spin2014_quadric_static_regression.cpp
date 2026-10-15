// Fresh squared-height numerical gate; not a complete material/support test.
#include "spin2014/QuadricStaticBalance.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool value,const char* why) {if(!value)throw std::runtime_error(why);}
template<class F> void rejected(F&& call) {
    bool seen=false;try {call();}catch(const std::invalid_argument&) {seen=true;}
    require(seen,"invalid static square-height state was accepted");
}
void freshSquaredHeight() {
    spin2014::QuadricObjectiveSettings objective;
    objective.objective=spin2014::Spin2014Objective::StaticLowering;
    objective.inertiaWeight=0;objective.centerWeight=1;
    Eigen::Matrix<double,10,1> s;
    // Mass 2, COM (0,0,1), positive central second moments. The compatible
    // retained halfspace is 1-z>0, not a height-minimizer pass flag.
    s<<2,0,0,2,1,1,3,0,0,0;
    const Eigen::Vector4d plane=Eigen::Vector4d(1,0,0,-1).normalized();
    const auto checked=spin2014::verifySpin2014StaticPlane(s,plane,1,objective);
    require(checked.squaredHeightConditionsVerified,"fresh squared-height KKT not checked");
    require(std::abs(checked.evaluation.objective-1)<1e-14,"height was not squared");
    require(std::abs(checked.evaluation.objectiveGradient[0]+1)<1e-14 &&
            std::abs(checked.evaluation.objectiveGradient[3]-1)<1e-14,
            "fresh objective derivatives do not match hand calculation");
    require(checked.multiplier>0 && checked.stationarityResidual<1e-14,
            "squared-height plane needs its own positive multiplier and residual");

    // Changing declared units changes both the squared objective and multiplier;
    // mass already includes density and must not enter a second time.
    objective.solverMassUnitInDeclaredUnits=2;
    objective.solverLengthUnitInDeclaredUnits=3;
    const auto scaled=spin2014::verifySpin2014StaticPlane(s,plane,1,objective);
    require(std::abs(scaled.evaluation.objective-9)<1e-13,"declared length not squared");
    require(std::abs(scaled.multiplier/checked.multiplier-9)<1e-12,
            "unit conversion failed on squared-height multiplier");

    const auto wrongPlane=spin2014::verifySpin2014StaticPlane(s,
        Eigen::Vector4d(.5,0,0,-1).normalized(),1,objective);
    require(!wrongPlane.squaredHeightConditionsVerified,
            "a candidate flag or horizontal feasibility cannot replace fresh stationarity");
    auto offaxis=s;offaxis[1]=.1;
    require(!spin2014::verifySpin2014StaticPlane(offaxis,plane,1,objective)
                .squaredHeightConditionsVerified,"horizontal constraint ignored");

    auto zero=s;zero[3]=0;
    const auto atSupport=spin2014::verifySpin2014StaticPlane(zero,plane,1,objective);
    require(atSupport.degenerateObjective && !atSupport.squaredHeightConditionsVerified,
            "zero height cannot borrow ordinary positive-multiplier optimality");
    objective.centerWeight=0;
    const auto noWeight=spin2014::verifySpin2014StaticPlane(s,plane,1,objective);
    require(noWeight.degenerateObjective && !noWeight.squaredHeightConditionsVerified,
            "zero weight is an explicit degenerate objective, not normal convergence");
    objective.centerWeight=1;
    auto bad=s;bad[0]=0;
    rejected([&] {spin2014::verifySpin2014StaticPlane(bad,plane,1,objective);});
    bad=s;bad[3]=-2;
    rejected([&] {spin2014::verifySpin2014StaticPlane(bad,plane,1,objective);});
}
spin::TriangleMesh box(const Eigen::Vector3d& low,const Eigen::Vector3d& high) {
    spin::TriangleMesh out;out.vertices.resize(8,3);out.faces.resize(12,3);
    const double a=low.x(),b=high.x(),c=low.y(),d=high.y(),e=low.z(),f=high.z();
    out.vertices<<a,c,e,b,c,e,b,d,e,a,d,e,a,c,f,b,c,f,b,d,f,a,d,f;
    out.faces<<0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;
    return out;
}
void realFullDomainIsRechecked() {
    spin::SpinDesignRequest request;
    request.mode=spin::SpinMode::StaticBalance;
    request.outerSurface=box(Eigen::Vector3d(0,0,0),Eigen::Vector3d(1,1,4));
    request.designDomain=box(Eigen::Vector3d(.1,.1,.1),Eigen::Vector3d(.9,.9,.9));
    request.targetFrame.origin=Eigen::Vector3d(.5,.5,0);
    request.materialDensity=3;
    spin2014::QuadricObjectiveSettings objective;
    objective.objective=spin2014::Spin2014Objective::StaticLowering;
    objective.inertiaWeight=0;
    const auto run=spin2014::solveSpin2014StaticQuadricHollowing(request,objective,2,
        Eigen::Vector4d(2,0,0,-1));
    require(run.candidate.evaluation.materialRegion==spin::PlaneMaterialRegion::StrictFull,
            "test must use the actual full-domain Boolean state");
    require(std::abs(run.candidate.evaluation.totalMomentsTarget.mass()-8.512)<1e-12,
            "actual shell/design densities were not applied once");
    const double expectedHeight=16.256/8.512;
    require(std::abs(run.squaredHeight.evaluation.objective-expectedHeight*expectedHeight)<1e-12,
            "full-domain candidate did not evaluate the original squared objective");
    require(run.squaredHeight.squaredHeightConditionsVerified,
            "full-domain squared-height optimum requires a freshly computed sign witness");
    require(run.squaredHeight.boundaryWitnessVerified && run.squaredHeight.boundaryMultiplier>0,
            "fresh full-domain witness/multiplier missing");
    require(std::abs(run.squaredHeight.boundaryGradientMinimum-
        2*expectedHeight*(expectedHeight-.9)/8.512)<1e-12,
        "full-domain witness uses the old height gradient instead of squared height");
    require((run.squaredHeight.boundaryWitnessCoefficients-run.candidate.coefficients).norm()>.1,
            "new witness replaced or borrowed the actual accepted plane");
    // This local gate is not a nontrivial optimized-case success.
    require(run.candidate.acceptedSteps==0,"fixture unexpectedly changed its initial full material");
    auto below=request;below.targetFrame.origin.z()=.01;
    rejected([&] {spin2014::solveSpin2014StaticQuadricHollowing(below,objective,2,
        Eigen::Vector4d(2,0,0,-1));});
    auto offaxis=request;offaxis.targetFrame.origin.x()=.6;
    const auto failed=spin2014::solveSpin2014StaticQuadricHollowing(offaxis,objective,2,
        Eigen::Vector4d(2,0,0,-1));
    require(!failed.squaredHeight.boundaryWitnessVerified &&
            !failed.squaredHeight.squaredHeightConditionsVerified,
            "full-domain sign witness ignored horizontal mass constraints");
}
}
int main() {
    try {freshSquaredHeight();realFullDomainIsRechecked();std::cout<<"[PASS] fresh static squared-height gate\n";return EXIT_SUCCESS;}
    catch(const std::exception& error) {std::cerr<<"[FAIL] "<<error.what()<<'\n';return EXIT_FAILURE;}
}
