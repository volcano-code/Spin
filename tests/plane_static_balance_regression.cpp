#include "spin/PlaneStaticBalanceSolver.h"
#include "quadric_certified_empty_cases.h"
#include <cmath>
#include <iomanip>
#include <iostream>

namespace {
using certified_empty_tests::require;

spin::TriangleMesh box(const Eigen::Vector3d& lo,const Eigen::Vector3d& hi) {
    spin::TriangleMesh m;
    m.vertices.resize(8,3);
    m.vertices<<lo.x(),lo.y(),lo.z(), hi.x(),lo.y(),lo.z(),
        hi.x(),hi.y(),lo.z(), lo.x(),hi.y(),lo.z(),
        lo.x(),lo.y(),hi.z(), hi.x(),lo.y(),hi.z(),
        hi.x(),hi.y(),hi.z(), lo.x(),hi.y(),hi.z();
    m.faces.resize(12,3);
    m.faces<<0,2,1,0,3,2,4,5,6,4,6,7,0,1,5,0,5,4,
        3,7,6,3,6,2,0,4,7,0,7,3,1,2,6,1,6,5;
    return m;
}

spin::SpinDesignRequest offCenterDomain() {
    auto r=certified_empty_tests::request();
    r.designDomain=box(Eigen::Vector3d(-1,-1,-1),Eigen::Vector3d(2,1,1));
    r.outerSurface=box(Eigen::Vector3d(-2,-2,-2),Eigen::Vector3d(3,2,2));
    r.fixedShellMoments=physics::RawMassMoments10();
    r.materialDensity=1.0;
    r.minimumWallThickness=0.0;
    return r;
}

spin::SpinDesignRequest centeredShellDomain() {
    auto r=certified_empty_tests::request();
    r.designDomain=box(Eigen::Vector3d(-1,-1,-1),Eigen::Vector3d(1,1,1));
    r.outerSurface=box(Eigen::Vector3d(-2,-2,-2),Eigen::Vector3d(2,2,2));
    r.materialDensity=1.0;
    r.minimumWallThickness=0.0;
    MassMoments shell;
    shell.mass=56.0; // [-2,2]^3 minus [-1,1]^3, unit density.
    shell.firstMoment.setZero();
    shell.secondMoment=Eigen::Matrix3d::Identity()*(248.0/3.0);
    r.fixedShellMoments=physics::RawMassMoments10::fromMassMoments(shell);
    return r;
}

void lowerAndParetoPlaneFamilies() {
    auto request=centeredShellDomain();
    request.staticObjective=spin::StaticObjective::LowerCenterOfMass;
    spin::PlaneStaticBalanceSettings settings;
    settings.maximumIterations=40;
    settings.constraintTolerance=2e-10;
    const double expected=(-60.0+std::sqrt(3584.0))/4.0;
    const auto lower=spin::solveStaticBalanceByPlaneFamily(
        request,Eigen::Vector4d(0.0,0.0,0.0,-1.0),settings);
    require(lower.status==spin::PlaneStaticBalanceStatus::Converged,
        "lower-center-of-mass plane-family solve did not converge");
    require(lower.paperOptimalityConditionsVerified,
        "lower-center-of-mass solution lacks the paper KKT check");
    require(std::abs(-lower.coefficients[0]/lower.coefficients[3]*4.0-expected)<1e-8,
        "lower-center-of-mass plane does not match the analytic box solution");
    require(lower.normalizedConstraints.norm()<=settings.constraintTolerance,
        "lower-center-of-mass horizontal balance constraints exceed tolerance");
    require(lower.acceptedSteps>0,
        "lower-center-of-mass solve was mistaken for an initial solution");
    const auto lowerMeshMoments=physics::RawMassMoments10::fromMassMoments(
        computeSurfaceMassMoments(lower.evaluation.variableMaterialBoundaryTarget.vertices,
                                  lower.evaluation.variableMaterialBoundaryTarget.faces,
                                  request.materialDensity));
    require((lowerMeshMoments.values()-lower.evaluation.variableMomentsTarget.values()).norm()<1e-10,
        "lower-center-of-mass output mesh failed independent ten-moment integration");

    auto massRequest=offCenterDomain();
    massRequest.staticObjective=spin::StaticObjective::ParetoTradeoff;
    massRequest.staticParetoWeight=1.0;
    const auto massEndpoint=spin::solveStaticBalanceByPlaneFamily(
        massRequest,Eigen::Vector4d(1.0,-5.0,0.0,0.0),settings);
    require(massEndpoint.status==spin::PlaneStaticBalanceStatus::Converged &&
            massEndpoint.paperOptimalityConditionsVerified,
        "Pareto maximum-mass endpoint did not use its exact vertical-plane limit");
    require(std::abs(massEndpoint.coefficients[3])<1e-10,
        "Pareto maximum-mass endpoint did not recover a vertical plane");

    request.staticObjective=spin::StaticObjective::ParetoTradeoff;
    request.staticParetoWeight=0.0;
    const auto lowerEndpoint=spin::solveStaticBalanceByPlaneFamily(
        request,Eigen::Vector4d(0.0,0.0,0.0,-1.0),settings);
    require(lowerEndpoint.status==spin::PlaneStaticBalanceStatus::Converged &&
            lowerEndpoint.paperOptimalityConditionsVerified,
        "Pareto lowering endpoint did not match the dedicated lowering objective");
    require((lowerEndpoint.coefficients-lower.coefficients).norm()<1e-9,
        "Pareto zero-weight endpoint differs from lower-center-of-mass solve");

    request.staticParetoWeight=0.01;
    const auto pareto=spin::solveStaticBalanceByPlaneFamily(
        request,lower.coefficients,settings);
    require(pareto.status==spin::PlaneStaticBalanceStatus::Converged &&
            pareto.paperOptimalityConditionsVerified,
        "interior static Pareto plane did not satisfy its KKT system");
    require(pareto.normalizedConstraints.norm()<=settings.constraintTolerance &&
            pareto.stationarityResidual<=settings.constraintTolerance,
        "interior Pareto endpoint failed independent balance or stationarity checks");

    spin::PlaneStaticBalanceSettings capped=settings;capped.maximumIterations=0;
    const auto limited=spin::solveStaticBalanceByPlaneFamily(
        request,lower.coefficients,capped);
    require(limited.status==spin::PlaneStaticBalanceStatus::IterationLimit &&
            !limited.paperOptimalityConditionsVerified,
        "plane-family iteration limit was promoted to KKT success");
    bool wrongOrientationRejected=false;
    try {
        spin::solveStaticBalanceByPlaneFamily(
            request,Eigen::Vector4d(0.0,0.0,0.0,1.0),settings);
    } catch(const std::invalid_argument&) { wrongOrientationRejected=true; }
    require(wrongOrientationRejected,"wrong below-plane orientation was silently accepted");
    std::cout<<"STATIC_FAMILY_RESULT lower="<<lower.coefficients.transpose()
        <<" lower-mass="<<lower.evaluation.variableMomentsTarget.mass()
        <<" pareto="<<pareto.coefficients.transpose()
        <<" constraints="<<pareto.normalizedConstraints.transpose()
        <<" stationarity="<<pareto.stationarityResidual<<"\n";
}
}

int main() {
    try {
        std::cout<<std::setprecision(17);
        const auto request=offCenterDomain();
        spin::PlaneStaticBalanceSettings settings;
        settings.maximumIterations=40;
        settings.constraintTolerance=1e-12;
        const auto result=spin::solveStaticBalanceByPlane(
            request,Eigen::Vector2d(-2.0,0.0),settings);
        require(result.status==spin::PlaneStaticBalanceStatus::Converged,
            "static plane solve did not converge");
        require(result.paperOptimalityConditionsVerified,
            "feasible plane threshold did not satisfy paper optimality conditions");
        require(std::abs(result.lambda.x()+1.0)<1e-8 && std::abs(result.lambda.y())<1e-8,
            "plane multipliers differ from analytic box solution");
        require(std::abs(result.variableMoments.mass()-8.0)<1e-9,
            "maximal filled halfspace has incorrect analytic mass");
        require(result.normalizedConstraints.norm()<=settings.constraintTolerance,
            "fixed-input normalized static constraints exceed tolerance");
        require(result.constraints.norm()<=settings.constraintTolerance*result.constraintScale,
            "raw constraints disagree with their fixed-input normalized gate");
        require(result.iterations>0 && result.acceptedSteps>0,
            "nontrivial plane optimization was mistaken for an initial solution");
        require(result.evaluation.coefficientJacobian.has_value(),
            "terminal plane differential unavailable");
        int acceptedTraceRows=0;
        for(const auto& row:result.trace)acceptedTraceRows+=row.accepted;
        require(result.trace.size()>static_cast<std::size_t>(result.iterations+1) &&
                acceptedTraceRows==result.acceptedSteps,
            "plane solve trace does not retain rejected and accepted line-search evaluations");
        for(const auto& row:result.trace)
            require(std::isfinite(row.residualNorm) && std::isfinite(row.variableMass),
                "plane solver trace contains nonfinite numbers");

        spin::PlaneStaticBalanceSettings capped=settings;
        capped.maximumIterations=0;
        const auto noRun=spin::solveStaticBalanceByPlane(
            request,Eigen::Vector2d(-2.0,0.0),capped);
        require(noRun.status==spin::PlaneStaticBalanceStatus::IterationLimit &&
                !noRun.paperOptimalityConditionsVerified,
            "iteration-limit result was promoted to a certified optimum");

        auto impossible=offCenterDomain();
        impossible.targetFrame.origin.x()=4.0;
        const auto notSolved=spin::solveStaticBalanceByPlane(
            impossible,Eigen::Vector2d(-2.0,0.0),settings);
        require(notSolved.status!=spin::PlaneStaticBalanceStatus::Converged &&
                !notSolved.paperOptimalityConditionsVerified,
            "nonconvergence was mislabeled as a verified optimum");

        auto tiny=offCenterDomain();
        tiny.designDomain.vertices*=1e-3;
        tiny.outerSurface.vertices*=1e-3;
        const auto scaled=spin::solveStaticBalanceByPlane(
            tiny,Eigen::Vector2d(-2000.0,0.0),settings);
        require(scaled.status==spin::PlaneStaticBalanceStatus::Converged &&
                scaled.paperOptimalityConditionsVerified,
            "uniformly scaled static problem failed normalized convergence");
        require(std::abs(scaled.lambda.x()+1000.0)<1e-5 &&
                std::abs(scaled.lambda.y())<1e-5,
            "small-unit input was incorrectly accepted before reaching its balancing plane");
        require(scaled.normalizedConstraints.norm()<=settings.constraintTolerance,
            "scaled problem did not meet the normalized constraint gate");
        lowerAndParetoPlaneFamilies();
        std::cout<<"[PASS] static plane balance restricted optimization and fail-closed status\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<"[FAIL] "<<error.what()<<'\n';
        return 1;
    }
}
