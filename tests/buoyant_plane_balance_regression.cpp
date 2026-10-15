#include "spin/BuoyantPlaneBalanceSolver.h"
#include "quadric_certified_empty_cases.h"

#include <cmath>
#include <iostream>

namespace {
using certified_empty_tests::require;

spin::TriangleMesh box(const Eigen::Vector3d& lo, const Eigen::Vector3d& hi) {
    spin::TriangleMesh mesh;
    mesh.vertices.resize(8, 3);
    mesh.vertices << lo.x(),lo.y(),lo.z(), hi.x(),lo.y(),lo.z(),
        hi.x(),hi.y(),lo.z(), lo.x(),hi.y(),lo.z(),
        lo.x(),lo.y(),hi.z(), hi.x(),lo.y(),hi.z(),
        hi.x(),hi.y(),hi.z(), lo.x(),hi.y(),hi.z();
    mesh.faces.resize(12, 3);
    mesh.faces << 0,2,1, 0,3,2, 4,5,6, 4,6,7,
        0,1,5, 0,5,4, 3,7,6, 3,6,2,
        0,4,7, 0,7,3, 1,2,6, 1,6,5;
    return mesh;
}

spin::SpinDesignRequest request() {
    auto result = certified_empty_tests::request();
    result.designDomain = box(Eigen::Vector3d(-1,-1,-1), Eigen::Vector3d(2,1,1));
    result.outerSurface = box(Eigen::Vector3d(-2,-2,-2), Eigen::Vector3d(3,2,2));
    result.fixedShellMoments = physics::RawMassMoments10();
    result.materialDensity = 1.0;
    result.minimumWallThickness = 0.0;
    result.mode = spin::SpinMode::BuoyantStability;
    result.buoyancy.liquidDensity = 1.0;
    result.buoyancy.displacedLiquidVolume = 6.0;
    result.buoyancy.centerOfBuoyancyTarget = Eigen::Vector3d(0.5, 0.0, -0.5);
    result.buoyancy.waterlineArea = 6.0;
    result.buoyancy.waterlineCentralSecondMoment =
        (Eigen::Matrix2d() << 4.5, 0.0, 0.0, 2.0).finished();
    result.buoyancy.fixedShellVolume = 0.0;
    return result;
}
}

int main() {
    try {
        const auto input = request();
        spin::BuoyantPlaneBalanceSettings settings;
        settings.maximumIterations = 40;
        settings.constraintTolerance = 1e-11;
        const auto result = spin::solveBuoyantStabilityByPlane(
            input, Eigen::Vector3d(0.1, 0.2, -0.1), settings);
        require(result.status == spin::BuoyantPlaneBalanceStatus::Converged,
            "buoyant plane solve did not converge from a non-solution start");
        require(result.equilibriumConstraintsVerified,
            "converged buoyant plane did not meet all three constraints with a nonsingular reduced system");
        require(result.iterations > 0 && result.acceptedSteps > 0,
            "nontrivial buoyant solve was misreported as an initial solution");
        require(std::abs(result.planeCoefficients.x()) < 1e-9 &&
                std::abs(result.planeCoefficients.y()) < 1e-9 &&
                std::abs(result.planeCoefficients.z()) < 1e-9 &&
                std::abs(result.planeCoefficients.w() + 1.0) < 1e-14,
            "analytic half-box buoyant interface was not recovered");
        require(std::abs(result.totalMoments.values()[0] - 6.0) < 1e-9 &&
                std::abs(result.totalMoments.values()[1] - 3.0) < 1e-9 &&
                std::abs(result.totalMoments.values()[2]) < 1e-9,
            "final total moments do not satisfy mass and horizontal buoyancy balance");
        require(result.normalizedConstraints.norm() <= settings.constraintTolerance,
            "fixed-input normalized buoyant constraints exceed tolerance");
        require(result.equalityMultipliers.size() == 3 &&
                std::abs(result.equalityMultipliers(0)) < 1e-12 &&
                std::abs(result.equalityMultipliers(1)) < 1e-12 &&
                std::abs(result.equalityMultipliers(2) - 1.0 / 12.0) < 1e-10 &&
                std::abs(result.multiplierMu - 1.0 / 6.0) < 1e-10,
            "reduced plane coefficients were not converted to the published Eq. 7 multipliers");
        require(result.theoremOneNumericalAudit.has_value() &&
                result.theoremOneNumericalAudit->status ==
                    spin::FasterTheoremConditionStatus::SatisfiedNumerically &&
                result.eq7NumericalConditionsVerified,
            "converged plane did not pass the independent finite-precision Eq. 7 audit");
        require(result.theoremOneNumericalAudit->unitNormError <= 1e-12 &&
                result.theoremOneNumericalAudit->constraintResidual <= 1e-10 &&
                result.theoremOneNumericalAudit->stationarityResidual <= 1e-12 &&
                result.theoremOneNumericalAudit->proofScope ==
                    "NUMERICAL_EQ7_AND_EQUALITY_CONDITIONS_ONLY_NOT_A_THEOREM_PROOF",
            "Eq. 7 audit did not expose bounded residuals with its restricted proof scope");
        require(result.stabilityMargin > 0.0 && result.stableAtPrescribedWaterline,
            "independent metacenter test did not classify the analytic body as stable");
        require(std::abs(result.metacenterHeight + 1.0 / 6.0) < 1e-12,
            "buoyant solver metacenter differs from the independently specified waterline");

        auto tinyInput = input;
        tinyInput.outerSurface.vertices *= 1e-3;
        tinyInput.designDomain.vertices *= 1e-3;
        tinyInput.buoyancy.displacedLiquidVolume *= 1e-9;
        tinyInput.buoyancy.centerOfBuoyancyTarget *= 1e-3;
        tinyInput.buoyancy.waterlineArea *= 1e-6;
        tinyInput.buoyancy.waterlineCentralSecondMoment *= 1e-12;
        const auto tiny = spin::solveBuoyantStabilityByPlane(
            tinyInput, Eigen::Vector3d(1e-4, 0.2, -0.1), settings);
        require(tiny.status == spin::BuoyantPlaneBalanceStatus::Converged &&
                tiny.equilibriumConstraintsVerified &&
                tiny.eq7NumericalConditionsVerified &&
                tiny.theoremOneNumericalAudit.has_value() &&
                tiny.theoremOneNumericalAudit->status ==
                    spin::FasterTheoremConditionStatus::SatisfiedNumerically &&
                tiny.normalizedConstraints.norm() <= settings.constraintTolerance &&
                std::abs(tiny.planeCoefficients.x()) < 1e-12,
            "fixed-input normalized buoyant solve and Eq. 7 audit failed under uniform unit scaling");

        auto unstableInput = input;
        unstableInput.buoyancy.centerOfBuoyancyTarget.z() = -1.0;
        const auto unstable = spin::solveBuoyantStabilityByPlane(
            unstableInput, Eigen::Vector3d::Zero(), settings);
        require(unstable.status == spin::BuoyantPlaneBalanceStatus::Converged &&
                unstable.equilibriumConstraintsVerified &&
                unstable.eq7NumericalConditionsVerified &&
                unstable.theoremOneNumericalAudit.has_value() &&
                unstable.theoremOneNumericalAudit->status ==
                    spin::FasterTheoremConditionStatus::SatisfiedNumerically &&
                !unstable.stableAtPrescribedWaterline && unstable.stabilityMargin < 0.0,
            "unstable metacenter comparison was conflated with optimizer failure");

        auto capped = settings;
        capped.maximumIterations = 0;
        const auto stopped = spin::solveBuoyantStabilityByPlane(
            input, Eigen::Vector3d(0.1, 0.2, -0.1), capped);
        require(stopped.status == spin::BuoyantPlaneBalanceStatus::IterationLimit &&
                !stopped.equilibriumConstraintsVerified &&
                !stopped.eq7NumericalConditionsVerified &&
                !stopped.theoremOneNumericalAudit.has_value() &&
                !stopped.stableAtPrescribedWaterline,
            "iteration limit was promoted to buoyant feasibility or stability");
        std::cout << "[PASS] prescribed-waterline plane solve and finite-precision Eq. 7 audit passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
