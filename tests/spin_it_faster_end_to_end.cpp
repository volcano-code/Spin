#include "spin/QuadricHollowingSolver.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>


namespace
{


using Quadric = spin::QuadricCoefficients10;
using Vector10 = Quadric::Vector;


spin::TriangleMesh cube(double minimum, double maximum)
{
    spin::TriangleMesh mesh;
    mesh.vertices.resize(8, 3);
    mesh.vertices <<
        minimum, minimum, minimum,
        maximum, minimum, minimum,
        maximum, maximum, minimum,
        minimum, maximum, minimum,
        minimum, minimum, maximum,
        maximum, minimum, maximum,
        maximum, maximum, maximum,
        minimum, maximum, maximum;
    mesh.faces.resize(12, 3);
    mesh.faces <<
        0, 2, 1, 0, 3, 2,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5;
    return mesh;
}


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=spin-it-faster-end-to-end "
              << message << '\n';
    return EXIT_FAILURE;
}


spin::SpinDesignRequest baseRequest(spin::SpinMode mode)
{
    spin::SpinDesignRequest request;
    request.outerSurface = cube(-5.0, 5.0);
    request.designDomain = cube(-1.0, 1.0);
    request.minimumWallThickness = 4.0;
    request.mode = mode;
    request.settings.maximumRestarts = 1;
    return request;
}


physics::RawMassMoments10 rawMoments(
    double mass,
    const Eigen::Vector3d& first,
    const Eigen::Vector3d& diagonalSecond
)
{
    MassMoments moments;
    moments.mass = mass;
    moments.firstMoment = first;
    moments.secondMoment = diagonalSecond.asDiagonal();
    return physics::RawMassMoments10::fromMassMoments(moments);
}


spin::SpinDesignRequest requestWithTarget(
    spin::SpinDesignRequest request,
    const Quadric& quadric,
    const physics::RawMassMoments10& target
)
{
    const spin::QuadricMomentEvaluation variable =
        spin::evaluateQuadricMaterialMoments(request, quadric);
    request.fixedShellMoments = target - variable.variableMomentsTarget;
    spin::validateSpinDesignRequest(request);
    return request;
}


bool solveAndCheck(
    const char* name,
    const spin::SpinDesignRequest& request,
    const Quadric& expected,
    spin::QuadricSurfaceType expectedType,
    std::string* outputHash
)
{
    const spin::QuadricMomentEvaluation expectedEvaluation =
        spin::evaluateQuadricMaterialMoments(request, expected);
    if (expectedEvaluation.quadricTriangulation.type != expectedType)
    {
        std::cerr << "unexpected interface type for " << name << '\n';
        return false;
    }
    const spin::QuadricHollowingResult result = spin::solveQuadricHollowing(
        request,
        spin::QuadricMomentEvaluatorSettings(),
        {expected}
    );
    const spin::QuadricHollowingResult repeated = spin::solveQuadricHollowing(
        request,
        spin::QuadricMomentEvaluatorSettings(),
        {expected}
    );
    if (
        result.status != spin::SpinSolverStatus::Converged ||
        result.normalizedEquationResidual > 1e-12 ||
        result.physicalRelativeResidual > 1e-6 ||
        result.materialBoundary.faces.rows() == 0 ||
        result.cavityBoundary.faces.rows() == 0 ||
        result.geometryHash.size() != 64 || result.traceHash.size() != 64 ||
        result.outputHash.size() != 64 ||
        result.outputHash != repeated.outputHash ||
        !result.quadricCoefficients.isApprox(expected.values(), 1e-13) ||
        result.outerSurface.vertices != request.outerSurface.vertices ||
        result.outerSurface.faces != request.outerSurface.faces
    )
    {
        std::cerr << name << " status="
                  << spin::spinSolverStatusToken(result.status)
                  << " residual=" << result.normalizedEquationResidual
                  << " physical=" << result.physicalRelativeResidual
                  << " faces=" << result.materialBoundary.faces.rows()
                  << " cavity_faces=" << result.cavityBoundary.faces.rows()
                  << " geometry_hash=" << result.geometryHash.size()
                  << " output_hash=" << result.outputHash.size()
                  << " repeat_equal="
                  << (result.outputHash == repeated.outputHash)
                  << " coefficient_error="
                  << (result.quadricCoefficients - expected.values()).norm()
                  << " outer_v_equal="
                  << (result.outerSurface.vertices == request.outerSurface.vertices)
                  << " outer_f_equal="
                  << (result.outerSurface.faces == request.outerSurface.faces)
                  << " diagnostic=" << result.diagnostic << '\n';
        for (const auto& event : result.events)
        {
            if (event.detail.find("BOOLEAN_ARRANGEMENT_PWN_FAILURE") ==
                std::string::npos)
                continue;
            std::cerr << std::setprecision(17)
                      << "[DIAGNOSTIC] failed_boolean_evaluation="
                      << event.evaluationIndex
                      << " restart_evaluation=" << event.restartEvaluationIndex
                      << " step=" << event.step.value_or(-1.0)
                      << " q=";
            if (event.quadricCoefficients)
            {
                for (double value : *event.quadricCoefficients)
                    std::cerr << value << ',';
            }
            else
            {
                std::cerr << "missing";
            }
            std::cerr << " detail=" << event.detail << '\n';
        }
        return false;
    }
    if (request.mode == spin::SpinMode::BuoyantStability)
    {
        if (
            !std::isfinite(result.buoyantMetacenterHeight) ||
            !std::isfinite(result.buoyantStabilityMargin) ||
            result.buoyantStabilityMargin <= 0.0
        )
        {
            std::cerr << name << " invalid buoyancy audit: metacenter="
                      << result.buoyantMetacenterHeight << " margin="
                      << result.buoyantStabilityMargin << '\n';
            return false;
        }
    }
    *outputHash = result.outputHash;
    return true;
}


} // namespace


int main()
{
    std::vector<std::string> hashes;

    Vector10 values = Vector10::Zero();
    values(0) = 0.371;
    values(1) = 0.8;
    values(2) = 0.123;
    const Quadric staticPlane = Quadric::fromUnnormalized(values);
    spin::SpinDesignRequest staticRequest =
        baseRequest(spin::SpinMode::StaticBalance);
    const spin::QuadricMomentEvaluation staticVariable =
        spin::evaluateQuadricMaterialMoments(staticRequest, staticPlane);
    staticRequest.fixedShellMoments = rawMoments(
        10.0,
        -staticVariable.variableMomentsTarget.toMassMoments().firstMoment,
        Eigen::Vector3d::Constant(20.0)
    );
    hashes.emplace_back();
    if (!solveAndCheck(
            "static",
            staticRequest,
            staticPlane,
            spin::QuadricSurfaceType::Plane,
            &hashes.back()
        ))
    {
        return fail("static planar fixture failed");
    }

    // Exercise the full design-domain -> Boolean -> mass-moment -> Newton
    // path from a deterministic seed that is not the known feasible solution.
    // This is a synthetic algorithm-integration check, not a paper case.
    Vector10 offTargetValues = staticPlane.values();
    offTargetValues(0) += 0.025;
    offTargetValues(1) -= 0.02;
    const Quadric offTargetSeed =
        Quadric::fromUnnormalized(offTargetValues);
    if (offTargetSeed.hash() == staticPlane.hash())
    {
        return fail("off-target seed unexpectedly equals known solution");
    }
    const spin::QuadricHollowingResult offTargetResult =
        spin::solveQuadricHollowing(
            staticRequest,
            spin::QuadricMomentEvaluatorSettings(),
            {offTargetSeed}
        );
    const spin::QuadricHollowingResult offTargetRepeat =
        spin::solveQuadricHollowing(
            staticRequest,
            spin::QuadricMomentEvaluatorSettings(),
            {offTargetSeed}
        );
    if (
        offTargetResult.status != spin::SpinSolverStatus::Converged ||
        offTargetResult.normalizedEquationResidual > 1e-12 ||
        offTargetResult.physicalRelativeResidual > 1e-6 ||
        offTargetResult.trace.size() < 2 ||
        (offTargetResult.quadricCoefficients - offTargetSeed.values()).norm() <
            1e-8 ||
        offTargetResult.outputHash != offTargetRepeat.outputHash ||
        offTargetResult.traceHash != offTargetRepeat.traceHash ||
        offTargetResult.materialBoundary.faces.rows() == 0 ||
        offTargetResult.cavityBoundary.faces.rows() == 0
    )
    {
        return fail(
            "off-target integrated Newton solve failed: status=" +
            std::string(spin::spinSolverStatusToken(offTargetResult.status)) +
            " residual=" +
            std::to_string(offTargetResult.normalizedEquationResidual) +
            " physical=" +
            std::to_string(offTargetResult.physicalRelativeResidual) +
            " trace=" + std::to_string(offTargetResult.trace.size())
        );
    }

    values.setZero();
    values(0) = 0.2;
    values(1) = 0.12;
    values(2) = -0.08;
    values(3) = -1.0;
    const Quadric staticLowPlane = Quadric::fromUnnormalized(values);
    spin::SpinDesignRequest staticLowRequest =
        baseRequest(spin::SpinMode::StaticBalance);
    staticLowRequest.staticObjective =
        spin::StaticObjective::LowerCenterOfMass;
    const spin::QuadricMomentEvaluation staticLowVariable =
        spin::evaluateQuadricMaterialMoments(
            staticLowRequest,
            staticLowPlane
        );
    staticLowRequest.fixedShellMoments = rawMoments(
        10.0 - staticLowVariable.variableMomentsTarget.mass(),
        Eigen::Vector3d(0.0, 0.0, 2.0) -
            staticLowVariable.variableMomentsTarget.toMassMoments().firstMoment,
        Eigen::Vector3d::Constant(20.0)
    );
    hashes.emplace_back();
    if (!solveAndCheck(
            "static-lower-center",
            staticLowRequest,
            staticLowPlane,
            spin::QuadricSurfaceType::Plane,
            &hashes.back()
        ))
    {
        return fail("static lower-center planar fixture failed");
    }

    values.setZero();
    values(0) = 0.2;
    values(1) = 0.2;
    values(2) = -0.1;
    values(3) = -1.0;
    const Quadric buoyantPlane = Quadric::fromUnnormalized(values);
    spin::SpinDesignRequest buoyantRequest =
        baseRequest(spin::SpinMode::BuoyantStability);
    buoyantRequest.buoyancy.liquidDensity = 1.0;
    buoyantRequest.buoyancy.displacedLiquidVolume = 20.0;
    buoyantRequest.buoyancy.centerOfBuoyancyTarget =
        Eigen::Vector3d(0.02, -0.03, 0.0);
    buoyantRequest.buoyancy.waterlineArea = 4.0;
    buoyantRequest.buoyancy.waterlineCentralSecondMoment =
        2.0 * Eigen::Matrix2d::Identity();
    const spin::QuadricMomentEvaluation buoyantVariable =
        spin::evaluateQuadricMaterialMoments(
            buoyantRequest,
            buoyantPlane
        );
    const physics::RawMassMoments10 buoyantTarget = rawMoments(
        20.0,
        Eigen::Vector3d(0.4, -0.6, 0.0),
        Eigen::Vector3d::Constant(100.0)
    );
    buoyantRequest.fixedShellMoments =
        buoyantTarget - buoyantVariable.variableMomentsTarget;
    buoyantRequest.buoyancy.fixedShellVolume =
        buoyantRequest.fixedShellMoments.mass();
    hashes.emplace_back();
    if (!solveAndCheck(
            "buoyant",
            buoyantRequest,
            buoyantPlane,
            spin::QuadricSurfaceType::Plane,
            &hashes.back()
        ))
    {
        return fail("buoyant planar fixture failed");
    }

    values.setZero();
    values(1) = 0.8;
    values(4) = 1.0;
    values(5) = 1.0;
    const Quadric fidgetCircle = Quadric::fromUnnormalized(values);
    spin::SpinDesignRequest fidgetRequest =
        baseRequest(spin::SpinMode::FidgetSpinner);
    const spin::QuadricMomentEvaluation fidgetVariable =
        spin::evaluateQuadricMaterialMoments(fidgetRequest, fidgetCircle);
    fidgetRequest.fixedShellMoments = rawMoments(
        10.0,
        -fidgetVariable.variableMomentsTarget.toMassMoments().firstMoment,
        Eigen::Vector3d::Constant(20.0)
    );
    hashes.emplace_back();
    if (!solveAndCheck(
            "fidget",
            fidgetRequest,
            fidgetCircle,
            spin::QuadricSurfaceType::EllipticCylinder,
            &hashes.back()
        ))
    {
        return fail("fidget circular-void fixture failed");
    }

    values.setZero();
    values(3) = 1.0;
    values(4) = 1.0;
    values(5) = 1.0;
    values(6) = -1.0;
    const Quadric yoyoHyperboloid = Quadric::fromUnnormalized(values);
    spin::SpinDesignRequest yoyoRequest = baseRequest(spin::SpinMode::YoYo);
    yoyoRequest = requestWithTarget(
        yoyoRequest,
        yoyoHyperboloid,
        rawMoments(
            20.0,
            Eigen::Vector3d::Zero(),
            Eigen::Vector3d(100.0, 100.0, 200.0)
        )
    );
    Vector10 offTargetYoyoValues = yoyoHyperboloid.values();
    offTargetYoyoValues(0) += 0.015;
    offTargetYoyoValues(6) -= 0.025;
    const Quadric offTargetYoyoSeed =
        Quadric::fromUnnormalized(offTargetYoyoValues);
    const spin::QuadricHollowingResult offTargetYoyo =
        spin::solveQuadricHollowing(
            yoyoRequest,
            spin::QuadricMomentEvaluatorSettings(),
            {offTargetYoyoSeed}
        );
    const spin::QuadricMomentEvaluation baseEvaluation =
        spin::evaluateQuadricMaterialMoments(yoyoRequest, offTargetYoyoSeed);
    const std::vector<double> derivativeSteps = {1e-4, 1e-5, 1e-6};
    std::vector<double> derivativeErrors;
    std::vector<int> topologyChangeCounts;
    for (const double derivativeStep : derivativeSteps)
    {
        Eigen::Matrix<double, 10, 10> finiteDifferenceJacobian;
        Eigen::Matrix<double, 10, 10> projectedAnalyticJacobian;
        int topologyChanges = 0;
        for (int coordinate = 0; coordinate < 10; ++coordinate)
        {
            Eigen::Matrix<double, 10, 1> basis =
                Eigen::Matrix<double, 10, 1>::Zero();
            basis(coordinate) = 1.0;
            const Eigen::Matrix<double, 10, 1> direction =
                offTargetYoyoSeed.projectToTangent(basis);
            const auto plus = spin::evaluateQuadricMaterialMoments(
                yoyoRequest,
                offTargetYoyoSeed.geodesicUpdate(direction, derivativeStep)
            );
            const auto minus = spin::evaluateQuadricMaterialMoments(
                yoyoRequest,
                offTargetYoyoSeed.geodesicUpdate(-direction, derivativeStep)
            );
            finiteDifferenceJacobian.col(coordinate) =
                (plus.variableMomentsTarget.values() -
                 minus.variableMomentsTarget.values()) /
                (2.0 * derivativeStep);
            projectedAnalyticJacobian.col(coordinate) =
                baseEvaluation.momentJacobian * direction;
            topologyChanges += plus.topologyHash != minus.topologyHash ? 1 : 0;
        }
        derivativeErrors.push_back(
            (projectedAnalyticJacobian - finiteDifferenceJacobian).norm() /
            std::max(1.0, finiteDifferenceJacobian.norm())
        );
        topologyChangeCounts.push_back(topologyChanges);
    }
    const bool derivativeMismatch = std::any_of(
        derivativeErrors.begin(), derivativeErrors.end(),
        [](double error) { return error > 0.2; }
    );
    const bool topologyTransitionObserved = std::any_of(
        topologyChangeCounts.begin(), topologyChangeCounts.end(),
        [](int count) { return count > 0; }
    );
    if (
        offTargetYoyo.status != spin::SpinSolverStatus::LineSearchExhausted ||
        offTargetYoyo.trace.size() != 1 || !derivativeMismatch ||
        !topologyTransitionObserved ||
        baseEvaluation.derivativeDiagnostic.find(
            "not the discrete geometry derivative"
        ) == std::string::npos
    )
    {
        return fail(
            "hyperboloid derivative-mismatch negative control changed: status=" +
            std::string(spin::spinSolverStatusToken(offTargetYoyo.status)) +
            " derivative_mismatch=" + std::to_string(derivativeMismatch) +
            " topology_transition=" +
            std::to_string(topologyTransitionObserved)
        );
    }
    std::cerr << "[BLOCKER] scenario=synthetic-yoyo-hyperboloid-off-target"
              << " status=" << spin::spinSolverStatusToken(offTargetYoyo.status)
              << " residual=" << offTargetYoyo.normalizedEquationResidual
              << " derivative_error_1e4=" << derivativeErrors[0]
              << " derivative_error_1e5=" << derivativeErrors[1]
              << " derivative_error_1e6=" << derivativeErrors[2]
              << " topology_changed_directions_1e4=" << topologyChangeCounts[0]
              << " topology_changed_directions_1e5=" << topologyChangeCounts[1]
              << " topology_changed_directions_1e6=" << topologyChangeCounts[2]
              << " classification=LEGACY_CONTINUOUS_DERIVATIVE_NOT_DISCRETE_BOOLEAN_DERIVATIVE\n";

    spin::QuadricMomentEvaluatorSettings discreteSettings;
    discreteSettings.derivativeMethod =
        spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
    const spin::QuadricHollowingResult discreteOffTargetYoyo =
        spin::solveQuadricHollowing(
            yoyoRequest, discreteSettings, {offTargetYoyoSeed});
    if (discreteOffTargetYoyo.status != spin::SpinSolverStatus::ResourceLimitReached ||
        discreteOffTargetYoyo.diagnostic.find("degenerate_faces=2") == std::string::npos ||
        discreteOffTargetYoyo.materialBoundary.faces.rows() != 0)
    {
        return fail("tangent-contact negative control did not fail closed: status=" +
            std::string(spin::spinSolverStatusToken(discreteOffTargetYoyo.status)) +
            " diagnostic=" + discreteOffTargetYoyo.diagnostic);
    }
    std::cout << "[PASS] discrete Boolean-to-Newton tangency failure remains fail-closed status="
              << spin::spinSolverStatusToken(discreteOffTargetYoyo.status)
              << " iterations=" << discreteOffTargetYoyo.trace.size() << '\n';

    spin::SpinDesignRequest nonTangentYoyoRequest =
        baseRequest(spin::SpinMode::YoYo);
    nonTangentYoyoRequest.designDomain = cube(-0.8, 0.8);
    nonTangentYoyoRequest = requestWithTarget(
        nonTangentYoyoRequest,
        yoyoHyperboloid,
        rawMoments(20.0, Eigen::Vector3d::Zero(),
                   Eigen::Vector3d(100.0, 100.0, 200.0)));
    const auto discreteNonTangentYoyo = spin::solveQuadricHollowing(
        nonTangentYoyoRequest, discreteSettings, {offTargetYoyoSeed});
    const auto discreteNonTangentYoyoRepeat = spin::solveQuadricHollowing(
        nonTangentYoyoRequest, discreteSettings, {offTargetYoyoSeed});
    bool everyEvaluatedMapIsBoundToItsDifferential = true;
    int evaluatedMapCount = 0;
    for (const auto& event : discreteNonTangentYoyo.events)
    {
        if (event.type != "evaluation" || event.outcome == "failed")
            continue;
        ++evaluatedMapCount;
        everyEvaluatedMapIsBoundToItsDifferential =
            everyEvaluatedMapIsBoundToItsDifferential &&
            event.detail.find("evaluation_map=") != std::string::npos &&
            event.detail.find("backend=") != std::string::npos &&
            event.detail.find("derivative=") != std::string::npos;
    }
    if (discreteNonTangentYoyo.status != spin::SpinSolverStatus::Converged ||
        discreteNonTangentYoyo.normalizedEquationResidual > 1e-12 ||
        discreteNonTangentYoyo.physicalRelativeResidual > 1e-6 ||
        discreteNonTangentYoyo.derivativeMethod !=
            spin::quadricMomentDerivativeMethodToken(discreteSettings.derivativeMethod) ||
        discreteNonTangentYoyo.trace.size() < 2 ||
        (discreteNonTangentYoyo.quadricCoefficients -
            offTargetYoyoSeed.values()).norm() < 1e-8 ||
        discreteNonTangentYoyo.materialBoundary.faces.rows() == 0 ||
        discreteNonTangentYoyo.cavityBoundary.faces.rows() == 0 ||
        discreteNonTangentYoyo.outputHash != discreteNonTangentYoyoRepeat.outputHash ||
        discreteNonTangentYoyo.traceHash != discreteNonTangentYoyoRepeat.traceHash ||
        discreteNonTangentYoyo.eventHash != discreteNonTangentYoyoRepeat.eventHash ||
        discreteNonTangentYoyo.geometryHash != discreteNonTangentYoyoRepeat.geometryHash)
    {
        return fail("non-tangent discrete Boolean-to-Newton integration failed: status=" +
            std::string(spin::spinSolverStatusToken(discreteNonTangentYoyo.status)) +
            " residual=" + std::to_string(discreteNonTangentYoyo.normalizedEquationResidual) +
            " physical=" + std::to_string(discreteNonTangentYoyo.physicalRelativeResidual) +
            " trace=" + std::to_string(discreteNonTangentYoyo.trace.size()) +
            " diagnostic=" + discreteNonTangentYoyo.diagnostic);
    }
    if (evaluatedMapCount == 0 || !everyEvaluatedMapIsBoundToItsDifferential)
    {
        return fail("strict Newton evaluation events do not bind the actual Boolean backend and derivative identity to the evaluated moment map");
    }
    std::cout << "[PASS] discrete Boolean-to-Newton yoyo status="
              << spin::spinSolverStatusToken(discreteNonTangentYoyo.status)
              << " iterations=" << discreteNonTangentYoyo.trace.size()
              << " residual=" << discreteNonTangentYoyo.normalizedEquationResidual
              << " derivative=" << discreteNonTangentYoyo.derivativeMethod << '\n';

    hashes.emplace_back();
    if (!solveAndCheck(
            "yo-yo",
            yoyoRequest,
            yoyoHyperboloid,
            spin::QuadricSurfaceType::Hyperboloid,
            &hashes.back()
        ))
    {
        return fail("yo-yo symmetric fixture failed");
    }

    values(3) = -1.0;
    const Quadric topHyperboloid = Quadric::fromUnnormalized(values);
    spin::SpinDesignRequest topRequest =
        baseRequest(spin::SpinMode::SpinningTop);
    topRequest.centerOfMassWeight = 0.6;
    topRequest = requestWithTarget(
        topRequest,
        topHyperboloid,
        rawMoments(
            20.0,
            Eigen::Vector3d::Zero(),
            Eigen::Vector3d(100.0, 100.0, 200.0)
        )
    );
    hashes.emplace_back();
    if (!solveAndCheck(
            "top",
            topRequest,
            topHyperboloid,
            spin::QuadricSurfaceType::Hyperboloid,
            &hashes.back()
        ))
    {
        return fail("spinning-top hyperboloid fixture failed");
    }

    std::cout << "[PASS] scenario=spin-it-faster-end-to-end"
              << " off_target_iterations=" << offTargetResult.trace.size()
              << " off_target_seed_delta="
              << (offTargetResult.quadricCoefficients -
                  offTargetSeed.values()).norm()
              << " off_target_residual="
              << offTargetResult.normalizedEquationResidual
              << " off_target_output_hash=" << offTargetResult.outputHash
              << " off_target_trace_hash=" << offTargetResult.traceHash
              << " off_target_hyperboloid_iterations="
              << offTargetYoyo.trace.size()
              << " off_target_hyperboloid_residual="
              << offTargetYoyo.normalizedEquationResidual
              << " off_target_hyperboloid_derivative_error="
              << derivativeErrors[1];
    for (const std::string& hash : hashes)
    {
        std::cout << ' ' << hash;
    }
    std::cout << '\n';
    return EXIT_SUCCESS;
}
